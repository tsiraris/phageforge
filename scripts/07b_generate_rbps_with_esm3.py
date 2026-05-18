"""Stage 07b: Generate family-constrained RBP candidates with local ESM3 or Forge-backed ESM3.

This script reads the stage 07 context JSON and generates mutated sequence candidates. 
It supports two primary generation backends:

    - Forge API: Calls the frontier EvolutionaryScale API to perform scaffold-constrained masked editing.
    - Local ESM3: Runs open-weights ESM3 models locally if the API isn't used.

Note: In this upgraded pure-ESM3 path, legacy ESM2 fallbacks are intentionally disabled to ensure structural and generative quality.

It generates multiple candidates and strictly logs the provenance (generation status, temperature, random seeds, specific mutations made) into a CSV file (all_generated_candidates.csv).
"""

from __future__ import annotations
import argparse
import os
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from phageforge.stage07_utils import (
    Stage07Regime,
    candidate_guidance_score,
    choose_hotspots,
    make_masked_prompt,
    mutation_list,
    mutation_penalty,
    parse_regimes_json,
    read_json,
    seed_everything,
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for Stage 07 ESM3 generation."""
    ap = argparse.ArgumentParser(description="Generate Stage 07 RBP candidates with ESM3.")                                                                                               # Initializes the argument parser with a description
    ap.add_argument("--context_json", type=str, required=True, help="Stage 07 context JSON produced by 07a_prepare_stage07_design_context.py.")                                           # Adds argument for the input context JSON path
    ap.add_argument("--out_csv", type=str, required=True, help="Where to write the generated-candidate CSV.")                                                                             # Adds argument for the output CSV path
    ap.add_argument("--n_samples", type=int, default=12, help="Total number of candidate rows to request across all regimes.")                                                            # Adds argument for total sequences to generate
    ap.add_argument("--temperature", type=float, default=0.7, help="Default generation temperature when explicit regimes are not provided.")                                              # Adds argument for generation randomness factor
    ap.add_argument("--top_k", type=int, default=5, help="Default top-k metadata value when explicit regimes are not provided.")                                                          # Adds argument for vocabulary limiting during gen
    ap.add_argument("--max_esm3_masked_positions", type=int, default=24, help="Maximum number of masked sequence positions per prompt.")                                                  # Adds argument for max simultaneous mutations
    ap.add_argument("--sampling_seed", type=int, default=42, help="Base random seed for hotspot choice and generation reproducibility.")                                                  # Adds argument for reproducible random seeding
    ap.add_argument("--esm3_backend", type=str, default="local", choices=["local", "forge"], help="Run local esm3-open weights or call the Forge API.")                                   # Adds argument to choose local model or remote API
    ap.add_argument("--esm3_model", type=str, default="esm3-open", help="Model label for local or forge generation.")                                                                     # Adds argument to specify the exact model name
    ap.add_argument("--esm3_num_steps", type=int, default=8, help="Number of ESM3 iterative unmasking steps.")                                                                            # Adds argument for generative unmasking iterations
    ap.add_argument("--esm3_error_fallback", type=str, default="none", choices=["none", "esm2"], help="Optional fallback; 'none' keeps the run pure ESM3.")                               # Adds argument for fallback logic if model fails
    ap.add_argument("--regimes_json", type=str, default=None, help="Optional JSON list describing multiple generation regimes.")                                                          # Adds argument for advanced multi-strategy configs
    ap.add_argument("--hotspot_strategy", type=str, default="mixed", choices=["even", "priority", "mixed"], help="Fallback hotspot selection strategy when regimes_json is omitted.")     # Adds argument for how to pick mutation sites
    ap.add_argument("--max_attempts_per_sample", type=int, default=3, help="How many guided ESM3 attempts to try before keeping the best candidate.")                                     # Adds argument for retry loops to maximize score
    return ap.parse_args()                                                                                                                                                                # Parses the command line and returns the namespace


def load_esm3_client(backend: str, model_name: str):
    """Instantiate a local or Forge ESM3 client using the public ESM python interface."""
    import esm
    import torch
    from esm.models.esm3 import ESM3

    if backend == "local":                                                                                                                        # Checks if the user requested a local model
        client = ESM3.from_pretrained(model_name)                                                                                                 # Loads the model weights from disk into memory
        if hasattr(client, "to"):                                                                                                                 # Checks if the model object supports device mapping
            client = client.to("cuda" if torch.cuda.is_available() else "cpu")                                                                    # Moves model to GPU if available, else leaves on CPU
        return client                                                                                                                             # Returns the ready-to-use local model

    token = os.environ.get("ESM_API_KEY") or os.environ.get("FORGE_API_TOKEN") or os.environ.get("ESM3_FORGE_TOKEN")                              # Attempts to pull auth tokens from env variables
    if not token:                                                                                                                                 # Checks if token retrieval failed
        raise RuntimeError("Forge generation was requested but no Forge token was found in ESM_API_KEY / FORGE_API_TOKEN / ESM3_FORGE_TOKEN.")    # Crashes the script if API requested without auth
    return esm.sdk.client(model_name, token=token)                                                                                                # Authenticates and returns the remote API client


def generate_sequence(client, prompt_sequence: str, num_steps: int, temperature: float, top_k: int, sample_seed: int) -> str:
    """Generate a sequence completion from an ESM3 prompt string."""
    import torch
    from esm.sdk.api import ESMProtein, GenerationConfig

    # Apply random seed to ensure reproducibility
    seed_everything(sample_seed)                                                                                               # Applies the seed to python's random and numpy
    torch.manual_seed(sample_seed)                                                                                             # Applies the seed to PyTorch CPU operations
    if torch.cuda.is_available():                                                                                              # Checks if GPU is being used
        torch.cuda.manual_seed_all(sample_seed)                                                                                # Applies the seed to PyTorch GPU operations

    # Generate the sequence and return it
    protein = ESMProtein(sequence=prompt_sequence)                                                                             # Wraps the raw string into an ESM protein object
    cfg = GenerationConfig(track="sequence", num_steps=num_steps, temperature=float(temperature))                              # Configures the generation track and parameters
    if hasattr(cfg, "top_k"):                                                                                                  # Checks which top_k attribute name the library uses
        cfg.top_k = int(top_k)                                                                                                 # Assigns the parameter to top_k
    elif hasattr(cfg, "topk"):                                                                                                 # Checks for alternative top_k attribute spelling
        cfg.topk = int(top_k)                                                                                                  # Assigns the parameter to topk
    out = client.generate(protein, cfg)                                                                                        # Executes the inference generation call
    sequence = getattr(out, "sequence", None)                                                                                  # Safely extracts the sequence property from output
    if not sequence:                                                                                                           # Validates that an actual string was returned
        raise RuntimeError("ESM3 returned no sequence output.")                                                                # Crashes if the model output was empty
    return str(sequence)                                                                                                       # Returns the generated sequence as a standard string


def _attempt_generation(client, seed_sequence: str, prompt_sequence: str, masked_positions: list[int], regime: Stage07Regime, args: argparse.Namespace, position_features: list[dict], sample_seed: int) -> tuple[str, str | None, int, float]:
    """ Run a few (max_attempts) guided ESM3 attempts and keep the best score successful candidate.
        Returns the best sequence, error, score, and seed."""
    
    best_sequence = ""                                                                                                                                   # Initializes tracking variable for top sequence
    best_error = None                                                                                                                                    # Initializes tracking variable for generation errors
    best_score = float("-inf")                                                                                                                           # Initializes top score tracker to lowest possible
    best_seed = sample_seed                                                                                                                              # Initializes tracker for the most successful seed

    for attempt_idx in range(max(1, int(args.max_attempts_per_sample))):                                                                                 # Loops over allowed retries ensuring at least 1
        attempt_seed = int(sample_seed + attempt_idx)                                                                                                    # Creates a unique, reproducible seed per attempt
        try:                                                                                                                                             # Starts a block to catch inference failures gracefully
            candidate_sequence = generate_sequence(                                                                                                      # Calls the helper function to run the model
                client=client,                                                                                                                           # Passes the instantiated client
                prompt_sequence=prompt_sequence,                                                                                                         # Passes the template with blanks
                num_steps=regime.num_steps,                                                                                                              # Passes how many steps the model should take
                temperature=regime.temperature,                                                                                                          # Passes randomness level
                top_k=regime.top_k,                                                                                                                      # Passes token restriction limit
                sample_seed=attempt_seed,                                                                                                                # Passes the attempt-specific seed
            )                                                                                                                                            # Closes generation call
            # Identifies which spots actually changed
            mutated_positions = [pos for pos in masked_positions if pos <= len(candidate_sequence) and candidate_sequence[pos - 1] != seed_sequence[pos - 1]] 
            # Every time a sequence is generated, it calculates a total_score, 
            # rewarding the model for making mutations in the correct functional spots (guidance_score) 
            guidance_score = candidate_guidance_score(seed_sequence, candidate_sequence, position_features, mutated_positions)                           # Computes a biological logic score for the variant
            # plus a small mathematical reward for making novel changes (novelty_bonus).
            novelty_bonus = 0.005 * len(mutated_positions)                                                                                               # Calculates a small mathematical reward for mutating
            total_score = guidance_score + novelty_bonus                                                                                                 # Combines the logic score and the mutation reward
            if total_score > best_score:                                                                                                                 # Checks if this attempt is the best one so far
                best_sequence = candidate_sequence                                                                                                       # Updates the tracking variable to the new leader
                best_error = None                                                                                                                        # Clears error memory since it succeeded
                best_score = total_score                                                                                                                 # Updates the highest score tracker
                best_seed = attempt_seed                                                                                                                 # Records which seed produced the winning sequence
        except Exception as exc:                                                                                                                         # Catches any crashes during the generate sequence
            best_error = str(exc)                                                                                                                        # Saves the error message for logging if it's the last try

    return best_sequence, best_error, best_seed, best_score                                                                                              # Returns a tuple of all top candidate details


def build_rows(seed_sequence: str, target_host: str, family_label: str, regimes: list[Stage07Regime], context: dict, client, args: argparse.Namespace) -> list[dict]:
    """ Generate candidate rows across regimes (temperatures, topks, num_steps, hotspot strategy), deduplicate sequences, and retain provenance. """
    hotspot_positions = context["editable_region"].get("hotspot_positions", [])                                                                                          # Extracts allowed edit sites from context
    hotspot_weights = context["editable_region"].get("hotspot_priority_weights", {})                                                                                     # Extracts custom probabilities for those sites
    position_features = context["editable_region"].get("position_features", [])                                                                                          # Extracts feature data used for scoring later
    per_regime = max(1, (args.n_samples + len(regimes) - 1) // len(regimes))                                                                                             # Calculates even distribution of samples per strategy

    rows = []                                                                                                                                                            # Initializes master list for all output rows
    seen_sequences = set()                                                                                                                                               # Initializes a set to quickly check for duplicates
    sample_id = 0   
    # Loops through the different generation "regimes"
    for regime_rank, regime in enumerate(regimes):                                                                                                                       # Iterates through generation strategies with ranking
        for local_idx in tqdm(range(per_regime), total=per_regime):                                                                                                      # Iterates to produce N samples per strategy with UI bar
            if len(rows) >= args.n_samples:                                                                                                                              # Checks if global sample quota has been met
                break                                                                                                                                                    # Escapes the loop early if quota is met
            base_seed = int(args.sampling_seed + 1000 * regime_rank + 100 * local_idx)                                                                                   # Computes highly unique seed based on loop state
            # For each regime, it calculates which hotspots to hide
            masked_positions = choose_hotspots(                                                                                                                          # Calls utility to select which amino acids to hide
                hotspots_1based=hotspot_positions,                                                                                                                       # Provides total pool of editable sites
                priority_weights=hotspot_weights,                                                                                                                        # Provides the bias weighting
                max_positions=regime.max_masked_positions,                                                                                                               # Caps the maximum number of edits allowed
                strategy=regime.hotspot_strategy,                                                                                                                        # Uses the strategy defined by current regime
                sample_seed=base_seed,                                                                                                                                   # Uses the unique loop seed for reproducibility
            ) 
            # Creates the prompt string sequence with masked positions for ESM3
            prompt_sequence = make_masked_prompt(seed_sequence, masked_positions)                                                                                        # Replaces chosen sites in string with mask characters
            try:                                                                                                                                                         # Begins fault-tolerant generation block
                # and sends it to _attempt_generation to return the best sequence of a few attempts
                candidate_sequence, best_error, chosen_seed, guidance_score = _attempt_generation(                                                                       # Calls multi-try model function
                    client=client,                                                                                                                                       # Passes the ESM3 connection
                    seed_sequence=seed_sequence,                                                                                                                         # Passes wildtype sequence
                    prompt_sequence=prompt_sequence,                                                                                                                     # Passes masked template string
                    masked_positions=masked_positions,                                                                                                                   # Passes list of hidden indices
                    regime=regime,                                                                                                                                       # Passes current generation settings
                    args=args,                                                                                                                                           # Passes global config arguments
                    position_features=position_features,                                                                                                                 # Passes scoring data
                    sample_seed=base_seed,                                                                                                                               # Passes the baseline loop seed
                )                                                                                                                                                        # Closes multi-try model function call
                if not candidate_sequence:                                                                                                                               # Verifies something was actually returned
                    raise RuntimeError(best_error or "ESM3 generation did not return a candidate sequence.")                                                             # Forces failure path if sequence is empty
                status = "ok"                                                                                                                                            # Marks status as successful for CSV logging
                error_text = ""                                                                                                                                          # Clears error column string
                generator_mode = f"esm3_{args.esm3_backend}:{args.esm3_model}"                                                                                           # Formats the model provenance string
            except Exception as exc:                                                                                                                                     # Catches any error that surfaced from the try block
                if args.esm3_error_fallback != "none":                                                                                                                   # Checks if a fallback architecture is permitted
                    raise RuntimeError("ESM2 fallback is intentionally not implemented in the upgraded pure-ESM3 Stage 07 path.") from exc                               # Refuses to run fallback, preserving pipeline integrity
                candidate_sequence = ""                                                                                                                                  # Defaults to empty string on failure
                chosen_seed = base_seed                                                                                                                                  # Logs baseline seed regardless of failure
                guidance_score = float("-inf")                                                                                                                           # Defaults score to worst possible outcome
                status = "error"                                                                                                                                         # Marks status as failure for CSV logging
                error_text = str(exc)                                                                                                                                    # Records exact error string for CSV logging
                generator_mode = f"esm3_{args.esm3_backend}:{args.esm3_model}"                                                                                           # Formats the model provenance string

            # If sequence was successfully made but already exists
            if status == "ok" and candidate_sequence in seen_sequences:                                                                                                  
                continue                                                                                                                                                 # Skips appending this to CSV to maintain diversity
            if status == "ok":                                                                                                                                           # Confirms candidate is both successful and novel
                seen_sequences.add(candidate_sequence)                                                                                                                   # Adds to uniqueness tracker to block future duplicates

            # Builds a massive dictionary for the current sequence that contains all the metadata needed
            row = {                                                                                                                                                      # Begins constructing dictionary for current sequence
                "sample_id": sample_id,                                                                                                                                  # Saves the iterative ID
                "generator_mode": generator_mode,                                                                                                                        # Saves the model hardware/architecture info
                "generation_regime": regime.name,                                                                                                                        # Saves the strategy name that generated this
                "generation_regime_rank": regime_rank,                                                                                                                   # Saves the hierarchy rank of the strategy
                "generation_status": status,                                                                                                                             # Saves success or fail status
                "generation_error": error_text or None,                                                                                                                  # Saves error message or null
                "target_host": target_host,                                                                                                                              # Saves the intended bacterial target
                "family_id": family_label,                                                                                                                               # Saves the protein family classification
                "seed_sequence": seed_sequence,                                                                                                                          # Saves the parent baseline sequence
                "candidate_sequence": candidate_sequence if candidate_sequence else seed_sequence,                                                                       # Saves the new variant, reverting to parent on error
                "editable_hotspots": ",".join(str(pos) for pos in masked_positions),                                                                                     # Saves comma-separated list of selected mutation targets
                "editable_hotspot_count": int(len(masked_positions)),                                                                                                    # Saves count of mutation targets
                "mutation_positions": ";".join(mutation_list(seed_sequence, candidate_sequence if candidate_sequence else seed_sequence)),                               # Saves explicitly computed mutation differences
                "mutation_penalty": int(mutation_penalty(seed_sequence, candidate_sequence if candidate_sequence else seed_sequence)),                                   # Saves penalty score for excessive changes
                "esm3_prompt_sequence": prompt_sequence,                                                                                                                 # Saves the exact input given to ESM3
                "esm3_hotspot_positions": ",".join(str(pos) for pos in masked_positions),                                                                                # Redundantly logs targets specific to ESM3 config
                "esm3_temperature": float(regime.temperature),                                                                                                           # Saves the randomness factor used
                "esm3_top_k": int(regime.top_k),                                                                                                                         # Saves the vocabulary limit used
                "esm3_sampling_seed": int(chosen_seed),                                                                                                                  # Saves the precise seed that yielded the result
                "esm3_model": args.esm3_model,                                                                                                                           # Saves the explicit model version name
                "esm3_num_steps": int(regime.num_steps),                                                                                                                 # Saves the number of generation steps used
                "hotspot_strategy": regime.hotspot_strategy,                                                                                                             # Saves the strategy used to select mutation sites
                "guided_mutation_score": float(guidance_score) if np.isfinite(guidance_score) else None,                                                                 # Saves the biological logic score if valid
                "used_esm3_api": bool(args.esm3_backend == "forge"),                                                                                                     # Flags if a remote API was charged
                "used_local_generator": bool(args.esm3_backend == "local"),                                                                                              # Flags if inference was local
                "used_local_esm3": bool(args.esm3_backend == "local"),                                                                                                   # Redundantly flags local usage
                "used_esm2_fallback": False,                                                                                                                             # Hardcodes False since fallback is disabled
            }                                                                                                                                                            # Closes the dictionary construction
            rows.append(row)                                                                                                                                             # Appends the completed dictionary to the master list
            sample_id += 1                                                                                                                                               # Increments the ID for the next sequence
    return rows                                                                                                                                                          # Returns the full list of generated rows


def main() -> None:
    # Parse command-line arguments, load context, and build generation strategies
    args = parse_args()                                                                                                                                                                        # Retrieves parsed command-line arguments
    context = read_json(args.context_json)                                                                                                                                                     # Loads the stage 07 design context from the JSON file
    # Parses or builds the generation strategies/regimes (temperatures, topks, num_steps, hotspot strategy)
    regimes = parse_regimes_json(args.regimes_json, args.temperature, args.top_k, args.max_esm3_masked_positions, args.esm3_num_steps)                                                         
    if args.regimes_json is None:                                                                                                                                                              # Checks if specific strategies were omitted
        for regime in regimes:                                                                                                                                                                 # Iterates through the default generated strategies
            regime.hotspot_strategy = args.hotspot_strategy if regime.name == "balanced" else regime.hotspot_strategy                                                                          # Overrides default hotspot strategies with the user arg
    client = load_esm3_client(args.esm3_backend, args.esm3_model)                                                                                                                              # Instantiates the ESM3 client based on backend

    # Load seed sequence, target host, and family label
    seed_sequence = str(context["selected_seed"]["seed_sequence"])                                                                                                                             # Extracts the baseline sequence string from context
    target_host = str(context["target_host"])                                                                                                                                                  # Extracts the target bacterial species from context
    family_label = str(context["family_context"].get("family_product_majority", "receptor-binding protein"))                                                                                   # Extracts the functional family label from context
    # Generates all candidate sequences based on inputs: seed_sequence, target_host, family_label, regimes, context.json
    # rows is a massive dictionary containing not only the sequences for each regime but also all the metadata needed
    rows = build_rows(seed_sequence, target_host, family_label, regimes, context, client, args)                                                                                                
    if not rows:                                                                                                                                                                               # Validates that at least one sequence was generated
        raise RuntimeError("ESM3 generation did not produce any rows.")                                                                                                                        # Crashes if generation completely failed

    # Write results to CSV
    out_path = Path(args.out_csv)                                                                                                                                                              # Converts the output string path to a Path object
    out_path.parent.mkdir(parents=True, exist_ok=True)                                                                                                                                         # Creates necessary parent directories for the output file
    pd.DataFrame(rows).sort_values(["generation_status", "guided_mutation_score", "mutation_penalty", "sample_id"], ascending=[True, False, True, True]).to_csv(out_path, index=False)         # Converts to DataFrame, sorts by quality, and writes to CSV
    print(f"Wrote: {out_path}")                                                                                                                                                                # Prints confirmation of the written file path
    print(f"used_backend: {args.esm3_backend}")                                                                                                                                                # Prints which backend architecture was used
    print(f"rows_written: {len(rows)}")                                                                                                                                                        # Prints the total number of successful candidates saved


if __name__ == "__main__":
    main()