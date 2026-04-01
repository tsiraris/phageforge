"""
====================================================================================================
Stage 07b: Practical ESM3-aware masked editing over the validated Stage 07 design window.
====================================================================================================

This script reads the Stage 07 context, chooses a generator backend, generates candidates, and writes them to CSV.

It supports three generation backends:

1. ESM3 Forge API generation (with a Forge API key).
2. A fallback local masked conditional generator checkpoint (`07b_train_conditional_generator.py`).
3. ESM2 masked language model fallback (when neither of the above is used).

The script will write the generated candidates to a CSV file.
"""

from __future__ import annotations                                               # Enable postponed annotation evaluation for cleaner typing.
import argparse                                                                  # Parse command-line arguments.
import json                                                                      # Read the Stage 07 context JSON file.
import os                                                                        # Read the ESM_API_KEY environment variable.
import random                                                                    # Seed Python's RNG for more reproducible candidate generation.
from pathlib import Path                                                         # Handle filesystem paths cleanly and portably.
import pandas as pd                                                              # Store generated candidates in a DataFrame and write CSV output.
import torch                                                                     # Run the local generator checkpoint and ESM2 fallback inference.
from transformers import AutoTokenizer, EsmForMaskedLM                           # Use an ESM2 masked language model as the final fallback editor.
from phageforge.generation.dataset import VOCAB                                  # Reuse the shared generator vocabulary size.
from phageforge.generation.generator_model import ConditionalMaskedRBPGenerator  # Load the practical local Stage 07 generator.
from phageforge.generation.guidance import (                                     # Reuse shared helper utilities from the generation package.
    decode_ids,
    encode_seed,
    hotspot_masked_ids,
    mutation_burden,
    sample_hotspot_edits,
)


def parse_args() -> argparse.Namespace:
    """Return command-line arguments for Stage 07 candidate generation."""
    ap = argparse.ArgumentParser(description="Generate Stage 07 candidate RBP edits.")                    # Create the generation CLI parser.
    ap.add_argument("--context_json", type=str, required=True, help="Context JSON from 07a_prepare_stage07_design_context.py.")  # Point to the prepared Stage 07 context.
    ap.add_argument("--out_csv", type=str, required=True, help="Where to write generated candidates.")   # Point to the CSV that will hold sampled candidates.
    ap.add_argument("--n_samples", type=int, default=128, help="How many candidates to sample.")         # Control how many candidates to generate.
    ap.add_argument("--temperature", type=float, default=0.7, help="Sampling temperature for token selection.")  # Control exploration during sampling.
    ap.add_argument("--top_k", type=int, default=5, help="Top-k token sampling inside hotspot positions for local/ESM2 fallback.")  # Restrict local token sampling to top-k choices.
    ap.add_argument("--generator_checkpoint", type=str, default="", help="Optional local checkpoint from 07b_train_conditional_generator.py.")  # Enable the local generator mode.
    ap.add_argument("--esm_model", type=str, default="facebook/esm2_t33_650M_UR50D", help="Fallback ESM2 masked-LM editor model.")  # Choose the fallback ESM2 masked-LM model.
    ap.add_argument("--use_esm3_api", action="store_true", help="Use the ESM3 Forge API instead of the local generator or ESM2 fallback.")      # Turn on real ESM3 remote generation explicitly.
    ap.add_argument("--esm3_model", type=str, default="esm3-medium-2024-08", help="Forge ESM3 model name to request.")                           # Choose the Forge model variant you have access to.
    ap.add_argument("--esm3_num_steps", type=int, default=8, help="Number of iterative unmasking steps for ESM3 generation.")                   # Control how many iterative decode steps ESM3 runs.
    ap.add_argument("--sampling_seed", type=int, default=42, help="Random seed used for reproducible local / ESM2 sampling and logged as provenance.")  # Expose the seed so repeated runs can be reproduced more easily.
    return ap.parse_args()                                                                                      # Parse the CLI and return the resulting namespace.


def load_context(path: Path) -> dict:
    """Load the Stage 07 context JSON, containing all the information needed for candidate generation."""
    with open(path, "r", encoding="utf-8") as handle:                                                           # Open the JSON file using UTF-8.
        return json.load(handle)                                                                                # Parse and return the context dictionary.


def resolve_hotspots(context: dict, seed_seq: str) -> list[int]:
    """
    Return the editable protein positions as a list of 0-based residue indices.

    If the context JSON has hotspot positions, it returns them directly as a list.
    Otherwise, it uses the window start and end to expand the window into a list of explicit positions.
    """
    hotspot_positions = list(context["editable_region"].get("hotspot_positions", []))                          # Read any explicit hotspot list from the context JSON.
    if hotspot_positions:                                                                                       # Use the explicit hotspot list when present.
        return [int(x) for x in hotspot_positions]                                                              # Normalize hotspot positions to integers.
    start = int(context["editable_region"]["window_start"])                                                     # Read the editable window start when hotspots are missing.
    end = int(context["editable_region"]["window_end"])                                                         # Read the editable window end when hotspots are missing.
    return list(range(start, min(end, len(seed_seq))))                                                          # Expand the window into explicit residue positions.


def build_esm3_masked_prompt(seed_seq: str, hotspot_positions: list[int]) -> str:
    """
    Convert a seed amino-acid sequence into an ESM3 partial prompt.

    ESM3 expects masked sequence positions to be written as underscores.
    Example:
        seed_seq = "ACDEFG"
        hotspot_positions = [2, 4]
        output = "AC_E_G"

    This lets ESM3 keep the unmasked scaffold fixed while it fills only the editable sites.
    """
    chars = list(str(seed_seq).strip().upper())                                                                 # Convert the sequence string into a mutable list of residues.
    for pos in hotspot_positions:                                                                               # Visit each editable residue position.
        if 0 <= pos < len(chars):                                                                               # Ignore malformed positions outside the sequence.
            chars[pos] = "_"                                                                                    # Use ESM3's partial-sequence mask token for that residue.
    return "".join(chars)                                                                                       # Rejoin the list into one masked sequence prompt.


def set_generation_seed(seed: int) -> None:
    """Seed Python and PyTorch RNGs so local stochastic generation is easier to reproduce."""
    random.seed(seed)                                                                                           # Seed Python's built-in random module.
    torch.manual_seed(seed)                                                                                     # Seed PyTorch on CPU.
    if torch.cuda.is_available():                                                                               # Also seed CUDA when a GPU is present.
        torch.cuda.manual_seed_all(seed)                                                                        # Seed all CUDA devices for consistency.


def describe_mutations(seed_seq: str, candidate_seq: str) -> str:
    """Return a compact human-readable mutation summary as a string.
    For example: `3:E→Y;4:F→W` corresponds to E replaced by Y at position 3, F replaced by W at position 4."""
    edits = []                                                                                                  # Collect one textual mutation description per changed residue.
    max_len = min(len(seed_seq), len(candidate_seq))                                                            # Compare only over the aligned shared prefix length.
    for pos in range(max_len):                                                                                  # Visit each aligned residue position in order.
        src = seed_seq[pos]                                                                                     # Read the residue from the original seed.
        dst = candidate_seq[pos]                                                                                # Read the residue from the generated candidate.
        if src != dst:                                                                                          # Keep only true amino-acid substitutions.
            edits.append(f"{pos}:{src}→{dst}")                                                                  # Record one 0-based mutation event in compact notation.
    if len(candidate_seq) > max_len:                                                                            # Record any appended tail beyond the shared prefix.
        edits.append(f"append:{candidate_seq[max_len:]}")                                                       # Describe extra generated suffix residues when present.
    if len(seed_seq) > max_len:                                                                                 # Record any deleted tail if the candidate is shorter.
        edits.append(f"truncate:{seed_seq[max_len:]}")                                                          # Describe missing suffix residues when present.
    return ";".join(edits) if edits else "WT"                                                                   # Return `WT` when no residues changed.

# This helper builds one consistent output row, independent of the generation backend used.
def build_generation_record(
    *,
    sample_id: int,
    generator_mode: str,
    target_host: str,
    family_id: str,
    seed_seq: str,
    candidate_seq: str,
    hotspot_positions: list[int],
    masked_prompt: str,
    temperature: float,
    top_k: int | None,
    sampling_seed: int,
    esm3_model: str | None = None,
    esm3_num_steps: int | None = None,
) -> dict:
    """Build one fully annotated candidate row with reusable provenance fields."""
    return {
        "sample_id": sample_id,                                                                              # Store the integer sample index.
        "generator_mode": generator_mode,                                                                    # Record the high-level generation backend used.
        "target_host": target_host,                                                                          # Store the intended target host label.
        "family_id": family_id,                                                                              # Store the scaffold family label used for design context.
        "seed_sequence": seed_seq,                                                                           # Keep the original seed sequence for comparison.
        "candidate_sequence": candidate_seq,                                                                 # Store the generated candidate amino-acid sequence.
        "editable_hotspots": ",".join(str(x) for x in hotspot_positions),                                   # Store editable positions in a human-readable comma-separated form.
        "editable_hotspot_count": int(len(hotspot_positions)),                                               # Store how many sites were allowed to mutate.
        "mutation_positions": describe_mutations(seed_seq, candidate_seq),                                   # Store an explicit mutation summary for interpretation and debugging.
        "mutation_penalty": mutation_burden(seed_seq, candidate_seq),                                        # Store the mutation burden relative to the seed.
        "esm3_prompt_sequence": masked_prompt,                                                               # Persist the masked partial prompt that was actually generated from.
        "esm3_hotspot_positions": ",".join(str(x) for x in hotspot_positions),                              # Duplicate hotspots under an ESM3-specific provenance label for downstream filtering.
        "esm3_temperature": float(temperature),                                                              # Log the generation temperature for reproducibility.
        "esm3_top_k": int(top_k) if top_k is not None else pd.NA,                                            # Log top-k when applicable; leave missing for ESM3 API mode.
        "esm3_sampling_seed": int(sampling_seed),                                                            # Log the random seed used for local stochastic generation.
        "esm3_model": esm3_model if esm3_model is not None else pd.NA,                                       # Log the Forge model name only when ESM3 API mode was used.
        "esm3_num_steps": int(esm3_num_steps) if esm3_num_steps is not None else pd.NA,                     # Log Forge iterative decode steps only when ESM3 API mode was used.
    }


def generate_with_esm3_api(
    context: dict,
    model_name: str,
    n_samples: int,
    temperature: float,
    num_steps: int,
    sampling_seed: int,
) -> pd.DataFrame:
    """
    Generate candidates with the remote ESM3 Forge API.

    Workflow:
    1. Read the validated seed scaffold from the Stage 07 context.
    2. Replace only editable hotspot residues with underscores.
    3. Send the partially specified sequence prompt to ESM3 Forge.
    4. Ask ESM3 to iteratively fill the sequence track.
    5. Save each generated completion as one candidate.

    Notes:
    - This mode uses the real ESM3 client and requires ESM_API_KEY in the environment.
    - We keep host/family metadata in the output rows even though this script's practical ESM3 prompt is sequence-only.
    - Downstream scoring still uses the host-transfer and family-validity stack.
    """
    # Import ESM3 Forge pieces so the script still runs when only fallback modes are used.
    try:
        from esm.sdk import client as esm_client                                                          # Import the official Forge client constructor only when needed.
        from esm.sdk.api import ESMProtein, GenerationConfig                                              # Import the protein prompt container and generation config only when needed.
    except Exception as exc:                                                                              # Catch import problems and raise a helpful message.
        raise ImportError(
            "ESM3 API mode requested, but the `esm` package is not available. "
            "Install it with `pip install esm`."
        ) from exc

    token = os.environ.get("ESM_API_KEY", "").strip()                                                     # Read the Forge API token from the environment.
    if not token:                                                                                          # Refuse to continue if the token is missing.
        raise RuntimeError(
            "ESM_API_KEY is not set. In Windows CMD use: set ESM_API_KEY=your_key_here"
        )

    # Create the remote Forge client exactly once, then reuse it for all requested samples.
    model = esm_client(model_name, token=token)                                                           # Create a Forge-backed ESM3 inference client.

    # Read the validated seed and determine which positions are allowed to change.
    seed_seq = context["selected_seed"]["seed_sequence"]                                                  # Read the amino-acid seed sequence from the context JSON.
    target_host = context["target_host"]                                                                   # Read the requested target host from the context JSON.
    family_id = context["family_context"].get("family_product_majority") or "unknown_family"             # Recover the scaffold family label when available.
    hotspot_positions = resolve_hotspots(context, seed_seq)                                               # Determine which residue positions are editable.
    masked_prompt = build_esm3_masked_prompt(seed_seq, hotspot_positions)                                 # Convert the seed into an ESM3 partial masked prompt.
    rows = []                                                                                             # Collect one dictionary per generated candidate.

    # Sample multiple completions from the same masked scaffold prompt.
    for sample_id in range(n_samples):                                                                    # Loop once per requested ESM3 sample.
        protein = ESMProtein(sequence=masked_prompt)                                                      # Build the ESM3 prompt object from the masked partial sequence.
        generated = model.generate(                                                                       # Ask Forge ESM3 to fill the masked sequence track.
            protein,
            GenerationConfig(
                track="sequence",                                                                         # Tell ESM3 to generate amino-acid sequence tokens.
                num_steps=int(num_steps),                                                                  # Control the number of iterative unmasking steps.
                temperature=float(temperature),                                                            # Control generation diversity through temperature.
            ),
        )
        candidate_seq = str(generated.sequence).replace(" ", "")                                           # Extract the generated sequence and remove any spacing artifacts.
        rows.append(                                                                                       # Append one fully annotated provenance row.
            build_generation_record(
                sample_id=sample_id,
                generator_mode=f"esm3_api:{model_name}",
                target_host=target_host,
                family_id=family_id,
                seed_seq=seed_seq,
                candidate_seq=candidate_seq,
                hotspot_positions=hotspot_positions,
                masked_prompt=masked_prompt,
                temperature=temperature,
                top_k=None,
                sampling_seed=sampling_seed,
                esm3_model=model_name,
                esm3_num_steps=num_steps,
            )
        )

    return pd.DataFrame(rows)                                                                               # Convert all candidate dictionaries into a DataFrame.


# Local fallback candidate generator (if the ESM3 model is not available/selected): a custom Transformer model trained locally.
def generate_with_local_generator(
    context: dict,
    checkpoint_path: Path,
    n_samples: int,
    temperature: float,
    top_k: int,
    sampling_seed: int,
) -> pd.DataFrame:
    """Generate candidates with the local conditional masked generator checkpoint (custom Transformer model trained in 07b_train_conditional_generator.py)."""
    # Load the generator checkpoint from disk, build the model, and switch to evaluation mode.
    ckpt = torch.load(checkpoint_path, map_location="cpu")                            # Load the generator checkpoint from disk onto CPU first.
    host_to_id = ckpt["host_to_id"]                                                   # Recover the host conditioning vocabulary.
    family_to_id = ckpt["family_to_id"]                                               # Recover the family conditioning vocabulary.
    device = "cuda" if torch.cuda.is_available() else "cpu"                           # Prefer GPU when available.
    model = ConditionalMaskedRBPGenerator(                                            # Rebuild the model architecture that matches the checkpoint.
        vocab_size=len(VOCAB),                                                        # Use the shared generator vocabulary size.
        n_hosts=len(host_to_id),                                                      # Set the number of host labels from the checkpoint.
        n_families=len(family_to_id),                                                 # Set the number of family labels from the checkpoint.
    ).to(device)                                                                      # Move the model to the selected device.
    model.load_state_dict(ckpt["model"])                                              # Load the trained model weights into the architecture.
    model.eval()                                                                      # Switch the model to evaluation mode.

    # Build the masked sequence with BOS/EOS markers and replace editable positions with mask tokens.
    seed_seq = context["selected_seed"]["seed_sequence"]                              # Read the amino-acid seed sequence from the context JSON.
    family_id = context["family_context"].get("family_product_majority") or "unknown_family"  # Prefer the majority family label from Stage 06.
    if family_id not in family_to_id:                                                 # Handle the case where the family label is absent from the checkpoint vocabulary.
        family_id = next(iter(family_to_id.keys()))                                   # Fall back to the first known family to keep the script runnable.
    target_host = context["target_host"]                                              # Read the requested target host from the context JSON.
    if target_host not in host_to_id:                                                 # Ensure the target host exists in the checkpoint vocabulary.
        raise ValueError(f"Target host {target_host} not found in generator host_to_id mapping.")

    hotspot_positions = resolve_hotspots(context, seed_seq)                           # Determine which residue positions are editable.
    base_ids = encode_seed(seed_seq)                                                  # Tokenize the seed sequence with BOS/EOS markers.
    masked_ids = hotspot_masked_ids(base_ids, hotspot_positions)                      # Replace editable positions with mask tokens.
    masked_prompt = build_esm3_masked_prompt(seed_seq, hotspot_positions)             # Build a human-readable masked prompt for provenance logging.
    rows = []                                                                         # Collect generated candidate rows in a list of dictionaries.

    host_id = torch.tensor([host_to_id[target_host]], dtype=torch.long, device=device)   # Build the single-example host conditioning tensor.
    fam_id = torch.tensor([family_to_id[family_id]], dtype=torch.long, device=device)    # Build the single-example family conditioning tensor.

    # Sample multiple edited candidates by repeatedly filling the same hotspot mask pattern.
    for sample_id in range(n_samples):                                                # Loop once per requested generated candidate.
        input_ids = torch.tensor([masked_ids], dtype=torch.long, device=device)       # Build a batch of one masked seed sequence.
        attention_mask = (input_ids != 0).long()                                      # Mark all non-PAD tokens as valid attention positions.
        with torch.no_grad():                                                         # Disable gradient tracking because this is inference only.
            logits = model(                                                           # Predict token distributions at every sequence position.
                input_ids=input_ids,                                                  # Provide the masked token IDs.
                attention_mask=attention_mask,                                        # Provide the valid-token attention mask.
                host_ids=host_id,                                                     # Provide the host condition label.
                family_ids=fam_id,                                                    # Provide the family condition label.
            )[0]                                                                      # Remove the batch dimension because batch size is one.
        sampled_ids = sample_hotspot_edits(                                           # Sample edits only at hotspot positions.
            logits=logits,                                                            # Use the predicted token logits.
            base_ids=base_ids,                                                        # Start sampling from the original seed token sequence.
            hotspot_positions=hotspot_positions,                                      # Restrict sampling to the editable positions.
            temperature=temperature,                                                  # Apply the requested temperature.
            top_k=top_k,                                                              # Apply the requested top-k restriction.
        )
        candidate_seq = decode_ids(sampled_ids)                                       # Convert sampled token IDs back to an amino-acid sequence.
        rows.append(                                                                  # Append one fully annotated provenance row.
            build_generation_record(
                sample_id=sample_id,
                generator_mode="local_conditional_generator",
                target_host=target_host,
                family_id=family_id,
                seed_seq=seed_seq,
                candidate_seq=candidate_seq,
                hotspot_positions=hotspot_positions,
                masked_prompt=masked_prompt,
                temperature=temperature,
                top_k=top_k,
                sampling_seed=sampling_seed,
            )
        )

    return pd.DataFrame(rows)                                                         # Convert all candidate dictionaries into a DataFrame.



def generate_with_esm_masked_lm(
    context: dict,
    model_name: str,
    n_samples: int,
    temperature: float,
    top_k: int,
    sampling_seed: int,
) -> pd.DataFrame:
    """Generate candidates with an ESM masked language model as a fallback editor."""
    device = "cuda" if torch.cuda.is_available() else "cpu"                           # Prefer GPU when available for faster ESM inference.
    tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)        # Load the ESM tokenizer for raw amino-acid strings.
    model = EsmForMaskedLM.from_pretrained(model_name).to(device).eval()              # Load the masked language model and switch it to eval mode.

    seed_seq = context["selected_seed"]["seed_sequence"]                              # Read the seed amino-acid sequence from the context JSON.
    target_host = context["target_host"]                                              # Read the requested target host for provenance logging.
    family_id = context["family_context"].get("family_product_majority") or "unknown_family"  # Recover the family label when available.
    hotspot_positions = resolve_hotspots(context, seed_seq)                           # Determine which residue positions should be edited.
    masked_prompt = build_esm3_masked_prompt(seed_seq, hotspot_positions)             # Build a human-readable masked prompt for provenance logging.
    rows = []                                                                         # Collect generated candidate rows.

    # Build a base string sequence where editable residues are replaced by tokenizer mask tokens.
    chars = list(seed_seq)                                                            # Convert the seed string into a mutable list of residues.
    for pos in hotspot_positions:                                                     # Visit each editable residue position.
        if 0 <= pos < len(chars):                                                     # Keep only positions that fall inside the sequence length.
            chars[pos] = tokenizer.mask_token                                         # Replace the editable residue with the ESM mask token string.
    masked_sequence = "".join(chars)                                                  # Join the list back into one masked protein string.

    # Sample one edited candidate at a time from the masked-LM logits.
    for sample_id in range(n_samples):                                                # Loop once per requested candidate.
        # Tokenize the masked seed sequence, compute the logits, and find the masked token positions.
        toks = tokenizer(masked_sequence, return_tensors="pt").to(device)             # Tokenize the masked sequence and move it to the chosen device.
        with torch.no_grad():                                                         # Disable gradient tracking during inference.
            logits = model(**toks).logits[0]                                          # Compute token logits and remove the batch dimension.
        sampled_ids = toks["input_ids"][0].clone()                                    # Start from the masked token IDs already produced by the tokenizer.

        mask_positions = (toks["input_ids"][0] == tokenizer.mask_token_id).nonzero(as_tuple=False).view(-1)  # Find every masked token position.

        # Sample one replacement token for each masked position.
        for mask_pos in mask_positions.tolist():                                      # Sample one replacement token for each masked position.
            row = logits[mask_pos] / max(float(temperature), 1e-6)                    # Apply temperature scaling to the current masked-token logits.
            probs = torch.softmax(row, dim=-1)                                        # Convert logits into probabilities.
            if 0 < top_k < probs.numel():                                             # Restrict sampling to top-k tokens when requested.
                top_vals, top_idx = torch.topk(probs, k=top_k)                        # Keep the top-k token candidates and their probabilities.
                top_vals = top_vals / top_vals.sum()                                  # Renormalize the truncated distribution.
                chosen = top_idx[torch.multinomial(top_vals, num_samples=1)[0]].item()# Sample one token from the top-k distribution.
            else:                                                                     # Otherwise sample from the full vocabulary distribution.
                chosen = torch.multinomial(probs, num_samples=1)[0].item()            # Draw one token ID from the full probability vector.
            sampled_ids[mask_pos] = int(chosen)                                       # Replace the mask token with the sampled token ID.

        # Decode the generated candidate sequence, append it to the output dataframe and compute its mutation burden.
        candidate_seq = tokenizer.decode(sampled_ids, skip_special_tokens=True).replace(" ", "")  # Decode the candidate and remove spacing artifacts.
        rows.append(                                                                  # Append one fully annotated provenance row.
            build_generation_record(
                sample_id=sample_id,
                generator_mode=f"esm_masked_lm:{model_name}",
                target_host=target_host,
                family_id=family_id,
                seed_seq=seed_seq,
                candidate_seq=candidate_seq,
                hotspot_positions=hotspot_positions,
                masked_prompt=masked_prompt,
                temperature=temperature,
                top_k=top_k,
                sampling_seed=sampling_seed,
            )
        )

    return pd.DataFrame(rows)                                                         # Convert the generated rows into a DataFrame.


# Main script: load the context, choose one generation backend, and write one candidate table with provenance.
def main() -> None:
    # Parse arguments, load the Stage 07 context, and fix the random seed before any stochastic generation happens.
    args = parse_args()                                                               # Parse command-line arguments.
    context = load_context(Path(args.context_json))                                   # Load the Stage 07 context JSON from disk.
    set_generation_seed(args.sampling_seed)                                           # Seed Python and PyTorch RNGs so stochastic generation is easier to reproduce.

    # Route generation to the requested backend, preferring ESM3 Forge when the user explicitly asks for it.
    if args.use_esm3_api:                                                             # Use real ESM3 generation when the API flag was supplied.
        generated = generate_with_esm3_api(                                           # Sample candidates with the remote Forge ESM3 model.
            context=context,
            model_name=args.esm3_model,
            n_samples=args.n_samples,
            temperature=args.temperature,
            num_steps=args.esm3_num_steps,
            sampling_seed=args.sampling_seed,
        )
    elif args.generator_checkpoint:                                                   # Otherwise prefer the local conditional generator if a checkpoint was provided.
        generated = generate_with_local_generator(                                    # Sample candidates with the local conditional generator.
            context=context,
            checkpoint_path=Path(args.generator_checkpoint),
            n_samples=args.n_samples,
            temperature=args.temperature,
            top_k=args.top_k,
            sampling_seed=args.sampling_seed,
        )
    else:                                                                             # Fall back to a plain ESM2 masked language model when no other backend is selected.
        generated = generate_with_esm_masked_lm(                                      # Sample candidates with the ESM2 masked-LM fallback editor.
            context=context,
            model_name=args.esm_model,
            n_samples=args.n_samples,
            temperature=args.temperature,
            top_k=args.top_k,
            sampling_seed=args.sampling_seed,
        )

    # Write the final candidate table, now including explicit generation provenance fields.
    out_path = Path(args.out_csv)                                                     # Convert the requested output path into a Path object.
    out_path.parent.mkdir(parents=True, exist_ok=True)                                # Create the output directory if needed.
    generated.to_csv(out_path, index=False)                                           # Write the generated candidate table to CSV.
    print(f"Wrote: {out_path}")                                                       # Print the written CSV path for quick confirmation.


if __name__ == "__main__":                                                            # Standard Python entrypoint guard.
    main()                                                                            # Execute the generation CLI.
