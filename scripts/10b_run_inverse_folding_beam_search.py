#!/usr/bin/env python
"""Stage 10b: Run structure-conditioned local redesign with inverse-folding scoring.

This script is the core methodological step of Stage 10. Instead of searching mainly in
sequence space and only approximating structure with proxies, it explicitly scores each
candidate against the fixed seed scaffold with an inverse-folding model.

Operationally it does the following:
- starts from the validated seed sequence,
- expands only a compact set of allowed local substitutions,
- scores every proposal for target-host compatibility,
- scores every proposal for backbone compatibility with ESM-IF1,
- and keeps only the strongest beam of candidates round by round.

The result is a candidate table that is genuinely structure-conditioned rather than merely
sequence-first with downstream structural heuristics.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from phageforge.stage10_utils import (
    Stage10Candidate,
    apply_mutation,
    choose_top_substitutions,
    composite_stage10_score,
    compute_family_centroid,
    evaluate_candidate_table,
    greedy_diverse_subset,
    load_embedding_backend,
    load_inverse_folding_model,
    load_inverse_folding_structure,
    mutation_list,
    read_json,
    seed_everything,
    write_json,
)
from phageforge.stage09_utils import load_target_predictor


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 10 inverse-folding beam search."""
    ap = argparse.ArgumentParser(description="Run Stage 10 structure-conditioned inverse-folding search.")                                                                                          # Initialize the CLI argument parser with a description
    ap.add_argument("--stage10_context_json", type=str, required=True, help="Stage 10 context JSON produced by 10a_prepare_stage10_structure_context.py.")                                          # Require the strict Stage 10 inverse-folding blueprint
    ap.add_argument("--predictor_model", type=str, required=True, help="Serialized trained host predictor used for target-host probability scoring.")                                               # Require the serialized Logistic Regression host predictor model
    ap.add_argument("--label_classes_json", type=str, required=True, help="JSON file containing the predictor label order.")                                                                        # Require the JSON mapping the predictor output neurons to specific bacteria names
    ap.add_argument("--embedding_model", type=str, default="facebook/esm2_t33_650M_UR50D", help="Embedding backbone used by the trained predictor and family centroid scoring.")                    # Define the specific ESM-2 language model required for grammatical embeddings
    ap.add_argument("--beam_width", type=int, default=24, help="Number of candidates retained after each redesign round.")                                                                          # Define the absolute maximum number of parent candidates that survive each generational loop
    ap.add_argument("--rounds", type=int, default=4, help="Number of sequential redesign rounds to perform.")                                                                                       # Define the total number of iterative mutation steps the search will execute
    ap.add_argument("--proposals_per_parent", type=int, default=8, help="Maximum number of single-site proposal branches created from each parent.")                                                # Define the maximum number of child sequences spawned from a single parent candidate
    ap.add_argument("--substitutions_per_position", type=int, default=3, help="Maximum number of amino-acid substitutions tried at each editable position.")                                        # Restrict the AI to only testing the top-N best amino acids at any specific coordinate
    ap.add_argument("--if_chain_id", type=str, default="A", help="Chain identifier used when loading the seed scaffold for inverse folding.")                                                       # Define the specific biological chain to extract from the PDB file for structural grading
    ap.add_argument("--if_device", type=str, default="cuda", help="Device used to load and run the inverse-folding model.")                                                                         # Command the heavy Inverse-Folding Graph Neural Network to load onto specific hardware
    ap.add_argument("--batch_size", type=int, default=4, help="Batch size used for embedding candidates with the predictor backbone.")                                                              # Define the inference batch size to prevent VRAM Out-Of-Memory crashes
    ap.add_argument("--out_csv", type=str, required=True, help="Where to write the full Stage 10 search table.")                                                                                    # Require the destination path for the final search log CSV
    ap.add_argument("--out_json", type=str, required=True, help="Where to write the compact Stage 10 search summary JSON.")                                                                         # Require the destination path for the run metadata JSON
    ap.add_argument("--seed", type=int, default=42, help="Random seed used for deterministic search ordering.")                                                                                     # Define the random seed to guarantee strict scientific reproducibility
    return ap.parse_args()                                                                                                                                                                          # Parse the terminal arguments and return the Namespace object



def build_round_children(
    parents: list[Stage10Candidate],
    proposal_rows: list[dict],
    seed_sequence: str,
    max_mutations: int,
    proposals_per_parent: int,
    substitutions_per_position: int,
) -> list[Stage10Candidate]:
    """
    Expands the Beam Search by spawning a new generation of structure-aware "child" candidates.

    This function represents one round of the generative "step" of the Inverse-Folding engine. It strictly 
    controls how parents mutate into children by enforcing geographic limits, mutation budgets, 
    and biological substitution priorities.

    Specifically:
    1. It sorts the allowed edit rules (`proposal_rows`) prioritizing spots with high biological 
       reward and low evolutionary risk.
    2. For every parent sequence in the beam, it checks if the parent has already hit the 
       `max_mutations` ceiling. If so, it freezes that parent.
    3. It iterates over the sorted edit rules. It prevents the engine from mutating the exact 
       same position twice in one sequence.
    4. At an approved position, it calculates a blended priority score for all allowed amino acids 
       (weighing target-host success heavily). It takes the top `substitutions_per_position` (e.g., top 3).
    5. It physically alters the parent string, logs the mutation history into a new `Stage10Candidate` 
       object, and stops generating when it hits the `proposals_per_parent` cap.

    Example:
        children = build_round_children(beam, rules, "MKA", 4, 8, 3)
        Returns a list of ~192 new Stage10Candidate objects ready for physical scoring.
    """
    children: list[Stage10Candidate] = []                                                                                                                                                           # Initialize an empty list to stockpile the newly generated sequence objects
    # Sort the allowed edit rules (`proposal_rows`) so high biological reward mutations are prioritized
    proposal_rows = sorted(                                                                                                                                                                         # Begin sorting the allowed substitution rulebook to ensure the AI attacks high-value targets first
        proposal_rows,                                                                                                                                                                              # Pass the raw list of dictionaries defining the edit space
        key=lambda row: (float(row.get("functional_weight", 0.0)), -float(row.get("conservation_penalty", 0.0)), int(row["position"])),                                                             # Sort hierarchically: highest biological necessity, lowest evolutionary risk, then positional order
        reverse=True,                                                                                                                                                                               # Enforce descending order so the absolute best proposals sit at index zero
    )                                                                                                                                                                                               # Close the sorting block

    # For each parent in the current beam
    for parent in parents:                                                                                                                                                                          # Iterate sequentially through the elite candidates that survived the previous round
        # Check if this parent has already hit the structural mutation budget, and if so, break
        if len(parent.mutated_positions) >= max_mutations:                                                                                                                                          # Verify if this specific parent has already hit the absolute ceiling of the structural mutation budget
            continue                                                                                                                                                                                # Immediately freeze this parent and prevent it from spawning further mutations
        
        # Extract the integer coordinates of all edits this parent has already undergone
        parent_mutated = set(int(p) for p in parent.mutated_positions)                                                                                                                              
        emitted = 0                                                                                                                                                                                 # Initialize a counter tracking how many child branches this parent has successfully spawned
        # Iterate through the sorted substitution rules
        for row in proposal_rows:                                                                                                                                                                   
            pos = int(row["position"])                                                                                                                                                              # Extract the biological 1-based index for the currently evaluated rule
            if pos in parent_mutated:                                                                                                                                                               # Check if the parent sequence has already been mutated at this exact geometric coordinate
                continue                                                                                                                                                                            # Forbid double-editing the same spot and skip to the next available position

            allowed = list(row.get("allowed_aas", []))                                                                                                                                              # Extract the strictly curated menu of legally permitted amino acids for this spot
            target_preference = dict(row.get("target_preference", {}))                                                                                                                              # Extract the host-specific probability weights
            family_preference = dict(row.get("family_preference", {}))                                                                                                                              # Extract the evolutionary family probability weights
            functional_weight = float(row.get("functional_weight", 0.0))                                                                                                                            # Extract the raw biological necessity score

            aa_scores = []                                                                                                                                                                          # Initialize an empty list to rank the allowed amino acids
            # Iterate through all the allowed to be substituted amino acids
            for aa in allowed:                                                                                                                                                                      # Iterate through every legally permitted amino acid option
                # Calculate a substitution priority score, heavily prioritizing target-host infectivity metrics
                score = 0.55 * float(target_preference.get(aa, 0.0)) + 0.30 * float(family_preference.get(aa, 0.0)) + 0.15 * functional_weight                                                      
                aa_scores.append((aa, score))                                                                                                                                                       # Bundle the character and its calculated score into a tuple
            # Sort the amino acid options descending by score and then alphabetically
            aa_scores.sort(key=lambda item: (item[1], item[0]), reverse=True)                                                                                                                       # Sort the amino acid options descending by score, using alphabetical order for deterministic tie-breaking

            # Iterate over only the absolute best N amino acids to prevent combinatorial explosion
            for aa, _ in aa_scores[:substitutions_per_position]:                                                                                                                                    
                # Apply the substitution
                next_seq = apply_mutation(parent.candidate_sequence, pos, aa)                                                                                                                       # Physically generate the new child sequence by executing the string replacement
                if next_seq == parent.candidate_sequence:                                                                                                                                           # Double-check if the proposed "mutation" resulted in the exact same sequence
                    continue                                                                                                                                                                        # Discard the redundant sequence to save downstream compute
                # Log the mutation
                mutations = mutation_list(seed_sequence, next_seq)                                                                                                                                  # Generate the complete array of human-readable mutation tags relative to the original seed
                # Log the child
                child = Stage10Candidate(                                                                                                                                                           # Instantiate a new structured tracking object for the child
                    candidate_sequence=next_seq,                                                                                                                                                    # Inject the newly mutated amino acid text
                    parent_sequence=parent.parent_sequence,                                                                                                                                         # Log the text of the originating parent sequence
                    mutations=mutations,                                                                                                                                                            # Inject the array of human-readable mutation tags
                    mutated_positions=sorted({*parent_mutated, pos}),                                                                                                                               # Mathematically union the old and new edit positions and sort them
                    proposal_trace=list(parent.proposal_trace) + [f"{pos}:{row['seed_aa']}→{aa}"],                                                                                                  # Append the exact generative action taken to the unbroken audit trail
                    round_index=int(parent.round_index) + 1,                                                                                                                                        # Increment the generational counter
                )                                                                                                                                                                                   # Close the object instantiation
                children.append(child)                                                                                                                                                              # Officially add the newly minted child to the round's survivor pool
                emitted += 1                                                                                                                                                                        # Increment the counter tracking the parent's total spawned branches
                if emitted >= proposals_per_parent:                                                                                                                                                 # Check if the parent has exhausted its allowed generative branching quota
                    break                                                                                                                                                                           # Immediately terminate further mutation proposals for this specific parent
            if emitted >= proposals_per_parent:                                                                                                                                                     # Redundantly check the quota at the outer loop level to ensure strict compliance
                break                                                                                                                                                                               # Terminate the position-scanning loop and move to the next parent in the beam
    return children                                                                                                                                                                                 # Return the massive batch of new child sequences to the main orchestrator



def main() -> None:
    # Read the Stage 10 redesign context and initialize the exact structure-conditioned components required for search.
    args = parse_args()                                                                                                                                                                             # Execute argument parsing and secure user parameters into the args namespace
    seed_everything(args.seed)                                                                                                                                                                      # Lock the numpy and python random number generators to guarantee reproducible science
    context = read_json(args.stage10_context_json)                                                                                                                                                  # Load the strict Stage 10 mutation blueprint dictionary from disk

    # Recover the fixed scaffold backbone and the associated native sequence used by the inverse-folding model.
    seed_pdb_path = Path(context["seed_pdb_path"])                                                                                                                                                  # Extract and format the absolute filepath to the physical 3D wild-type scaffold
    coords, native_sequence = load_inverse_folding_structure(seed_pdb_path, chain_id=args.if_chain_id)                                                                                              # Parse the PDB file, extracting the raw atomic coordinate matrix and the original sequence
    seed_sequence = str(context["selected_seed"]["seed_sequence"])                                                                                                                                  # Extract the reference wild-type sequence from the metadata for comparison
    if len(native_sequence) != len(seed_sequence):                                                                                                                                                  # Check if the sequence extracted from the 3D file mismatches the 1D metadata length
        print("[WARN] Inverse-folding native sequence length differs from the selected seed; Stage 10 will still score the provided seed sequence on the fixed backbone.")                          # Warn the user of potential indexing misalignments between the PDB and the blueprint

    # Load the inverse-folding model once so the beam search can repeatedly score local redesign candidates against the fixed scaffold.
    _, if_model, if_alphabet = load_inverse_folding_model(device=args.if_device)                                                                                                                    # Load the massive ESM-IF1 Graph Neural Network and freeze it onto the requested compute hardware

    # Load the predictor model and host mappings, and the ESM-2 embedding backbone once so Stage 10 can reuse them across all search rounds.
    predictor, label_classes = load_target_predictor(args.predictor_model, args.label_classes_json)                                                                                                 # Load the trained Logistic Regression model and its bacterial classification mapping
    torch_emb, tokenizer_emb, model_emb, emb_device = load_embedding_backend(args.embedding_model)                                                                                                  # Load the massive ESM-2 sequence language model and freeze it onto the compute hardware

    # Build the family centroid under the same embedding backbone used by the trained host predictor.
    family_rows = pd.DataFrame(context.get("family_context", {}).get("family_rows", []))                                                                                                            # Reconstruct the dataframe containing thousands of natural evolutionary cousin sequences
    family_sequences = family_rows.get("aa_sequence", pd.Series(dtype=str)).astype(str).tolist() if not family_rows.empty else []                                                                   # Extract just the raw text strings of those sequences into a flat list
    family_centroid = compute_family_centroid(                                                                                                                                                      # Begin calculating the central point of the evolutionary manifold
        family_sequences=family_sequences,                                                                                                                                                          # Pass the raw text list of family sequences
        embedding_model=args.embedding_model,                                                                                                                                                       # Pass the specific embedding model architecture
        batch_size=args.batch_size,                                                                                                                                                                 # Pass the safe VRAM chunking parameter
        torch=torch_emb,                                                                                                                                                                            # Pass the pre-loaded PyTorch library
        tokenizer=tokenizer_emb,                                                                                                                                                                    # Pass the pre-loaded sequence tokenizer
        model=model_emb,                                                                                                                                                                            # Pass the pre-loaded ESM-2 neural network
        device=emb_device,                                                                                                                                                                          # Pass the active compute device map
    ) if family_sequences else None                                                                                                                                                                 # Safely bypass the calculation entirely if no evolutionary sequences exist

    # Extract the redesign constraints that Stage 10 will obey during beam expansion.
    editable_region = dict(context["editable_region"])                                                                                                                                              # Isolate the specific nested dictionary containing all the spatial mutation rules
    proposal_rows = list(editable_region.get("proposal_rows", []))                                                                                                                                  # Extract the massively detailed list of allowed substitutions per-position
    max_mutations = int(editable_region["max_mutations"])                                                                                                                                           # Secure the absolute upper bound of sequence edits permitted per candidate
    target_host = str(context["target_host"])                                                                                                                                                       # Extract the precise name of the bacterial species the AI must optimize for

    # Initialize the active beam with the untouched seed sequence so every redesign trajectory remains explicitly traceable.
    beam = [                                                                                                                                                                                        # Open the active tracking list representing the current "Beam" of surviving candidates
        Stage10Candidate(                                                                                                                                                                           # Instantiate the very first tracking object
            candidate_sequence=seed_sequence,                                                                                                                                                       # Set the current state to the untouched wild-type sequence
            parent_sequence=seed_sequence,                                                                                                                                                          # Set the parent state to the untouched wild-type sequence
            mutations=[],                                                                                                                                                                           # Initialize an empty array as zero edits have occurred
            mutated_positions=[],                                                                                                                                                                   # Initialize an empty array as zero positions have been altered
            proposal_trace=[],                                                                                                                                                                      # Initialize an empty audit log
            round_index=0,                                                                                                                                                                          # Tag this initial sequence explicitly as Round Zero (Generation 0)
        )                                                                                                                                                                                           # Close the object instantiation
    ]                                                                                                                                                                                               # Close the beam list initialization
    # Initialize a master list to collect the massive DataFrames generated after each scoring loop
    all_round_frames: list[pd.DataFrame] = []                                                                                                                                                       

    # Expand the beam one localized edit round at a time and rescore each candidate using both target and scaffold-conditioned objectives.
    for round_idx in range(1, int(args.rounds) + 1):                                                                                                                                                # Initiate the primary loop executing the predetermined number of evolutionary steps
        children = build_round_children(                                                                                                                                                            # Call the generative function to spawn mutated variations of the current beam survivors
            parents=beam,                                                                                                                                                                           # Feed the active pool of elite parents into the generator
            proposal_rows=proposal_rows,                                                                                                                                                            # Feed the strict biological rulebook controlling what edits are legal
            seed_sequence=seed_sequence,                                                                                                                                                            # Pass the original wild-type string for mutation tracking
            max_mutations=max_mutations,                                                                                                                                                            # Enforce the absolute structural mutation ceiling
            proposals_per_parent=args.proposals_per_parent,                                                                                                                                         # Enforce the limit on how many branches a single parent can spawn
            substitutions_per_position=args.substitutions_per_position,                                                                                                                             # Restrict how many different amino acids can be tested at a single coordinate
        )                                                                                                                                                                                           # Close the generator execution
        if not children:                                                                                                                                                                            # Check if the generator completely exhausted all legal mutation pathways
            break                                                                                                                                                                                   # Terminate the massive search loop immediately to prevent processing empty datasets

        # Deduplicate raw sequence proposals before expensive scoring so each candidate is evaluated only once per round.
        dedup: dict[str, Stage10Candidate] = {}                                                                                                                                                     # Initialize an empty dictionary utilizing raw strings as keys to enforce uniqueness
        for item in children:                                                                                                                                                                       # Iterate through the massive batch of newly spawned child candidates
            dedup.setdefault(item.candidate_sequence, item)                                                                                                                                         # Add the candidate to the dictionary, effortlessly overwriting/ignoring exact sequence clones
        children = list(dedup.values())                                                                                                                                                             # Extract the cleanly deduplicated candidate objects back into a flat list

        # Score all candidate sequences for the current round under the target predictor and the inverse-folding model on the fixed seed scaffold.
        score_frame = evaluate_candidate_table(                                                                                                                                                     # Call the ultimate multi-modal scoring gauntlet
            sequences=[item.candidate_sequence for item in children],                                                                                                                               # Pass the raw text sequences of all unique children
            target_host=target_host,                                                                                                                                                                # Pass the designated bacterial target name
            predictor_model_path=args.predictor_model,                                                                                                                                              # Provide the disk path to the host-prediction AI
            predictor_label_classes_path=args.label_classes_json,                                                                                                                                   # Provide the disk path to the classification mappings
            embedding_model=args.embedding_model,                                                                                                                                                   # Provide the name of the grammatical embedding backbone
            family_centroid=family_centroid,                                                                                                                                                        # Pass the mathematical center of the evolutionary manifold
            coords=coords,                                                                                                                                                                          # Pass the absolutely critical physical 3D atomic coordinates of the scaffold
            if_model=if_model,                                                                                                                                                                      # Pass the loaded Inverse-Folding Graph Neural Network
            if_alphabet=if_alphabet,                                                                                                                                                                # Pass the specialized token dictionary for the 3D network
            batch_size=args.batch_size,                                                                                                                                                             # Enforce the safe VRAM chunking parameter
            predictor=predictor,                                                                                                                                                                    # Pass the pre-loaded host prediction model
            label_classes=label_classes,                                                                                                                                                            # Pass the pre-loaded classification mappings
            torch=torch_emb,                                                                                                                                                                        # Pass the pre-loaded PyTorch library
            tokenizer=tokenizer_emb,                                                                                                                                                                # Pass the pre-loaded grammatical tokenizer
            model=model_emb,                                                                                                                                                                        # Pass the pre-loaded grammatical ESM-2 network
            device=emb_device,                                                                                                                                                                      # Pass the active compute device map
        )                                                                                                                                                                                           # Close the scoring gauntlet execution, returning a populated DataFrame

        # Attach trajectory-level metadata so the resulting search table remains auditable and easy to rank later.
        meta_frame = pd.DataFrame(                                                                                                                                                                  # Begin constructing a secondary DataFrame purely for tracking data
            {                                                                                                                                                                                       # Open the dictionary mapping columns to data lists
                "candidate_sequence": [item.candidate_sequence for item in children],                                                                                                               # Inject the primary key: the candidate sequence text
                "mutation_count": [len(item.mutated_positions) for item in children],                                                                                                               # Inject the calculated total edit volume
                "mutated_positions": [";".join(map(str, item.mutated_positions)) for item in children],                                                                                             # Format and inject the list of altered biological coordinates
                "mutation_text": [";".join(item.mutations) for item in children],                                                                                                                   # Format and inject the human-readable substitution tags
                "proposal_trace": [";".join(item.proposal_trace) for item in children],                                                                                                             # Format and inject the chronological audit log of generative actions
                "round_index": [int(item.round_index) for item in children],                                                                                                                        # Inject the generational tier number
                "seed_identity": [sum(a == b for a, b in zip(seed_sequence, item.candidate_sequence)) / len(seed_sequence) for item in children],                                                   # Compute and inject the exact fractional sequence identity compared to the original wild-type
            }                                                                                                                                                                                       # Close the dictionary mapping
        )                                                                                                                                                                                           # Close the tracking DataFrame construction
        round_frame = score_frame.merge(meta_frame, on="candidate_sequence", how="inner")                                                                                                           # Horizontally join the heavy scoring data with the audit tracking data using the sequence as the key

        # Add a new "stage10_composite_score" column combining structure-conditioned and target-conditioned terms.
        round_frame["stage10_composite_score"] = composite_stage10_score(                                                                                                                           # Execute the final overarching fitness calculus
            target_probability=round_frame["target_probability"].to_numpy(dtype=np.float32),                                                                                                        # Inject the predictive infectivity likelihood array
            if1_log_likelihood=round_frame["if1_log_likelihood"].to_numpy(dtype=np.float32),                                                                                                        # Inject the critical physical 3D structural compatibility array
            family_cosine=round_frame["family_cosine"].to_numpy(dtype=np.float32),                                                                                                                  # Inject the evolutionary realism array
            seed_identity=round_frame["seed_identity"].to_numpy(dtype=np.float32),                                                                                                                  # Inject the scaffold preservation metrics
            mutation_count=round_frame["mutation_count"].to_numpy(dtype=np.float32),                                                                                                                # Inject the edit burden penalty array
        )                                                                                                                                                                                           # Close the composite fitness calculation
        # Attach the explicit filepath of the anchor scaffold to every row, and sort the completely evaluated round ("round_frame") by fitness, strip old indices, and vault it into the master list ("all_round_frames").
        round_frame["seed_pdb_path"] = str(seed_pdb_path)                                                                                                                                           
        all_round_frames.append(round_frame.sort_values("stage10_composite_score", ascending=False).reset_index(drop=True))                                                                         
        
        # Embed the current round candidates with the same backbone to be able to apply the diverse subset selection later based on them.
        from phageforge.stage10_utils import embed_sequences_with_backend                                                                                                                           # Locally import the specific high-performance embedding tool
        cand_emb = embed_sequences_with_backend(                                                                                                                                                    # Execute the heavy transformer inference to map the survivors into spatial coordinates
            round_frame["candidate_sequence"].astype(str).tolist(),                                                                                                                                 # Feed the raw text sequences of the evaluated candidates
            torch_emb,                                                                                                                                                                              # Pass the active PyTorch library
            tokenizer_emb,                                                                                                                                                                          # Pass the active tokenizer
            model_emb,                                                                                                                                                                              # Pass the active ESM-2 network
            device=emb_device,                                                                                                                                                                      # Ensure execution remains on the correct hardware
            batch_size=args.batch_size,                                                                                                                                                             # Enforce the safe memory chunking limits
        )                                                                                                                                                                                           # Close the embedding execution
        # Extract the highly diverse, top-scoring sub-panel of candidates and keep the raw text sequences of the top-k winners in a fast lookup set.
        keep_idx = greedy_diverse_subset(cand_emb, round_frame["stage10_composite_score"].to_numpy(dtype=np.float32), top_k=args.beam_width)                                                        # Execute spatial repulsion math, extracting the indices of a highly diverse, top-scoring sub-panel
        keep_sequences = set(round_frame.iloc[keep_idx]["candidate_sequence"].astype(str).tolist())                                                                                                 # Extract the raw text sequences of those elite winners and cast them to a fast lookup set

        # Carry only these strongest diverse candidates (as Stage 10 objects) into the next beam round.
        next_beam: list[Stage10Candidate] = []                                                                                                                                                      # Initialize an empty roster for the upcoming generation
        for item in children:                                                                                                                                                                       # Iterate through the massive initial pool of all generated children
            if item.candidate_sequence in keep_sequences:                                                                                                                                           # Check if the specific child survived the brutal multi-objective and spatial diversity gauntlet
                next_beam.append(item)                                                                                                                                                              # Promote the surviving child to become a parent in the next round
        beam = next_beam                                                                                                                                                                            # Officially overwrite the active beam with the new generation

    # If the search successfully completed at least one functional round, vertically stack all generational DataFrames into one historical log of the entire search
    if all_round_frames:                                                                                                                                                                            # Check if the search successfully completed at least one functional round
        search_df = pd.concat(all_round_frames, ignore_index=True)                                                                                                                                  
    else:                                                                                                                                                                                           # Engage fallback logic if the entire search failed immediately
        # If the entire search failed, construct a completely empty DataFrame with the correct column structure
        search_df = pd.DataFrame(columns=[                                                                                                                                                          # Construct a completely empty, but correctly structured DataFrame
            "candidate_sequence",                                                                                                                                                                   # Define text sequence column
            "target_probability",                                                                                                                                                                   # Define host infectivity column
            "if1_log_likelihood",                                                                                                                                                                   # Define inverse-folding 3D physics column
            "family_cosine",                                                                                                                                                                        # Define evolutionary safety column
            "mutation_count",                                                                                                                                                                       # Define total edit volume column
            "mutated_positions",                                                                                                                                                                    # Define altered biological coordinates column
            "mutation_text",                                                                                                                                                                        # Define human-readable substitution tags column
            "proposal_trace",                                                                                                                                                                       # Define chronological audit log column
            "round_index",                                                                                                                                                                          # Define generational tier column
            "seed_identity",                                                                                                                                                                        # Define wild-type preservation percentage column
            "stage10_composite_score",                                                                                                                                                              # Define the ultimate ranking metric column
            "seed_pdb_path",                                                                                                                                                                        # Define the physical 3D anchor traceability column
        ])                                                                                                                                                                                          # Close the empty DataFrame structure definition

    # Sort the final search table by the Stage 10 score and write both the detailed CSV and the compact JSON summary.
    search_df = search_df.sort_values(["stage10_composite_score", "target_probability", "if1_log_likelihood"], ascending=False).reset_index(drop=True)                                              # Perform the ultimate hierarchical sort: overall fitness, then infectivity, then structural integrity
    out_csv = Path(args.out_csv)                                                                                                                                                                    # Convert the user-provided destination string into a robust Path object
    out_csv.parent.mkdir(parents=True, exist_ok=True)                                                                                                                                               # Ensure the entire directory tree exists on the filesystem
    search_df.to_csv(out_csv, index=False)                                                                                                                                                          # Dump the massive historical search log to disk without numerical indices

    summary = {                                                                                                                                                                                     # Begin constructing the compact run manifest
        "stage": "10b",                                                                                                                                                                             # Tag the artifact with its specific originating pipeline stage
        "stage10_context_json": str(args.stage10_context_json),                                                                                                                                     # Record the precise physical constraints used during the run
        "search_rows": int(len(search_df)),                                                                                                                                                         # Log the total volume of generated and evaluated hypotheses
        "rounds_completed": int(search_df["round_index"].max()) if len(search_df) else 0,                                                                                                           # Identify how deeply the generational loop penetrated before terminating
        "top_candidate_sequence": str(search_df.iloc[0]["candidate_sequence"]) if len(search_df) else "",                                                                                           # Extract the text string of the absolute #1 ranked protein design
        "top_stage10_score": float(search_df.iloc[0]["stage10_composite_score"]) if len(search_df) else float("nan"),                                                                               # Extract the highest achieved composite fitness metric
        "top_target_probability": float(search_df.iloc[0]["target_probability"]) if len(search_df) else float("nan"),                                                                               # Extract the highest achieved host infectivity likelihood
        "top_if1_log_likelihood": float(search_df.iloc[0]["if1_log_likelihood"]) if len(search_df) else float("nan"),                                                                               # Extract the highest achieved 3D physics stability score
        "out_csv": str(out_csv),                                                                                                                                                                    # Record the permanent location of the heavy data log
    }                                                                                                                                                                                               # Close the dictionary mapping
    write_json(summary, args.out_json)                                                                                                                                                              # Serialize and dump the compact manifest to disk
    print(f"Wrote: {args.out_csv}")                                                                                                                                                                 # Output terminal confirmation for the heavy CSV
    print(f"Wrote: {args.out_json}")                                                                                                                                                                # Output terminal confirmation for the JSON manifest


if __name__ == "__main__":                                                                                                                                                                          # Standard programmatic execution validation guard check
    main()                                                                                                                                                                                          # Invoke the primary orchestration function