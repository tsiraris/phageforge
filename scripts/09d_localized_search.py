#!/usr/bin/env python
"""Stage 09d: Run structure-aware localized search around the selected seed.

This search stage deliberately replaces broad free-form generation with a narrow edit-space search.
It proposes only seed-local substitutions, scores candidates with the existing target model plus the
Stage 09 structural surrogate, and keeps a beam of sequences that balance retargeting and scaffold
preservation.
"""

from __future__ import annotations

import argparse
import math
import random
from pathlib import Path

import numpy as np
import pandas as pd

from phageforge.stage07_utils import cosine_similarity
from phageforge.stage09_utils import (
    EditProposal,
    build_basic_sequence_features,
    build_edit_proposals_from_context,
    embed_sequences,
    greedy_diverse_pick,
    load_target_predictor,
    maybe_load_surrogate,
    mutation_list,
    parse_mutation_positions,
    predict_target_probability,
    read_json,
    sequence_identity,
    substitution_priority,
    surrogate_structural_risk,
    write_json,
)

def parse_args() -> argparse.Namespace:
    """ Initializes and configures the command-line interface for the Stage 09 Search Engine. """
    ap = argparse.ArgumentParser(description="Run Stage 09 structure-aware localized search.")                                                                                                # Initialize the argument parser for the command-line interface
    ap.add_argument("--context_json", type=str, required=True, help="Stage 07 context JSON produced by 07a_prepare_stage07_design_context.py.")                                               # Define the required context JSON file argument
    ap.add_argument("--edit_space_json", type=str, required=True, help="Stage 09 edit-space JSON produced by 09a_define_edit_space.py.")                                                      # Define the required edit-space JSON file argument
    ap.add_argument("--strict_csv", type=str, default=None, help="Optional strict-bank CSV that can sharpen substitution proposals if the edit-space JSON needs rebuilding.")                 # Define the optional strict CSV file argument
    ap.add_argument("--predictor_model", type=str, required=True, help="Host predictor joblib file, typically results/broad/linear_probe/seed_42/model.joblib.")                              # Define the required predictor model joblib file argument
    ap.add_argument("--label_classes_json", type=str, required=True, help="JSON file storing the predictor label order.")                                                                     # Define the required label classes JSON file argument
    ap.add_argument("--surrogate_model", type=str, default=None, help="Optional Stage 09 surrogate joblib bundle produced by 09c_train_structure_surrogate.py.")                              # Define the optional surrogate model joblib file argument
    ap.add_argument("--out_csv", type=str, required=True, help="Where to write the full Stage 09 localized-search candidate table.")                                                          # Define the required output CSV file argument
    ap.add_argument("--out_json", type=str, required=True, help="Where to write the Stage 09 run metadata JSON.")                                                                             # Define the required output JSON file argument
    ap.add_argument("--esm_model", type=str, default="facebook/esm2_t33_650M_UR50D", help="Embedding model used for target scoring and manifold similarity.")                                 # Define the ESM embedding model argument
    ap.add_argument("--batch_size", type=int, default=4, help="Batch size used by the embedding model during proposal scoring.")                                                              # Define the batch size argument for embedding calculations
    ap.add_argument("--max_aa", type=int, default=2048, help="Maximum sequence length passed into the ESM embedding model.")                                                                  # Define the maximum sequence length argument
    ap.add_argument("--rounds", type=int, default=4, help="How many beam-search rounds to run.")                                                                                              # Define the number of beam-search rounds argument
    # Strict maximum number of candidates allowed to survive and advance to the next round.
    ap.add_argument("--beam_width", type=int, default=24, help="How many candidates to keep after each round.")                                                                               # Define the beam width argument per round
    # The number of mutated children spawned from each parent per round.
    ap.add_argument("--proposals_per_parent", type=int, default=18, help="How many local edits to propose per parent candidate in each round.")                                               # Define the number of proposals per parent argument
    ap.add_argument("--max_mutations", type=int, default=8, help="Maximum allowed mutation count for any candidate.")                                                                        # Define the maximum allowed mutations budget argument
    ap.add_argument("--seed", type=int, default=42, help="Random seed for deterministic search reproducibility.")                                                                             # Define the random seed argument for reproducibility
    ap.add_argument("--w_target", type=float, default=0.45, help="Weight of the target-host probability term in the combined Stage 09 score.")                                                # Define the weight for the target probability score
    ap.add_argument("--w_family", type=float, default=0.15, help="Weight of the family-centroid cosine term in the combined Stage 09 score.")                                                 # Define the weight for the family cosine score
    ap.add_argument("--w_seed", type=float, default=0.15, help="Weight of the seed-similarity term in the combined Stage 09 score.")                                                          # Define the weight for the seed cosine score
    ap.add_argument("--w_guidance", type=float, default=0.10, help="Weight of the local substitution-prior term in the combined Stage 09 score.")                                             # Define the weight for the guidance score
    ap.add_argument("--w_surrogate", type=float, default=0.15, help="Weight of the structural-risk penalty in the combined Stage 09 score.")                                                  # Define the weight for the surrogate structural risk penalty
    return ap.parse_args()                                                                                                                                                                    # Parse and return the provided command-line arguments


def load_proposal_table(edit_space: dict, context: dict, strict_df: pd.DataFrame | None) -> list[EditProposal]:
    """
    Parses and serializes the authorized mutation rules for the current redesign project.

    This function attempts to load pre-calculated substitution rules from the Stage 09a 
    edit-space artifact. If the artifact is empty or missing data, it falls back to 
    extracting family and host preferences from the raw databases.

    It iterates through the JSON "proposal_rows," casting raw dictionary values back into 
    the strictly typed `EditProposal` dataclass. This ensures that downstream search 
    logic can access positions, frequencies, and weights using object-oriented dot notation.

    Example:
        proposals = load_proposal_table(edit_space_dict, context_dict, None)
        Returns: [EditProposal(position=42, allowed_aas=['S', 'T'], ...), ...]
    """
    proposal_rows = edit_space.get("edit_space", {}).get("proposal_rows", [])                                                       # Extract proposal rows from the edit_space dictionary safely
    if proposal_rows:                                                                                                               # Check if the extracted proposal rows list is not empty
        proposals: list[EditProposal] = []                                                                                          # Initialize an empty list to hold the parsed EditProposal objects
        for row in proposal_rows:                                                                                                   # Iterate over each row in the extracted proposal rows
            proposals.append(                                                                                                       # Append a newly created EditProposal object to the list
                EditProposal(                                                                                                       # Instantiate an EditProposal object with mapped fields
                    position=int(row["position"]),                                                                                  # Extract and cast the mutation position to an integer
                    seed_aa=str(row["seed_aa"]),                                                                                    # Extract and cast the original seed amino acid to a string
                    allowed_aas=[str(x) for x in row.get("allowed_aas", [])],                                                       # Extract and cast the allowed amino acids list to strings
                    target_preference={str(k): float(v) for k, v in row.get("target_preference", {}).items()},                      # Extract and cast the target preference mapping dictionary
                    family_preference={str(k): float(v) for k, v in row.get("family_preference", {}).items()},                      # Extract and cast the family preference mapping dictionary
                    functional_weight=float(row.get("functional_weight", 0.0)),                                                     # Extract and cast the functional weight score to a float
                    conservation_penalty=float(row.get("conservation_penalty", 0.0)),                                               # Extract and cast the conservation penalty score to a float
                    region_name=str(row.get("region_name", "ungrouped")),                                                           # Extract and cast the region name to a string with a fallback
                )                                                                                                                   # Close the EditProposal instantiation
            )                                                                                                                       # Close the list append call
        return proposals                                                                                                            # Return the populated list of EditProposal objects
    return build_edit_proposals_from_context(context=context, strict_df=strict_df)                                                  # Fallback to building proposals from context if rows were empty


def sequence_from_parent(parent_sequence: str, position: int, new_aa: str) -> str:
    """
    Generates a child sequence by applying a single amino acid substitution.

    Since strings are immutable in Python, this function converts the parent sequence 
    into a list of characters, translates the 1-based biological position into a 
    0-indexed list offset, swaps the character, and joins the list back into a string.

    Example:
        child = sequence_from_parent("MKA", 2, "G")
        Returns: "MGA"
    """
    chars = list(parent_sequence)                   # Convert the parent sequence string into a mutable list of characters
    chars[int(position) - 1] = str(new_aa)          # Replace the character at the specified 0-indexed position
    return "".join(chars)                           # Join the character list back into a string and return the mutated sequence


def local_guidance_score(proposal_map: dict[int, EditProposal], seed_sequence: str, candidate_sequence: str) -> float:
    """
    Evaluates how closely a candidate adheres to the pre-approved biological "rulebook."

    This score acts as a "loyalty" check, ensuring the generator hasn't hallucinated 
    amino acids that violate evolutionary or host-specific frequencies.

    It identifies all mutation positions by comparing the candidate to the seed. For 
    each position, it looks up the `EditProposal`. It computes the average 
    `substitution_priority` (target vs. family likelihood) for all edits made. If an 
    edit was made at a position not in the map, it applies a heavy -1.0 penalty.

    Example:
        If the candidate mutated position 42 to 'S', and 'S' was the #1 target preference:
        score = local_guidance_score(pmap, seed, cand)
        Returns: e.g., 0.92
    """
    mutations = parse_mutation_positions(mutation_list(seed_sequence, candidate_sequence))                                          # Extract mutation positions between the seed and candidate sequence
    if not mutations:                                                                                                               # Check if there are no mutations present between the sequences
        return -1.0                                                                                                                 # Return a default negative score since it's identical
    scores = []                                                                                                                     # Initialize an empty list to collect individual mutation scores
    for pos in mutations:                                                                                                           # Iterate through each parsed mutation position
        item = proposal_map.get(int(pos))                                                                                           # Retrieve the corresponding EditProposal for the current position
        if item is None:                                                                                                            # Check if no proposal exists for this mutation position
            scores.append(-1.0)                                                                                                     # Append a penalty score for making an unguided mutation
            continue                                                                                                                # Skip the rest of the loop for this unguided position
        aa = candidate_sequence[int(pos) - 1]                                                                                       # Extract the actual mutated amino acid from the candidate
        scores.append(substitution_priority(item, aa))                                                                              # Calculate and append the mapped priority score for this substitution
    return float(np.mean(scores)) if scores else -1.0                                                                               # Return the average of collected scores, or a default penalty if empty


def build_round_proposals(parents: list[dict], proposal_rows: list[EditProposal], seed_sequence: str, max_mutations: int, proposals_per_parent: int, rng: random.Random) -> list[dict]:
    """
    Expands the Beam Search by spawning a new generation of single-edit "child" candidates.

    For every parent sequence in the beam, it scans the `EditProposal` table. It 
    evaluates every possible substitution at every allowed position. It skips 
    identical edits and enforces the `max_mutations` budget. Potential edits are 
    ranked by biological priority, shuffled for stochastic diversity, and the top-N 
    are used to create new candidate records for scoring.

    Example:
        children = build_round_proposals(beam_survivors, proposals, seed_seq, 8, 18, rng)
        Returns a list of ~432 children (24 parents * 18 proposals each).
    """
    proposal_map = {int(item.position): item for item in proposal_rows}                                                             # Create a dictionary mapping positions to proposals for quick O(1) lookups
    rows: list[dict] = []                                                                                                           # Initialize an empty list to hold the newly proposed candidate dictionaries

    for parent in parents:                                                                                                          # Iterate over each parent candidate supplied from the current beam
        parent_sequence = str(parent["candidate_sequence"])                                                                         # Extract the parent's sequence string representation
        # Extract a set of integers representing already mutated positions on the parent
        parent_positions = set(int(x) for x in parent.get("mutated_positions", []))                                                 

        # Rank mutation opportunities by structural locality and proposal priority so the search edits high-value sites first.
        ranked_sites = []                                                                                                           # Initialize an empty list to evaluate and rank potential mutation sites
        # Iterate through all structurally allowed proposal definitions
        for item in proposal_rows:                                                                                                  
            current_aa = parent_sequence[int(item.position) - 1]                                                                    # Get the amino acid currently residing at the proposal's specific position
            for aa in item.allowed_aas:                                                                                             # Iterate through all globally allowed substitutions for this position
                # Skip identical edits
                if aa == current_aa:                                                                                                # Check if the proposed new amino acid is identical to the current one
                    continue                                                                                                        # Skip this specific substitution as it would cause no net change
                # Calculate the total cumulative mutations if this edit is actually applied
                child_mut_count = len(parent_positions | {int(item.position)})                                                      
                # If the edit exceeds the maximum mutation budget, skip it
                if child_mut_count > int(max_mutations):                                                                            # Check if the potential edit strictly exceeds the maximum allowed mutations budget
                    continue                                                                                                        # Skip this substitution to strictly enforce the maximum mutation budget
                # Record the substitution priority score, the integer position, and the new amino acid
                ranked_sites.append((substitution_priority(item, aa), int(item.position), aa))                                      
        # Sort the collected substitution opportunities by priority (descending) and position (ascending)
        ranked_sites.sort(key=lambda x: (x[0], -x[1]), reverse=True)                                                                
        if not ranked_sites:                                                                                                        # Check if absolutely no valid mutation sites were found for this parent
            continue                                                                                                                # Move immediately to the next parent if there are no valid edit opportunities

        # Keep the pool of ranked sites to twice the target proposal count, then randomly shuffle, re-sort by priority score and select the top N proposals
        ranked_sites = ranked_sites[: max(1, int(proposals_per_parent) * 2)]                                                        # Truncate the ranked list to a manageable pool strictly twice the target proposal size
        rng.shuffle(ranked_sites)                                                                                                   # Randomly shuffle the pool to organically introduce necessary diversity in search exploration
        # When the final sorted() command runs on this shuffled list: The mutations with the absolute highest scores still rise to the top.
        # But among the mutations with tied scores, the winner is now chosen randomly because the "original relative order" was randomized by the shuffle.
        ranked_sites = sorted(ranked_sites, key=lambda x: x[0], reverse=True)[: int(proposals_per_parent)]                          # Re-sort by priority score and cleanly slice off exactly the top N desired proposals

        # Iterate over the finalized top proposals assigned for this particular parent, applying the substitutions to the parent sequence, 
        # and create a new candidate record (dict) to the batch rows list for each proposal applied.
        for local_rank, (_, position, aa) in enumerate(ranked_sites, start=1):                                                      # Iterate over the finalized top proposals assigned for this particular parent
            candidate_sequence = sequence_from_parent(parent_sequence=parent_sequence, position=position, new_aa=aa)                # Generate the new mutated sequence by applying the selected edit to the parent
            rows.append(                                                                                                            # Append the newly formulated candidate dictionary to the batch rows list
                {                                                                                                                   # Open the candidate dictionary definition
                    "candidate_sequence": candidate_sequence,                                                                       # Store the newly generated sequence text
                    "parent_sequence": parent_sequence,                                                                             # Record the originating parent sequence text
                    "parent_id": str(parent["candidate_id"]),                                                                       # Record the originating parent's unique string identifier
                    "proposal_position": int(position),                                                                             # Record the integer position that was explicitly mutated
                    "proposal_aa": str(aa),                                                                                         # Record the new amino acid character that was introduced
                    "proposal_rank_within_parent": int(local_rank),                                                                 # Record the generated rank of this proposal specific to the parent
                    "round_index": int(parent["round_index"]) + 1,                                                                  # Increment the round index identifier for the child candidate
                    "mutations": mutation_list(seed_sequence, candidate_sequence),                                                  # Compute and store the full cumulative list of mutations relative to the original seed
                }                                                                                                                   # Close the candidate dictionary definition
            )                                                                                                                       # Close the list append call
    return rows                                                                                                                     # Return the fully completed batch of newly proposed candidates


def score_candidates(candidate_rows: list[dict], seed_sequence: str, target_host: str, predictor_model, label_classes: list[str], surrogate_bundle: dict | None, proposal_map: dict[int, EditProposal], family_centroid: np.ndarray, seed_embedding: np.ndarray, args: argparse.Namespace) -> pd.DataFrame:
    """
    Executes the multi-objective Stage 09 grading function for a batch of candidates.

    This function represents the mathematical core of the search, calculating target 
    probability, manifold similarity, and predicted structural risk.

    1. It deduplicates sequences to minimize expensive GPU calls.
    2. It generates ESM-2 embeddings for the entire batch.
    3. It runs the Host Predictor (Logistic Regression) and computes Cosine Similarities 
       to the seed and family centroid.
    4. It calculates sequence pathology metrics (entropy, homopolymers) and queries 
       the Structural Surrogate (Random Forest) for a risk estimate.
    5. It blends these metrics into the final `stage09_score` using user-defined weights.

    Example:
        scored_df = score_candidates(child_rows, seed, "Klebsiella", model, labels, surrogate, pmap, centroid, seed_emb, args)
        Returns a DataFrame where rows have a calculated "stage09_score" (e.g., 0.785), and other already calculated metrics such as "mutation_count", "family_cosine_similarity", "guidance_score", etc.
    """
    if not candidate_rows:                                                                                                          # Check if the incoming batch of candidate rows is completely empty
        return pd.DataFrame()                                                                                                       # Return an empty DataFrame immediately to gracefully avoid downstream errors

    # Deduplicate identical sequences first so expensive embedding/model evaluation is done only once per candidate sequence.
    tmp_df = pd.DataFrame(candidate_rows).drop_duplicates(subset=["candidate_sequence"]).reset_index(drop=True)                      # Create a DataFrame and aggressively remove duplicate sequences to save precious compute
    embeddings = embed_sequences(tmp_df["candidate_sequence"].astype(str).tolist(), model_name=args.esm_model, batch_size=args.batch_size, max_length=args.max_aa)             # Generate heavy ESM language model embeddings for the unique candidate sequences

    # Score candidates with the existing host predictor, then add seed/family similarity terms from the same embedding batch.
    target_prob = predict_target_probability(model=predictor_model, label_classes=label_classes, target_host=target_host, embeddings=embeddings)                             # Run the host predictor model to score the biological target probability
    seed_cos = np.asarray([cosine_similarity(row, seed_embedding) for row in embeddings], dtype=np.float32)                                                                 # Calculate cosine similarity distances between all candidates and the original seed
    family_cos = np.asarray([cosine_similarity(row, family_centroid) for row in embeddings], dtype=np.float32) if family_centroid.size else np.zeros(len(tmp_df), dtype=np.float32)  # Calculate cosine similarity to the family centroid, or smartly default to an array of zeros

    # Compute lightweight sequence features and ask the structural surrogate to estimate risk before any expensive folding step.
    # Create a quick lookup set containing all designated editable positions
    editable_positions = set(proposal_map)                                                                                          
    feature_rows = []                                                                                                               # Initialize an empty list to systematically hold computed sequence feature dictionaries
    guidance_scores = []                                                                                                            # Initialize an empty list to systematically hold calculated local guidance scores
    # Iterate sequentially over each unique candidate sequence string
    for seq in tmp_df["candidate_sequence"].astype(str):                                                                            
        # Compute standard biological metrics like sequence complexity and mutation counts
        feature_row = build_basic_sequence_features(seed_sequence=seed_sequence, candidate_sequence=seq, editable_positions=editable_positions)                              
        # Initialize a zero placeholder score reserved for downstream pipelines
        feature_row["final_multimodal_rank_score"] = 0.0                                                                            # Initialize a zero placeholder score reserved for downstream pipelines
        feature_row["strict_manifold_score"] = 0.0                                                                                  # Initialize a zero placeholder strict structural manifold score
        feature_row["structure_score"] = 0.0                                                                                        # Initialize a zero placeholder general structure score
        feature_row["target_score"] = 0.0                                                                                           # Initialize a zero placeholder general target alignment score
        feature_row["family_cosine"] = 0.0                                                                                          # Initialize a zero placeholder raw family cosine score
        feature_row["seed_cosine"] = 0.0                                                                                            # Initialize a zero placeholder raw seed cosine score
        feature_row["target_anchor_cosine"] = 0.0                                                                                   # Initialize a zero placeholder target anchor cosine score
        # Append the feature row to the collection, and compute the local guidance score for each candidate
        feature_rows.append(feature_row)                                                                                            # Append the completely populated feature row to the collection list
        guidance_scores.append(local_guidance_score(proposal_map=proposal_map, seed_sequence=seed_sequence, candidate_sequence=seq))                                         # Calculate and append the local guidance score reflecting edit compliance for this sequence
    # Convert the raw list of feature dictionaries into a DataFrame, and attach the computed family cosine similarities,
    # seed cosine similarities, and target prediction probabilities to the features
    feature_df = pd.DataFrame(feature_rows)                                                                                         # Convert the raw list of feature dictionaries into a structured DataFrame
    feature_df["family_cosine"] = family_cos                                                                                        # Explicitly attach the computed family cosine similarities to the features
    feature_df["seed_cosine"] = seed_cos                                                                                            # Explicitly attach the computed seed cosine similarities to the features
    feature_df["target_score"] = target_prob                                                                                        # Explicitly attach the computed target prediction probabilities to the features
    # Run the structural surrogate model to estimate folding risk, 
    # If no Random Forrest bundle is provided, falls back to a mathematical heuristic equation that penalizes sequence pathologies (like low complexity or extreme drift).
    structural_risk, pred_plddt, pred_rmsd = surrogate_structural_risk(bundle=surrogate_bundle, feature_frame=feature_df)           # Execute the structural surrogate model to swiftly estimate folding risk metrics

    # Combine all terms into one Stage 09 search score that explicitly rewards target gain while penalizing structural risk.
    tmp_df["target_probability"] = target_prob                                                                                      # Add the target probability column to the main candidate tracking DataFrame
    tmp_df["seed_cosine"] = seed_cos                                                                                                # Add the seed cosine similarity column to the main candidate tracking DataFrame
    tmp_df["family_cosine"] = family_cos                                                                                            # Add the family cosine similarity column to the main candidate tracking DataFrame
    tmp_df["guidance_score"] = np.asarray(guidance_scores, dtype=np.float32)                                                        # Add the cleanly computed guidance scores as a 32-bit float array column
    tmp_df["predicted_structural_risk"] = structural_risk                                                                           # Add the surrogate's computed structural risk score column
    tmp_df["predicted_mean_plddt"] = pred_plddt                                                                                     # Add the surrogate's predicted mean pLDDT confidence score column
    tmp_df["predicted_rmsd"] = pred_rmsd                                                                                            # Add the surrogate's predicted structural RMSD score column
    tmp_df["sequence_identity"] = feature_df["seed_identity"].to_numpy(dtype=np.float32)                                            # Add the seed identity percentage cleanly extracted from the feature set
    tmp_df["mutation_count"] = feature_df["mutation_count"].to_numpy(dtype=int)                                                     # Add the exact sequence mutation count cleanly extracted from the feature set
    tmp_df["mutation_span"] = feature_df["mutation_span"].to_numpy(dtype=int)                                                       # Add the overall mutation span width cleanly extracted from the feature set
    tmp_df["low_complexity_fraction"] = feature_df["low_complexity_fraction"].to_numpy(dtype=np.float32)                            # Add the specific fraction of low complexity regions extracted from the feature set
    tmp_df["outside_editable_fraction"] = feature_df["outside_editable_fraction"].to_numpy(dtype=np.float32)                        # Add the exact fraction of mutations strictly outside allowed bounds from the feature set
    tmp_df["stage09_score"] = (                                                                                                     # Start calculating the comprehensive combined multi-objective Stage 09 score
        args.w_target * tmp_df["target_probability"]                                                                                # Add the user-weighted target probability term to the total
        + args.w_family * tmp_df["family_cosine"]                                                                                   # Add the user-weighted family cosine alignment term to the total
        + args.w_seed * tmp_df["seed_cosine"]                                                                                       # Add the user-weighted seed cosine alignment term to the total
        + args.w_guidance * tmp_df["guidance_score"]                                                                                # Add the user-weighted edit guidance conformity score term to the total
        - args.w_surrogate * tmp_df["predicted_structural_risk"]                                                                    # Subtract the user-weighted structural risk penalty term from the total
    )                                                                                                                               # Close the Stage 09 score mathematical calculation
    return tmp_df                                                                                                                   # Return the fully populated, scored, and annotated DataFrame

def main() -> None:
    # Read the Stage 07 context.json, Stage 09 edit-space definition (edit_space.json), strict CSV, pre-trained host predictor (model.joblib from 03b), and the structural surrogate bundle (surrogate_model.joblib from 09c).
    args = parse_args()                                                                                                                                                                                     # Parse and store the command-line arguments securely
    rng = random.Random(args.seed)                                                                                                                                                                          # Initialize a seeded random number generator for absolute reproducibility
    context = read_json(args.context_json)                                                                                                                                                                  # Read and strictly parse the design context JSON configuration file
    edit_space = read_json(args.edit_space_json)                                                                                                                                                            # Read and strictly parse the explicit edit-space definition JSON file
    strict_df = pd.read_csv(args.strict_csv) if args.strict_csv else None                                                                                                                                   # Conditionally load the strict CSV dataframe safely if it was provided
    predictor_model, label_classes = load_target_predictor(model_path=args.predictor_model, label_classes_path=args.label_classes_json)                                                                     # Load the pre-trained host predictor model alongside its classification labels
    surrogate_bundle = maybe_load_surrogate(args.surrogate_model)                                                                                                                                           # Attempt to intelligently load the structural surrogate model bundle into memory

    # Recover the proposal table (list of EditableProposal objects)
    proposal_rows = load_proposal_table(edit_space=edit_space, context=context, strict_df=strict_df)                                                                                                        # Load all discrete mutation proposals dynamically based on edit space and context
    if not proposal_rows:                                                                                                                                                                                   # Verify that we have actionable programmatic proposals to actively work with
        raise ValueError("No proposal rows were available for Stage 09 localized search.")                                                                                                                  # Immediately abort execution if absolutely no proposal rows are available
    # Filter available proposals to tightly include only explicitly editable positions
    proposal_map = {int(item.position): item for item in proposal_rows if int(item.position) in set(edit_space["edit_space"]["editable_positions"])}                                                        
    # Sort and sequentially store the finalized allowed and constrained proposal rows
    proposal_rows = [proposal_map[pos] for pos in sorted(proposal_map)]                                                                                                                                     
    # Extract the "Follow-up" seed sequence, target host, and family centroid from the context, and calculate the dense baseline embedding for the seed
    seed_sequence = str(context["selected_seed"]["seed_sequence"])                                                                                                                                          # Extract the authoritative root seed sequence string out of the context
    target_host = str(context["target_host"])                                                                                                                                                               # Extract the designated desired target biological host string out of the context
    seed_embedding = embed_sequences([seed_sequence], model_name=args.esm_model, batch_size=1, max_length=args.max_aa)[0]                                                                                   # Generate the dense baseline embedding representation exclusively for the seed
    family_centroid = np.asarray(context.get("family_context", {}).get("family_centroid", []), dtype=np.float32)                                                                                            # Extract the family centroid vector safely, converting it natively to a numpy array

    # Initialize the beam with the "Follow-up" seed so the first round explores only tiny local edits.
    beam = [                                                                                                                                                                                                # Initialize the master beam search list with the unmodified seed state
        {                                                                                                                                                                                                   # Open the initial origin candidate dictionary block
            "candidate_id": "seed_0",                                                                                                                                                                       # Set the explicit tracking ID for the unmutated origin seed
            "candidate_sequence": seed_sequence,                                                                                                                                                            # Populate the raw root sequence directly into the search beam
            "mutations": [],                                                                                                                                                                                # Initialize an explicitly empty mutation list for the untouched root
            "mutated_positions": [],                                                                                                                                                                        # Initialize an explicitly empty mutation positions list for the untouched root
            "round_index": 0,                                                                                                                                                                               # Set the authoritative starting search round index cleanly to 0
        }                                                                                                                                                                                                   # Close the initial origin candidate dictionary block
    ]                                                                                                                                                                                                       # Close the master beam list initialization
    all_rows: list[pd.DataFrame] = []                                                                                                                                                                       # Initialize a dynamic list to safely accumulate robust dataframes from all rounds

    # Expand the beam round by round, score proposals in batches, and keep the strongest localized candidates only.
    for round_idx in range(1, int(args.rounds) + 1):                                                                                                                                                        # Loop sequentially through each expressly requested beam search round index
        # Generate a new generation of #args.proposals_per_parent single-edit "child" candidates (list of candidate dictionaries), sorted primarily by which amino acid substitution has the highest substitution_priority score, calculated from:
        # 50% of the frequency of the amino acid in the target host bacteria (target_preference), 30% of the frequency in the broad evolutionary family (family_preference), and general biological importance of mutating that position (functional_weight) (20%)
        candidate_rows = build_round_proposals(                                                                                                                                                             # Generate localized single-site edit proposals specifically for the current beam
            parents=beam,                                                                                                                                                                                   # Pass the current active beam efficiently as the parent generating pool
            proposal_rows=proposal_rows,                                                                                                                                                                    # Provide the rigidly allowed mutation proposal lookup table
            seed_sequence=seed_sequence,                                                                                                                                                                    # Provide the original authoritative seed sequence anchor reference
            max_mutations=args.max_mutations,                                                                                                                                                               # Pass the uncompromising maximum mutation budget numeric constraint
            proposals_per_parent=args.proposals_per_parent,                                                                                                                                                 # Dictate precisely how many discrete edits to branch off each parent
            rng=rng,                                                                                                                                                                                        # Pass the seeded random number generator cleanly for replicability
        )                                                                                                                                                                                                   # Close the build round proposal localized function call
        # Evaluate and thoroughly score the newly proposed batch of candidates based on host-binding probability, family & seed cosine similarity, guidance score, and predicted structural risk
        scored_df = score_candidates(                                                                                                                                                                       
            candidate_rows=candidate_rows,                                                                                                                                                                  # Supply the raw generated candidate dictionaries effectively for scoring
            seed_sequence=seed_sequence,                                                                                                                                                                    # Supply the reference root seed sequence efficiently for comparisons
            target_host=target_host,                                                                                                                                                                        # Supply the target host label precisely for predictive alignment scoring
            predictor_model=predictor_model,                                                                                                                                                                # Supply the heavily loaded target predictor model robustly into memory
            label_classes=label_classes,                                                                                                                                                                    # Supply the specific target predictor classes natively for mapping
            surrogate_bundle=surrogate_bundle,                                                                                                                                                              # Supply the dynamic structural surrogate bundle functionally for risk checks
            proposal_map=proposal_map,                                                                                                                                                                      # Supply the static lookup map securely for edit constraint guidance
            family_centroid=family_centroid,                                                                                                                                                                # Supply the vector manifold reference centroid cleanly for calculations
            seed_embedding=seed_embedding,                                                                                                                                                                  # Supply the root seed anchor embedding actively for similarity tracking
            args=args,                                                                                                                                                                                      # Pass the securely parsed command-line arguments explicitly for weights
        )                                                                                                                                                                                                   # Close the heavy candidate batch scoring function securely
        if scored_df.empty:                                                                                                                                                                                 # Check if the intensive scoring phase alarmingly resulted in an empty dataframe
            break                                                                                                                                                                                           # Gracefully exit the beam search loop entirely if absolutely no candidates remain

        # Apply hard prefilters that keep the search inside the structure-preserving design neighborhood.
        scored_df = scored_df.loc[                                                                                                                                                                          # Apply rigid biological sequence constraints to stringently filter out malformations
            (scored_df["outside_editable_fraction"] <= 0.0)                                                                                                                                                 # Enforce relentlessly that 0% of mutations wander strictly outside the editable zones
            & (scored_df["low_complexity_fraction"] <= 0.20)                                                                                                                                                # Enforce strictly that general sequence complexity organically remains robust
            & (scored_df["mutation_count"] <= int(args.max_mutations))                                                                                                                                      # Enforce identically the strict maximum mutation count ceiling budget parameter
            & (scored_df["sequence_identity"] >= 0.92)                                                                                                                                                      # Enforce uncompromisingly a strict minimum 92% structural identity directly to the seed
        ].copy()                                                                                                                                                                                            # Safely create a disconnected working copy of the fully filtered dataframe
        if scored_df.empty:                                                                                                                                                                                 # Check definitively if all surviving candidates were aggressively eliminated by filters
            break                                                                                                                                                                                           # Immediately exit the beam search loop proactively if no valid candidates survive

        # Embed the remaining candidates and use a diversity-aware pick so the beam does not collapse onto near-duplicate local optima too early.
        embeddings = embed_sequences(scored_df["candidate_sequence"].astype(str).tolist(), model_name=args.esm_model, batch_size=args.batch_size, max_length=args.max_aa)                                   # Re-embed the heavily filtered candidate sequences exclusively to assess spatial diversity
        chosen_idx, diversity_penalty = greedy_diverse_pick(embeddings=embeddings, scores=scored_df["stage09_score"].to_numpy(dtype=np.float32), top_k=args.beam_width, penalty_weight=0.20)                # Select top elite candidates explicitly while strictly penalizing localized structural clustering
        # Add the final diversity penalties and unique IDs of the survivors to the dataframe.
        beam_df = scored_df.iloc[chosen_idx].copy().reset_index(drop=True)                                                                                                                                  # Filter the dataframe securely right down to the final diverse beam algorithmic selection
        beam_df["diversity_penalty"] = diversity_penalty[chosen_idx] if len(diversity_penalty) else 0.0                                                                                                     # Record the strictly applied diversity spatial penalties natively for forensic transparency
        beam_df["candidate_id"] = [f"r{round_idx}_cand_{i+1}" for i in range(len(beam_df))]                                                                                                                 # Generate dynamically unique textual tracking IDs securely for this round's survivors
        # Append the round's elite survivor dataframe cleanly to the master accumulator list.
        all_rows.append(beam_df)                                                                                                                                                                            # Append the round's elite survivor dataframe cleanly to the master accumulator list

        # Convert dataframe rows iteratively into pure dictionaries simply to seed the next round
        beam = beam_df[["candidate_id", "candidate_sequence", "mutations", "round_index"]].to_dict(orient="records")                                                                                        
        # For each of the surviving beam candidates, which are going to be the new parent candidates for the next round,
        # re-parse their mutated positions and assign the correct round index to their dictionary.
        for item in beam:                                                                                                                                                                                   # Iterate strictly over each newly selected parent actively to systematically prep them
            item["mutated_positions"] = parse_mutation_positions(item.get("mutations", []))                                                                                                                 # Safely re-parse and cleanly cache the numeric mutated positions actively back in
            item["round_index"] = round_idx                                                                                                                                                                 # Assign the correct current round chronological index rigidly to the surviving parent

    # Concatenate all scored rounds into one candidate table, sorted by stage09_score, pLDDT, and RMSD
    if not all_rows:                                                                                                                                                                                        # Check critically if the master aggregated list of all rounds fundamentally is entirely empty
        raise ValueError("Stage 09 localized search did not produce any surviving candidates. Relax the hard filters or reduce the mutation budget.")                                                       # Violently raise an error openly indicating the complete collapse of the guided search
    out_df = pd.concat(all_rows, ignore_index=True)                                                                                                                                                         # Swiftly merge all valid collected round dataframes fully into a single master tabular output
    out_df = out_df.sort_values(["stage09_score", "predicted_mean_plddt", "predicted_rmsd"], ascending=[False, False, True]).reset_index(drop=True)                                                         # Sort candidates securely globally primarily by score, secondary pLDDT, and tertiary RMSD
    out_df["sample_id"] = np.arange(1, len(out_df) + 1)                                                                                                                                                     # Assign safely a strictly sequential global numeric rank ID plainly to every candidate
    out_df["generation_regime"] = "stage09_localized_search"                                                                                                                                                # Hardcode permanently the explicit overarching programmatic generation pipeline source tag
    out_df["mutation_positions"] = out_df["mutations"].apply(lambda muts: ";".join(muts) if isinstance(muts, list) else str(muts))                                                                          # Format robustly the raw python mutation list strings strictly into a much cleaner column
    out_df["editable_hotspots"] = ",".join(str(int(x)) for x in edit_space["edit_space"]["editable_positions"])                                                                                             # Record textually the globally statically defined natively editable hotspots rigidly in the dataframe

    out_path = Path(args.out_csv)                                                                                                                                                                           # Swiftly create a robust pathlib Path object properly for the target CSV destination output file
    out_path.parent.mkdir(parents=True, exist_ok=True)                                                                                                                                                      # Relentlessly ensure the entire parent directory hierarchical structure strictly exists before saving
    out_df.to_csv(out_path, index=False)                                                                                                                                                                    # Powerfully export the complete calculated overarching candidate table safely to the CSV file

    # Write run metadata for the downstream prefilter stage
    meta = {                                                                                                                                                                                                # Powerfully open the master operational metadata tracking foundational dictionary block
        "stage": "09d",                                                                                                                                                                                     # Record stably the strict programmatic operating stage structural identifier string
        "context_json": str(args.context_json),                                                                                                                                                             # Record cleanly the specified explicit input configuration context JSON lookup path
        "edit_space_json": str(args.edit_space_json),                                                                                                                                                       # Record natively the strictly specified input localized edit-space definition JSON path
        "predictor_model": str(args.predictor_model),                                                                                                                                                       # Record natively the rigidly explicitly designated overarching input predictor model path
        "surrogate_model": str(args.surrogate_model) if args.surrogate_model else None,                                                                                                                     # Record safely the conditionally provided external input surrogate active model path
        "target_host": target_host,                                                                                                                                                                         # Record structurally the utilized primary biological target host optimization label
        "n_candidates": int(len(out_df)),                                                                                                                                                                   # Record natively the final absolute total structural number of valid candidates generated
        "rounds_completed": int(max(out_df["round_index"])) if not out_df.empty else 0,                                                                                                                     # Record explicitly strictly precisely how many search rounds definitively completed successfully
        "score_weights": {                                                                                                                                                                                  # Open deeply the secondary foundational explicit mathematical score weights sub-dictionary
            "w_target": float(args.w_target),                                                                                                                                                               # Record safely the final heavily configured user primary target algorithmic weight multiplier
            "w_family": float(args.w_family),                                                                                                                                                               # Record natively the final carefully calibrated secondary family centroid cosine metric weight
            "w_seed": float(args.w_seed),                                                                                                                                                                   # Record strictly the final completely applied primary seed cosine distance measurement weight
            "w_guidance": float(args.w_guidance),                                                                                                                                                           # Record securely the final rigorously assigned strict edit guidance programmatic matrix weight
            "w_surrogate": float(args.w_surrogate),                                                                                                                                                         # Record natively the explicitly utilized final surrogate structural risk safety penalty weight
        },                                                                                                                                                                                                  # Close heavily the secondary explicit mathematical score weights tracking sub-dictionary
    }                                                                                                                                                                                                       # Close cleanly the overarching fully operational foundational metadata robust dictionary block
    write_json(meta, args.out_json)                                                                                                                                                                         # Commandingly save the finalized comprehensively configured metadata dictionary cleanly to a JSON file
    print(f"Wrote: {out_path}")                                                                                                                                                                             # Visually print a highly transparent clean execution confirmation message explicitly for the CSV file
    print(f"Wrote: {args.out_json}")                                                                                                                                                                        # Visually print a structurally secure overarching execution confirmation message strictly for the JSON file
    print(out_df[["sample_id", "round_index", "stage09_score", "target_probability", "predicted_mean_plddt", "predicted_rmsd", "sequence_identity", "mutation_count"]].head(10).to_string(index=False))     # Commandingly print immediately directly the absolute top 10 best elite candidates' critical metric stats


if __name__ == "__main__":                                                                                                                                                                                  # Standard primary foundational python global programmatic execution validation guard check
    main()                                                                                                                                                                                                  # Definitively cleanly definitively strictly invoke the primary main function to actively start the script