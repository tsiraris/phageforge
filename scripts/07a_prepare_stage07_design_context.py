"""
Stage 07a: Build the compact design context for the final Acinetobacter retargeting step.

This script reads the outputs from previous stages (the validated seed sequence and host targets). 
It then defines the "editable window" by extracting functional hotspots (positions that tolerate mutations and affect binding). 
It writes all of this into a reusable JSON file (stage07_context.json). 
This ensures the generative model knows exactly what it is allowed to process.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from phageforge.stage07_utils import (
    build_position_feature_table,
    build_structured_windows,
    read_json,
    write_json,
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 07 context builder."""
    ap = argparse.ArgumentParser(description="Prepare the Stage 07 design context JSON.")                                           # Initializes the argument parser with a description of the script's purpose.
    ap.add_argument("--phaseA_plan_json", type=str, required=True, help="Plan JSON from 06a_select_phaseA_family.py.")              # Adds argument for the Phase A planning JSON, which contains family info and hotspots.
    ap.add_argument("--phase06c_followup_summary_json", type=str, required=True, help="Follow-up seed JSON from 06c...")            # Adds argument for the Phase 06c JSON, which defines the chosen candidate seed.
    ap.add_argument("--strict_csv", type=str, required=True, help="Strict processed RBP dataset used to recover family...")         # Adds argument for the filtered dataset CSV to extract structural/family contexts.
    ap.add_argument("--target_host", type=str, required=True, help="Target host genus for Stage 07, such as Acinetobacter.")        # Adds argument for the target host genus (e.g., Acinetobacter) to guide retargeting.
    ap.add_argument("--output_json", type=str, required=True, help="Where to write the Stage 07 context JSON.")                     # Adds argument for the destination path where the final context JSON will be saved.
    return ap.parse_args()                                                                                                          # Parses the provided command-line arguments and returns them as a Namespace object.


def main() -> None:
    # Read the upstream planning artifacts and the strict family dataset.
    args = parse_args()                                                                                             # Executes the argument parsing function and stores the result.
    plan = read_json(args.phaseA_plan_json)                                                                         # Loads the Phase A plan JSON into a Python dictionary.
    followup = read_json(args.phase06c_followup_summary_json)                                                       # Loads the 06c follow-up JSON containing seed choices into a dictionary.
    strict_df = pd.read_csv(args.strict_csv)                                                                        # Reads the strict RBP dataset into a Pandas DataFrame for context retrieval.

    # Recover the full aa seed sequence of the selected seed row from the previous ladder step to propagate into Stage 07.
    source_top_candidates_csv = Path(followup["source_top_candidates_csv"])                                         # Extracts and wraps the file path for the source candidates CSV from the followup dict.
    if not source_top_candidates_csv.exists():                                                                      # Checks if the extracted candidate CSV file path actually exists on the filesystem.
        raise FileNotFoundError(f"Missing file: {source_top_candidates_csv}")                                       # Raises an error stopping execution if the candidate CSV file is missing.
    source_df = pd.read_csv(source_top_candidates_csv)                                                              # Reads the candidate CSV into a DataFrame to access sequence and metadata.
    chosen_id = str(followup["chosen_candidate_id"])                                                                # Extracts the specific candidate ID chosen in the previous stage as a string.
    selected_rows = source_df.loc[source_df["candidate_id"].astype(str) == chosen_id].reset_index(drop=True)        # Filters the DataFrame for rows matching the chosen ID and resets the index cleanly.
    if selected_rows.empty:                                                                                         # Checks if the filtered DataFrame is empty (meaning the chosen ID was not found).
        raise ValueError(f"Could not find chosen candidate_id={chosen_id} in {source_top_candidates_csv}")          # Raises a ValueError if the chosen candidate ID is missing from the source DataFrame.
    selected_row = selected_rows.iloc[0]                                                                            # Selects the first (and supposedly only) matching row representing our target seed protein.
    seed_sequence = str(selected_row["aa_sequence"])                                                                # Extracts the full amino acid sequence of the chosen seed protein as a string.

    # Recover family sequences from the strict dataset so the family context remains grounded in real proteins.
    family_member_ids = list(plan.get("family_member_ids", []))                                                                              # Retrieves the list of protein IDs belonging to the same family from the plan dictionary.
    family_rows = strict_df.loc[strict_df["protein_id"].astype(str).isin([str(x) for x in family_member_ids])].copy().reset_index(drop=True) # Filters the strict dataset for family member IDs, copies it, and resets the index.
    family_sequences = family_rows.get("aa_sequence", pd.Series(dtype=str)).astype(str).tolist()                                             # Extracts the amino acid sequences for all matching family members into a Python list.

    # Build a functionally informed editable region (functional mutation hotspots) using family variability and the target host position priors.
    hotspots_0based = [int(x) for x in plan.get("mutation_window_positions_0based", [])]                                                    # Retrieves 0-indexed mutation hotspots from the plan and casts them to integers.
    hotspots_1based = [pos + 1 for pos in hotspots_0based]                                                                                  # Converts the 0-indexed hotspot positions to 1-indexed format for biological standard.
    position_features = build_position_feature_table(  # Calls a utility function to calculate functional importance features for each position.
        seed_sequence=seed_sequence,                                                                                                        # Passes the target seed amino acid sequence to the feature builder.
        family_sequences=family_sequences,                                                                                                  # Passes the list of family sequences to compute positional variability and conservation.
        hotspots_1based=hotspots_1based,                                                                                                    # Passes the 1-indexed hotspot locations defining the permitted mutation window.
        target_position_priors=plan.get("target_position_priors", {}).get(args.target_host, []),                                            # Passes any existing positional biases or priors specific to the target host.
    )                                                                                                                                       # Closes the build_position_feature_table function call.
    # Rank the positions by their functional importance, assign a probability to each, select a max number of them to mask simultaneously.
    ranked_hotspots = [int(row["position"]) for row in position_features]                                                                   # Extracts and casts the ranked position indices from the computed features table.
    priority_weights = {str(int(row["position"])): float(row["functional_weight"]) for row in position_features}                            # Creates a dictionary mapping position strings to their computed functional weights.
    default_mask_count = max(24, int(plan.get("hotspot_min_count", 24)))                                                                    # Determines the maximum number of residues to mask simultaneously, defaulting to at least 24.
    # Create continuous regions of hotspots to mask, around the top_k most important positions, with a window size equal to the default_mask_count.
    structured_windows = build_structured_windows(position_features, seed_sequence=seed_sequence, window_size=default_mask_count, top_k=3)  
    # Define the overall editable region based on the min and max active hotspot positions, to lock out the rest of the sequence and save computational memory.
    active_positions = ranked_hotspots or hotspots_1based                                                                                   # Selects ranked hotspots if available, falling back to the raw 1-indexed hotspots otherwise.
    window_start = min(active_positions) if active_positions else 1                                                                         # Finds the lowest position index to define the start of the overall editable region.
    window_end = max(active_positions) + 1 if active_positions else len(seed_sequence) + 1                                                  # Finds the highest position index (+1) to define the end of the editable region.

    # Build the context dictionary that will become the stage07_context.json file.
    context = {                                                                                                                                     # Initializes the master dictionary that will become the stage07_context.json file.
        "stage": "07",                                                                                                                              # Hardcodes the pipeline stage identifier for tracking.
        "target_host": args.target_host,                                                                                                            # Records the targeted bacterial genus requested via command line.
        "canonical_seed": plan.get("canonical_seed", {}),                                                                                           # Copies the canonical seed metadata inherited from the phase A plan.
        "selected_seed": {                                                                                                                          # Opens a sub-dictionary to define the exact seed protein chosen for this generation run.
            "seed_rank": 0,                                                                                                                         # Sets an explicit baseline rank for the seed (useful when appending generated candidates).
            "seed_protein_id": str(selected_row.get("candidate_id", chosen_id)),                                                                    # Records the unique identifier of the selected seed protein.
            "seed_identifier_hint": str(selected_row.get("seed_protein_id", plan.get("canonical_seed", {}).get("seed_protein_id", ""))),            # Provides a fallback or original identifier hint for cross-referencing.
            "seed_source_kind": "06c_followup_summary_json",                                                                                        # Documents the specific pipeline stage/file that generated this seed choice.
            "seed_source_desc": f"candidate_id={chosen_id}",                                                                                        # Creates a human-readable description string of the seed's origin.
            "virus_accession": str(selected_row.get("candidate_id", selected_row.get("virus_accession", ""))),                                      # Extracts the viral accession number associated with the seed protein if available.
            "source_host": str(selected_row.get("source_host", plan.get("canonical_seed", {}).get("source_host", ""))),                             # Records the original bacterial host that this seed protein naturally infects.
            "seed_sequence": seed_sequence,                                                                                                         # Stores the exact amino acid sequence to be used as the generative starting point.
            "sequence_length": int(len(seed_sequence)),                                                                                             # Calculates and stores the total length of the seed amino acid sequence.
        },                                                                                                                                          # Closes the selected_seed sub-dictionary.
        "family_context": {                                                                                                                         # Opens a sub-dictionary containing aggregated structural and evolutionary family data.
            "family_member_count": int(len(family_rows)),                                                                                           # Counts and records the total number of valid family members found.
            "family_cosine_floor": float(plan.get("family_summary", {}).get("family_cosine_floor", 0.995)),                                         # Records the minimum cosine similarity threshold used to define this family cluster.
            "family_product_majority": str(family_rows["product"].mode().iat[0] if len(family_rows) else plan.get("canonical_seed", {}).get("product", "receptor-binding protein")), # Determines the most common protein product annotation among family members.
            "family_member_ids": family_member_ids,                                                                                                 # Stores the list of all protein IDs that make up the evolutionary family context.
            "family_centroid": plan.get("family_centroid", []),                                                                                     # Copies the embedding vector representing the center of the family cluster.
            "family_rows": family_rows.to_dict(orient="records"),                                                                                   # Serializes the entire subset of family DataFrames into a list of dictionaries.
        },                                                                                                                                          # Closes the family_context sub-dictionary.
        "target_context": {                                                                                                                         # Opens a sub-dictionary defining the properties of the intended new host target.
            "target_centroid": plan.get("target_reference_centroids", {}).get(args.target_host, []),                                                # Retrieves the central embedding vector representing known RBPs for the target host.
            "target_reference_count": int(len(plan.get("target_reference_rows", {}).get(args.target_host, []))) if isinstance(plan.get("target_reference_rows", {}).get(args.target_host, []), list) else int(plan.get("target_reference_rows", {}).get(args.target_host, 0)), # Safely calculates the number of reference sequences used to define the target host.
        },                                                                                                                                          # Closes the target_context sub-dictionary.
        "editable_region": {                                                                                                                        # Opens a sub-dictionary specifying exactly where the generative model is allowed to mutate.
            "window_start": int(window_start),                                                                                                      # Defines the 1-indexed start position of the broad editable window on the protein.
            "window_end": int(window_end),                                                                                                          # Defines the 1-indexed end position of the broad editable window on the protein.
            "hotspot_positions": ranked_hotspots,                                                                                                   # Stores the list of specific 1-indexed positions permitted for mutation, ranked by importance.
            "hotspot_priority_weights": priority_weights,                                                                                           # Stores the mapping of position indices to their computed evolutionary/functional weights.
            "position_features": position_features,                                                                                                 # Serializes the detailed tabular features calculated for every allowed position.
            "structured_windows": structured_windows,                                                                                               # Stores the chunked, contiguous regions of hotspots formulated for masked language modeling.
            "target_position_priors": plan.get("target_position_priors", {}).get(args.target_host, []),                                             # Copies any host-specific positional biases that guide the mutation strategy.
            "default_max_masked_positions": int(default_mask_count),                                                                                # Records the maximum number of residues that should be unmasked in a single generative pass.
        },                                                                                                                                          # Closes the editable_region sub-dictionary.
        "upstream_artifacts": {                                                                                                                     # Opens a sub-dictionary to maintain strict provenance of the files used in this run.
            "phaseA_plan_json": str(args.phaseA_plan_json),                                                                                         # Records the file path of the Phase A plan used.
            "phase06c_followup_summary_json": str(args.phase06c_followup_summary_json),                                                             # Records the file path of the Phase 06c followup summary used.
            "source_top_candidates_csv": str(source_top_candidates_csv),                                                                            # Records the file path of the upstream CSV where the seed sequence was actually extracted.
            "strict_csv": str(args.strict_csv),                                                                                                     # Records the file path of the raw dataset used to build the family context.
        },                                                                                                                                          # Closes the upstream_artifacts sub-dictionary.
    }                                                                                                                                               # Closes the master context dictionary.

    write_json(context, args.output_json)                                                                                           # Writes the compiled context dictionary to the filesystem as a formatted JSON file.
    print(f"Wrote: {args.output_json}")                                                                                             # Prints a confirmation message to standard output indicating successful creation.


if __name__ == "__main__":                                                                                                          # Standard Python idiom to check if the script is being executed directly rather than imported.
    main()                                                                                                                          # Calls the main execution function to run the context preparation logic.