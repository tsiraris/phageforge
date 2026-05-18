#!/usr/bin/env python
"""Stage 09b: Build a structural-surrogate training table from Stage 07 and Stage 08 artifacts.

This script merges the Stage 07 ranked candidate table with Stage 08 structural outcomes and
creates a compact modeling table for structural-risk estimation. It is intentionally permissive:
if only a few structurally labeled candidates exist, the output is still useful as a rule-based
calibration table and as a record of which sequence patterns already failed.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import json
import pandas as pd
from phageforge.stage09_utils import build_basic_sequence_features, read_json


FEATURE_COLUMNS = [                                                                    # Define the expected sequence and scoring feature columns
    "mutation_count",
    "mutation_span",
    "normalized_entropy",
    "max_single_residue_fraction",
    "longest_homopolymer_run",
    "low_complexity_fraction",
    "outside_editable_fraction",
    "seed_identity",
    "final_multimodal_rank_score",
    "strict_manifold_score",
    "structure_score",
    "target_score",
    "family_cosine",
    "seed_cosine",
    "target_anchor_cosine",
]



def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for building the Stage 09 surrogate dataset."""
    ap = argparse.ArgumentParser(description="Build the Stage 09 structural-surrogate feature table.")                         # Initialize the CLI argument parser with description
    ap.add_argument("--context_json", type=str, required=True, help="Stage 07 context JSON so sequence metrics can...")        # Define required context JSON argument
    ap.add_argument("--ranked_csv", type=str, required=True, help="Stage 07 ranked candidate CSV.")                            # Define required ranked CSV argument
    ap.add_argument("--structural_csv", type=str, required=True, help="Stage 08 structural summary CSV.")                      # Define required structural CSV argument
    ap.add_argument("--out_csv", type=str, required=True, help="Where to write the Stage 09 surrogate dataset CSV.")           # Define required output CSV argument
    ap.add_argument("--out_json", type=str, required=True, help="Where to write a short JSON summary of the dataset coverage.")# Define required output JSON argument
    return ap.parse_args()                                                                                                     # Execute parsing and return the populated arguments namespace



def main() -> None:
    # Read the Stage 07 ranking table, the Stage 08 structural outcomes, and the selected seed sequence from context.
    args = parse_args()                                                                                                      # Call the argument parser to get inputs
    context = read_json(args.context_json)                                                                                   # Load context JSON file into a Python dictionary
    ranked_df = pd.read_csv(args.ranked_csv)                                                                                 # Load ranked candidates CSV into a pandas DataFrame
    structural_df = pd.read_csv(args.structural_csv)                                                                         # Load structural outcomes CSV into a pandas DataFrame
    seed_sequence = str(context["selected_seed"]["seed_sequence"])                                                           # Extract the seed sequence as a string from the dictionary
    editable_positions = set(int(x) for x in context["editable_region"].get("hotspot_positions", []))                        # Extract hotspot string positions and cast to a set of integers

    # Select certain columns from the ranked candidates table (07e), merge them with the (physically folded proteins) structural_csv of 08a, so every labeled Stage 08 candidate carries its upstream metadata.
    merge_cols = [c for c in ["sample_id", "generation_regime", "final_multimodal_rank_score", "target_score"] if c in ranked_df.columns and c in structural_df.columns] # Intersect available columns for merging
    merged = structural_df.merge(ranked_df.drop(columns=["candidate_sequence"], errors="ignore"), on=merge_cols, how="left")                                             # Merge dfs, ignoring original sequence column to avoid duplication
    if "candidate_sequence" not in merged.columns and "candidate_sequence" in structural_df.columns:                                                                     # Check if we dropped sequence entirely but have it in structural
        merged["candidate_sequence"] = structural_df["candidate_sequence"]                                                                                               # Restore sequence column from structural df
    if "candidate_sequence" not in merged.columns:                                                                                                                       # Check if sequence column is still missing
        fallback = ranked_df[["sample_id", "candidate_sequence"]].drop_duplicates("sample_id")                                                                           # Extract an unambiguous sample_id to sequence map from ranked df
        merged = merged.merge(fallback, on="sample_id", how="left")                                                                                                      # Merge the fallback to ensure sequences exist

    # Compute compact sequence diagnostics so the surrogate can learn which edit patterns tend to fail structurally.
    rows = []                                                                                                                                                            # Initialize list to hold final dataset rows
    # Iterate over each row in the merged dataframe converted to a dictionary
    for row in merged.to_dict(orient="records"):                                                                                                                         
        # Retrieve the candidate sequence as a string, falling back to an empty string safely
        candidate_sequence = str(row.get("candidate_sequence", "") or "")                                                                                                # Retrieve sequence as string, falling back to empty string safely
        # For every row in the merged dictionary, compare the candidate sequence to the original wild-type seed and mathematically quantify its biological grammar (e.g., entropy, repetitiveness, mutation spread
        feature_row = build_basic_sequence_features(seed_sequence=seed_sequence, candidate_sequence=candidate_sequence, editable_positions=editable_positions)           # Compute biological sequence differences and features
        # Merge sequence features with structural metadata dict
        feature_row.update(                                                                                                                                              
            {
                "sample_id": row.get("sample_id"),                                                                                                                       # Attach sample identifier
                "generation_regime": row.get("generation_regime", "na"),                                                                                                 # Attach generative model origin string, default to "na"
                "candidate_sequence": candidate_sequence,                                                                                                                # Attach the raw candidate sequence string
                "stage08_pass": bool(row.get("stage08_pass", False)),                                                                                                    # Attach strict boolean pass/fail status
                "esmfold_mean_plddt": float(row.get("esmfold_mean_plddt", float("nan"))),                                                                                # Attach overall predicted confidence score as float
                "mutation_site_mean_plddt": float(row.get("mutation_site_mean_plddt", float("nan"))),                                                                    # Attach local predicted confidence score as float
                "rmsd_to_selected_seed": float(row.get("rmsd_to_selected_seed", float("nan"))),                                                                          # Attach spatial deviation from seed as float
                "stage08_decision_reason": str(row.get("stage08_decision_reason", "")),                                                                                  # Attach textual reason for pass/fail decision
            }
        )
        # Iterate through the global requirement columns 
        for col in FEATURE_COLUMNS:                                                                                                                                      # Iterate through the global requirement columns
            # If the feature has already been populated, skip to the next iteration
            if col in feature_row:                                                                                                                                       # Check if the feature has already been populated
                continue                                                                                                                                                 # Skip to the next iteration if so
            # If the feature has not been populated, extract its raw value, converting NaNs and nulls safely to 0.0 float
            feature_row[col] = float(row.get(col, 0.0)) if pd.notna(row.get(col, 0.0)) else 0.0                                                                          # Extract raw value, converting NaNs and nulls safely to 0.0 float
        rows.append(feature_row)                                                                                                                                         # Push the fully populated dictionary row to the array

    # Write the modeling table and a concise coverage summary that says how much labeled structural supervision is available.
    out_df = pd.DataFrame(rows)                                                                                          # Convert the final list of dictionaries back into a pandas DataFrame
    out_path = Path(args.out_csv)                                                                                        # Wrap the output csv string argument into a Path object
    out_path.parent.mkdir(parents=True, exist_ok=True)                                                                   # Create the directory tree for the file, suppressing exists errors
    out_df.to_csv(out_path, index=False)                                                                                 # Export the dataframe to a CSV file, ignoring the integer index

    summary = {                                                                                                          # Initialize dictionary to hold dataset overview statistics
        "n_rows": int(len(out_df)),                                                                                      # Calculate integer length of total rows
        "n_pass": int(out_df["stage08_pass"].sum()) if not out_df.empty else 0,                                          # Sum passed entries if not empty, otherwise set 0
        "n_fail": int((~out_df["stage08_pass"]).sum()) if not out_df.empty else 0,                                       # Sum failed entries if not empty, otherwise set 0
        "feature_columns": FEATURE_COLUMNS,                                                                              # Append list of active feature columns used
        "source_ranked_csv": str(args.ranked_csv),                                                                       # Log the explicit ranked CSV filepath origin
        "source_structural_csv": str(args.structural_csv),                                                               # Log the explicit structural CSV filepath origin
    }
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)                                                        # Create the directory tree for the JSON file
    Path(args.out_json).write_text(json.dumps(summary, indent=2))                                                        # Serialize the dictionary to a nicely-formatted JSON string and save
    print(f"Wrote: {out_path}")                                                                                          # Output the CSV path written to stdout
    print(f"Wrote: {args.out_json}")                                                                                     # Output the JSON path written to stdout
    print(f"rows: {len(out_df)} | passes: {summary['n_pass']} | fails: {summary['n_fail']}")                             # Output basic descriptive dataset statistics to stdout


if __name__ == "__main__":                      # Standard boilerplate check to ensure we only run when executed directly
    main()                                      # Execute the main entrypoint loop