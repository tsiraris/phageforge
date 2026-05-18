#!/usr/bin/env python
"""Stage 09e: Prefilter Stage 09 localized-search candidates before expensive structural validation.

This script keeps a compact top panel that already satisfies the main Stage 09 design goals:
strong target score, acceptable surrogate structural risk, limited seed drift, and sequence diversity.
The output CSV is designed to be directly consumable by the Stage 09 validation wrapper.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from phageforge.stage07_utils import embed_sequences
from phageforge.stage09_utils import greedy_diverse_pick, read_json


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 09 structural prefilter."""
    ap = argparse.ArgumentParser(description="Prefilter Stage 09 candidates before expensive structural validation.")                                 # Initialize argument parser with a description
    ap.add_argument("--search_csv", type=str, required=True, help="Full Stage 09 search CSV produced by 09d_localized_search.py.")                    # Define required argument for the input search results CSV
    ap.add_argument("--search_meta_json", type=str, required=True, help="Run metadata JSON produced by 09d_localized_search.py.")                     # Define required argument for the input metadata JSON
    ap.add_argument("--out_csv", type=str, required=True, help="Where to write the Stage 09 prefitered candidate CSV.")                               # Define required argument for the output CSV path
    ap.add_argument("--top_k", type=int, default=12, help="How many candidates to keep after structural prefiltering.")                               # Define optional argument for number of candidates to output
    ap.add_argument("--max_structural_risk", type=float, default=0.55, help="Maximum allowed surrogate structural risk.")                             # Define optional argument for the structural risk threshold
    ap.add_argument("--min_predicted_plddt", type=float, default=55.0, help="Minimum predicted mean pLDDT before expensive validation.")              # Define optional argument for the minimum allowed pLDDT score
    ap.add_argument("--max_predicted_rmsd", type=float, default=4.5, help="Maximum predicted RMSD to the seed before expensive validation.")          # Define optional argument for the maximum allowed RMSD drift
    ap.add_argument("--min_sequence_identity", type=float, default=0.93, help="Minimum allowed seed identity for Stage 09 candidates.")               # Define optional argument for the minimum sequence identity
    ap.add_argument("--esm_model", type=str, default="facebook/esm2_t12_35M_UR50D", help="Embedding model used for diversity-aware panel selection.") # Define optional argument for the ESM language model name
    ap.add_argument("--batch_size", type=int, default=4, help="Batch size used for diversity embeddings.")                                            # Define optional argument for the embedding inference batch size
    ap.add_argument("--max_aa", type=int, default=2048, help="Maximum sequence length passed into the diversity embedding model.")                    # Define optional argument for the maximum sequence length limit
    return ap.parse_args()                                                                                                                            # Parse and return the populated arguments namespace


def main() -> None:
    # Read the Localized Search (09d) candidate table and the metadata describing the search that produced it.
    args = parse_args()                                           # Parse the command-line arguments into the args variable
    search_df = pd.read_csv(args.search_csv)                      # Load the candidate sequences and metrics into a Pandas DataFrame
    search_meta = read_json(args.search_meta_json)                # Load the search metadata dictionary from the associated JSON file

    # Apply the main hard structural prefilters before diversity selection so obviously risky candidates are dropped early.
    keep_df = search_df.loc[                                                                                                                            # Start filtering the DataFrame rows based on multiple conditions
        (search_df["predicted_structural_risk"] <= float(args.max_structural_risk))                                                                     # Keep rows where structural risk is at or below the maximum allowed
        & (search_df["predicted_mean_plddt"] >= float(args.min_predicted_plddt))                                                                        # Keep rows where the predicted pLDDT score is at or above the minimum
        & (search_df["predicted_rmsd"] <= float(args.max_predicted_rmsd))                                                                               # Keep rows where the predicted RMSD is at or below the maximum allowed
        & (search_df["sequence_identity"] >= float(args.min_sequence_identity))                                                                         # Keep rows maintaining the minimum sequence identity to the original seed
        & (search_df["outside_editable_fraction"] <= 0.0)                                                                                               # Keep rows that only contain mutations within the allowed editable regions
    ].copy()                                                                                                                                            # Create an independent copy of the filtered DataFrame subset
    if keep_df.empty:                                                                                                                                   # Check if all candidates were removed by the filtering steps
        raise ValueError("No Stage 09 candidates survived the structural prefilter. Relax the thresholds or inspect the search outputs.")               # Raise an error if the resulting filtered DataFrame is empty

    # Embedd the candidate sequences and produce a diverse panel so the expensive structural validator sees several distinct local optima rather than near duplicates.
    embeddings = embed_sequences(keep_df["candidate_sequence"].astype(str).tolist(), model_name=args.esm_model, batch_size=args.batch_size, max_length=args.max_aa)                 # Generate ESM embeddings for the filtered candidate sequences
    chosen_idx, diversity_penalty = greedy_diverse_pick(embeddings=embeddings, scores=keep_df["stage09_score"].to_numpy(dtype=np.float32), top_k=args.top_k, penalty_weight=0.25)   # Select a diverse top-k subset using the embeddings and target scores
    out_df = keep_df.iloc[chosen_idx].copy().reset_index(drop=True)                                                                                                                 # Create the final DataFrame with the selected diverse candidates, resetting indices
    # For each candidate, record an applied diversity penalty, a numerical ranking from 1 to top_k, a tag indicating they passed the prefilter, and the score into the final ranking score column 
    out_df["prefilter_diversity_penalty"] = diversity_penalty[chosen_idx] if len(diversity_penalty) else 0.0                                                                        # Record the applied diversity penalty for each chosen candidate
    out_df["prefilter_rank"] = np.arange(1, len(out_df) + 1)                                                                                                                        # Assign a numerical ranking to each candidate from 1 to top_k
    out_df["prefilter_reason"] = "kept_by_stage09_structural_prefilter"                                                                                                             # Tag the rows indicating they successfully passed this prefilter stage
    out_df["final_multimodal_rank_score"] = out_df["stage09_score"]                                                                                                                 # Copy the stage 09 score into the final ranking score column

    # Write the prefiltered panel in a validator-compatible format that preserves candidate_sequence and mutation annotations.
    out_path = Path(args.out_csv)                                                                                                                                                   # Create a Path object for the output CSV file destination
    out_path.parent.mkdir(parents=True, exist_ok=True)                                                                                                                              # Ensure the destination directory exists, creating parents if needed
    out_df.to_csv(out_path, index=False)                                                                                                                                            # Save the final prefiltered candidates to the output CSV file without indices
    print(f"Wrote: {out_path}")                                                                                                                                                     # Print a confirmation message to standard output with the file path
    print(f"source_target_host: {search_meta.get('target_host', 'na')}")                                                                                                            # Print the target host extracted from the search metadata
    print(out_df[["prefilter_rank", "sample_id", "stage09_score", "target_probability", "predicted_mean_plddt", "predicted_rmsd", "sequence_identity"]].to_string(index=False))     # Print a summary table of the selected candidates and their key metrics


if __name__ == "__main__": # Entry point check to ensure the script is being run directly
    main()                 # Call the main function to execute the full prefiltering process