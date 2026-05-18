#!/usr/bin/env python
"""Stage 10c: Prefilter structure-conditioned candidates before expensive full validation.

This script is intentionally conservative. It keeps candidates only if they remain strong
under the new Stage 10 composite score and then applies a diversity-aware selection step.
The purpose is not to replace decisive structural validation, but to ensure the expensive
validator only sees the most scaffold-compatible and target-relevant candidates.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from phageforge.stage07_utils import embed_sequences
from phageforge.stage10_utils import greedy_diverse_subset, read_json, write_json


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 10 prefilter step."""
    ap = argparse.ArgumentParser(description="Prefilter and diversify Stage 10 redesign candidates.")                                                                 # Initialize the argument parser with a description
    ap.add_argument("--stage10_context_json", type=str, required=True, help="Stage 10 context JSON used to recover the mutation budget and target settings.")         # Define the required context JSON file argument
    ap.add_argument("--search_csv", type=str, required=True, help="Full search CSV written by 10b_run_inverse_folding_beam_search.py.")                               # Define the required search CSV file argument
    ap.add_argument("--embedding_model", type=str, default="facebook/esm2_t33_650M_UR50D", help="Embedding backbone used for diversity-aware filtering.")             # Define the default ESM embedding model to use
    ap.add_argument("--batch_size", type=int, default=4, help="Batch size used for candidate embedding during diversity filtering.")                                  # Define the batch size for embedding calculations
    ap.add_argument("--top_k", type=int, default=10, help="Number of candidates to keep in the main prefilter panel.")                                                # Define the extended panel top_k limit
    ap.add_argument("--top_k_final", type=int, default=3, help="Number of candidates to keep in the compact final validation panel.")                                 # Define the elite panel top_k limit
    ap.add_argument("--out_topk_csv", type=str, required=True, help="Where to write the top-k prefilter CSV.")                                                        # Define the required output top-k CSV file path
    ap.add_argument("--out_topk_final_csv", type=str, required=True, help="Where to write the compact final validation CSV.")                                         # Define the required output elite CSV file path
    ap.add_argument("--out_json", type=str, required=True, help="Where to write the compact prefilter summary JSON.")                                                 # Define the required output metadata JSON path
    return ap.parse_args()                                                                                                                                            # Parse the provided command-line arguments and return them


def main() -> None:
    # Read the search table and the Stage 10 redesign budget so only mutation counts inside the planned envelope are retained.
    args = parse_args()                                                                                                                                               # Extract all command-line configurations
    context = read_json(args.stage10_context_json)                                                                                                                    # Load the master structural constraints dict into memory
    search_df = pd.read_csv(args.search_csv)                                                                                                                          # Load the full, unfiltered Beam Search results table

    editable_region = dict(context["editable_region"])                                                                                                                # Isolate the specific mutation-rule dictionary
    min_mut = int(editable_region["min_mutations"])                                                                                                                   # Extract the minimum required mutation count limit
    max_mut = int(editable_region["max_mutations"])                                                                                                                   # Extract the maximum allowed mutation count limit

    # Apply a first-pass conservative filter that prioritizes the new structure-conditioned score and valid local mutation budgets.
    keep = search_df.copy()                                                                                                                                           # Create a safe, independent copy of the search dataframe
    if "mutation_count" in keep.columns:                                                                                                                              # Check if the dataframe contains the mutation tracking column
        keep = keep.loc[(keep["mutation_count"] >= min_mut) & (keep["mutation_count"] <= max_mut)].copy()                                                             # Mask and retain only rows strictly obeying the physical mutation budget
    if keep.empty:                                                                                                                                                    # Check if the budget filter completely erased all candidate data
        raise ValueError("Stage 10 prefilter retained zero candidates after applying the mutation-budget filter.")                                                    # Abort execution to prevent silent downstream failures

    # Remove exact duplicate sequences and keep only the best-scoring representative for each candidate sequence.
    keep = keep.sort_values(["stage10_composite_score", "target_probability", "if1_log_likelihood"], ascending=False)                                                 # Sort candidates hierarchically by the inverse-folding performance metrics
    keep = keep.drop_duplicates(subset=["candidate_sequence"], keep="first").reset_index(drop=True)                                                                   # Drop exact string duplicates, keeping the first (highest scoring) instance

    # Embed the surviving candidates and apply a diversity-aware selection, penalizing sequences that are geometrically similar, keeping a highly diverse top-k pool.
    sequences = keep["candidate_sequence"].astype(str).tolist()                                                                                                       # Extract the remaining unique sequence strings into a flat Python list
    embeddings = embed_sequences(sequences, model_name=args.embedding_model, batch_size=args.batch_size)                                                              # Run the heavy ESM-2 transformer to generate 3D embedding vectors for every candidate
    topk_idx = greedy_diverse_subset(embeddings, keep["stage10_composite_score"].to_numpy(dtype=np.float32), top_k=args.top_k)                                        # Execute the spatial clustering penalty algorithm to draft the top 10 most diverse leaders
    topk_df = keep.iloc[topk_idx].sort_values(["stage10_composite_score", "target_probability", "if1_log_likelihood"], ascending=False).reset_index(drop=True)        # Slice the dataframe to retain only those leaders, and re-sort them

    # Re-run the diversity-aware selection on the top-k pool to produce the compact final panel that is going to be used for expensive full validation.
    topk_embeddings = embed_sequences(topk_df["candidate_sequence"].astype(str).tolist(), model_name=args.embedding_model, batch_size=args.batch_size)                # Re-embed specifically the Top-10 panel to prepare for the secondary elite draft
    final_idx = greedy_diverse_subset(topk_embeddings, topk_df["stage10_composite_score"].to_numpy(dtype=np.float32), top_k=args.top_k_final)                         # Execute the spatial algorithm again to isolate the absolute top 3 most distinct theories
    final_df = topk_df.iloc[final_idx].sort_values(["stage10_composite_score", "target_probability", "if1_log_likelihood"], ascending=False).reset_index(drop=True)   # Slice the dataframe to retain only the top 3, and re-sort them

    # Add a compact 1-index sample identifier so the downstream validator and final report can reference candidates consistently.
    topk_df = topk_df.copy()                                                                                                                                          # Prevent pandas SettingWithCopy warnings by cloning the top-10 frame
    final_df = final_df.copy()                                                                                                                                        # Prevent pandas SettingWithCopy warnings by cloning the top-3 frame
    topk_df.insert(0, "sample_id", np.arange(1, len(topk_df) + 1))                                                                                                    # Programmatically inject a clean 1-based integer index column at the start of the frame
    final_df.insert(0, "sample_id", np.arange(1, len(final_df) + 1))                                                                                                  # Programmatically inject a clean 1-based integer index column at the start of the frame

    # Write the full top-k panel csv, the final compact validation panel csv, and a small JSON summary for notebooks and reports.
    out_topk = Path(args.out_topk_csv)                                                                                                                                # Construct the file path object for the top-10 CSV export
    out_topk.parent.mkdir(parents=True, exist_ok=True)                                                                                                                # Ensure the directory structure exists, creating it if necessary
    out_final = Path(args.out_topk_final_csv)                                                                                                                         # Construct the file path object for the top-3 CSV export
    out_final.parent.mkdir(parents=True, exist_ok=True)                                                                                                               # Ensure the directory structure exists, creating it if necessary
    topk_df.to_csv(out_topk, index=False)                                                                                                                             # Dump the full extended panel to disk
    final_df.to_csv(out_final, index=False)                                                                                                                           # Dump the elite validation panel to disk

    summary = {                                                                                                                                                       # Initiate the assembly of the metadata tracking dictionary
        "stage": "10c",                                                                                                                                               # Log the current pipeline stage identifier
        "input_search_csv": str(args.search_csv),                                                                                                                     # Log the absolute path of the originating search file
        "top_k_rows": int(len(topk_df)),                                                                                                                              # Log the exact number of candidates successfully placed in the extended panel
        "final_rows": int(len(final_df)),                                                                                                                             # Log the exact number of candidates successfully placed in the elite panel
        "best_stage10_score": float(topk_df.iloc[0]["stage10_composite_score"]) if len(topk_df) else float("nan"),                                                    # Log the absolute highest fitness score, or NaN if unavailable
        "best_target_probability": float(topk_df.iloc[0]["target_probability"]) if len(topk_df) else float("nan"),                                                    # Log the absolute highest host probability, or NaN if unavailable
        "best_if1_log_likelihood": float(topk_df.iloc[0]["if1_log_likelihood"]) if len(topk_df) else float("nan"),                                                    # Log the absolute highest thermodynamic stability score, or NaN if unavailable
        "out_topk_csv": str(out_topk),                                                                                                                                # Log the filepath where the top-10 panel was saved
        "out_topk_final_csv": str(out_final),                                                                                                                         # Log the filepath where the top-3 panel was saved
    }                                                                                                                                                                 # Close the dictionary instantiation
    write_json(summary, args.out_json)                                                                                                                                # Serialize the summary tracking data and write it directly to disk
    print(f"Wrote: {out_topk}")                                                                                                                                       # Print standard execution confirmation for the top-10 CSV
    print(f"Wrote: {out_final}")                                                                                                                                      # Print standard execution confirmation for the top-3 CSV
    print(f"Wrote: {args.out_json}")                                                                                                                                  # Print standard execution confirmation for the metadata JSON


if __name__ == "__main__":                                                                                                                                            # Standard Python check to isolate script execution from modular imports
    main()                                                                                                                                                            # Call the primary runtime logic