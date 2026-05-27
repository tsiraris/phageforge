#!/usr/bin/env python
"""Stage 11c: Prefilter Stage 11 search candidates before expensive validation.

This script is intentionally conservative. It enforces the Stage 11 mutation
budget, deduplicates exact sequence matches, and applies two passes of
embedding-space diversity reranking — first over the full search pool to pick
a top-K extended panel, then over that panel to pick the compact top-K-final
validation panel.

The two-pass pattern is the same one used in Stage 10c. Re-embedding the
top-K panel before the second pass guarantees the second-pass cosine
similarities are calibrated against the panel's actual variance rather than
the full pool's, which prevents the second pass from being driven by far
outliers in the larger search log.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from phageforge.stage11_utils import (
    EXIT_INFERENCE_ERROR,
    EXIT_INPUT_ERROR,
    EXIT_OK,
    embed_sequences,
    greedy_diverse_subset,
    read_json,
    seed_everything,
    write_json,
)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Prefilter and diversify Stage 11 redesign candidates.")                                               # Initialize the command-line argument parser to process terminal inputs
    ap.add_argument("--stage11_context_json", type=str, required=True, help="Stage 11 context JSON (used to recover the mutation budget).")         # Define the required context JSON file argument specifying mutational limits
    ap.add_argument("--search_csv", type=str, required=True, help="Full Stage 11 search CSV produced by 11b.")                                      # Define the required search CSV file argument containing the raw candidates
    ap.add_argument("--out_topk_csv", type=str, required=True, help="Where to write the top-K prefilter CSV.")                                      # Define the required output argument for the extended Top-K prefiltered CSV
    ap.add_argument("--out_topk_final_csv", type=str, required=True, help="Where to write the compact final validation CSV.")                       # Define the required output argument for the elite Top-K final validation CSV
    ap.add_argument("--out_json", type=str, required=True, help="Where to write the prefilter summary JSON.")                                       # Define the required output argument for the JSON tracking metadata
    ap.add_argument("--top_k", type=int, default=10, help="Number of candidates in the extended prefilter panel.")                                  # Define the configurable integer setting the extended panel target capacity
    ap.add_argument("--top_k_final", type=int, default=3, help="Number of candidates in the compact validation panel.")                             # Define the configurable integer setting the elite panel target capacity
    ap.add_argument("--embedding_model", type=str, default="facebook/esm2_t33_650M_UR50D", help="Embedding backbone used for diversity reranking.") # Define the specific HuggingFace transformer model ID used for vector embedding
    ap.add_argument("--diversity_penalty_weight", type=float, default=0.20, help="Weighting for the diversity penalty term used inside the diversity reranking.")                          # Define the configurable weight for the diversity penalty term
    ap.add_argument("--batch_size", type=int, default=4, help="Batch size for candidate embedding.")                                                # Define the batch size limit to protect hardware from VRAM exhaustion
    ap.add_argument("--seed", type=int, default=42, help="Random seed for deterministic ordering.")                                                 # Define the random seed to ensure exact reproducibility across multiple runs
    return ap.parse_args()                                                                                                                          # Parse the supplied arguments from the terminal and return the namespace


def main() -> None:
    args = parse_args()                                                                                                                             # Call the parser to map terminal inputs to accessible Python variables
    seed_everything(args.seed)                                                                                                                      # Force all pseudo-random number generators to follow a deterministic sequence
    # Load the Stage 11 contextual blueprint, the search log, and the mutation budget
    context_path = Path(args.stage11_context_json)                                                                                                  # Cast the context JSON filepath string into a robust Path object
    search_path = Path(args.search_csv)                                                                                                             # Cast the search CSV filepath string into a robust Path object
    if not context_path.exists():                                                                                                                   # Dynamically evaluate if the required context configuration file actually exists
        print(f"[ERROR] Missing context JSON: {context_path}", file=sys.stderr)                                                                     # Transmit a clear missing file error message to the standard error stream
        sys.exit(EXIT_INPUT_ERROR)                                                                                                                  # Terminate the application immediately with a standardized input error code
    if not search_path.exists():                                                                                                                    # Dynamically evaluate if the required search generation log actually exists
        print(f"[ERROR] Missing search CSV: {search_path}", file=sys.stderr)                                                                        # Transmit a clear missing file error message to the standard error stream
        sys.exit(EXIT_INPUT_ERROR)                                                                                                                  # Terminate the application immediately with a standardized input error code

    context = read_json(context_path)                                                                                                               # Deserialize the master Stage 11 contextual blueprint into a Python dictionary
    search_df = pd.read_csv(search_path)                                                                                                            # Load the comprehensive search candidate matrix into a pandas DataFrame

    if search_df.empty:                                                                                                                             # Check if the loaded search dataframe contains absolutely zero generated sequences
        print("[ERROR] Search CSV is empty — Stage 11b produced no candidates.", file=sys.stderr)                                                   # Emit a fatal error indicating the upstream generation engine completely failed
        sys.exit(EXIT_INPUT_ERROR)                                                                                                                  # Terminate the application immediately with a standardized input error code

    editable_region = dict(context.get("editable_region", {}))                                                                                      # Extract the mutational constraints dictionary from the master context blueprint
    min_mut = int(editable_region.get("min_mutations", 1))                                                                                          # Safely extract the absolute minimum required substitutions per candidate
    max_mut = int(editable_region.get("max_mutations", 4))                                                                                          # Safely extract the absolute maximum allowable substitutions per candidate
    
    # Apply the mutation budget filter to the search log, and verify that at least one candidate remains
    keep = search_df.copy()                                                                                                                         # Create an explicit memory copy of the dataframe to prevent SettingWithCopy warnings
    if "mutation_count" in keep.columns:                                                                                                            # Verify the existence of the mutation count telemetry column before filtering
        before = len(keep)                                                                                                                          # Record the total candidate volume prior to applying the budget constraints
        keep = keep.loc[(keep["mutation_count"] >= min_mut) & (keep["mutation_count"] <= max_mut)].copy()                                           # Definitively slice the dataframe retaining only sequences within the exact edit limits
        print(f"[INFO] Mutation budget filter [{min_mut}, {max_mut}] retained {len(keep)}/{before} candidates.", flush=True)                        # Log the exact attrition rate resulting from the mutational bounds check
    if keep.empty:                                                                                                                                  # Monitor system health to intercept complete candidate annihilation post-filtering
        print("[ERROR] Stage 11 prefilter retained zero candidates after the mutation-budget filter.", file=sys.stderr)                             # Emit a fatal error explaining that no candidates survived the strict budget rules
        sys.exit(EXIT_INPUT_ERROR)                                                                                                                  # Terminate the application immediately with a standardized input error code

    # Sort surviving candidates by the composite score, apply deduplication to the search log, and verify that at least one candidate remains
    sort_cols = [c for c in ["stage10_composite_score", "stage11_composite_score", "target_probability", "if1_log_likelihood"] if c in keep.columns]# Assemble an ordered list of available optimization metrics prioritizing the composite score
    if not sort_cols:                                                                                                                               # Check if the dataframe is completely devoid of recognized mathematical ranking columns
        print("[ERROR] Search CSV has no composite/target/IF1 score columns — cannot rank.", file=sys.stderr)                                       # Emit a fatal error indicating structural sorting is mathematically impossible
        sys.exit(EXIT_INPUT_ERROR)                                                                                                                  # Terminate the application immediately with a standardized input error code
    keep = keep.sort_values(sort_cols, ascending=False)                                                                                             # Hierarchically structure the survivors strictly optimizing for aggregate theoretical fitness
    keep = keep.drop_duplicates(subset=["candidate_sequence"], keep="first").reset_index(drop=True)                                                 # Completely eradicate repetitive structural sequence clones preserving only the best instance
    print(f"[INFO] After exact-sequence dedup: {len(keep)} unique candidates.", flush=True)                                                         # Log the surviving distinct sequence count confirming successful deduplication

    if len(keep) < int(args.top_k_final):                                                                                                           # Assess if the surviving pool is smaller than the mandated Oracle submission limit
        print(f"[WARN] Only {len(keep)} unique candidates available, less than --top_k_final={int(args.top_k_final)}. Continuing.", flush=True)     # Emit a non-fatal operational warning informing the user of the reduced panel size

    # Embed the remaining candidates into a high-dimensional vector space
    score_col = sort_cols[0]                                                                                                                        # Isolate the highest priority optimization metric available in the active dataframe
    print(f"[INFO] Embedding {len(keep)} sequences with '{args.embedding_model}' for first-pass diversity rerank…", flush=True)                     # Announce the initiation of the computationally intensive spatial embedding projection
    try:                                                                                                                                            # Wrap the heavy tensor operations in a protective error handling block
        embeddings = embed_sequences(keep["candidate_sequence"].astype(str).tolist(), model_name=args.embedding_model, batch_size=args.batch_size)  # Execute parallel transformation generating high-dimensional vector representations
    except Exception as exc:                                                                                                                        # Catch any hardware, memory, or library errors emitted during the embedding process
        print(f"[ERROR] First-pass embedding failed: {exc}", file=sys.stderr)                                                                       # Surface the specific transformer exception directly to the standard error channel
        sys.exit(EXIT_INFERENCE_ERROR)                                                                                                              # Terminate the application utilizing the specific inference failure exit code

    # Select the top-k candidates that are diverse enough from each other
    topk_idx = greedy_diverse_subset(embeddings, keep[score_col].to_numpy(dtype=np.float32), top_k=int(args.top_k), diversity_penalty_weight=0.20)                                 # Assess Euclidean penalties discarding parallel evolutionary paths to prevent mode collapse
    topk_df = keep.iloc[topk_idx].sort_values(sort_cols, ascending=False).reset_index(drop=True)                                                    # Siphon selected indices rendering the extended highly distinct secondary tracking list

    # Re-embeds the Top-k panel in isolation
    print(f"[INFO] Re-embedding the top-{len(topk_df)} panel for second-pass diversity rerank…", flush=True)                                        # Announce the initiation of the extreme-precision secondary geometric differentiation
    try:                                                                                                                                            # Wrap the secondary heavy tensor operations in a protective error handling block
        topk_embeddings = embed_sequences(topk_df["candidate_sequence"].astype(str).tolist(), model_name=args.embedding_model, batch_size=args.batch_size) # Re-project strictly the top-K panel to calibrate cosine variance exactly to the subset
    except Exception as exc:                                                                                                                        # Catch any hardware, memory, or library errors emitted during the secondary embedding
        print(f"[ERROR] Second-pass embedding failed: {exc}", file=sys.stderr)                                                                      # Surface the specific transformer exception directly to the standard error channel
        sys.exit(EXIT_INFERENCE_ERROR)                                                                                                              # Terminate the application utilizing the specific inference failure exit code

    # Runs the greedy_diverse_subset a second time to select the absolute elite panel (default: Top 3).
    final_idx = greedy_diverse_subset(topk_embeddings, topk_df[score_col].to_numpy(dtype=np.float32), top_k=int(args.top_k_final), diversity_penalty_weight=0.20)                  # Run exact secondary geometric differentiation extracting the absolute pinnacle panel
    final_df = topk_df.iloc[final_idx].sort_values(sort_cols, ascending=False).reset_index(drop=True)                                               # Score-sort the final maximum optimization panel ready for pure Oracle submission

    # Injects a 1-indexed sample_id column into the finalized Top 10 and Top 3 DataFrames, for compatibility with the unmodified Stage 08a structural validator.
    topk_df = topk_df.copy()                                                                                                                        # Instatiate explicit memory copy of the extended panel to prevent pandas slice warnings
    final_df = final_df.copy()                                                                                                                      # Instatiate explicit memory copy of the elite panel to prevent pandas slice warnings
    topk_df["sample_id"] = np.arange(1, len(topk_df) + 1)                                                                                           # Inject clean 1-indexed integers ensuring seamless legacy validation execution downstream
    final_df["sample_id"] = np.arange(1, len(final_df) + 1)                                                                                         # Inject identical logic specifically addressing the 3D Oracle submission panel
    
    # Write both the extended (1st diversity) and elite (2nd diversity) panels to compatible CSV files
    out_topk = Path(args.out_topk_csv)                                                                                                              # Resolve the target output path for the extended diversity candidate list
    out_final = Path(args.out_topk_final_csv)                                                                                                       # Resolve the target output path for the elite Oracle candidate list
    out_topk.parent.mkdir(parents=True, exist_ok=True)                                                                                              # Explicitly instantiate local directory hierarchies guaranteeing file write capabilities
    out_final.parent.mkdir(parents=True, exist_ok=True)                                                                                             # Explicitly instantiate local directory hierarchies guaranteeing file write capabilities
    topk_df.to_csv(out_topk, index=False)                                                                                                           # Emit the complete extended artifact to designated physical drive storage
    final_df.to_csv(out_final, index=False)                                                                                                         # Emit the explicit 3D Oracle target artifact physically protecting internal execution
    print(f"[OK] Wrote: {out_topk}", flush=True)                                                                                                    # Send basic user operational notification confirming extended CSV completion
    print(f"[OK] Wrote: {out_final}", flush=True)                                                                                                   # Send basic user operational notification confirming elite CSV completion

    # Generate a summary of the selected panels and write it to a JSON file
    summary = {                                                                                                                                     # Generate summary logic array packaging vital statistical metadata parameters
        "stage": "11c",                                                                                                                             # Anchor identity tagging the universal execution phase identifier explicitly
        "input_search_csv": str(search_path.resolve()),                                                                                             # Lock the absolute physical pathway resolving the raw generation history origin
        "stage11_context_json": str(context_path.resolve()),                                                                                        # Lock the absolute physical pathway resolving the overarching structural blueprint
        "top_k_rows": int(len(topk_df)),                                                                                                            # Store precise integer reflecting expanded selection pool total volume
        "final_rows": int(len(final_df)),                                                                                                           # Store exact constrained list sent identically to ESMFold fasttrack processing
        "score_column_used": score_col,                                                                                                             # Document the specific mathematical column utilized to arbitrate the sorting hierarchy
        "diversity_penalty_weight": 0.20,                                                                                                           # Record the immutable geometric penalty weight enforcing distinct structural theories
        "best_stage10_score": float(topk_df.iloc[0][score_col]) if len(topk_df) else float("nan"),                                                  # Record optimum normalized composite performance protecting against null value failures
        "best_target_probability": float(topk_df.iloc[0]["target_probability"]) if (len(topk_df) and "target_probability" in topk_df.columns) else float("nan"), # Record the highest target infectivity score preserving data integrity safely
        "best_if1_log_likelihood": float(topk_df.iloc[0]["if1_log_likelihood"]) if (len(topk_df) and "if1_log_likelihood" in topk_df.columns) else float("nan"), # Record the highest physical log-likelihood preserving data integrity safely
        "final_sample_ids": [int(x) for x in final_df["sample_id"].tolist()],                                                                       # Extract a clean list of the specific integer IDs that were forwarded to validation
        "out_topk_csv": str(out_topk.resolve()),                                                                                                    # Document the finalized operational location of the extended prefilter output
        "out_topk_final_csv": str(out_final.resolve()),                                                                                             # Document the finalized operational location of the elite prefilter output
    }
    write_json(summary, args.out_json)                                                                                                              # Perform explicit physical output command passing string JSON directly to memory
    print(f"[OK] Wrote: {args.out_json}", flush=True)                                                                                               # Send basic user operational notification confirming metadata JSON completion
    sys.exit(EXIT_OK)                                                                                                                               # Declare flawless architectural exit explicitly cleanly


if __name__ == "__main__":                                                                                                                          # Prevent unintended script initiation by assessing underlying execution triggers
    main()                                                                                                                                          # Execute root sequence triggering entirely completely defined logic chains