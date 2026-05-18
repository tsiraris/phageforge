#!/usr/bin/env python
"""Stage 09f: Validate the top Stage 09 candidates with the fixed Stage 08 structural validator.

This wrapper adapts the Stage 09 prefilter outputs into the validator-compatible file layout and
then launches the already-fixed 08a_structural_fasttrack_validation.py script. It keeps Stage 09
scientifically coherent by reusing the exact same structural decision layer that falsified Stage 07.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd



def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 09 structural-validation wrapper."""
    ap = argparse.ArgumentParser(description="Validate Stage 09 candidates with the fixed Stage 08 fast-track validator.")                      # Initialize the argument parser with a description
    ap.add_argument("--prefilter_csv", type=str, required=True, help="Prefiltered Stage 09 candidate CSV produced by 09e_structural_prefilter.py.") # Add required argument for prefiltered input CSV path
    ap.add_argument("--context_json", type=str, required=True, help="Stage 07 context JSON needed by 08a_structural_fasttrack_validation.py.")  # Add required argument for context JSON path
    ap.add_argument("--out_dir", type=str, required=True, help="Where to write the Stage 09 structural-validation outputs.")                    # Add required argument for the output directory
    ap.add_argument("--top_k", type=int, default=10, help="How many prefitered Stage 09 candidates to validate structurally.")                  # Add optional argument to define top-k candidates limit
    ap.add_argument("--device", type=str, default="auto", help="Device argument forwarded to 08a_structural_fasttrack_validation.py.")          # Add optional argument for compute device mapping
    ap.add_argument("--chunk_size", type=int, default=128, help="Chunk size forwarded to the ESMFold validator.")                               # Add optional argument for ESMFold processing chunk size
    ap.add_argument("--num_recycles", type=int, default=1, help="Number of ESMFold recycles forwarded to the validator.")                       # Add optional argument for ESMFold recycle iterations
    ap.add_argument("--resume", action="store_true", help="Reuse existing Stage 09 seed/candidate PDBs when present.")                          # Add boolean flag to skip re-computing existing PDBs
    return ap.parse_args()                                                                                                                      # Parse the provided command-line arguments and return them


def main() -> None:
    # Read the Stage 09 prefiltered candidate table and keep only the requested top-k panel for structural validation.
    args = parse_args()                                                                # Retrieve the parsed command-line arguments
    out_dir = Path(args.out_dir)                                                       # Convert the string output directory path to a Path object
    out_dir.mkdir(parents=True, exist_ok=True)                                         # Create the output directory, ignoring if it already exists
    panel_df = pd.read_csv(args.prefilter_csv).head(args.top_k).copy()                 # Read CSV, slice the top-k rows, and copy to a new DataFrame
    if panel_df.empty:                                                                 # Check if the filtered DataFrame is completely empty
        raise ValueError("The Stage 09 prefilter CSV is empty, so there is nothing to validate.") # Abort execution with an error if no data exists

    # Write validator-compatible ranked and validated CSV files so the fixed Stage 08 script can be reused unchanged.
    ranked_csv = out_dir / "stage09_ranked_for_validation.csv"                         # Define the file path for the ranked candidates CSV
    validated_csv = out_dir / f"stage09_validated_top{len(panel_df)}.csv"              # Define the file path for the validated candidates CSV
    panel_df.to_csv(ranked_csv, index=False)                                           # Save the panel DataFrame to the ranked CSV without row indices
    panel_df.to_csv(validated_csv, index=False)                                        # Save the same DataFrame to the validated CSV as a starting point

    # Launch the fixed Stage 08 structural validator from the same repository so Stage 09 uses the same falsification layer.
    stage08_script = Path(__file__).resolve().parent / "08a_structural_fasttrack_validation.py" # Resolve the absolute path to the Stage 08 validator script
    cmd = [                                                                                     # Start defining the list of strings for the subprocess command
        sys.executable,                                                                         # Add the path to the current Python executable
        str(stage08_script),                                                                    # Add the path to the external script to execute
        "--validated_csv",                                                                      # Pass the argument flag for validated CSV
        str(validated_csv),                                                                     # Pass the actual path of the validated CSV file
        "--ranked_csv",                                                                         # Pass the argument flag for ranked CSV
        str(ranked_csv),                                                                        # Pass the actual path of the ranked CSV file
        "--context_json",                                                                       # Pass the argument flag for context JSON
        str(args.context_json),                                                                 # Pass the actual path of the context JSON file
        "--out_dir",                                                                            # Pass the argument flag for output directory
        str(out_dir),                                                                           # Pass the actual path of the output directory
        "--top_k",                                                                              # Pass the argument flag for top_k limit
        str(len(panel_df)),                                                                     # Pass the dynamic length of the dataframe as top_k
        "--device",                                                                             # Pass the argument flag for compute device
        str(args.device),                                                                       # Pass the user-specified device selection
        "--chunk_size",                                                                         # Pass the argument flag for chunk size
        str(args.chunk_size),                                                                   # Pass the user-specified chunk size
        "--num_recycles",                                                                       # Pass the argument flag for number of recycles
        str(args.num_recycles),                                                                 # Pass the user-specified recycle count
    ]                                                                                           # Close the list of command arguments
    if args.resume:                                                                             # Check if the user enabled the resume flag
        cmd.append("--resume")                                                                  # Append the resume argument flag to the command list

    print("Running:")                                                                           # Print a notification that script execution is starting
    print(" ".join(cmd))                                                                        # Print the fully assembled command as a readable string
    subprocess.run(cmd, check=True)                                                             # Execute the external command, raising an exception if it fails


if __name__ == "__main__":      # Ensure this block executes only if the script is run directly, not imported
    main()                      # Invoke the main execution function