#!/usr/bin/env python
"""Stage 10d: Validate the final Stage 10 candidate panel with the fixed structural validator.

This script is a thin orchestration wrapper around the already-corrected Stage 08 validator.
It exists so Stage 10 stays operationally consistent with the earlier project stages while
keeping the final heavy structural check delegated to the established validation script.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from phageforge.stage10_utils import write_json


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 10 validation wrapper."""
    ap = argparse.ArgumentParser(description="Run the fixed structural validator on Stage 10 candidate panels.")                                  # Initialize the argument parser with a description
    ap.add_argument("--validated_csv", type=str, required=True, help="Top-k or top-3 Stage 10 candidate CSV produced by 10c_prefilter...")        # Add required argument for prefiltered input CSV path
    ap.add_argument("--ranked_csv", type=str, required=True, help="The broader Stage 10 search CSV used as the ranked reference table...")        # Add required argument for broader search CSV path
    ap.add_argument("--context_json", type=str, required=True, help="Stage 10 context JSON so the validator can recover the selected seed...")    # Add required argument for context JSON path
    ap.add_argument("--validator_script", type=str, default="scripts/08a_structural_fasttrack_validation.py", help="Path to the corrected...")    # Add optional argument for the validator script path
    ap.add_argument("--out_dir", type=str, required=True, help="Output directory where the structural validation artifacts will be written.")     # Add required argument for the output directory
    ap.add_argument("--top_k", type=int, default=3, help="Number of candidates to validate.")                                                     # Add optional argument to define top-k candidates limit
    ap.add_argument("--device", type=str, default="cuda", help="Device used by the structural validator.")                                        # Add optional argument for compute device mapping
    ap.add_argument("--chunk_size", type=int, default=128, help="Chunk size passed through to the structural validator.")                         # Add optional argument for ESMFold processing chunk size
    ap.add_argument("--num_recycles", type=int, default=1, help="Number of recycles passed through to the structural validator.")                 # Add optional argument for ESMFold recycle iterations
    ap.add_argument("--resume", action="store_true", help="Reuse existing PDB files when the validator supports resume mode.")                    # Add boolean flag to skip re-computing existing PDBs
    ap.add_argument("--out_json", type=str, required=True, help="Where to write the compact validation-launch summary JSON.")                     # Add required argument for the output summary JSON
    return ap.parse_args()                                                                                                                        # Parse the provided command-line arguments and return them


def main() -> None:
    # Build the exact subprocess command that reuses the already-corrected heavy structural validator from Stage 08.
    args = parse_args()                                                                                                                           # Retrieve the parsed command-line arguments
    validator = Path(args.validator_script)                                                                                                       # Convert the validator script string path to a Path object
    if not validator.exists():                                                                                                                    # Check if the validation script actually exists on disk
        raise FileNotFoundError(f"Missing validator script: {validator}")                                                                         # Abort execution if the validation script is not found

    # Start defining the list of strings for the subprocess command
    cmd = [                                                                                                                                       # Start defining the list of strings for the subprocess command
        sys.executable,                                                                                                                           # Add the path to the current Python executable
        str(validator),                                                                                                                           # Add the path to the external script to execute
        "--validated_csv",                                                                                                                        # Pass the argument flag for validated CSV
        str(args.validated_csv),                                                                                                                  # Pass the actual path of the validated CSV file
        "--ranked_csv",                                                                                                                           # Pass the argument flag for ranked CSV
        str(args.ranked_csv),                                                                                                                     # Pass the actual path of the ranked CSV file
        "--context_json",                                                                                                                         # Pass the argument flag for context JSON
        str(args.context_json),                                                                                                                   # Pass the actual path of the context JSON file
        "--out_dir",                                                                                                                              # Pass the argument flag for output directory
        str(args.out_dir),                                                                                                                        # Pass the actual path of the output directory
        "--top_k",                                                                                                                                # Pass the argument flag for top_k limit
        str(args.top_k),                                                                                                                          # Pass the user-specified top_k value
        "--device",                                                                                                                               # Pass the argument flag for compute device
        str(args.device),                                                                                                                         # Pass the user-specified device selection
        "--chunk_size",                                                                                                                           # Pass the argument flag for chunk size
        str(args.chunk_size),                                                                                                                     # Pass the user-specified chunk size
        "--num_recycles",                                                                                                                         # Pass the argument flag for number of recycles
        str(args.num_recycles),                                                                                                                   # Pass the user-specified recycle count
    ]                                                                                                                                             # Close the list of command arguments
    if args.resume:                                                                                                                               # Check if the user enabled the resume flag
        cmd.append("--resume")                                                                                                                    # Append the resume argument flag to the command list

    # Launch the validator directly so the terminal and notebook logs remain identical to the established Stage 08 execution flow.
    subprocess.run(cmd, check=True)                                                                                                               # Execute the external command, raising an exception if it fails

    # Write a tiny summary JSON so downstream reporting steps can record exactly how the heavy validation was invoked.
    summary = {                                                                                                                                   # Initialize a dictionary to store validation launch metadata
        "stage": "10d",                                                                                                                           # Tag the artifact with its specific pipeline stage
        "validated_csv": str(args.validated_csv),                                                                                                 # Record the validated CSV path used
        "ranked_csv": str(args.ranked_csv),                                                                                                       # Record the ranked CSV path used
        "context_json": str(args.context_json),                                                                                                   # Record the context JSON path used
        "validator_script": str(validator),                                                                                                       # Record the path of the validator script executed
        "out_dir": str(args.out_dir),                                                                                                             # Record the output directory path
        "top_k": int(args.top_k),                                                                                                                 # Record the number of top candidates requested
        "device": str(args.device),                                                                                                               # Record the compute device used
        "chunk_size": int(args.chunk_size),                                                                                                       # Record the chunk size parameter used
        "num_recycles": int(args.num_recycles),                                                                                                   # Record the number of recycles used
    }                                                                                                                                             # Close the summary dictionary
    write_json(summary, args.out_json)                                                                                                            # Serialize and write the summary dictionary to a JSON file
    print(f"Wrote: {args.out_json}")                                                                                                              # Print a confirmation message to standard output


if __name__ == "__main__":                                                                                                                        # Ensure this block executes only if the script is run directly
    main()                                                                                                                                        # Invoke the main execution function