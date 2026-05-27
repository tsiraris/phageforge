#!/usr/bin/env python
"""Stage 11d: Validate the Stage 11 panel with the unmodified Stage 08a validator.

This script is a thin orchestration wrapper around the already-corrected Stage 08 validator.
It exists so Stage 11 stays operationally consistent with the earlier project stages while
keeping the final heavy structural check delegated to the established validation script.

The subprocess isolation is important: it reclaims VRAM between the embedding
backbone (Stage 11b) and ESMFold (08a), preventing OOM crashes on SageMaker boxes.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from phageforge.stage11_utils import (
    EXIT_INFERENCE_ERROR,
    EXIT_INPUT_ERROR,
    EXIT_OK,
    write_json,
)


def parse_args() -> argparse.Namespace:
    """
    Parse CLI arguments for the Stage 11 validation wrapper.
    """
    ap = argparse.ArgumentParser(description="Run the fixed Stage 08a structural validator on Stage 11 candidate panels.")  # CLI description explaining wrapper purpose
    ap.add_argument("--validated_csv", type=str, required=True, help="Top-K Stage 11 candidate CSV produced by 11c.")       # Define the path to the compact elite sequence panel
    ap.add_argument("--ranked_csv", type=str, required=True, help="Full Stage 11 search CSV.")                              # Define the path to the broader search log for metadata lookup
    ap.add_argument("--context_json", type=str, required=True, help="Stage 11 context JSON.")                               # Define the blueprint providing structural anchor metadata
    ap.add_argument("--validator_script", type=str, default="scripts/08a_structural_fasttrack_validation.py", help="Path to the unmodified Stage 08a validator script.")               # Set the default path to the 08a structural oracle script
    ap.add_argument("--out_dir", type=str, required=True, help="Directory for structural artifacts.")                       # Define the target directory for PDB and summary outputs
    ap.add_argument("--out_json", type=str, required=True, help="Path for the compact wrapper-launch summary JSON.")        # Define the path to save the execution manifest JSON
    ap.add_argument("--top_k", type=int, default=3, help="Number of candidates to validate.")                               # Set the number of sequences to escalate to structural folding
    ap.add_argument("--device", type=str, default="cuda", help="Device passed to the validator.")                           # Assign the GPU target for the structural oracle
    ap.add_argument("--chunk_size", type=int, default=128, help="ESMFold chunk size.")                                      # Set the memory footprint chunking limit
    ap.add_argument("--num_recycles", type=int, default=1, help="ESMFold recycles.")                                        # Set the spatial coordinate refinement iterations
    ap.add_argument("--resume", action="store_true", help="Reuse PDBs.")                                                    # Flag to enable structural caching reuse behavior
    return ap.parse_args()                                                                                                  # Return the parsed terminal arguments object


def main() -> None:
    args = parse_args()                                                                                                    # Process user-provided arguments into a local namespace
    # Resolve the absolute filesystem path for the structural validator and the input data panel
    validator = Path(args.validator_script).resolve()                                                                      # Resolve the absolute filesystem path for the structural validator
    if not validator.exists():                                                                                             # Check if the requested script exists on the current host
        print(f"[ERROR] Missing validator script: {validator}", file=sys.stderr)                                           # Report pathing errors to standard error output
        sys.exit(EXIT_INPUT_ERROR)                                                                                         # Abort the process with an input error exit code
    validated_csv = Path(args.validated_csv).resolve()                                                                     # Convert string path to validated CSV into absolute Path object
    if not validated_csv.exists():                                                                                         # Validate existence of the input data panel
        print(f"[ERROR] Missing validated CSV: {validated_csv}", file=sys.stderr)                                          # Log input errors directly to the system console
        sys.exit(EXIT_INPUT_ERROR)                                                                                         # Signal an input-related failure termination
    ranked_csv = Path(args.ranked_csv).resolve()                                                                           # Resolve absolute path for the ranked history dataframe
    if not ranked_csv.exists():                                                                                            # Verify that the search reference history exists
        print(f"[ERROR] Missing ranked CSV: {ranked_csv}", file=sys.stderr)                                                # Log missing search data pathing error
        sys.exit(EXIT_INPUT_ERROR)                                                                                         # Abort with standard input failure exit code
    context_json = Path(args.context_json).resolve()                                                                       # Resolve absolute path for the structural context file
    if not context_json.exists():                                                                                          # Verify the existence of the biological anchor JSON
        print(f"[ERROR] Missing context JSON: {context_json}", file=sys.stderr)                                            # Log pathing error preventing structural anchor recovery
        sys.exit(EXIT_INPUT_ERROR)                                                                                         # Abort with standard input failure exit code

    out_dir = Path(args.out_dir).resolve()                                                                                 # Finalize the target root directory for 3D coordinate storage
    out_dir.mkdir(parents=True, exist_ok=True)                                                                             # Instantiate the output directory tree on the filesystem

    # Initiate the subprocess command array
    cmd = [                                                                                                                # Initiate the subprocess command array construction
        sys.executable,                                                                                                    # Identify the current running Python executable
        str(validator),                                                                                                    # Pass the validator script as the primary argument
        "--validated_csv", str(validated_csv),                                                                             # Provide the filtered sequence panel input
        "--ranked_csv", str(ranked_csv),                                                                                   # Provide the broader reference sequence table
        "--context_json", str(context_json),                                                                               # Provide the structural context anchor
        "--out_dir", str(out_dir),                                                                                         # Provide the target storage root directory
        "--top_k", str(int(args.top_k)),                                                                                   # Provide the structural validation depth limit
        "--device", str(args.device),                                                                                      # Provide hardware accelerator instructions
        "--chunk_size", str(int(args.chunk_size)),                                                                         # Provide memory protection chunking parameters
        "--num_recycles", str(int(args.num_recycles)),                                                                     # Provide spatial coordinate refinement cycles
    ]                                                                                                                      # Terminate the primary command argument list
    # If the structural caching mode is enabled, indicate that to the validator
    if args.resume:                                                                                                        # Assess if structural caching mode is enabled by operator
        cmd.append("--resume")                                                                                             # Append the caching flag to the system execution list

    # Monitor the structural validation process during execution
    print(f"[INFO] Launching Stage 08a validator: {' '.join(cmd)}", flush=True)                                            # Output the full system command for transparent logging
    start_ts = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")                             # Timestamp the start time for archival purposes

    try:                                                                                                                   # Wrap subprocess invocation protecting the pipeline from external runtime failures
        result = subprocess.run(cmd, check=False)                                                                          # Dispatch the structural validator process to an isolated hardware environment
    except FileNotFoundError as exc:                                                                                       # Handle runtime failures regarding missing binary files
        print(f"[ERROR] Could not launch validator: {exc}", file=sys.stderr)                                               # Surface the error message through standard streams
        sys.exit(EXIT_INFERENCE_ERROR)                                                                                     # Abort with the specific inference error exit status

    end_ts = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")                               # Timestamp the conclusion of the structural evaluation process

    # If the structural validator completed successfully, log the run metadata into a JSON metadata file                                                                                      # Assess if the structural validator completed successfully
    summary = {                                                                                                            # Begin documenting the execution metadata locally
        "stage": "11d",                                                                                                    # Define the stage identification metadata
        "validated_csv": str(validated_csv),                                                                               # Pin the provenance of the input data
        "ranked_csv": str(ranked_csv),                                                                                     # Pin the provenance of the reference data
        "context_json": str(context_json),                                                                                 # Pin the provenance of the structural context
        "validator_script": str(validator),                                                                                # Pin the path of the evaluation oracle
        "out_dir": str(out_dir),                                                                                           # Pin the location of generated 3D atomic files
        "top_k": int(args.top_k),                                                                                          # Log the validation depth constraint used
        "device": str(args.device),                                                                                        # Log the hardware platform utilized
        "chunk_size": int(args.chunk_size),                                                                                # Log the memory optimization setting
        "num_recycles": int(args.num_recycles),                                                                            # Log the geometry refinement count
        "resume": bool(args.resume),                                                                                       # Log the cache usage state
        "validator_exit_code": int(result.returncode),                                                                     # Store the oracle's process completion code
        "start_timestamp_utc": start_ts,                                                                                   # Archive the chronological start event
        "end_timestamp_utc": end_ts,                                                                                       # Archive the chronological end event
        "expected_outputs": {                                                                                              # Map expected standard output artifacts for clarity
            "summary_csv": str(out_dir / "stage08_structural_fasttrack_summary.csv"),                                      # Identify the path to the physical validation table
            "summary_json": str(out_dir / "stage08_structural_fasttrack_summary.json"),                                    # Identify the path to the JSON validation results
            "report_md": str(out_dir / "stage08_structural_fasttrack_report.md"),                                          # Identify the path to the physical validation markdown report
            "pdb_dir": str(out_dir / "pdbs"),                                                                              # Identify the path to the physical coordinate file directory
        },                                                                                                                 # Terminate output mapping
    }                                                                                                                      # Close telemetry record dictionary
    write_json(summary, args.out_json)                                                                                     # Flush the launch manifest record to the persistent disk storage
    print(f"[OK] Wrote: {args.out_json}", flush=True)                                                                      # Confirm successful archival of the process summary

    if result.returncode != 0:                                                                                             # Assess if the structural oracle reported failure
        print(f"[ERROR] Stage 08a validator exited with code {result.returncode}.", file=sys.stderr)                       # Notify console of oracle execution failure
        sys.exit(EXIT_INFERENCE_ERROR)                                                                                     # Terminate process identifying inference error status
    sys.exit(EXIT_OK)                                                                                                      # Signal standard execution success following validator pass


if __name__ == "__main__":                                                                                                 # Guard ensuring execution happens only on direct call
    main()                                                                                                                 # Run the validation orchestration logic