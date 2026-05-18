#!/usr/bin/env python
"""Stage 09g: Summarize the Stage 09 structure-aware redesign run.

This script compares the Stage 09 prefilter and structural-validation outputs against the Stage 08
baseline so the user can tell whether the redesign strategy is genuinely improving structural success.
"""

from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd



def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 09 reporting script."""
    ap = argparse.ArgumentParser(description="Build the Stage 09 final comparison report.")                                                                                   # Initialize the argument parser with a description
    ap.add_argument("--search_csv", type=str, required=True, help="Stage 09 full search CSV produced by 09d_localized_search.py.")                                          # Add the required search_csv argument
    ap.add_argument("--prefilter_csv", type=str, required=True, help="Stage 09 prefitered candidate CSV produced by 09e_structural_prefilter.py.")                          # Add the required prefilter_csv argument
    ap.add_argument("--validation_csv", type=str, required=True, help="Stage 09 structural-validation summary CSV produced by 09f_validate_stage09_candidates.py / 08a_structural_fasttrack_validation.py.")  # Add the required validation_csv argument
    ap.add_argument("--baseline_stage08_csv", type=str, default=None, help="Optional baseline Stage 08 structural summary CSV for direct comparison.")                      # Add the optional baseline_stage08_csv argument
    ap.add_argument("--out_dir", type=str, required=True, help="Directory where the Stage 09 summary artifacts will be written.")                                           # Add the required out_dir argument
    return ap.parse_args()                                                                                                                                                    # Parse the provided command-line arguments and return them



def safe_mean(series: pd.Series) -> float:
    """
    Calculates the arithmetic mean of a pandas Series safely, correctly handling 
    empty, completely null, or non-numeric datasets without raising exceptions.
    This might happen for example if a protein fails to fold entirely, with its pLDDT or RMSD columns might 
    be empty, or upstream processes might inject text like "error" or "timeout").
    
    Example:
        Standard case (valid floats)
        safe_mean(pd.Series([85.0, 90.0, 95.0]))
        90.0
        
        Corrupted data (strings and None are coerced to NaN and ignored)
        safe_mean(pd.Series([70.0, "folding_error", 80.0, None]))
        75.0
        
        Failure case (no valid numbers exist)
        safe_mean(pd.Series(["timeout", None, "OOM_error"]))
        nan
    """
    values = pd.to_numeric(series, errors="coerce")  # Convert series items to numeric types, forcing unparseable items to NaN
    if values.dropna().empty:                        # Check if dropping all NaNs leaves the sequence completely empty
        return math.nan                              # Return mathematical NaN if no numeric data is available
    return float(values.mean())                      # Calculate the mean of the numeric sequence and return as a float



def main() -> None:
    # Read the Stage 09 search, prefilter, and structural-validation tables that together define the redesigned pipeline outcome.
    args = parse_args()                                                                                         # Call the parse_args function to get the CLI arguments
    search_df = pd.read_csv(args.search_csv)                                                                    # Load the full search CSV file into a pandas DataFrame
    prefilter_df = pd.read_csv(args.prefilter_csv)                                                              # Load the prefiltered candidate CSV file into a pandas DataFrame
    validation_df = pd.read_csv(args.validation_csv)                                                            # Load the validation summary CSV file into a pandas DataFrame
    baseline_df = pd.read_csv(args.baseline_stage08_csv) if args.baseline_stage08_csv else None                 # Load the baseline CSV if an argument was provided, otherwise set to None

    # Compute the main Stage 09 summary metrics, with a strong emphasis on top-panel structural success rather than search-only scores.
    summary = {                                                                                                                                           # Initialize a dictionary to store Stage 09 summary metrics
        "n_search_candidates": int(len(search_df)),                                                                                                       # Record the total number of rows in the search dataframe
        "n_prefilter_candidates": int(len(prefilter_df)),                                                                                                 # Record the total number of rows in the prefilter dataframe
        "n_validated_candidates": int(len(validation_df)),                                                                                                # Record the total number of rows in the validation dataframe
        "stage09_pass_count": int(validation_df["stage08_pass"].sum()) if "stage08_pass" in validation_df.columns else 0,                                 # Sum up the boolean "pass" values if the column exists
        "stage09_pass_rate": float(validation_df["stage08_pass"].mean()) if "stage08_pass" in validation_df.columns and len(validation_df) else 0.0,      # Calculate the proportion of passing candidates
        "stage09_mean_plddt": safe_mean(validation_df.get("esmfold_mean_plddt", pd.Series(dtype=float))),                                                 # Extract and calculate the safe mean for overall ESMFold pLDDT
        "stage09_mean_mutation_site_plddt": safe_mean(validation_df.get("mutation_site_mean_plddt", pd.Series(dtype=float))),                             # Extract and calculate the safe mean for mutation site pLDDT
        "stage09_mean_rmsd": safe_mean(validation_df.get("rmsd_to_selected_seed", pd.Series(dtype=float))),                                               # Extract and calculate the safe mean for RMSD to the seed scaffold
        "stage09_mean_target_probability": safe_mean(prefilter_df.get("target_probability", pd.Series(dtype=float))),                                     # Extract and calculate the safe mean for the targeted probability
    }                                                                                                                                                     # End summary dictionary initialization

    # Optionally compare Stage 09 against the Stage 08 baseline so the redesign step can claim a real structural improvement.
    baseline_summary = None                                                                                                                         # Initialize baseline summary dictionary reference to None
    if baseline_df is not None and not baseline_df.empty:                                                                                           # Check to ensure baseline dataframe exists and is not empty
        baseline_summary = {                                                                                                                        # Initialize the baseline dictionary for direct comparison
            "baseline_n_candidates": int(len(baseline_df)),                                                                                         # Record total number of rows from the baseline dataframe
            "baseline_pass_count": int(baseline_df["stage08_pass"].sum()) if "stage08_pass" in baseline_df.columns else 0,                          # Sum the boolean "pass" values for baseline candidates
            "baseline_pass_rate": float(baseline_df["stage08_pass"].mean()) if "stage08_pass" in baseline_df.columns else 0.0,                      # Calculate baseline pass rate
            "baseline_mean_plddt": safe_mean(baseline_df.get("esmfold_mean_plddt", pd.Series(dtype=float))),                                        # Extract and calculate the safe mean for baseline ESMFold pLDDT
            "baseline_mean_mutation_site_plddt": safe_mean(baseline_df.get("mutation_site_mean_plddt", pd.Series(dtype=float))),                    # Extract and calculate the safe mean for baseline mutation site pLDDT
            "baseline_mean_rmsd": safe_mean(baseline_df.get("rmsd_to_selected_seed", pd.Series(dtype=float))),                                      # Extract and calculate the safe mean for baseline RMSD to seed scaffold
        }                                                                                                                                           # End baseline dictionary initialization

    # Write the final candidate table, a compact JSON summary, and a markdown report for project Stage 09 closeout.
    out_dir = Path(args.out_dir)                                                                                                                                                   # Create a Path object for the target output directory
    out_dir.mkdir(parents=True, exist_ok=True)                                                                                                                                     # Make the directory structure if it doesn't already exist
    validation_df.sort_values(["stage08_pass", "esmfold_mean_plddt", "rmsd_to_selected_seed"], ascending=[False, False, True]).to_csv(out_dir / "stage09_final_candidate_table.csv", index=False) # Sort candidates by pass, descending pLDDT, and ascending RMSD, then export to CSV

    summary_path = out_dir / "stage09_summary.json"                                                                                     # Define the destination path for the JSON summary output
    summary_path.write_text(json.dumps({"summary": summary, "baseline": baseline_summary}, indent=2))                                   # Dump both summary dicts into formatted JSON string and write to file

    lines = [                                                                                                                                                                                                                                           # Start a list of markdown lines to build the text report
        "# Stage 09 structure-aware redesign report",                                                                                                                                                                                                   # Add the main header title
        "",                                                                                                                                                                                                                                             # Add an empty string for paragraph spacing
        "## Stage 09 summary",                                                                                                                                                                                                                          # Add a subheader for the stage 09 section
        "",                                                                                                                                                                                                                                             # Add an empty string for paragraph spacing
        f"- search candidates generated: **{summary['n_search_candidates']}**",                                                                                                                                                                         # Inject the search candidates count into a bullet point
        f"- prefilter survivors: **{summary['n_prefilter_candidates']}**",                                                                                                                                                                              # Inject the prefilter survivor count into a bullet point
        f"- structurally validated candidates: **{summary['n_validated_candidates']}**",                                                                                                                                                                # Inject the validated candidates count into a bullet point
        f"- structural pass count: **{summary['stage09_pass_count']}**",                                                                                                                                                                                # Inject the total structural pass count into a bullet point
        f"- structural pass rate: **{summary['stage09_pass_rate']:.3f}**",                                                                                                                                                                              # Inject the structural pass rate into a bullet point formatting it to 3 decimals
        f"- mean ESMFold pLDDT: **{summary['stage09_mean_plddt']:.3f}**" if not math.isnan(summary['stage09_mean_plddt']) else "- mean ESMFold pLDDT: **nan**",                                                                                         # Inject the mean ESMFold pLDDT gracefully falling back to nan text
        f"- mean mutation-site pLDDT: **{summary['stage09_mean_mutation_site_plddt']:.3f}**" if not math.isnan(summary['stage09_mean_mutation_site_plddt']) else "- mean mutation-site pLDDT: **nan**",                                                 # Inject the mean mutation-site pLDDT gracefully falling back to nan text
        f"- mean RMSD to seed: **{summary['stage09_mean_rmsd']:.3f} Å**" if not math.isnan(summary['stage09_mean_rmsd']) else "- mean RMSD to seed: **nan**",                                                                                           # Inject the mean RMSD gracefully falling back to nan text
        "",                                                                                                                                                                                                                                             # Add an empty string for paragraph spacing
        "## Top validated candidates",                                                                                                                                                                                                                  # Add a subheader for the candidates table
        "",                                                                                                                                                                                                                                             # Add an empty string for paragraph spacing
        "| rank | sample_id | stage09_score | mean_pLDDT | mut_mean_pLDDT | RMSD_to_seed | pass | reason |",                                                                                                                                            # Define the markdown table headers
        "|---:|---:|---:|---:|---:|---:|:---:|---|",                                                                                                                                                                                                    # Define the markdown table text alignments
    ]                                                                                                                                                                                                                                                   # End line list initialization
    # Iterate over the top 10 best-scoring dataframe rows, giving each a rank starting at 1, prioritizing stage08_pass and dynamically pull dataframe attributes into markdown format
    for rank, (_, row) in enumerate(validation_df.sort_values(["stage08_pass", "esmfold_mean_plddt", "rmsd_to_selected_seed"], ascending=[False, False, True]).head(10).iterrows(), start=1):                                                                                                                                                       # Iterate over the top 10 best-scoring dataframe rows giving each a rank starting at 1
        lines.append(                                                                                                                                                                                                                                                                                                                               # Append a newly formatted markdown table row to the lines list
            f"| {rank} | {int(row['sample_id'])} | {float(row.get('final_multimodal_rank_score', row.get('stage09_score', 0.0))):.6f} | {float(row.get('esmfold_mean_plddt', float('nan'))):.2f} | {float(row.get('mutation_site_mean_plddt', float('nan'))):.2f} | {float(row.get('rmsd_to_selected_seed', float('nan'))):.3f} | {bool(row.get('stage08_pass', False))} | {row.get('stage08_decision_reason', '')} |"  # Dynamically pull dataframe attributes into markdown format
        )                                                                                                                                                                                                                                                                                                                                           # End the append operation

    if baseline_summary is not None:                                                                                                                                                                                                                              # Check if the baseline stats successfully computed
        lines.extend(                                                                                                                                                                                                                                             # Add multiple strings at once to the end of the markdown list
            [                                                                                                                                                                                                                                                     # Start a list block of baseline-related strings
                "",                                                                                                                                                                                                                                               # Add an empty string for paragraph spacing
                "## Comparison to failed Stage 08 baseline",                                                                                                                                                                                                      # Add a subheader for the baseline data comparison section
                "",                                                                                                                                                                                                                                               # Add an empty string for paragraph spacing
                f"- baseline validated candidates: **{baseline_summary['baseline_n_candidates']}**",                                                                                                                                                              # Inject the baseline validated candidates count into a bullet point
                f"- baseline pass rate: **{baseline_summary['baseline_pass_rate']:.3f}**",                                                                                                                                                                        # Inject the baseline pass rate into a bullet point
                f"- baseline mean ESMFold pLDDT: **{baseline_summary['baseline_mean_plddt']:.3f}**" if not math.isnan(baseline_summary['baseline_mean_plddt']) else "- baseline mean ESMFold pLDDT: **nan**",                                                     # Inject baseline mean pLDDT into a bullet point falling back to nan
                f"- baseline mean mutation-site pLDDT: **{baseline_summary['baseline_mean_mutation_site_plddt']:.3f}**" if not math.isnan(baseline_summary['baseline_mean_mutation_site_plddt']) else "- baseline mean mutation-site pLDDT: **nan**",             # Inject baseline mutation pLDDT into a bullet point falling back to nan
                f"- baseline mean RMSD to seed: **{baseline_summary['baseline_mean_rmsd']:.3f} Å**" if not math.isnan(baseline_summary['baseline_mean_rmsd']) else "- baseline mean RMSD to seed: **nan**",                                                       # Inject baseline mean RMSD into a bullet point falling back to nan
            ]                                                                                                                                                                                                                                                     # End the baseline block list
        )                                                                                                                                                                                                                                                         # Complete the lines.extend operation

    (out_dir / "stage09_report.md").write_text("\n".join(lines))                                 # Concatenate the final markdown lines using newline characters and write to file

    print(f"Wrote: {out_dir / 'stage09_final_candidate_table.csv'}")                             # Print terminal confirmation for the candidates table
    print(f"Wrote: {summary_path}")                                                              # Print terminal confirmation for the JSON summary path
    print(f"Wrote: {out_dir / 'stage09_report.md'}")                                             # Print terminal confirmation for the text report path


if __name__ == "__main__":                                                                       # Check if the python script is executed as the main runner
    main()                                                                                       # Call the entrypoint logic