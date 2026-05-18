#!/usr/bin/env python
"""Stage 10e: Build the final Stage 10 report and comparison summary.

This script packages the key Stage 10 artifacts into one compact scientific report.
It emphasizes the exact shift made in this phase:
- Stage 09 still relied on sequence-first local search with structure proxies,
- Stage 10 moved the structure signal directly into candidate generation with inverse folding.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from phageforge.stage10_utils import read_json, write_json


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 10 reporting step."""
    ap = argparse.ArgumentParser(description="Build the final Stage 10 report.")                                                                             # Initialize the CLI argument parser with a description of the reporting script
    ap.add_argument("--stage10_context_json", type=str, required=True, help="Stage 10 context JSON produced by 10a_prepare_stage10_structure_context.py.") # Define the required argument for the Stage 10 structural blueprint JSON path
    ap.add_argument("--search_csv", type=str, required=True, help="Search CSV produced by 10b_run_inverse_folding_beam_search.py.")                        # Define the required argument for the full inverse-folding search log CSV path
    ap.add_argument("--prefilter_csv", type=str, required=True, help="Top-k Stage 10 prefilter CSV produced by 10c_prefilter_stage10_candidates.py.")      # Define the required argument for the prefiltered top candidate CSV path
    ap.add_argument("--validation_csv", type=str, default=None, help="Optional structural validation summary CSV produced by the Stage 08 validator.")     # Define the optional argument for the final ESMFold structural validation CSV path
    ap.add_argument("--baseline_validation_csv", type=str, default=None, help="Optional earlier Stage 08 or Stage 09 validation CSV used for comparison.") # Define the optional argument for historical baseline data to enable direct A/B comparison
    ap.add_argument("--out_dir", type=str, required=True, help="Directory where the Stage 10 markdown report and summary JSON will be written.")           # Define the required argument specifying where all report artifacts will be saved
    return ap.parse_args()                                                                                                                                 # Parse the provided terminal arguments and return them as a Namespace object


def summarize_validation(validation_df: pd.DataFrame) -> dict:
    """Summarize the decisive structural validation outcomes into one compact dictionary."""
    if validation_df.empty:                                                                                                                                # Check if the validation dataframe is completely empty (e.g., validation was skipped)
        return {                                                                                                                                           # Return a safely zeroed-out dictionary to prevent metric calculation crashes
            "rows": 0,                                                                                                                                     # Indicate that zero candidates were validated
            "pass_count": 0,                                                                                                                               # Indicate that zero candidates passed the structural gauntlet
            "best_mean_plddt": None,                                                                                                                       # Provide a null placeholder since no thermodynamic confidence data exists
            "best_mutation_site_mean_plddt": None,                                                                                                         # Provide a null placeholder since no local mutation confidence data exists
            "best_rmsd": None,                                                                                                                             # Provide a null placeholder since no geometric drift data exists
        }                                                                                                                                                  # Close the fallback dictionary
    pass_col = "stage08_pass" if "stage08_pass" in validation_df.columns else None                                                                         # Safely identify the boolean pass/fail column, returning None if it is missing
    return {                                                                                                                                               # Return a populated dictionary containing peak structural discovery metrics
        "rows": int(len(validation_df)),                                                                                                                   # Count and record the exact number of rows present in the validation dataframe
        "pass_count": int(validation_df[pass_col].astype(bool).sum()) if pass_col else 0,                                                                  # Calculate the total number of perfectly passing candidates, defaulting to 0 if the column is missing
        "best_mean_plddt": float(validation_df["esmfold_mean_plddt"].max()) if "esmfold_mean_plddt" in validation_df.columns else None,                    # Extract the absolute highest global structural confidence score achieved
        "best_mutation_site_mean_plddt": float(validation_df["mutation_site_mean_plddt"].max()) if "mutation_site_mean_plddt" in validation_df.columns else None, # Extract the highest local mutation confidence score achieved
        "best_rmsd": float(validation_df["rmsd_to_selected_seed"].min()) if "rmsd_to_selected_seed" in validation_df.columns else None,                    # Extract the absolute lowest (best) geometric drift RMSD value achieved
    }                                                                                                                                                      # Close the metric aggregation dictionary


def main() -> None:
    # Read the Stage 10 context and the central candidate tables generated during structure-conditioned redesign.
    args = parse_args()                                                                                                                                    # Execute the argument parser to retrieve user-defined paths and settings
    out_dir = Path(args.out_dir)                                                                                                                           # Convert the raw output directory string into a robust Python Path object
    out_dir.mkdir(parents=True, exist_ok=True)                                                                                                             # Safely create the destination directory tree if it does not already exist

    context = read_json(args.stage10_context_json)                                                                                                         # Load and parse the master Stage 10 constraint and context JSON file into a dictionary
    search_df = pd.read_csv(args.search_csv)                                                                                                               # Read the comprehensive inverse-folding search history into a pandas DataFrame
    prefilter_df = pd.read_csv(args.prefilter_csv)                                                                                                         # Read the diversity-prefiltered candidate list into a pandas DataFrame
    validation_df = pd.read_csv(args.validation_csv) if args.validation_csv else pd.DataFrame()                                                            # Read the structural validation CSV if provided, otherwise initialize an empty DataFrame
    baseline_df = pd.read_csv(args.baseline_validation_csv) if args.baseline_validation_csv else pd.DataFrame()                                            # Read the historical baseline CSV if provided, otherwise initialize an empty DataFrame

    # Summarize the main Stage 10 redesign table so the report records the best structure-conditioned candidate statistics clearly.
    search_summary = {                                                                                                                                     # Initialize a dictionary to aggregate peak performance metrics from the generative search phase
        "rows": int(len(search_df)),                                                                                                                       # Record the total number of hypotheses evaluated during the beam search
        "best_stage10_score": float(search_df["stage10_composite_score"].max()) if len(search_df) else None,                                               # Extract the absolute highest composite fitness score achieved during search
        "best_target_probability": float(search_df["target_probability"].max()) if len(search_df) else None,                                               # Extract the highest target-host infectivity probability achieved
        "best_if1_log_likelihood": float(search_df["if1_log_likelihood"].max()) if len(search_df) else None,                                               # Extract the highest ESM-IF1 3D structural compatibility score achieved
        "best_seed_identity": float(search_df["seed_identity"].max()) if len(search_df) else None,                                                         # Extract the highest sequence identity retained among the generated candidates
    }                                                                                                                                                      # Close the search summary dictionary
    prefilter_summary = {                                                                                                                                  # Initialize a dictionary to summarize the outcomes of the prefilter stage
        "rows": int(len(prefilter_df)),                                                                                                                    # Record the exact number of candidates that survived the diversity prefilter
        "top_sample_ids": prefilter_df["sample_id"].tolist() if "sample_id" in prefilter_df.columns else [],                                               # Extract a clean list of the specific integer IDs that were forwarded to validation
    }                                                                                                                                                      # Close the prefilter summary dictionary
    validation_summary = summarize_validation(validation_df)                                                                                               # Call the helper function to calculate peak structural metrics for the Stage 10 panel
    baseline_summary = summarize_validation(baseline_df)                                                                                                   # Call the helper function to calculate peak structural metrics for the historical baseline

    # Build the compact JSON artifact so the notebook and any downstream manuscript/reporting step can consume one canonical summary.
    summary = {                                                                                                                                            # Compile the ultimate master summary dictionary representing the entire Stage 10 lifecycle
        "stage": "10e",                                                                                                                                    # Tag the metadata payload with the specific pipeline stage identifier
        "target_host": str(context["target_host"]),                                                                                                        # Record the specific bacteria the AI was trying to target
        "seed_protein_id": str(context["selected_seed"].get("seed_protein_id", "")),                                                                       # Record the specific wild-type protein ID used as the chassis
        "seed_pdb_path": str(context.get("seed_pdb_path", "")),                                                                                            # Record the absolute path to the physical 3D scaffold file used for inverse folding
        "search_summary": search_summary,                                                                                                                  # Embed the previously calculated search metric dictionary
        "prefilter_summary": prefilter_summary,                                                                                                            # Embed the previously calculated prefilter metric dictionary
        "validation_summary": validation_summary,                                                                                                          # Embed the previously calculated structural validation metric dictionary
        "baseline_summary": baseline_summary,                                                                                                              # Embed the previously calculated historical baseline comparison dictionary
    }                                                                                                                                                      # Close the master summary dictionary
    write_json(summary, out_dir / "stage10_report_summary.json")                                                                                           # Serialize the master summary dictionary to disk as a cleanly formatted JSON file

    # Write the human-readable markdown report that explains why Stage 10 was introduced and what it achieved relative to earlier phases.
    lines = []                                                                                                                                             # Initialize an empty list to systematically build the Markdown report lines
    lines.append("# Stage 10 structure-conditioned redesign report")                                                                                       # Append the primary document title header
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    lines.append("## Why Stage 10 exists")                                                                                                                 # Append a secondary header explaining the project's architectural shift
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    lines.append("Stage 08 and Stage 09 showed that sequence-first redesign, even when tightened with structure-aware proxies, still failed decisive full structural validation. Stage 10 therefore moved the structure signal upstream and redesigned candidates directly against a fixed seed scaffold with an inverse-folding objective.") # Append the core scientific rationale for the Inverse-Folding pivot
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    lines.append("## Core redesign setup")                                                                                                                 # Append a secondary header detailing the generative constraints
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    lines.append(f"- target host: **{context['target_host']}**")                                                                                           # Inject the target bacteria strain into a formatted Markdown bullet point
    lines.append(f"- selected seed id: **{context['selected_seed'].get('seed_protein_id', '')}**")                                                         # Inject the base scaffold ID into a formatted Markdown bullet point
    lines.append(f"- seed scaffold: **{context.get('seed_pdb_path', '')}**")                                                                               # Inject the physical path of the 3D anchor into a formatted Markdown bullet point
    lines.append(f"- hard editable positions: **{context['editable_region'].get('hard_positions', [])}**")                                                 # Inject the primary permitted mutation indices into a formatted Markdown bullet point
    lines.append(f"- soft editable positions: **{context['editable_region'].get('soft_positions', [])}**")                                                 # Inject the secondary fallback mutation indices into a formatted Markdown bullet point
    lines.append(f"- mutation budget: **{context['editable_region'].get('min_mutations')} to {context['editable_region'].get('max_mutations')}**")         # Inject the strict upper and lower mutation limits into a bullet point
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    lines.append("## Search summary")                                                                                                                      # Append a secondary header for generative search results
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    lines.append(f"- scored candidates: **{search_summary['rows']}**")                                                                                     # Inject the total evaluated sequence count into a formatted Markdown bullet point
    lines.append(f"- best Stage 10 composite score: **{search_summary['best_stage10_score']}**")                                                           # Inject the peak mathematical fitness score into a formatted Markdown bullet point
    lines.append(f"- best target probability: **{search_summary['best_target_probability']}**")                                                            # Inject the peak infectivity prediction into a formatted Markdown bullet point
    lines.append(f"- best inverse-folding log-likelihood: **{search_summary['best_if1_log_likelihood']}**")                                                # Inject the peak 3D compatibility score into a formatted Markdown bullet point
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    lines.append("## Prefilter summary")                                                                                                                   # Append a secondary header for diversity prefiltering results
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    lines.append(f"- prefilter panel size: **{prefilter_summary['rows']}**")                                                                               # Inject the count of diverse survivors into a formatted Markdown bullet point
    lines.append(f"- sample ids carried forward: **{prefilter_summary['top_sample_ids']}**")                                                               # Inject the exact identifiers of the forwarded candidates into a bullet point
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    lines.append("## Structural validation summary")                                                                                                       # Append a secondary header for the ESMFold oracle results
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    if validation_df.empty:                                                                                                                                # Branch logic: check if the structural validation was skipped or aborted
        lines.append("Structural validation CSV was not provided to the report builder. The Stage 10 report therefore ends at the prefilter stage.")       # Append a graceful fallback message explaining the absence of structural data
    else:                                                                                                                                                  # Branch logic: execute if valid structural data exists
        lines.append(f"- validated rows: **{validation_summary['rows']}**")                                                                                # Inject the count of physically modeled candidates into a bullet point
        lines.append(f"- pass count: **{validation_summary['pass_count']}**")                                                                              # Inject the exact number of physically viable, passing candidates into a bullet point
        lines.append(f"- best mean pLDDT: **{validation_summary['best_mean_plddt']}**")                                                                    # Inject the peak global folding confidence achieved into a bullet point
        lines.append(f"- best mutation-site mean pLDDT: **{validation_summary['best_mutation_site_mean_plddt']}**")                                        # Inject the peak local mutation confidence achieved into a bullet point
        lines.append(f"- best RMSD to seed: **{validation_summary['best_rmsd']} Å**")                                                                      # Inject the lowest (best) geometric deviation from the seed into a bullet point
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    lines.append("## Comparison to earlier structural validation")                                                                                         # Append a secondary header to establish the scientific A/B test
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    if baseline_df.empty:                                                                                                                                  # Branch logic: check if historical baseline data was provided
        lines.append("No baseline validation CSV was provided for comparison.")                                                                            # Append a graceful fallback message indicating no comparative baseline exists
    else:                                                                                                                                                  # Branch logic: execute if historical baseline data is present
        lines.append(f"- baseline validated rows: **{baseline_summary['rows']}**")                                                                         # Inject the baseline candidate count into a comparative bullet point
        lines.append(f"- baseline pass count: **{baseline_summary['pass_count']}**")                                                                       # Inject the baseline success count into a comparative bullet point
        lines.append(f"- baseline best mean pLDDT: **{baseline_summary['best_mean_plddt']}**")                                                             # Inject the baseline peak global confidence into a comparative bullet point
        lines.append(f"- baseline best mutation-site mean pLDDT: **{baseline_summary['best_mutation_site_mean_plddt']}**")                                 # Inject the baseline peak local confidence into a comparative bullet point
        lines.append(f"- baseline best RMSD to seed: **{baseline_summary['best_rmsd']} Å**")                                                               # Inject the baseline lowest drift into a comparative bullet point
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    lines.append("## Interpretation")                                                                                                                      # Append the final secondary header for scientific conclusions
    lines.append("")                                                                                                                                       # Append a blank line for Markdown visual spacing
    lines.append("Stage 10 is the first phase in which the scaffold itself becomes part of candidate generation rather than only downstream evaluation. This makes it the correct methodological successor to Stage 08 and Stage 09, regardless of whether the final heavy structural validation fully succeeds on the first attempt.") # Append the definitive project doctrine securing Inverse-Folding as the correct scientific path
    (out_dir / "stage10_report.md").write_text("\n".join(lines), encoding="utf-8")                                                                         # Compile all list items into a massive string with newlines and save it to disk as a Markdown file
    print(f"Wrote: {out_dir / 'stage10_report_summary.json'}")                                                                                             # Output a confirmation to the terminal that the JSON summary was saved successfully
    print(f"Wrote: {out_dir / 'stage10_report.md'}")                                                                                                       # Output a confirmation to the terminal that the human-readable Markdown report was saved successfully


if __name__ == "__main__":                                                                                                                                 # Check if the python script is executed directly from the command line
    main()                                                                                                                                                 # Invoke the primary orchestration function