#!/usr/bin/env python
"""
Stage 08b Closeout Packager.

This script serves as the final data aggregation and packaging step (Stage 08b) 
in the phage receptor-binding protein (RBP) retargeting pipeline. 

It takes the sequence data, initial multimodal rankings (Stage 07), and structural 
validation metrics (Stage 08) across multiple top-candidate CSV files, and combines 
them into a single consolidated table. It then deduplicates the data, applies a 
final sorting based on viability and structural integrity (pLDDT, RMSD), and assigns 
a final tier-based decision label to each candidate. Finally, it exports this data 
into a comprehensive handoff bundle containing a CSV table, a FASTA file of the 
protein sequences, a JSON summary, and a Markdown case study report.
"""

from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()                                      # Initialize the command-line argument parser object
    ap.add_argument("--ranked_csv", type=Path, required=True)           # Define the required argument for the ranked CSV path
    ap.add_argument("--validated_top3_csv", type=Path, required=True)   # Define the required argument for the validated top 3 CSV
    ap.add_argument("--validated_top5_csv", type=Path, required=True)   # Define the required argument for the validated top 5 CSV
    ap.add_argument("--stage08_top3_csv", type=Path, required=True)     # Define the required argument for the Stage 08 top 3 CSV
    ap.add_argument("--stage08_top5_csv", type=Path, required=True)     # Define the required argument for the Stage 08 top 5 CSV
    ap.add_argument("--context_json", type=Path, required=True)         # Define the required argument for the context JSON path
    ap.add_argument("--out_dir", type=Path, required=True)              # Define the required argument for the output directory path
    return ap.parse_args()                                              # Parse and return the provided command-line arguments


def write_fasta(df: pd.DataFrame, fasta_path: Path) -> None:
    lines: list[str] = []                                               # Initialize an empty list to store the formatted FASTA lines
    for _, row in df.iterrows():                                        # Iterate over each row (candidate) in the provided DataFrame
        header = (                                                      # Begin constructing the FASTA header string
            f">sample{int(row['sample_id'])}|rank={int(row['final_rank'])}|" # Embed the sample ID and the final computed rank
            f"pass={bool(row['stage08_pass'])}|score={float(row['final_multimodal_rank_score']):.6f}" # Embed the pass status and score
        )                                                               # Finish constructing the FASTA header string
        lines.append(header)                                            # Append the finalized header string to the lines list
        lines.append(str(row['candidate_sequence']))                    # Append the actual candidate protein sequence string
    fasta_path.write_text("\n".join(lines) + "\n")                      # Join all lines with newline characters and write to file


def decision_label(pass_flag: bool, rank: int) -> str:
    if pass_flag and rank <= 2:                                         # Check if the candidate passed validation and is in the top 2
        return "primary"                                                # Return 'primary' label for these top-tier candidates
    if pass_flag:                                                       # Check if the candidate passed validation but ranked > 2
        return "backup"                                                 # Return 'backup' label for these secondary candidates
    return "near_pass_or_fail"                                          # Return 'near_pass_or_fail' for candidates that failed validation


def main() -> None:
    args = parse_args()                                                                 # Call parse_args to retrieve parsed command-line arguments
    args.out_dir.mkdir(parents=True, exist_ok=True)                                     # Create the output directory, ignoring if it already exists

    ranked = pd.read_csv(args.ranked_csv)                                               # Load the primarily ranked candidates CSV into a DataFrame
    valid3 = pd.read_csv(args.validated_top3_csv)                                       # Load the validated top 3 candidates CSV into a DataFrame
    valid5 = pd.read_csv(args.validated_top5_csv)                                       # Load the validated top 5 candidates CSV into a DataFrame
    s08_3 = pd.read_csv(args.stage08_top3_csv)                                          # Load the Stage 08 top 3 structural metrics into a DataFrame
    s08_5 = pd.read_csv(args.stage08_top5_csv)                                          # Load the Stage 08 top 5 structural metrics into a DataFrame
    context = json.loads(args.context_json.read_text())                                 # Read and decode the project context variables from JSON

    stage08 = pd.concat([s08_3, s08_5], ignore_index=True)                              # Concatenate the Stage 08 top 3 and top 5 structural dataframes
    stage08 = stage08.sort_values(["sample_id", "stage08_structural_rank"]).drop_duplicates("sample_id", keep="first") # Sort by rank and drop duplicate IDs

    valid = pd.concat([valid3, valid5], ignore_index=True)                              # Concatenate the validated top 3 and top 5 dataframes
    valid = valid.sort_values(["sample_id", "validated_rank"]).drop_duplicates("sample_id", keep="first") # Sort by validated rank and drop duplicate IDs

    keep_cols = [                                                                       # Define a list of essential columns to retain from ranked df
        "sample_id", "generation_regime", "candidate_sequence", "mutation_positions",   # Include sequence identifiers and mutation tracking columns
        "final_multimodal_rank_score", "target_score", "strict_manifold_score", "structure_score", # Include various scores calculated in earlier stages
        "rank_raw", "rank_diverse"                                                      # Include the raw and diversity-based ranking positions
    ]                                                                                   # Finish defining the list of essential columns
    final_df = ranked[keep_cols].merge(valid[["sample_id", "validated_rank"]], on="sample_id", how="inner") # Merge ranked data with validation ranks
    final_df = final_df.merge(                                                          # Begin merging the resulting dataframe with structural metrics
        stage08[[                                                                       # Select only the specific columns needed from the stage08 dataframe
            "sample_id", "stage08_structural_rank", "esmfold_mean_plddt", "esmfold_ptm",# Include structural ranks and global ESMFold confidence scores
            "mutation_site_mean_plddt", "mutation_site_confidence_ge70_fraction",       # Include localized confidence scores around the mutated sites
            "rmsd_to_selected_seed", "stage08_pass", "stage08_decision_reason", "candidate_pdb" # Include RMSD drift, pass flags, and structural paths
        ]],                                                                             # Finish selecting the specific columns from stage08
        on="sample_id",                                                                 # Specify that the merge should align on the 'sample_id' key
        how="left",                                                                     # Perform a left join so no rows from final_df are dropped
    )                                                                                   # Finish merging the structural metrics into final_df
    final_df = final_df.sort_values(                                                    # Begin sorting the fully assembled final dataframe
        ["stage08_pass", "final_multimodal_rank_score", "esmfold_mean_plddt", "rmsd_to_selected_seed"], # Sort hierarchically by pass, score, pLDDT, RMSD
        ascending=[False, False, False, True],                                          # Set sorting orders: descending for best scores, ascending for lowest RMSD
    ).reset_index(drop=True)                                                            # Apply the sorting and reset the dataframe index purely sequentially
    final_df["final_rank"] = np.arange(1, len(final_df) + 1)                            # Generate and assign sequential integers as the absolute final rank
    final_df["final_decision"] = [decision_label(bool(p), int(r)) for p, r in zip(final_df["stage08_pass"], final_df["final_rank"])] # Compute and map decision labels

    csv_path = args.out_dir / "final_candidate_table.csv"                               # Define the full output file path for the final CSV table
    fasta_path = args.out_dir / "final_top_candidates.fasta"                            # Define the full output file path for the final FASTA file
    md_path = args.out_dir / "final_case_study.md"                                      # Define the full output file path for the Markdown case study
    json_path = args.out_dir / "final_closeout_summary.json"                            # Define the full output file path for the JSON summary

    final_df.to_csv(csv_path, index=False)                                              # Export the assembled final dataframe to the defined CSV path
    write_fasta(final_df, fasta_path)                                                   # Invoke the write_fasta helper to generate the FASTA file

    selected_seed = context["selected_seed"]                                            # Extract the selected seed dictionary from the context JSON
    pass_count = int(final_df["stage08_pass"].fillna(False).sum())                      # Calculate the total number of candidates that passed Stage 08
    primary = final_df[final_df["final_decision"] == "primary"]                         # Filter the dataframe to isolate only 'primary' categorized candidates

    summary = {                                                                         # Begin constructing the summary dictionary for JSON export
        "selected_seed": {                                                              # Create a nested dictionary to store properties of the selected seed
            "seed_protein_id": selected_seed.get("seed_protein_id"),                    # Retrieve the selected seed's protein ID
            "virus_accession": selected_seed.get("virus_accession"),                    # Retrieve the selected seed's originating virus accession
            "source_host": selected_seed.get("source_host"),                            # Retrieve the selected seed's natural host bacteria
            "sequence_length": selected_seed.get("sequence_length"),                    # Retrieve the total sequence length of the selected seed
        },                                                                              # End the nested selected seed dictionary
        "num_final_candidates": int(len(final_df)),                                     # Record the total amount of final candidates processed
        "num_stage08_pass": pass_count,                                                 # Record the final count of candidates that successfully passed
        "primary_sample_ids": [int(x) for x in primary["sample_id"].tolist()],          # Extract and record a clean list of 'primary' sample IDs
        "final_table": str(csv_path),                                                   # Record the string file path for the generated CSV
        "final_fasta": str(fasta_path),                                                 # Record the string file path for the generated FASTA
    }                                                                                   # End constructing the summary dictionary
    json_path.write_text(json.dumps(summary, indent=2))                                 # Serialize the dictionary to a formatted JSON string and write to file

    md = [                                                                              # Begin compiling a list of strings representing Markdown document lines
        "# Final closeout case study",                                                  # Insert the primary Header 1 for the document
        "",                                                                             # Insert a blank line to maintain proper Markdown formatting
        "## Project framing",                                                           # Insert a Header 2 for the project overview section
        "A validity-aware, scaffold-constrained phage RBP retargeting workflow was closed out with local ESM3 generation, Stage 07 multimodal reranking, and Stage 08 structural fast-track validation.", # Insert the project overview paragraph text
        "",                                                                             # Insert a blank line to maintain proper Markdown formatting
        "## Selected seed",                                                             # Insert a Header 2 for the selected seed details
        f"- protein_id: **{selected_seed.get('seed_protein_id', 'unknown')}**",         # Format and append a bullet point for the seed protein ID
        f"- source host: **{selected_seed.get('source_host', 'unknown')}**",            # Format and append a bullet point for the seed source host
        f"- virus accession: **{selected_seed.get('virus_accession', 'unknown')}**",    # Format and append a bullet point for the seed virus accession
        f"- sequence length: **{selected_seed.get('sequence_length', 'unknown')}**",    # Format and append a bullet point for the seed sequence length
        "",                                                                             # Insert a blank line to maintain proper Markdown formatting
        "## Structural closeout summary",                                               # Insert a Header 2 for the summary metrics section
        f"- final candidate count: **{len(final_df)}**",                                # Format and append a bullet point for total candidates
        f"- Stage 08 passes: **{pass_count}**",                                         # Format and append a bullet point for passing candidates
        f"- primary candidates: **{', '.join(str(int(x)) for x in primary['sample_id']) if len(primary) else 'none'}**", # Format and append a bullet point for primary IDs
        "",                                                                             # Insert a blank line to maintain proper Markdown formatting
        "## Final ranking",                                                             # Insert a Header 2 for the final tabular ranking section
        "",                                                                             # Insert a blank line before starting the Markdown table
        "| final_rank | sample_id | decision | pass | final_score | mean_pLDDT | mut_mean_pLDDT | RMSD_to_seed | reason |", # Insert the column headers for the Markdown table
        "|---:|---:|---|:---:|---:|---:|---:|---:|---|",                                # Insert the column alignment dividers for the Markdown table
    ]                                                                                   # Finish compiling the static upper lines of the Markdown document
    for _, r in final_df.iterrows():                                                    # Iterate through every row in the sorted dataframe to build the table
        md.append(                                                                      # Append a newly constructed row string to the Markdown document list
            f"| {int(r['final_rank'])} | {int(r['sample_id'])} | {r['final_decision']} | {bool(r['stage08_pass'])} | " # Format rank, ID, decision, and pass boolean
            f"{float(r['final_multimodal_rank_score']):.6f} | {float(r['esmfold_mean_plddt']):.2f} | " # Format score and overall pLDDT structural metrics
            f"{float(r['mutation_site_mean_plddt']):.2f} | {float(r['rmsd_to_selected_seed']):.3f} | {r['stage08_decision_reason']} |" # Format localized pLDDT, RMSD, and rationale
        )                                                                               # Finish appending the data row string to the Markdown document list
    md.extend([                                                                         # Extend the Markdown list with the remaining trailing deliverables section
        "",                                                                             # Insert a blank line below the completed Markdown table
        "## Deliverables",                                                              # Insert a Header 2 for the deliverables references
        f"- final table: `{csv_path}`",                                                 # Format and append a code-block formatted link to the CSV
        f"- final FASTA: `{fasta_path}`",                                               # Format and append a code-block formatted link to the FASTA
        f"- summary JSON: `{json_path}`",                                               # Format and append a code-block formatted link to the JSON
    ])                                                                                  # Finish extending the Markdown document list
    md_path.write_text("\n".join(md))                                                   # Join all Markdown lines with line breaks and save to file

    print(f"Wrote: {csv_path}")                                                         # Print a console notification that the CSV was successfully written
    print(f"Wrote: {fasta_path}")                                                       # Print a console notification that the FASTA was successfully written
    print(f"Wrote: {json_path}")                                                        # Print a console notification that the JSON was successfully written
    print(f"Wrote: {md_path}")                                                          # Print a console notification that the Markdown was successfully written
    print("\nFinal table preview:")                                                     # Print a spacer and header for the console preview representation
    print(final_df[[                                                                    # Begin extracting the subset of columns to display in the console preview
        "final_rank", "sample_id", "final_decision", "stage08_pass",                    # Select rank, ID, decision, and boolean pass status
        "final_multimodal_rank_score", "esmfold_mean_plddt", "mutation_site_mean_plddt",# Select the primary scoring and local/global confidence metrics
        "rmsd_to_selected_seed", "stage08_decision_reason"                              # Select the measured RMSD drift and the explanation string
    ]].to_string(index=False))                                                          # Render the subset dataframe as a formatted string without indexes and print


if __name__ == "__main__":                                                              # Check if the python script is executed as the primary module
    main()                                                                              # Invoke the main execution function