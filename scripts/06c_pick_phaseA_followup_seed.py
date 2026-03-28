"""
Phase A / Step 06c: Choosing the follow-up seed
===================

This script chooses the best follow-up seed from the first host-ladder step and writes it to FASTA
so the next host-ladder step can start immediately.

The script is intentionally simple and transparent:
- read the first step's `top_candidates.csv`,
- choose the strongest candidate after deduplicating identical sequences,
- write a FASTA file and a small JSON summary,
- print the exact second-step command arguments you should use next.
"""

from __future__ import annotations
import argparse                                                                        # Parse command-line arguments.
import json                                                                            # Read the Phase A plan and write the follow-up metadata JSON.
from pathlib import Path                                                               # Work with filesystem paths cleanly.
import pandas as pd                                                                    # Read and rank top candidate tables.


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for choosing the Phase A follow-up seed."""
    ap = argparse.ArgumentParser(description="Pick the best follow-up seed for the next Phase A ladder step.")
    ap.add_argument("--phaseA_plan_json", type=str, required=True, help="Plan JSON produced by 06a_select_phaseA_family.py.")
    ap.add_argument("--top_candidates_csv", type=str, required=True, help="Top candidate CSV produced by the first Phase A optimization step.")
    ap.add_argument("--from_target", type=str, default="Enterobacter", help="Human-readable label for the source ladder step.")
    ap.add_argument("--next_target", type=str, default="Acinetobacter", help="Human-readable label for the next ladder step.")
    ap.add_argument("--output_dir", type=str, required=True, help="Directory where the follow-up FASTA and metadata will be written.")
    return ap.parse_args()


def main():
    # Read command-line arguments, create output directory, and read the Phase A plan.
    args = parse_args()                                                                  # Parse command-line arguments.
    output_dir = Path(args.output_dir)                                                   # Convert the output directory into a Path object.
    output_dir.mkdir(parents=True, exist_ok=True)                                        # Ensure the output directory exists before writing files.

    with open(args.phaseA_plan_json, "r", encoding="utf-8") as handle:                   # Read the Phase A plan so the output metadata retains project context.
        phaseA_plan = json.load(handle)

    # Load the first ladder step's final candidate table, sort by selection score, and choose the best follow-up seed as the follow-up seed.
    df = pd.read_csv(args.top_candidates_csv)                                            # Load the first ladder step's final candidate table.
    if df.empty:
        raise ValueError(f"Top candidate table is empty: {args.top_candidates_csv}")     # Fail loudly if the previous run produced no candidates.

    ranking_columns = [col for col in ["selection_score", "target_score", "family_cosine", "target_anchor_cosine", "n_mutations"] if col in df.columns]  # Keep only ranking columns that actually exist in the file.
    ascending = [False, False, False, False, True][: len(ranking_columns)]               # Rank descending on scores and ascending on mutation count.
    best_df = (
        df.sort_values(by=ranking_columns, ascending=ascending)                          # Sort by the strongest available ranking signal.
        .drop_duplicates(subset=["aa_sequence"])                                         # Collapse duplicate sequences so the follow-up seed is truly unique.
        .reset_index(drop=True)
    )
    best_row = best_df.iloc[0]                                                           # Select the strongest unique candidate as the follow-up seed.

    # Write the follow-up seed to FASTA and a small JSON summary.
    fasta_name = f"phaseA_followup_seed_{args.from_target}_to_{args.next_target}.fasta"  # Create a descriptive FASTA filename for the ladder transition.
    fasta_path = output_dir / fasta_name                                                 # Build the full FASTA output path.
    fasta_header = f"phaseA_{args.from_target}_best_for_{args.next_target}|candidate_id={best_row['candidate_id']}"  # Include the candidate ID in the FASTA header for traceability.
    with open(fasta_path, "w", encoding="utf-8") as handle:                              # Write the chosen sequence to FASTA.
        handle.write(f">{fasta_header}\n")
        handle.write(str(best_row["aa_sequence"]).strip() + "\n")

    summary = {
        "from_target": args.from_target,
        "next_target": args.next_target,
        "phaseA_plan_json": str(args.phaseA_plan_json),
        "source_top_candidates_csv": str(args.top_candidates_csv),
        "chosen_candidate_id": str(best_row["candidate_id"]),
        "chosen_mutations": str(best_row.get("mutations", "")),
        "chosen_target_score": float(best_row.get("target_score", 0.0)),
        "chosen_selection_score": float(best_row.get("selection_score", 0.0)),
        "chosen_family_cosine": float(best_row.get("family_cosine", 0.0)),
        "chosen_target_anchor_cosine": float(best_row.get("target_anchor_cosine", 0.0)),
        "chosen_n_mutations": int(best_row.get("n_mutations", 0)),
        "fasta_path": str(fasta_path),
        "canonical_seed_protein_id": phaseA_plan["canonical_seed"]["seed_protein_id"],
    }                                                                                    # Store a compact summary of why this candidate became the follow-up seed.
    with open(output_dir / "phaseA_followup_seed_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)                                             # Write the follow-up seed metadata JSON.

    print(f"Chosen follow-up candidate: {best_row['candidate_id']}")                     # Print the chosen candidate ID for quick confirmation.
    print(f"Saved FASTA: {fasta_path}")                                                  # Print the FASTA path you should pass into the next ladder step.
    print("\nNext command (core arguments only):")                                       # Print the exact key arguments needed for the next Phase A step.
    print(
        f"--seed_fasta {fasta_path} --seed_label phaseA_{args.from_target}_best --target_host {args.next_target}"
    )                                                                                    # Provide a copy-paste-ready argument fragment for the next run.


if __name__ == "__main__":
    main()
