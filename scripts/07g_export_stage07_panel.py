"""Stage 07g: export the top Stage 07 validation panel to CSV, FASTA, and JSON.

This script turns the final Stage 07 ranking table into a compact, validation-ready handoff bundle.
It keeps the diversity-aware validation panel when present, then writes CSV, FASTA, and JSON outputs
that can be consumed by downstream structural validation and reporting steps.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from phageforge.stage07_utils import write_fasta, write_json


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 07 panel exporter."""
    parser = argparse.ArgumentParser(description="Export the Stage 07 validation panel to compact deliverables.")
    parser.add_argument("--ranked_csv", type=str, required=True, help="Ranked candidate CSV from 07e_rank_multimodal_candidates.py.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory where FASTA, CSV, and JSON exports will be written.")
    parser.add_argument("--top_k", type=int, default=5, help="Number of panel members to export.")
    return parser.parse_args()  # Return the parsed CLI namespace.


def main() -> None:
    # Load the final Stage 07e ranked table and choose the intended validation panel.
    args = parse_args()                             # Read the command-line arguments for this export run.
    ranked_df = pd.read_csv(args.ranked_csv)        # Read the final Stage 07 ranking table.
    if "selected_for_panel" in ranked_df.columns and ranked_df["selected_for_panel"].any():  # Prefer the explicit diverse panel when available.
        panel_df = ranked_df.loc[ranked_df["selected_for_panel"]].sort_values(["rank_diverse", "rank_raw"]).head(args.top_k).copy()  # Keep the requested top-k panel rows.
    else:  # Fall back to the top-ranked rows if no explicit panel flag exists.
        panel_df = ranked_df.sort_values(["rank_diverse", "rank_raw"]).head(args.top_k).copy()  # Use the best-ranked rows directly.

    # Prepare the output paths for the three compact handoff artifacts.
    out_dir = Path(args.output_dir)                                         # Convert the output directory argument into a Path object.
    out_dir.mkdir(parents=True, exist_ok=True)                              # Create the output directory if it does not exist.
    csv_path = out_dir / f"top{len(panel_df)}_validation_panel.csv"         # Define the CSV export path.
    fasta_path = out_dir / f"top{len(panel_df)}_validation_panel.fasta"     # Define the FASTA export path.
    json_path = out_dir / f"top{len(panel_df)}_validation_panel.json"       # Define the JSON export path.

    # Write the CSV, FASTA, and JSON outputs used by downstream validation steps.
    panel_df.to_csv(csv_path, index=False)  # Persist the compact tabular panel.
    write_fasta(  # Write FASTA records with score-rich headers.
        [
            (
                f"stage07_sample{int(row['sample_id'])}|diverse_rank={int(row.get('rank_diverse', 0))}|score={row.get('final_multimodal_rank_score', 0.0):.6f}",
                str(row["candidate_sequence"]),
            )
            for _, row in panel_df.iterrows()
        ],
        fasta_path,
    )
    write_json(  # Write a compact JSON summary for easier programmatic downstream use.
        {
            "top_k": int(len(panel_df)),
            "csv_path": str(csv_path),
            "fasta_path": str(fasta_path),
            "candidates": panel_df[[
                "sample_id",
                "rank_diverse",
                "rank_raw",
                "final_multimodal_rank_score",
                "target_score",
                "strict_manifold_score",
                "structure_score",
                "mutation_positions",
                "candidate_sequence",
            ]].to_dict(orient="records"),
        },
        json_path,
    )
    print(f"Wrote: {csv_path}")     # Report the written CSV path.
    print(f"Wrote: {fasta_path}")   # Report the written FASTA path.
    print(f"Wrote: {json_path}")    # Report the written JSON path.


if __name__ == "__main__":
    main()  # Execute the CLI entrypoint when the script is run directly.
