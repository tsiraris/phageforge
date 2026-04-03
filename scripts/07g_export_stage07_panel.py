"""Stage 07g: Export the top Stage 07 validation panel to CSV, FASTA, and JSON."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from phageforge.stage07_utils import write_fasta, write_json


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 07 panel exporter."""
    ap = argparse.ArgumentParser(description="Export the Stage 07 validation panel to compact deliverables.")
    ap.add_argument("--ranked_csv", type=str, required=True, help="Ranked candidate CSV from 07e_rank_multimodal_candidates.py.")
    ap.add_argument("--output_dir", type=str, required=True, help="Directory where FASTA, CSV, and JSON exports will be written.")
    ap.add_argument("--top_k", type=int, default=5, help="Number of panel members to export.")
    return ap.parse_args()


def main() -> None:
    # Read the ranked table and choose the diversity-aware panel when available.
    args = parse_args()
    ranked_df = pd.read_csv(args.ranked_csv)
    if "selected_for_panel" in ranked_df.columns and ranked_df["selected_for_panel"].any():
        panel_df = ranked_df.loc[ranked_df["selected_for_panel"]].sort_values(["rank_diverse", "rank_raw"]).head(args.top_k).copy()
    else:
        panel_df = ranked_df.sort_values(["rank_diverse", "rank_raw"]).head(args.top_k).copy()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"top{len(panel_df)}_validation_panel.csv"
    fasta_path = out_dir / f"top{len(panel_df)}_validation_panel.fasta"
    json_path = out_dir / f"top{len(panel_df)}_validation_panel.json"

    # Write the compact CSV, FASTA, and JSON handoff artifacts for downstream validation.
    panel_df.to_csv(csv_path, index=False)
    write_fasta(
        [
            (
                f"stage07_sample{int(row['sample_id'])}|diverse_rank={int(row.get('rank_diverse', 0))}|score={row.get('final_multimodal_rank_score', 0.0):.6f}",
                str(row["candidate_sequence"]),
            )
            for _, row in panel_df.iterrows()
        ],
        fasta_path,
    )
    write_json(
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
    print(f"Wrote: {csv_path}")
    print(f"Wrote: {fasta_path}")
    print(f"Wrote: {json_path}")


if __name__ == "__main__":
    main()
