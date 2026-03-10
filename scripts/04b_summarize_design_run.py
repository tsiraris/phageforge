from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main():
    # Create an ArgumentParser object to run the same script with different parameters without editing code.
    ap = argparse.ArgumentParser()
    ap.add_argument("--run_dir", type=str, required=True)
    ap.add_argument("--top_k", type=int, default=10)
    args = ap.parse_args()

    # Load data
    run_dir = Path(args.run_dir)
    top_csv = run_dir / "top_candidates.csv"
    all_csv = run_dir / "all_candidates.csv"
    metadata_json = run_dir / "run_metadata.json"

    # Check that the files exist
    if not top_csv.exists():
        raise FileNotFoundError(f"Missing file: {top_csv}")
    if not all_csv.exists():
        raise FileNotFoundError(f"Missing file: {all_csv}")
    if not metadata_json.exists():
        raise FileNotFoundError(f"Missing file: {metadata_json}")

    # Load data and metadata
    top_df = pd.read_csv(top_csv)
    all_df = pd.read_csv(all_csv)
    metadata = json.loads(metadata_json.read_text())

    # Create summary files and preview
    pred_label_counts = (
        top_df["pred_label"]
        .value_counts()
        .rename_axis("pred_label")
        .reset_index(name="count")
    )
    
    # Count the number of mutations
    mutation_count_counts = (
        top_df["n_mutations"]
        .value_counts()
        .sort_index()
        .rename_axis("n_mutations")
        .reset_index(name="count")
    )
    
    # Compute summary statistics
    score_summary = pd.DataFrame(
        {
            "metric": ["min", "mean", "median", "max"],
            "target_score": [
                float(top_df["target_score"].min()),
                float(top_df["target_score"].mean()),
                float(top_df["target_score"].median()),
                float(top_df["target_score"].max()),
            ],
        }
    )

    # Save summary files
    pred_label_counts.to_csv(run_dir / "pred_label_counts.csv", index=False)
    mutation_count_counts.to_csv(run_dir / "mutation_count_counts.csv", index=False)
    score_summary.to_csv(run_dir / "target_score_summary.csv", index=False)

    # Save preview
    preview_cols = [
        c for c in [
            "candidate_id",
            "source_host",
            "target_host",
            "pred_label",
            "target_score",
            "mutations",
            "n_mutations",
        ] if c in top_df.columns
    ]
    top_preview = top_df[preview_cols].head(args.top_k)

    # Create the summary file
    lines = []
    lines.append("# Design Run Summary")
    lines.append("")
    lines.append("## Run metadata")
    lines.append("")
    lines.append(f"- Seed protein ID: `{metadata.get('seed_protein_id')}`")
    lines.append(f"- Virus accession: `{metadata.get('virus_accession')}`")
    lines.append(f"- Source host: `{metadata.get('source_host')}`")
    lines.append(f"- Target host: `{metadata.get('target_host')}`")
    lines.append(f"- ESM model: `{metadata.get('esm_model')}`")
    lines.append(f"- Rounds: `{metadata.get('rounds')}`")
    lines.append(f"- Candidates per round: `{metadata.get('candidates_per_round')}`")
    lines.append(f"- Min mutations: `{metadata.get('min_mutations', metadata.get('mutations_per_candidate'))}`")
    lines.append(f"- Max mutations: `{metadata.get('max_mutations', metadata.get('mutations_per_candidate'))}`")
    lines.append(f"- Top-K kept per round: `{metadata.get('keep_top_k')}`")
    lines.append(f"- Proposal top-K: `{metadata.get('proposal_top_k', 'N/A')}`")
    lines.append(f"- Seed: `{metadata.get('seed')}`")
    lines.append(f"- Total candidates evaluated: `{metadata.get('n_total_candidates')}`")
    lines.append(f"- Top candidates saved: `{metadata.get('n_top_candidates_saved')}`")
    lines.append("")

    lines.append("## Target score summary")
    lines.append("")
    lines.append(score_summary.to_markdown(index=False))
    lines.append("")

    lines.append("## Predicted label distribution in top candidates")
    lines.append("")
    lines.append(pred_label_counts.to_markdown(index=False))
    lines.append("")

    lines.append("## Mutation count distribution in top candidates")
    lines.append("")
    lines.append(mutation_count_counts.to_markdown(index=False))
    lines.append("")

    lines.append(f"## Top {args.top_k} candidates")
    lines.append("")
    lines.append(top_preview.to_markdown(index=False))
    lines.append("")

    # Save report
    report_path = run_dir / "design_run_report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"✅ Saved: {run_dir / 'pred_label_counts.csv'}")
    print(f"✅ Saved: {run_dir / 'mutation_count_counts.csv'}")
    print(f"✅ Saved: {run_dir / 'target_score_summary.csv'}")
    print(f"✅ Saved: {report_path}")


if __name__ == "__main__":
    main()