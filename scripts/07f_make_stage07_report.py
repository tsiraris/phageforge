"""Stage 07f: Write a compact markdown report for the Stage 07 generation and ranking run."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from phageforge.stage07_utils import read_json


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 07 markdown report."""
    ap = argparse.ArgumentParser(description="Create a compact markdown report for Stage 07.")
    ap.add_argument("--context_json", type=str, required=True, help="Stage 07 context JSON from 07a_prepare_stage07_design_context.py.")
    ap.add_argument("--ranked_csv", type=str, required=True, help="Ranked candidate CSV from 07e_rank_multimodal_candidates.py.")
    ap.add_argument("--report_md", type=str, required=True, help="Where to write the markdown report.")
    ap.add_argument("--top_k", type=int, default=5, help="How many top-ranked candidates to summarize in the report.")
    return ap.parse_args()


def main() -> None:
    # Read the Stage 07 context and ranked candidate table.
    args = parse_args()
    context = read_json(args.context_json)
    ranked_df = pd.read_csv(args.ranked_csv)
    top_df = ranked_df.sort_values(["rank_diverse", "rank_raw"]).head(args.top_k).copy()

    # Build a concise decision-oriented markdown report instead of only dumping raw JSON again.
    lines = [
        "# Stage 07 report",
        "",
        "## Run summary",
        f"- Target host: **{context['target_host']}**",
        f"- Selected seed: **{context['selected_seed']['seed_protein_id']}** from **{context['selected_seed']['source_host']}**",
        f"- Family members: **{context['family_context']['family_member_count']}**",
        f"- Editable hotspot count: **{len(context['editable_region'].get('hotspot_positions', []))}**",
        f"- Raw ranked candidates: **{len(ranked_df)}**",
        f"- Validation panel size: **{int(ranked_df.get('selected_for_panel', pd.Series(dtype=bool)).sum())}**",
        f"- Used local ESM3: **{bool(ranked_df.get('used_local_esm3', pd.Series([False])).any())}**",
        f"- Used Forge / API ESM3: **{bool(ranked_df.get('used_esm3_api', pd.Series([False])).any())}**",
        f"- Used ESM2 fallback: **{bool(ranked_df.get('used_esm2_fallback', pd.Series([False])).any())}**",
        "",
        "## Top candidates",
        "",
        "| diverse_rank | raw_rank | sample_id | regime | target_score | manifold | structure | final_score | selected_panel | mutations |",
        "|---:|---:|---:|---|---:|---:|---:|---:|:---:|---|",
    ]
    for _, row in top_df.iterrows():
        lines.append(
            f"| {int(row.get('rank_diverse', 0))} | {int(row.get('rank_raw', 0))} | {int(row['sample_id'])} | {row.get('generation_regime', 'balanced')} | "
            f"{row.get('target_score', 0.0):.6f} | {row.get('strict_manifold_score', 0.0):.6f} | {row.get('structure_score', 0.0):.6f} | "
            f"{row.get('final_multimodal_rank_score', 0.0):.6f} | {bool(row.get('selected_for_panel', False))} | {row.get('mutation_positions', '')} |"
        )

    lines.extend([
        "",
        "## Context JSON",
        "```json",
        Path(args.context_json).read_text(encoding='utf-8'),
        "```",
    ])

    out_path = Path(args.report_md)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
