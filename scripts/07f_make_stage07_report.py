"""Stage 07f: write a compact markdown report for a completed Stage 07 generation run.

This script reads the Stage 07 design context and the final diversity-aware ranking table,
then writes a short markdown report that summarizes the selected seed, the target host,
the size of the generated panel, and the top candidates that should move forward.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from phageforge.stage07_utils import read_json


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 07 markdown report builder."""
    parser = argparse.ArgumentParser(description="Create a compact markdown report for Stage 07.")
    parser.add_argument("--context_json", type=str, required=True, help="Stage 07 context JSON from 07a_prepare_stage07_design_context.py.")
    parser.add_argument("--ranked_csv", type=str, required=True, help="Ranked candidate CSV from 07e_rank_multimodal_candidates.py.")
    parser.add_argument("--report_md", type=str, required=True, help="Where to write the markdown report.")
    parser.add_argument("--top_k", type=int, default=5, help="How many top-ranked candidates to summarize in the report.")
    return parser.parse_args()  # Return the parsed CLI namespace.


def main() -> None:
    # Load the Stage 07 design context and the final ranking table.
    args = parse_args()                         # Read the command-line configuration.
    context = read_json(args.context_json)      # Load the Stage 07 context JSON.
    ranked_df = pd.read_csv(args.ranked_csv)    # Read the final diversity-aware ranking table.
    top_df = ranked_df.sort_values(["rank_diverse", "rank_raw"]).head(args.top_k).copy()  # Keep only the top rows for the compact summary section.

    # Build a concise markdown report that is easy to archive and share.
    lines = [  # Assemble the markdown report content in order.
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

    # Add one markdown table row per candidate in the compact top-k summary.
    for _, row in top_df.iterrows():    # Iterate over the shortlisted rows in final ranking order.
        lines.append(                   # Append a formatted markdown table row for one candidate.
            f"| {int(row.get('rank_diverse', 0))} | {int(row.get('rank_raw', 0))} | {int(row['sample_id'])} | {row.get('generation_regime', 'balanced')} | "
            f"{row.get('target_score', 0.0):.6f} | {row.get('strict_manifold_score', 0.0):.6f} | {row.get('structure_score', 0.0):.6f} | "
            f"{row.get('final_multimodal_rank_score', 0.0):.6f} | {bool(row.get('selected_for_panel', False))} | {row.get('mutation_positions', '')} |"
        )

    # Attach the full context JSON at the end for provenance-complete archiving.
    lines.extend([  # Add the literal context JSON block below the human-readable summary.
        "",
        "## Context JSON",
        "```json",
        Path(args.context_json).read_text(encoding="utf-8"),
        "```",
    ])

    # Write the markdown report to disk.
    out_path = Path(args.report_md)                             # Convert the report path into a Path object.
    out_path.parent.mkdir(parents=True, exist_ok=True)          # Create the parent directory if needed.
    out_path.write_text("\n".join(lines), encoding="utf-8")     # Persist the markdown report.
    print(f"Wrote: {out_path}")                                 # Report the generated report path.


if __name__ == "__main__":
    main()  # Execute the CLI entrypoint when the script is run directly.
