"""
==========================================================
07f: Generate a compact markdown report for Stage 07 runs.
==========================================================

This script writes a small human-readable report for one completed Stage 07 run.
It is intentionally lightweight to quickly summarize results without first
building a large notebook or PDF pipeline.
"""

from __future__ import annotations                                                    # Enable postponed annotation evaluation for cleaner typing.
import argparse                                                                       # Parse command-line arguments.
from pathlib import Path                                                              # Build report output paths robustly.
import pandas as pd                                                                   # Read the ranked candidate table and build the markdown summary.


def parse_args() -> argparse.Namespace:
    """Return command-line arguments for writing a compact Stage 07 markdown report."""
    ap = argparse.ArgumentParser(description="Write a markdown summary for a completed Stage 07 run.")  # Create the parser for the reporting stage.
    ap.add_argument("--context_json", type=str, required=True, help="Stage 07 context JSON from 07a.")  # Point to the context JSON used for the run.
    ap.add_argument("--ranked_csv", type=str, required=True, help="Final ranked CSV from 07e.")         # Point to the final ranked candidate table.
    ap.add_argument("--report_md", type=str, required=True, help="Where to write the markdown report.") # Point to the markdown output path.
    ap.add_argument("--top_k", type=int, default=20, help="How many top candidates to include in the report table.")  # Control the number of displayed top candidates.
    return ap.parse_args()                                                                              # Parse the CLI and return the resulting namespace.


def summarize_provenance(top: pd.DataFrame) -> list[str]:
    """Return short markdown bullet points summarizing which generator backends and ESM settings produced the shortlist."""
    lines = []                                                                                          # Collect the provenance summary lines before returning them.
    if "generator_mode" in top.columns:                                                                 # Summarize which generation backends appear in the shortlist.
        counts = top["generator_mode"].astype(str).value_counts()                                       # Count how often each backend appears in the top table.
        lines.append("- Generator modes: " + ", ".join(f"{k} ({v})" for k, v in counts.items()))        # Report backend counts in compact markdown form.
    if "esm3_model" in top.columns:                                                                     # Summarize which Forge model labels appear when present.
        models = [x for x in top["esm3_model"].dropna().astype(str).unique().tolist() if x != ""]       # Keep only non-empty model labels.
        if models:                                                                                      # Only emit the line when at least one explicit model label exists.
            lines.append("- ESM3 models observed: " + ", ".join(models))                                # Report distinct Forge model names used in the shortlist.
    if "esm3_temperature" in top.columns:                                                               # Summarize temperature settings when present.
        temps = top["esm3_temperature"].dropna().astype(float)                                          # Keep only non-missing temperatures.
        if len(temps) > 0:                                                                              # Only emit the line when at least one temperature value exists.
            lines.append(f"- Temperature range in shortlist: {temps.min():.3f} to {temps.max():.3f}")   # Report the observed temperature range.
    if "esm3_num_steps" in top.columns:                                                                 # Summarize Forge iterative decode steps when present.
        steps = top["esm3_num_steps"].dropna().astype(int).unique().tolist()                            # Gather distinct decode-step values.
        if steps:                                                                                       # Only emit the line when at least one explicit num_steps value exists.
            lines.append("- ESM3 num_steps values observed: " + ", ".join(str(x) for x in sorted(steps)))  # Report distinct num_steps settings.
    return lines                                                                                        # Return the assembled provenance summary bullets.


# Main entrypoint: Gather context and writes the markdown report.
def main() -> None:
    # Read the full context JSON and the final ranked candidate table.
    args = parse_args()                                                               # Parse command-line arguments.
    context_text = Path(args.context_json).read_text(encoding="utf-8")                # Read the full context JSON as text for direct inclusion in the report.
    ranked = pd.read_csv(args.ranked_csv)                                             # Read the final ranked candidate table from disk.
    # Keep only the top-k rows for the summary section.
    top = ranked.head(args.top_k).copy()                                              

    # Write the markdown report.
    lines = []                                                                        # Collect report lines before joining them into one markdown string.
    lines.append("# Stage 07 report")                                                 # Write the top-level report title.
    lines.append("")                                                                  # Add a blank line for markdown readability.
    lines.append("## Stage 07 context JSON")                                          # Add a section heading for the run context.
    lines.append("```json")                                                           # Start a fenced JSON code block.
    lines.append(context_text)                                                        # Insert the raw context JSON text.
    lines.append("```")                                                               # Close the fenced JSON block.
    lines.append("")                                                                  # Add a blank line between sections.
    lines.append("## Generation provenance summary")                                  # Add a section explaining how the shortlist was actually generated.
    prov_lines = summarize_provenance(top)                                            # Build concise provenance bullets from the top-k candidate table.
    if prov_lines:                                                                    # Only print the summary when provenance fields are present.
        lines.extend(prov_lines)                                                      # Append each provenance bullet line to the markdown report.
    else:                                                                             # Handle cases where provenance columns are absent.
        lines.append("- No explicit generation provenance columns were available.")   # Insert a fallback note instead of leaving the section empty.
    lines.append("")                                                                  # Add a blank line between sections.
    lines.append("## Top candidates")                                                 # Add a section heading for the top-ranked candidates.

    cols = [                                                                          # Choose a compact set of important ranking columns when present.
        c for c in [
            "candidate_sequence",
            "target_host",
            "generator_mode",
            "esm3_model",
            "esm3_temperature",
            "esm3_num_steps",
            "editable_hotspot_count",
            "mutation_positions",
            "target_score",
            "strict_manifold_score",
            "family_cosine",
            "target_anchor_cosine",
            "structure_score",
            "tissue_score",
            "mutation_penalty",
            "final_multimodal_rank_score",
        ]
        if c in top.columns
    ]
    if cols:                                                                          # Only render a table when at least one relevant column exists.
        lines.append(top[cols].to_markdown(index=False))                              # Convert the top-k candidate table into markdown format.
    else:                                                                             # Handle degenerate cases where the ranked CSV lacks expected columns.
        lines.append("_No ranking columns available._")                               # Insert a fallback note instead of an empty table.

    # Write the markdown report to disk.
    report_path = Path(args.report_md)                                                # Convert the output markdown path into a Path object.
    report_path.parent.mkdir(parents=True, exist_ok=True)                             # Create the report directory if needed.
    report_path.write_text("\n".join(lines), encoding="utf-8")                        # Join the report lines with newlines and write the markdown file.
    print(f"Wrote: {report_path}")                                                    # Print the report path for quick confirmation.


if __name__ == "__main__":                                                            # Standard Python entrypoint guard.
    main()                                                                            # Execute the report-writing CLI.
