#!/usr/bin/env python
"""Stage 11e: Build the Stage 11 comparative report.

Reads everything produced by 11a → 11d and emits:

- `stage11_report.md`: a human-readable summary with the headline outcome,
  Baseline Qualification result, search statistics, top-3 panel table, and a
  side-by-side comparison against historical Stage 08/09/10 baselines (each
  of which produced 0/N pass rates).
- `stage11_report_summary.json`: the same numbers in machine-readable form.
- `stage11_handoff.fasta`: the top-3 candidate sequences with metadata-rich
  headers following the Stage 07g convention.

The report explicitly contrasts this Stage 11 run against the
corrupted-chassis baselines so a reviewer can immediately see whether the
Baseline Qualification fix landed structurally viable redesigns or whether
the seed-chassis hypothesis is falsified.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from phageforge.stage11_utils import (
    EXIT_INPUT_ERROR,
    EXIT_OK,
    read_json,
    write_fasta,
    write_json,
)


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for the Stage 11 reporter."""
    ap = argparse.ArgumentParser(description="Build the final Stage 11 comparative report.")                                                                                # Create argument parser.
    ap.add_argument("--stage11_context_json", type=str, required=True, help="stage11_context.json produced by 11a.")                                                        # Add context file argument. Explain context input.
    ap.add_argument("--search_csv", type=str, required=True, help="Full Stage 11 search CSV produced by 11b.")                                                              # Add search CSV argument. Explain search input.
    ap.add_argument("--prefilter_csv", type=str, required=True, help="Top-K prefilter CSV produced by 11c.")                                                                # Add prefilter CSV argument. Explain prefilter input.
    ap.add_argument("--validation_csv", type=str, default=None, help="Stage 11 structural validation summary CSV produced by 11d (08a output).")                            # Add validation CSV argument. Explain validation input.
    ap.add_argument("--baseline_stage08_csv", type=str, default=None, help="Historical Stage 08 structural validation summary CSV (for comparison).")                       # Add S08 baseline argument. Explain S08 input.
    ap.add_argument("--baseline_stage09_csv", type=str, default=None, help="Historical Stage 09 structural validation summary CSV (for comparison).")                       # Add S09 baseline argument. Explain S09 input.
    ap.add_argument("--baseline_stage10_csv", type=str, default=None, help="Historical Stage 10 structural validation summary CSV (for comparison).")                       # Add S10 baseline argument. Explain S10 input.
    ap.add_argument("--baseline_qualification_json", type=str, default=None, help="Optional baseline_qualification.json (auto-resolved next to the context if omitted).")   # Add qualification JSON argument. Explain qualification input.
    ap.add_argument("--out_dir", type=str, required=True, help="Directory where the Stage 11 report artifacts will be written.")                                            # Add output directory argument. Explain output dir.
    return ap.parse_args()                                                                                                                                                  # Parse and return arguments.


def _safe_read_csv(path: str | Path | None) -> pd.DataFrame:
    """Reads a CSV file from a given path safely. Returns an empty DataFrame if the file is missing or unreadable, preventing execution interruption.

    Example:
        df = _safe_read_csv("results.csv")
        Returns loaded DataFrame, or an empty DataFrame if the path is invalid.
    """
    if not path:                                                                                                                # Check if path is None or empty.
        return pd.DataFrame()                                                                                                   # Return empty frame if no path.
    p = Path(path)                                                                                                              # Convert to Path object.
    if not p.exists():                                                                                                          # Verify file existence.
        print(f"[WARN] Skipping missing CSV: {p}", flush=True)                                                                  # Alert user if missing.
        return pd.DataFrame()                                                                                                   # Return empty frame.
    try:                                                                                                                        # Try to read the file.
        return pd.read_csv(p)                                                                                                   # Load CSV into DataFrame.
    except Exception as exc:                                                                                                    # Handle potential read errors.
        print(f"[WARN] Failed to read {p}: {exc}", flush=True)                                                                  # Alert user if read fails.
        return pd.DataFrame()                                                                                                   # Return empty frame on failure.


def summarize_validation(df: pd.DataFrame) -> dict:
    """Aggregates structural validation metrics (pLDDT, RMSD, pass rates) from a CSV into a dictionary for reporting. Useful for extracting key performance indicators.

    Example:
        metrics = summarize_validation(df)
        Returns: {'rows': 10, 'pass_count': 5, 'best_mean_plddt': 85.0, ...}
    """
    if df.empty:                                                                                                                # Handle empty input data.
        return {                                                                                                                # Initialize empty result structure.
            "rows": 0, "pass_count": 0,                                                                                         # Zero counts.
            "best_mean_plddt": None, "mean_mean_plddt": None,                                                                   # None values for metrics.
            "best_mutation_site_mean_plddt": None, "best_rmsd": None,                                                           # None values for metrics.
            "mean_rmsd": None, "best_mut_count": None,                                                                          # None values for metrics.
        }                                                                                                                       # Return zeroed dictionary.
    pass_col = "stage08_pass" if "stage08_pass" in df.columns else None                                                         # Identify pass column.
    plddt_col = "esmfold_mean_plddt" if "esmfold_mean_plddt" in df.columns else None                                            # Identify pLDDT column.
    rmsd_col = "rmsd_to_selected_seed" if "rmsd_to_selected_seed" in df.columns else None                                       # Identify RMSD column.
    return {                                                                                                                    # Build metric summary dictionary.
        "rows": int(len(df)),                                                                                                   # Calculate total row count.
        "pass_count": int(df[pass_col].astype(bool).sum()) if pass_col else 0,                                                  # Calculate pass count.
        "best_mean_plddt": float(df[plddt_col].max()) if plddt_col else None,                                                   # Get maximum pLDDT.
        "mean_mean_plddt": float(df[plddt_col].mean()) if plddt_col else None,                                                  # Get average pLDDT.
        "best_mutation_site_mean_plddt": (                                                                                      # Extract mutation-site pLDDT.
            float(df["mutation_site_mean_plddt"].max())                                                                         # Compute max mutation pLDDT.
            if "mutation_site_mean_plddt" in df.columns else None                                                               # Safety check.
        ),                                                                                                                      # Close expression.
        "best_rmsd": float(df[rmsd_col].min()) if rmsd_col else None,                                                           # Get minimum RMSD (best).
        "mean_rmsd": float(df[rmsd_col].mean()) if rmsd_col else None,                                                          # Get average RMSD.
        "best_mut_count": (                                                                                                     # Find mutation count of best entry.
            int(df.sort_values(plddt_col, ascending=False).iloc[0]["mutation_count"])                                           # Sort by pLDDT and get top.
            if (plddt_col and "mutation_count" in df.columns) else None                                                         # Ensure columns exist.
        ),                                                                                                                      # Close expression.
    }                                                                                                                           # Return populated metrics.


def _fmt(value) -> str:
    """Formats values for clean presentation in markdown, handling nulls, NaNs, and numerical precision.

    Example:
        fmt_val = _fmt(1234.5678) # returns "1235"
        fmt_val = _fmt(None)      # returns "—"
    """
    if value is None:                                                                                                           # Handle None types.
        return "—"                                                                                                              # Return dash placeholder.
    if isinstance(value, float):                                                                                                # Handle float values.
        if math.isnan(value):                                                                                                   # Check for NaN.
            return "—"                                                                                                          # Return dash placeholder.
        if abs(value) >= 1000:                                                                                                  # Check for large numbers.
            return f"{value:.0f}"                                                                                               # Format without decimals.
        return f"{value:.3f}"                                                                                                   # Format with 3 decimals.
    return str(value)                                                                                                           # Convert others to string.


def main() -> None:
    args = parse_args()                                                                                                         # Parse command line arguments.

    # ----- Resolve inputs -----
    context_path = Path(args.stage11_context_json)                                                                              # Convert context path.
    if not context_path.exists():                                                                                               # Validate context file existence.
        print(f"[ERROR] Missing Stage 11 context JSON: {context_path}", file=sys.stderr)                                        # Print error message to stderr.
        sys.exit(EXIT_INPUT_ERROR)                                                                                              # Exit with input error code.
    context = read_json(context_path)                                                                                           # Load context data.

    search_df = _safe_read_csv(args.search_csv)                                                                                 # Read full search log.
    prefilter_df = _safe_read_csv(args.prefilter_csv)                                                                           # Read prefilter panel.
    validation_df = _safe_read_csv(args.validation_csv)                                                                         # Read validation results.
    baseline_08 = _safe_read_csv(args.baseline_stage08_csv)                                                                     # Read S08 baseline.
    baseline_09 = _safe_read_csv(args.baseline_stage09_csv)                                                                     # Read S09 baseline.
    baseline_10 = _safe_read_csv(args.baseline_stage10_csv)                                                                     # Read S10 baseline.

    # ----- Baseline Qualification — pull from explicit JSON or co-located file -----
    if args.baseline_qualification_json:                                                                                        # Check for explicit JSON argument.
        bq_path = Path(args.baseline_qualification_json)                                                                        # Use provided path.
    else:                                                                                                                       # Resolve locally if missing.
        bq_path = context_path.parent / "baseline_qualification.json"                                                           # Path near context.
    if bq_path.exists():                                                                                                        # Verify existence.
        baseline_qualification = read_json(bq_path)                                                                             # Read explicit qualification.
    else:                                                                                                                       # Fallback to context data.
        baseline_qualification = dict(context.get("baseline_qualification", {}))                                                # Use context backup.

    # ----- Headline metadata from the context -----
    seed_protein_id = str(context["selected_seed"].get("seed_protein_id", ""))                                                  # Extract seed ID.
    seed_length = int(context["selected_seed"].get("sequence_length", 0))                                                       # Extract sequence length.
    target_host = str(context.get("target_host", ""))                                                                           # Extract target host.
    source_host = str(context.get("source_host", ""))                                                                           # Extract source host.
    family_member_count = int(context.get("family_context", {}).get("family_member_count", 0))                                  # Extract family count.
    target_member_count = int(context.get("target_context", {}).get("target_member_count", 0))                                  # Extract target count.
    editable_region = dict(context.get("editable_region", {}))                                                                  # Extract edit space.
    hard_positions = list(editable_region.get("hard_positions", []))                                                            # Extract hard positions.
    soft_positions = list(editable_region.get("soft_positions", []))                                                            # Extract soft positions.
    min_mutations = int(editable_region.get("min_mutations", 0))                                                                # Extract min mutations.
    max_mutations = int(editable_region.get("max_mutations", 0))                                                                # Extract max mutations.

    # ----- Search summary -----
    score_col = next(                                                                                                           # Find highest priority score column.
        (c for c in ["stage10_composite_score", "stage11_composite_score", "final_multimodal_rank_score"]                       # Check candidate columns.
         if c in search_df.columns),                                                                                            # Check existence in DF.
        None,                                                                                                                   # Return None if no column.
    )                                                                                                                           # Close iterator.
    search_summary = {                                                                                                          # Compile search metrics.
        "rows": int(len(search_df)),                                                                                            # Count total search candidates.
        "rounds_completed": (                                                                                                   # Determine search depth.
            int(search_df["round_index"].max()) if (len(search_df) and "round_index" in search_df.columns) else 0                # Get max round index.
        ),                                                                                                                      # Close expression.
        "best_composite_score": (                                                                                               # Identify best score.
            float(search_df[score_col].max()) if (score_col and len(search_df)) else float("nan")                               # Calculate max or NaN.
        ),                                                                                                                      # Close expression.
        "best_target_probability": (                                                                                            # Identify best target probability.
            float(search_df["target_probability"].max())                                                                        # Calculate max probability.
            if (len(search_df) and "target_probability" in search_df.columns) else float("nan")                                 # Safety check.
        ),                                                                                                                      # Close expression.
        "best_if1_log_likelihood": (                                                                                            # Identify best IF1 score.
            float(search_df["if1_log_likelihood"].max())                                                                        # Calculate max likelihood.
            if (len(search_df) and "if1_log_likelihood" in search_df.columns) else float("nan")                                 # Safety check.
        ),                                                                                                                      # Close expression.
        "best_seed_identity": (                                                                                                 # Identify best seed identity.
            float(search_df["seed_identity"].max())                                                                             # Calculate max identity.
            if (len(search_df) and "seed_identity" in search_df.columns) else float("nan")                                      # Safety check.
        ),                                                                                                                      # Close expression.
        "mutation_count_histogram": (                                                                                           # Create mutation distribution.
            {int(k): int(v) for k, v in search_df["mutation_count"].value_counts().items()}                                     # Map counts to frequency.
            if (len(search_df) and "mutation_count" in search_df.columns) else {}                                               # Safety check.
        ),                                                                                                                      # Close expression.
    }                                                                                                                           # End dictionary.

    # ----- Prefilter summary -----
    prefilter_summary = {                                                                                                       # Compile prefilter metrics.
        "rows": int(len(prefilter_df)),                                                                                         # Count prefilter rows.
        "top_sample_ids": (                                                                                                     # List top sample IDs.
            [int(x) for x in prefilter_df["sample_id"].tolist()]                                                                # Extract IDs.
            if (len(prefilter_df) and "sample_id" in prefilter_df.columns) else []                                              # Safety check.
        ),                                                                                                                      # Close list.
    }                                                                                                                           # End dictionary.

    # ----- Stage 11 + baseline validation summaries -----
    stage11_validation = summarize_validation(validation_df)                                                                    # Summarize this run.
    stage08_validation = summarize_validation(baseline_08)                                                                      # Summarize S08 baseline.
    stage09_validation = summarize_validation(baseline_09)                                                                      # Summarize S09 baseline.
    stage10_validation = summarize_validation(baseline_10)                                                                      # Summarize S10 baseline.

    # ----- Persist the machine-readable summary -----
    out_dir = Path(args.out_dir)                                                                                                # Ensure path object.
    out_dir.mkdir(parents=True, exist_ok=True)                                                                                  # Create directory.
    summary = {                                                                                                                 # Construct main summary object.
        "stage": "11e",                                                                                                         # Set stage version.
        "run_metadata": {                                                                                                       # Add metadata section.
            "seed_protein_id": seed_protein_id,                                                                                 # Insert seed ID.
            "seed_length": seed_length,                                                                                         # Insert sequence length.
            "source_host": source_host,                                                                                         # Insert source host.
            "target_host": target_host,                                                                                         # Insert target host.
            "family_member_count": family_member_count,                                                                         # Insert family count.
            "target_member_count": target_member_count,                                                                         # Insert target count.
            "hard_positions": [int(p) for p in hard_positions],                                                                 # Insert hard positions.
            "soft_positions": [int(p) for p in soft_positions],                                                                 # Insert soft positions.
            "mutation_budget": [min_mutations, max_mutations],                                                                  # Insert mutation budget.
        },                                                                                                                      # Close metadata.
        "baseline_qualification": baseline_qualification,                                                                       # Add qualification results.
        "search_summary": search_summary,                                                                                       # Add search summary.
        "prefilter_summary": prefilter_summary,                                                                                 # Add prefilter summary.
        "stage11_validation": stage11_validation,                                                                               # Add current validation.
        "comparative": {                                                                                                        # Add comparative metrics.
            "stage08": stage08_validation,                                                                                      # Add S08 stats.
            "stage09": stage09_validation,                                                                                      # Add S09 stats.
            "stage10": stage10_validation,                                                                                      # Add S10 stats.
        },                                                                                                                      # Close comparative section.
        "verdict": _stage11_verdict(stage11_validation),                                                                        # Generate and add verdict.
    }                                                                                                                           # End dictionary.
    summary_path = out_dir / "stage11_report_summary.json"                                                                      # Define path.
    write_json(summary, summary_path)                                                                                           # Write JSON report.
    print(f"[OK] Wrote: {summary_path}", flush=True)                                                                            # Notify user.

    # ----- Build the human-readable markdown report -----
    md = _render_markdown_report(                                                                                               # Generate markdown string.
        seed_protein_id=seed_protein_id, seed_length=seed_length,                                                               # Pass metadata args.
        source_host=source_host, target_host=target_host,                                                                       # Pass host args.
        family_member_count=family_member_count, target_member_count=target_member_count,                                       # Pass counts.
        hard_positions=hard_positions, soft_positions=soft_positions,                                                           # Pass positions.
        min_mutations=min_mutations, max_mutations=max_mutations,                                                               # Pass mutation constraints.
        baseline_qualification=baseline_qualification,                                                                          # Pass qualification.
        search_summary=search_summary, prefilter_summary=prefilter_summary,                                                     # Pass summaries.
        prefilter_df=prefilter_df, validation_df=validation_df,                                                                 # Pass dataframes.
        stage11_validation=stage11_validation,                                                                                  # Pass current val.
        stage08_validation=stage08_validation,                                                                                  # Pass S08 val.
        stage09_validation=stage09_validation,                                                                                  # Pass S09 val.
        stage10_validation=stage10_validation,                                                                                  # Pass S10 val.
        verdict=summary["verdict"],                                                                                             # Pass verdict.
    )                                                                                                                           # End function call.
    report_path = out_dir / "stage11_report.md"                                                                                 # Define path.
    report_path.write_text(md, encoding="utf-8")                                                                                # Write MD report.
    print(f"[OK] Wrote: {report_path}", flush=True)                                                                             # Notify user.

    # ----- Top-3 handoff FASTA -----
    handoff_path = out_dir / "stage11_handoff.fasta"                                                                            # Define path.
    handoff_records = _build_handoff_fasta_records(                                                                             # Generate fasta records.
        prefilter_df=prefilter_df, validation_df=validation_df,                                                                 # Pass DFs.
        target_host=target_host, seed_protein_id=seed_protein_id,                                                               # Pass metadata.
    )                                                                                                                           # End function call.
    write_fasta(handoff_records, handoff_path)                                                                                  # Write FASTA.
    print(f"[OK] Wrote: {handoff_path}", flush=True)                                                                            # Notify user.

    sys.exit(EXIT_OK)                                                                                                           # Exit cleanly.


def _stage11_verdict(stage11_validation: dict) -> dict:
    """Evaluates the success of the run based on validation metrics and returns a descriptive status dictionary.

    Example:
        outcome = _stage11_verdict({'pass_count': 1, 'rows': 10})
        Returns: {'category': 'scientifically_validated', ...}
    """
    rows = int(stage11_validation.get("rows", 0))                                                                               # Get validated count.
    pass_count = int(stage11_validation.get("pass_count", 0))                                                                   # Get pass count.
    if rows == 0:                                                                                                               # Handle no validation data.
        return {                                                                                                                # Build error structure.
            "category": "no_validation_data",                                                                                   # Set category.
            "description": "No Stage 11 validation CSV was provided to the reporter.",                                          # Describe issue.
            "pass_rate": None,                                                                                                  # No rate.
        }                                                                                                                       # End dictionary.
    pass_rate = pass_count / max(rows, 1)                                                                                       # Calculate rate.
    if pass_count >= 1:                                                                                                         # Check for pass.
        return {                                                                                                                # Build success structure.
            "category": "scientifically_validated",                                                                             # Set category.
            "description": (                                                                                                    # Formulate description.
                f"{pass_count}/{rows} Stage 11 candidates passed the six 08a hard "                                             # Count stats.
                "thresholds. The Baseline Qualification hypothesis is supported: "                                              # Theoretical success.
                "starting from a wild-type foldable seed produces structurally "                                                # Mechanism.
                "viable inverse-folded redesigns."                                                                              # Outcome.
            ),                                                                                                                  # Close description.
            "pass_rate": pass_rate,                                                                                             # Set rate.
        }                                                                                                                       # End dictionary.
    return {                                                                                                                    # Build failure structure.
        "category": "falsified_or_partial",                                                                                     # Set category.
        "description": (                                                                                                        # Formulate description.
            f"0/{rows} Stage 11 candidates passed the six 08a hard thresholds. "                                                # Count stats.
            "The Baseline Qualification hypothesis is not supported by this run. "                                              # Theoretical fail.
            "Consider: (a) trying a different seed protein_id, (b) widening the "                                               # Mitigation strategies.
            "mutation budget, (c) upgrading the structural oracle to AlphaFold-Multimer."                                        # Further ideas.
        ),                                                                                                                      # Close description.
        "pass_rate": pass_rate,                                                                                                 # Set rate.
    }                                                                                                                           # End dictionary.


def _render_markdown_report(
    *,
    seed_protein_id: str, seed_length: int,
    source_host: str, target_host: str,
    family_member_count: int, target_member_count: int,
    hard_positions: list, soft_positions: list,
    min_mutations: int, max_mutations: int,
    baseline_qualification: dict,
    search_summary: dict, prefilter_summary: dict,
    prefilter_df: pd.DataFrame, validation_df: pd.DataFrame,
    stage11_validation: dict,
    stage08_validation: dict,
    stage09_validation: dict,
    stage10_validation: dict,
    verdict: dict,
) -> str:
    """Assembles the final report contents into a structured markdown document for review.

    Example:
        md_text = _render_markdown_report(seed_protein_id="P1", ...)
        # Returns a string containing the complete markdown report.
    """
    lines: list[str] = []                                                                                                       # Initialize line buffer.
    lines.append("# Stage 11 — Universal Structure-Conditioned RBP Redesign Report")                                            # Add report title.
    lines.append("")                                                                                                            # Blank line.
    lines.append("## Why Stage 11 exists")                                                                                      # Add section header.
    lines.append("")                                                                                                            # Blank line.
    lines.append(                                                                                                               # Write project intro.
        "Stages 07–10 all failed ESMFold validation with mean pLDDT ≈ 0.22 across the panels. "                                 # Describe problem.
        "The Stage 10 post-mortem traced this to a corrupted Stage 06 seed chassis "                                           # Mention root cause.
        "(`round8_cand473` itself folds at pLDDT 0.22). Stage 11 mandates a Baseline "                                         # Mention chassis.
        "Qualification gate (default pLDDT ≥ 70.0 on the wild-type seed) before any "                                           # Mention fix.
        "redesign compute is spent, and recomputes every family/target/edit-space "                                             # Mention methodology.
        "artifact in-stage so no failed downstream artifact can re-enter the loop."                                             # Conclusion.
    )                                                                                                                           # End paragraph.
    lines.append("")                                                                                                            # Blank line.

    # Headline outcome
    lines.append("## Headline outcome")                                                                                         # Add section header.
    lines.append("")                                                                                                            # Blank line.
    lines.append(f"- **Verdict:** `{verdict['category']}`")                                                                     # Print category.
    lines.append(f"- **Interpretation:** {verdict['description']}")                                                             # Print description.
    if verdict.get("pass_rate") is not None:                                                                                    # Add pass rate if available.
        lines.append(f"- **Pass rate:** {verdict['pass_rate']:.0%}")                                                            # Format as percent.
    lines.append(f"- **Best Stage 11 mean pLDDT:** {_fmt(stage11_validation.get('best_mean_plddt'))}")                          # Print max pLDDT.
    lines.append(f"- **Best Stage 11 RMSD to seed:** {_fmt(stage11_validation.get('best_rmsd'))} Å")                            # Print min RMSD.
    lines.append("")                                                                                                            # Blank line.

    # Baseline Qualification
    lines.append("## Baseline Qualification")                                                                                   # Add section header.
    lines.append("")                                                                                                            # Blank line.
    bq_passed = bool(baseline_qualification.get("passed", False))                                                               # Check pass status.
    bq_actual = float(baseline_qualification.get("seed_mean_plddt", float("nan")))                                              # Get seed pLDDT.
    bq_threshold = float(baseline_qualification.get("min_plddt_threshold", float("nan")))                                       # Get threshold.
    lines.append(f"- **Gate result:** {'PASS' if bq_passed else 'FAIL'}")                                                       # Output gate status.
    lines.append(f"- **Seed mean pLDDT:** {bq_actual:.3f}" if math.isfinite(bq_actual) else "- **Seed mean pLDDT:** —")         # Format/display pLDDT.
    lines.append(f"- **Threshold:** {bq_threshold:.3f}" if math.isfinite(bq_threshold) else "- **Threshold:** —")               # Format/display threshold.
    lines.append(f"- **Reason:** {baseline_qualification.get('reason', '—')}")                                                  # Output reason.
    lines.append(                                                                                                               # Output comparison note.
        "- **For comparison:** the historical Stage 06 chassis (`round8_cand473`) folded "                                        # Describe historical context.
        "at pLDDT ≈ 0.22, which is the structural defect Stage 11 is built to refuse."                                          # Link context to current issue.
    )                                                                                                                           # End paragraph.
    lines.append("")                                                                                                            # Blank line.

    # Core setup
    lines.append("## Core redesign setup")                                                                                      # Add section header.
    lines.append("")                                                                                                            # Blank line.
    lines.append(f"- **Source host:** {source_host}")                                                                           # List source host.
    lines.append(f"- **Target host:** {target_host}")                                                                           # List target host.
    lines.append(f"- **Seed protein_id:** `{seed_protein_id}`")                                                                 # List seed ID.
    lines.append(f"- **Seed length:** {seed_length} AA")                                                                        # List length.
    lines.append(f"- **Family members used:** {family_member_count}")                                                           # List family count.
    lines.append(f"- **Target members used:** {target_member_count}")                                                           # List target count.
    lines.append(f"- **Hard editable positions:** {hard_positions}")                                                            # List hard positions.
    lines.append(f"- **Soft editable positions:** {soft_positions}")                                                            # List soft positions.
    lines.append(f"- **Mutation budget:** {min_mutations} … {max_mutations}")                                                   # List budget.
    lines.append("")                                                                                                            # Blank line.

    # Search statistics
    lines.append("## Search statistics")                                                                                        # Add section header.
    lines.append("")                                                                                                            # Blank line.
    lines.append(f"- **Unique candidates evaluated:** {search_summary.get('rows', 0)}")                                         # Total candidates.
    lines.append(f"- **Rounds completed:** {search_summary.get('rounds_completed', 0)}")                                        # Search depth.
    lines.append(f"- **Best composite score:** {_fmt(search_summary.get('best_composite_score'))}")                             # Max score.
    lines.append(f"- **Best target probability:** {_fmt(search_summary.get('best_target_probability'))}")                       # Max probability.
    lines.append(f"- **Best inverse-folding log-likelihood:** {_fmt(search_summary.get('best_if1_log_likelihood'))}")           # Max likelihood.
    lines.append(f"- **Best seed identity:** {_fmt(search_summary.get('best_seed_identity'))}")                                 # Max identity.
    hist = search_summary.get("mutation_count_histogram", {}) or {}                                                             # Get histogram data.
    if hist:                                                                                                                    # Only write if populated.
        lines.append("- **Mutation count histogram:** "                                                                         # Histogram heading.
                     + ", ".join(f"{k}→{hist[k]}" for k in sorted(hist)))                                                       # Flattened string.
    lines.append("")                                                                                                            # Blank line.

    # Top panel table
    lines.append("## Top-K compact panel")                                                                                      # Add section header.
    lines.append("")                                                                                                            # Blank line.
    if prefilter_df.empty:                                                                                                      # Check if prefilter empty.
        lines.append("Prefilter CSV not provided; panel table omitted.")                                                        # Error/Skip note.
    else:                                                                                                                       # Build table logic.
        cols_to_show = [                                                                                                        # Select columns.
            c for c in [                                                                                                        # Filter available columns.
                "sample_id", "mutation_count", "mutation_text",                                                                 # Standard columns.
                "target_probability", "if1_log_likelihood",                                                                     # Score columns.
                "family_cosine", "seed_identity", "stage10_composite_score",                                                    # Context columns.
            ] if c in prefilter_df.columns                                                                                      # Existence check.
        ]                                                                                                                       # End list.
        header = "| " + " | ".join(cols_to_show) + " |"                                                                         # Build Markdown header.
        sep = "|" + "|".join("---:" if c != "mutation_text" else "---" for c in cols_to_show) + "|"                             # Build Markdown separator.
        lines.append(header)                                                                                                    # Add header row.
        lines.append(sep)                                                                                                       # Add separator row.
        for _, row in prefilter_df.iterrows():                                                                                  # Iterate data.
            cells = []                                                                                                          # Clear row cells.
            for c in cols_to_show:                                                                                              # Iterate selected columns.
                val = row[c]                                                                                                    # Get value.
                if isinstance(val, float):                                                                                      # Check float.
                    cells.append(_fmt(val))                                                                                     # Format numbers.
                else:                                                                                                           # Handle string/int.
                    cells.append(str(val))                                                                                      # Cast to string.
            lines.append("| " + " | ".join(cells) + " |")                                                                       # Add row to report.
    lines.append("")                                                                                                            # Blank line.

    # Structural validation table
    lines.append("## Stage 11 structural validation")                                                                           # Add section header.
    lines.append("")                                                                                                            # Blank line.
    if validation_df.empty:                                                                                                     # Check validation availability.
        lines.append("Validation CSV was not provided; structural metrics are unavailable for this run.")                       # Error/Skip note.
    else:                                                                                                                       # Add metrics.
        lines.append(f"- **Validated rows:** {stage11_validation.get('rows', 0)}")                                              # Row count.
        lines.append(f"- **Pass count:** {stage11_validation.get('pass_count', 0)} "                                            # Pass info.
                     f"({verdict.get('pass_rate', 0):.0%} pass rate)" if verdict.get('pass_rate') is not None else "")          # Rate if applicable.
        lines.append(f"- **Best mean pLDDT:** {_fmt(stage11_validation.get('best_mean_plddt'))}")                               # Max pLDDT.
        lines.append(f"- **Mean mean pLDDT:** {_fmt(stage11_validation.get('mean_mean_plddt'))}")                               # Mean pLDDT.
        lines.append(f"- **Best mutation-site mean pLDDT:** {_fmt(stage11_validation.get('best_mutation_site_mean_plddt'))}")     # Max mut-site.
        lines.append(f"- **Best RMSD to seed:** {_fmt(stage11_validation.get('best_rmsd'))} Å")                                 # Min RMSD.
        lines.append(f"- **Mean RMSD to seed:** {_fmt(stage11_validation.get('mean_rmsd'))} Å")                                 # Mean RMSD.
    lines.append("")                                                                                                            # Blank line.

    # Comparative table
    lines.append("## Comparative table — Stage 08 / 09 / 10 / 11")                                                              # Add section header.
    lines.append("")                                                                                                            # Blank line.
    lines.append("| Stage | Rows | Pass | Best pLDDT | Mean pLDDT | Best RMSD (Å) | Mean RMSD (Å) | Best mut-site pLDDT |") # Create header.
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")                                                                    # Create divider.
    for tag, summ in (                                                                                                          # Iterate through stages.
        ("Stage 08", stage08_validation),                                                                                       # Add S08 row.
        ("Stage 09", stage09_validation),                                                                                       # Add S09 row.
        ("Stage 10", stage10_validation),                                                                                       # Add S10 row.
        ("Stage 11", stage11_validation),                                                                                       # Add S11 row.
    ):                                                                                                                          # Loop.
        lines.append(                                                                                                           # Generate table row.
            f"| {tag} | {summ.get('rows', 0)} | {summ.get('pass_count', 0)} | "                                                # Standard stats.
            f"{_fmt(summ.get('best_mean_plddt'))} | {_fmt(summ.get('mean_mean_plddt'))} | "                                      # pLDDT stats.
            f"{_fmt(summ.get('best_rmsd'))} | {_fmt(summ.get('mean_rmsd'))} | "                                                 # RMSD stats.
            f"{_fmt(summ.get('best_mutation_site_mean_plddt'))} |"                                                              # Mut-site stats.
        )                                                                                                                       # Append to lines.
    lines.append("")                                                                                                            # Blank line.

    # Conclusion
    lines.append("## Conclusion")                                                                                               # Add section header.
    lines.append("")                                                                                                            # Blank line.
    lines.append(verdict.get("description", ""))                                                                                # Write verdict text.
    lines.append("")                                                                                                            # Blank line.
    lines.append(                                                                                                               # Write conclusion note.
        "Stage 11's core architectural commitment is that the wild-type seed must be "                                          # Restate methodology.
        "structurally healthy before any redesign compute is spent. This run's "                                                # Justify why.
        "comparative table is the first chance in the project to isolate the seed-chassis "                                     # Explain novelty.
        "hypothesis from the inverse-folding engine — if Stage 11 lifts the pass rate "                                         # Mention success condition.
        "above zero while Stages 08–10 stayed at 0/N, the chassis hypothesis is supported."                                     # Closing logic.
    )                                                                                                                           # End paragraph.
    return "\n".join(lines)                                                                                                     # Join all lines.


def _build_handoff_fasta_records(
    *,
    prefilter_df: pd.DataFrame,
    validation_df: pd.DataFrame,
    target_host: str,
    seed_protein_id: str,
) -> list[tuple[str, str]]:
    """Generates a list of FASTA formatted (header, sequence) tuples for top candidates.

    Example:
        records = _build_handoff_fasta_records(df_pre, df_val, "host", "seed_id")
        # Returns: [(">header1", "SEQ1"), (">header2", "SEQ2"), ...]
    """
    if prefilter_df.empty:                                                                                                      # Handle empty DF.
        return []                                                                                                               # Return empty record list.
    sub = prefilter_df.copy()                                                                                                   # Create working DF.
    sub = sub.head(3)                                                                                                           # Limit to top 3.

    # Optionally merge validation metrics if the CSV has them.
    if not validation_df.empty and "sample_id" in validation_df.columns:                                                        # Check if validation exists.
        merge_cols = [c for c in ["esmfold_mean_plddt", "rmsd_to_selected_seed",                                                # Define merge columns.
                                 "mutation_site_mean_plddt", "stage08_pass"]                                                    # List potential columns.
                      if c in validation_df.columns]                                                                            # Verify columns exist.
        if merge_cols:                                                                                                          # Check if merge valid.
            sub = sub.merge(validation_df[["sample_id"] + merge_cols],                                                          # Perform merge.
                            on="sample_id", how="left")                                                                         # Left join on sample_id.

    records: list[tuple[str, str]] = []                                                                                         # Result list.
    for _, row in sub.iterrows():                                                                                               # Iterate rows.
        score = float(row.get("stage10_composite_score", row.get("stage11_composite_score", float("nan"))))                     # Get composite score.
        plddt = float(row.get("esmfold_mean_plddt", float("nan")))                                                              # Get pLDDT if exists.
        rmsd = float(row.get("rmsd_to_selected_seed", float("nan")))                                                            # Get RMSD if exists.
        passed = bool(row.get("stage08_pass", False)) if "stage08_pass" in row.index else None                                  # Get pass if exists.
        # Build the header (Stage 07g-style metadata-rich FASTA header).
        parts = [                                                                                                               # List metadata parts.
            f"stage11_sample{int(row['sample_id'])}",                                                                           # Identifier.
            f"target={target_host}",                                                                                            # Target info.
            f"seed={seed_protein_id}",                                                                                          # Seed info.
            f"score={score:.6f}" if math.isfinite(score) else "score=na",                                                       # Score value.
            f"plddt={plddt:.2f}" if math.isfinite(plddt) else "plddt=na",                                                       # pLDDT value.
            f"rmsd={rmsd:.3f}" if math.isfinite(rmsd) else "rmsd=na",                                                           # RMSD value.
            f"mutations={row.get('mutation_text', row.get('mutation_positions', ''))}",                                         # Mutations string.
        ]                                                                                                                       # End list.
        if passed is not None:                                                                                                  # Add pass info.
            parts.append(f"stage08_pass={passed}")                                                                              # Append to metadata.
        header = "|".join(parts)                                                                                                # Join header fields.
        sequence = str(row["candidate_sequence"])                                                                               # Get sequence string.
        records.append((header, sequence))                                                                                      # Add record.
    return records                                                                                                              # Return records list.


if __name__ == "__main__":                                                                                                      # Guard for CLI.
    main()                                                                                                                      # Run main.