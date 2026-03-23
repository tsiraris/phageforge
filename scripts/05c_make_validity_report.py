# This script generates the final presentation artifacts (plots and markdown documents) to summarize the validation phase.

from __future__ import annotations                                                      
import argparse                                                                         
from pathlib import Path                                                                
import matplotlib.pyplot as plt                                                         
import pandas as pd                                                                     


def make_score_vs_manifold_plot(master_df: pd.DataFrame, shortlist_df: pd.DataFrame, out_path: Path) -> None:  
    """
    This helper creates a scatter plot of target score versus strict manifold score.
    """
    plt.figure(figsize=(8, 6))                                                          # Create a fresh figure with a readable size.
    plt.scatter(master_df["target_score"], master_df["strict_manifold_score"], alpha=0.25, s=18)  # Plot all candidates as a light background cloud.
    plt.scatter(shortlist_df["target_score"], shortlist_df["strict_manifold_score"], alpha=0.95, s=40)  # Overlay shortlist candidates as the highlighted set.
    plt.xlabel("Target host score")                                                     # Label the x-axis.
    plt.ylabel("Strict manifold score")                                                 # Label the y-axis.
    plt.title("Target score vs. strict RBP-manifold score")                             # Add a descriptive figure title.
    plt.tight_layout()                                                                  # Reduce extra whitespace and prevent label clipping.
    plt.savefig(out_path, dpi=200)                                                      # Save the figure at publication-friendly resolution.
    plt.close()                                                                         # Close the figure to free memory.


def make_score_vs_pll_plot(master_df: pd.DataFrame, shortlist_df: pd.DataFrame, out_path: Path) -> None:  
    """
    This helper creates a scatter plot of target score versus local PLL delta.
    """
    plt.figure(figsize=(8, 6))                                                          # Create a fresh figure with a readable size.
    plt.scatter(master_df["target_score"], master_df["pll_delta_candidate_minus_seed"], alpha=0.25, s=18)  # Plot all candidates as a light background cloud.
    plt.scatter(shortlist_df["target_score"], shortlist_df["pll_delta_candidate_minus_seed"], alpha=0.95, s=40)  # Overlay shortlist candidates as the highlighted set.
    plt.xlabel("Target host score")                                                     # Label the x-axis.
    plt.ylabel("Mutation-local PLL delta (candidate - seed)")                           # Label the y-axis.
    plt.title("Target score vs. mutation-local naturalness")                            # Add a descriptive figure title.
    plt.tight_layout()                                                                  # Reduce extra whitespace and prevent label clipping.
    plt.savefig(out_path, dpi=200)                                                      # Save the figure at publication-friendly resolution.
    plt.close()                                                                         # Close the figure to free memory.


def write_markdown_report(master_df: pd.DataFrame, unique_df: pd.DataFrame, shortlist_df: pd.DataFrame, run_summary_df: pd.DataFrame, out_path: Path) -> None:  # Accept the deduplicated table so the report can distinguish rows from unique sequences.
    """
    This helper writes a compact markdown report using the computed tables.
    """
    n_total = len(master_df)                                                            # Count all candidate rows in the master table.
    n_validated = int(master_df["validated_pass"].sum())                                 # Count rows that pass the fast validity filters.
    validated_fraction = (n_validated / n_total) if n_total > 0 else 0.0                # Compute the validated fraction of the candidate pool.
    n_unique = len(unique_df)                                                           # Count the number of unique target-host sequence entries after deduplication.
    n_validated_unique = int(unique_df["validated_pass"].sum()) if "validated_pass" in unique_df.columns else 0  # Count validated unique sequences after deduplication.
    best_raw = master_df.sort_values("target_score", ascending=False).head(10)          # Extract the top raw candidates by target score.
    best_validated = shortlist_df.sort_values(["target_host", "rank_score"], ascending=[True, False]).head(20)  # Extract the strongest shortlisted candidates.

    lines = []                                                                          # Create a list that will accumulate markdown lines.
    lines.append("# PhageForge fast-close validity report")                             # Add the report title.
    lines.append("")                                                                    # Add a blank line for markdown readability.
    lines.append("## Executive summary")                                                # Add the executive-summary heading.
    lines.append("")                                                                    # Add a blank line for readability.
    lines.append(f"- Total candidate rows re-evaluated: **{n_total:,}**")                # Report how many candidate rows were re-evaluated.
    lines.append(f"- Candidates passing fast validity filters: **{n_validated:,}** ({validated_fraction:.1%})")  # Report the validated fraction.
    lines.append(f"- Total unique target-host sequences after deduplication: **{n_unique:,}**")  # Report how many unique target-host sequences remain after deduplication.
    lines.append(f"- Unique sequences passing fast validity filters: **{n_validated_unique:,}**")  # Report the validated count after deduplication.
    lines.append(f"- Final shortlist size: **{len(shortlist_df):,}**")                   # Report the final shortlist size.
    lines.append(f"- Unique sequences in final shortlist: **{shortlist_df['aa_sequence'].nunique():,}**")  # Make it explicit when the shortlist contains repeated versus unique sequences.
    lines.append("")                                                                    # Add a blank line before the interpretation text.
    lines.append("The fast-close analysis reframes the project from raw host-score optimization to **validity-aware RBP design**. The key result is that candidate selection should not rely on host score alone; instead, candidates are retained only when they also remain close to the strict RBP manifold, preserve an RBP-like nearest-neighbor family annotation, and avoid large local naturalness drops around mutated residues.")  # Add the main conclusion in plain language.
    lines.append("")                                                                    # Add a blank line before the run-summary section.
    lines.append("## Run-level summary")                                                # Add the run-summary heading.
    lines.append("")                                                                    # Add a blank line for readability.
    lines.append(run_summary_df.to_markdown(index=False))                                # Insert the run-level summary table directly into the markdown report.
    lines.append("")                                                                    # Add a blank line before the shortlist section.
    lines.append("## Final shortlist")                                                  # Add the shortlist heading.
    lines.append("")                                                                    # Add a blank line for readability.
    shortlist_cols = [                                                                   # Define the most important columns to display in the shortlist table.
        "target_host",                                                                   # Keep the target host so the shortlist is grouped meaningfully.
        "run_name",                                                                      # Keep the originating run name for traceability.
        "candidate_id",                                                                  # Keep the candidate ID for exact lookup.
        "actual_mutation_list" if "actual_mutation_list" in best_validated.columns else "mutations",  # Prefer the exact sequence-derived mutation list in the report.
        "target_score",                                                                  # Keep the target-host score.
        "strict_manifold_score",                                                         # Keep the strict manifold score.
        "pll_delta_candidate_minus_seed",                                                # Keep the local PLL delta.
        "strict_nn_product",                                                             # Keep the nearest strict annotation for family-retention context.
        "strict_nn_host",                                                                # Keep the nearest strict host for quick biological interpretation.
    ]                                                                                     # Finish the shortlist column selection.
    lines.append(best_validated[shortlist_cols].to_markdown(index=False))                # Insert the formatted shortlist table.
    lines.append("")                                                                    # Add a blank line before the raw top-score section.
    lines.append("## Highest raw target-score candidates")                              # Add the raw top-score heading.
    lines.append("")                                                                    # Add a blank line for readability.
    raw_cols = [                                                                          # Define the columns shown for the raw top-score view.
        "target_host",                                                                   # Keep the target host label.
        "run_name",                                                                      # Keep the run name for traceability.
        "candidate_id",                                                                  # Keep the candidate ID for exact lookup.
        "actual_mutation_list" if "actual_mutation_list" in best_raw.columns else "mutations",  # Prefer sequence-derived mutations in the raw top-score table too.
        "target_score",                                                                  # Keep the target-host score.
        "strict_manifold_score",                                                         # Keep the strict manifold score.
        "pll_delta_candidate_minus_seed",                                                # Keep the local PLL delta.
        "validated_pass",                                                                # Show whether the raw high-scoring candidate passed the fast filters.
    ]                                                                                     # Finish the raw-table column selection.
    lines.append(best_raw[raw_cols].to_markdown(index=False))                            # Insert the formatted raw top-score table.
    lines.append("")                                                                    # Add a blank line before the interpretation section.
    lines.append("## Interpretation")                                                   # Add the interpretation heading.
    lines.append("")                                                                    # Add a blank line for readability.
    lines.append("- **What worked:** the original optimization runs do produce candidates with improved target-host score.")  # State the first positive finding.
    lines.append("- **What needed fixing:** raw target score alone is not enough, because some high-scoring candidates drift toward lower manifold support or worse local mutation plausibility.")  # State the core limitation discovered by the analysis.
    lines.append("- **What closes the project credibly:** a final shortlist should be selected from the intersection of host score, RBP-manifold support, family retention, mutation-local naturalness, and **unique-sequence evidence after deduplication**.")  # Make the corrected deduplication lesson explicit in the final framing.
    lines.append("")                                                                    # Add a blank line before the next-step section.
    lines.append("## Immediate next step")                                              # Add the next-step heading.
    lines.append("")                                                                    # Add a blank line for readability.
    lines.append("Use `final_validated_shortlist.csv` as the final set for any optional downstream structural screen, wet-lab prioritization, or portfolio write-up. Only if this shortlist is too weak should you do a minimal constrained rerun with distal-region mutation masking and the same fast validity gates used here.")  # Give the user the practical next action.

    out_path.write_text("\n".join(lines), encoding="utf-8")                             # Write the accumulated markdown lines to disk.


def main() -> None:                                                                     
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Create final plots and a markdown report for the PhageForge fast-close analysis.")  
    parser.add_argument("--repo_root", type=str, default=".", help="Path to the phageforge repository root.")  
    parser.add_argument("--master_csv", type=str, default="results/analysis/validity/master_candidates_validated.csv", help="Master candidate table with validity metrics.")  
    parser.add_argument("--shortlist_csv", type=str, default="results/analysis/validity/final_validated_shortlist.csv", help="Final shortlist table produced by the ranking script.")  
    parser.add_argument("--run_summary_csv", type=str, default="results/analysis/validity/validity_summary_by_run.csv", help="Run-level summary table produced by step 1.")  
    parser.add_argument("--unique_csv", type=str, default="results/analysis/validity/master_candidates_unique_corrected.csv", help="Deduplicated target-host sequence table produced by step 1.")  # Read the deduplicated table from the merged 05 pipeline.
    parser.add_argument("--output_dir", type=str, default="results/analysis/validity", help="Directory where plots and the markdown report will be written.")  
    args = parser.parse_args()                                                          

    # Resolve paths
    repo_root = Path(args.repo_root).resolve()                                          # Resolve the repository root path.
    master_csv = (repo_root / args.master_csv).resolve()                                # Resolve the master-table path.
    shortlist_csv = (repo_root / args.shortlist_csv).resolve()                          # Resolve the shortlist-table path.
    run_summary_csv = (repo_root / args.run_summary_csv).resolve()                      # Resolve the run-summary-table path.
    unique_csv = (repo_root / args.unique_csv).resolve()                                # Resolve the deduplicated-table path.
    output_dir = (repo_root / args.output_dir).resolve()                                # Resolve the output directory path.
    output_dir.mkdir(parents=True, exist_ok=True)                                       # Create the output directory if needed.

    # Load dataframes
    master_df = pd.read_csv(master_csv)                                                 # Load the master candidate table.
    shortlist_df = pd.read_csv(shortlist_csv)                                           # Load the final shortlist table.
    run_summary_df = pd.read_csv(run_summary_csv)                                       # Load the run-level summary table.
    unique_df = pd.read_csv(unique_csv)                                                 # Load the deduplicated corrected table.

    # Create the validity report figures
    make_score_vs_manifold_plot(master_df=master_df, shortlist_df=shortlist_df, out_path=output_dir / "score_vs_manifold.png")  # Create the target-score vs manifold-score figure.
    make_score_vs_pll_plot(master_df=master_df, shortlist_df=shortlist_df, out_path=output_dir / "score_vs_naturalness.png")  # Create the target-score vs PLL-delta figure.
    write_markdown_report(master_df=master_df, unique_df=unique_df, shortlist_df=shortlist_df, run_summary_df=run_summary_df, out_path=output_dir / "final_report.md")  # Write the report using both row-level and deduplicated evidence.

    print("[INFO] Saved:")                                                               # Print a concise success summary.
    print(f"  - {output_dir / 'score_vs_manifold.png'}")                                 # Show the first figure path.
    print(f"  - {output_dir / 'score_vs_naturalness.png'}")                              # Show the second figure path.
    print(f"  - {output_dir / 'final_report.md'}")                                       # Show the markdown report path.


if __name__ == "__main__":                                                               # Run the CLI only when executed directly.
    main()                                                                                # Invoke the main function.
