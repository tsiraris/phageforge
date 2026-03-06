from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd


def fmt(mean_value: float, std_value: float) -> str:
    """ Format mean and std to 3 decimal places and return as string. """
    if pd.isna(std_value):
        return f"{mean_value:.3f}"
    return f"{mean_value:.3f} ± {std_value:.3f}"


def load_split_counts(results_root: Path, dataset_name: str, model_name: str, seed: int = 42):
    """ Load split counts for viruses and proteins and return as dataframes. """
    base = results_root / dataset_name / model_name / f"seed_{seed}"
    virus_csv = base / "split_counts_viruses.csv"
    protein_csv = base / "split_counts_proteins.csv"

    virus_df = pd.read_csv(virus_csv) if virus_csv.exists() else None
    protein_df = pd.read_csv(protein_csv) if protein_csv.exists() else None
    return virus_df, protein_df


def df_to_markdown(df: pd.DataFrame) -> str:
    """ Convert dataframe to markdown table. """
    return df.to_markdown(index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results_root", type=str, default="results")
    ap.add_argument("--output_path", type=str, default="results/benchmark_report.md")
    args = ap.parse_args()

    # Load summary and aggregate csv
    results_root = Path(args.results_root)
    summary_csv = results_root / "summary.csv"
    aggregate_csv = results_root / "summary_aggregate.csv"

    if not summary_csv.exists() or not aggregate_csv.exists():
        raise FileNotFoundError("Expected results/summary.csv and results/summary_aggregate.csv to exist. Run scripts/03c_run_benchmark.py first.")

    summary_df = pd.read_csv(summary_csv)
    aggregate_df = pd.read_csv(aggregate_csv)

    # Format summary and aggregate dataframes
    display_df = aggregate_df.copy()
    display_df["accuracy"] = display_df.apply(lambda r: fmt(r["accuracy_mean"], r["accuracy_std"]), axis=1)
    display_df["macro_f1"] = display_df.apply(lambda r: fmt(r["macro_f1_mean"], r["macro_f1_std"]), axis=1)
    display_df["weighted_f1"] = display_df.apply(lambda r: fmt(r["weighted_f1_mean"], r["weighted_f1_std"]), axis=1)

    # Keep only relevant columns: dataset_name, model_name, accuracy, macro_f1, weighted_f1
    display_df = display_df[
        ["dataset_name", "model_name", "accuracy", "macro_f1", "weighted_f1", "n_total", "n_train", "n_val", "n_test"]
    ]

    # Load split counts for viruses and proteins for broad and strict datasets
    broad_virus_df, broad_protein_df = load_split_counts(results_root, "broad", "linear_probe", seed=42)
    strict_virus_df, strict_protein_df = load_split_counts(results_root, "strict", "linear_probe", seed=42)

    # Get best run per dataset/model
    best_runs = (
        summary_df.sort_values(["macro_f1", "accuracy"], ascending=False)
        .groupby(["dataset_name", "model_name"], as_index=False)
        .first()
    )

    # Generate report
    lines = []
    lines.append("# PhageForge Benchmark Report")
    lines.append("")
    lines.append("This report compares:")
    lines.append("- broad vs strict dataset definitions")
    lines.append("- MLP vs linear probe baselines")
    lines.append("- 3 random seeds (42, 43, 44)")
    lines.append("")
    lines.append("## Aggregate benchmark table")
    lines.append("")
    lines.append(df_to_markdown(display_df))
    lines.append("")
    lines.append("## Best run per dataset/model")
    lines.append("")
    lines.append(df_to_markdown(best_runs[["dataset_name", "model_name", "seed", "accuracy", "macro_f1", "weighted_f1"]]))
    lines.append("")

    if broad_virus_df is not None:
        lines.append("## Broad dataset split counts (seed 42)")
        lines.append("")
        lines.append("### Virus counts")
        lines.append("")
        lines.append(df_to_markdown(broad_virus_df))
        lines.append("")
        lines.append("### Protein counts")
        lines.append("")
        lines.append(df_to_markdown(broad_protein_df))
        lines.append("")

    if strict_virus_df is not None:
        lines.append("## Strict dataset split counts (seed 42)")
        lines.append("")
        lines.append("### Virus counts")
        lines.append("")
        lines.append(df_to_markdown(strict_virus_df))
        lines.append("")
        lines.append("### Protein counts")
        lines.append("")
        lines.append(df_to_markdown(strict_protein_df))
        lines.append("")

    lines.append("## Interpretation")
    lines.append("")
    lines.append("- The broad dataset tests whether larger but noisier adsorption-module collections support host prediction.")
    lines.append("- The strict dataset tests whether high-confidence RBPs provide cleaner host-specific signal.")
    lines.append("- Comparing MLP vs linear probe helps distinguish representation quality from classifier instability.")
    lines.append("- Macro-F1 is the most important metric here because host genera are not equally represented.")
    lines.append("")

    output_path = Path(args.output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"✅ Saved: {output_path}")


if __name__ == "__main__":
    main()