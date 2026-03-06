from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--python_bin", type=str, default=sys.executable)
    ap.add_argument("--results_root", type=str, default="results")
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--batch_size", type=int, default=64)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--use_class_weights", action="store_true")
    args = ap.parse_args()

    results_root = Path(args.results_root)
    results_root.mkdir(parents=True, exist_ok=True)

    # Run benchmarks for both datasets (broad and strict)
    datasets = [
        {
            "dataset_name": "broad",
            "emb": "data/processed/esm2_embeddings.pt",
            "idx": "data/processed/esm2_embeddings_index.csv",
        },
        {
            "dataset_name": "strict",
            "emb": "data/processed/strict/esm2_embeddings.pt",
            "idx": "data/processed/strict/esm2_embeddings_index.csv",
        },
    ]

    seeds = [42, 43, 44]
    all_metrics = []

    # Run benchmarks for each dataset
    for ds in datasets:
        emb_path = Path(ds["emb"])
        idx_path = Path(ds["idx"])

        if not emb_path.exists() or not idx_path.exists():
            print(f"[WARN] Skipping dataset '{ds['dataset_name']}' because files are missing:")
            print(f"       emb: {emb_path}")
            print(f"       idx: {idx_path}")
            continue

        for seed in seeds:
            # Run MLP
            
            # Output directories
            mlp_out_dir = results_root / ds["dataset_name"] / "mlp" / f"seed_{seed}"
            mlp_out_dir.mkdir(parents=True, exist_ok=True)

            mlp_cmd = [
                args.python_bin,
                "scripts/03_train_phi_mlp.py",
                "--emb", ds["emb"],
                "--idx", ds["idx"],
                "--out_dir", str(mlp_out_dir),
                "--epochs", str(args.epochs),
                "--batch_size", str(args.batch_size),
                "--lr", str(args.lr),
                "--seed", str(seed),
                "--dataset_name", ds["dataset_name"],
            ]
            if args.use_class_weights:
                mlp_cmd.append("--use_class_weights")
            
            print("\n[RUNNING]", " ".join(mlp_cmd)) # Print the command that will be executed
            subprocess.run(mlp_cmd, check=True)     # Run the command

            mlp_metrics_path = mlp_out_dir / "metrics.json"                 # Get the path to the metrics file
            all_metrics.append(json.loads(mlp_metrics_path.read_text()))    # Append the metrics to the list

            # Run logistic regression probe
            lin_out_dir = results_root / ds["dataset_name"] / "linear_probe" / f"seed_{seed}"
            lin_out_dir.mkdir(parents=True, exist_ok=True)

            lin_cmd = [
                args.python_bin,
                "scripts/03b_linear_probe.py",
                "--emb", ds["emb"],
                "--idx", ds["idx"],
                "--seed", str(seed),
                "--dataset_name", ds["dataset_name"],
                "--out_dir", str(lin_out_dir),
            ]

            print("\n[RUNNING]", " ".join(lin_cmd))
            subprocess.run(lin_cmd, check=True)

            lin_metrics_path = lin_out_dir / "metrics.json"
            all_metrics.append(json.loads(lin_metrics_path.read_text()))
    
    # Check if any benchmark runs were completed
    if not all_metrics:
        print("[WARN] No benchmark runs were completed.")
        return

    # Save benchmark results to CSV
    summary_df = pd.DataFrame(all_metrics)
    summary_csv = results_root / "summary.csv"
    summary_df.to_csv(summary_csv, index=False)

    # Aggregate benchmark results
    aggregate_df = (
        summary_df.groupby(["dataset_name", "model_name"])
        .agg(
            accuracy_mean=("accuracy", "mean"),
            accuracy_std=("accuracy", "std"),
            macro_f1_mean=("macro_f1", "mean"),
            macro_f1_std=("macro_f1", "std"),
            weighted_f1_mean=("weighted_f1", "mean"),
            weighted_f1_std=("weighted_f1", "std"),
            n_total=("n_total", "mean"),
            n_train=("n_train", "mean"),
            n_val=("n_val", "mean"),
            n_test=("n_test", "mean"),
        )
        .reset_index()
    )

    # Save aggregate benchmark results
    aggregate_csv = results_root / "summary_aggregate.csv"
    aggregate_df.to_csv(aggregate_csv, index=False)

    # Save benchmark manifest (the parameters used for each benchmark run)
    manifest = {
        "results_root": str(results_root),
        "seeds": seeds,
        "datasets": datasets,
        "epochs": args.epochs,
        "batch_size": args.batch_size,
        "lr": args.lr,
        "use_class_weights": bool(args.use_class_weights),
    }
    (results_root / "benchmark_manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"\n✅ Saved: {summary_csv}")
    print(f"✅ Saved: {aggregate_csv}")
    print(f"✅ Saved: {results_root / 'benchmark_manifest.json'}")


if __name__ == "__main__":
    main()