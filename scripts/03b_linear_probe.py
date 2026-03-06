from __future__ import annotations
import argparse
import json
from pathlib import Path
import torch
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, f1_score, accuracy_score
from sklearn.preprocessing import LabelEncoder


def split_by_virus(df, train_frac=0.8, val_frac=0.1, seed=42):
    """ Group split by virus_accession per class to reduce leakage:
        - for each host_genus: split its unique viruses into train/val/test
        - assign all proteins of a virus to the same split
        Returns a split array where each element is either "train", "val" or "test"."""
    
    rng = np.random.default_rng(seed)
    split = np.full(len(df), "train", dtype=object)

    for genus, g in df.groupby("host_genus"):
        viruses = g["virus_accession"].unique().tolist()
        rng.shuffle(viruses)

        n = len(viruses)
        n_train = int(round(train_frac * n))
        n_val = int(round(val_frac * n))

        train_set = set(viruses[:n_train])
        val_set = set(viruses[n_train:n_train+n_val])
        test_set = set(viruses[n_train+n_val:])

        idx = g.index
        v = df.loc[idx, "virus_accession"]

        split[idx[v.isin(train_set)]] = "train"
        split[idx[v.isin(val_set)]] = "val"
        split[idx[v.isin(test_set)]] = "test"

    return split  # Return a split array where each element is either "train", "val" or "test" based on the spitting ratio


def save_split_tables(df: pd.DataFrame, out_dir: Path):
    """ Save virus-level and protein-level split summaries as CSV files. """
    virus_counts = (
        df.groupby(["host_genus", "split"])["virus_accession"]
        .nunique()
        .unstack(fill_value=0)
        .reset_index()
    )
    virus_counts.to_csv(out_dir / "split_counts_viruses.csv", index=False)

    protein_counts = (
        df.groupby(["host_genus", "split"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    protein_counts.to_csv(out_dir / "split_counts_proteins.csv", index=False)


def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("--emb", type=str)  # Embeddings
    ap.add_argument("--idx", type=str)  # Index dataframe
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--dataset_name", type=str, default="broad")
    ap.add_argument("--out_dir", type=str, required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    emb = torch.load(args.emb, map_location="cpu").numpy()
    df = pd.read_csv(args.idx)

    le = LabelEncoder()
    y = le.fit_transform(df["host_genus"])  # Scans host_genus column, finds all unique strings, sorts them alphabetically, assigns an integer to each and replaces every string in the dataframe with its corresponding integer.

    df["split"] = split_by_virus(df, seed=args.seed)        # Add a new column "split" to the dataframe

    print("\nViruses per genus per split:")
    print(
        df.groupby(["host_genus", "split"])["virus_accession"]
        .nunique()
        .unstack(fill_value=0)
    )
    print("\nProteins per genus per split:")
    print(df["host_genus"].value_counts())
    print(df[df["split"]=="test"]["host_genus"].value_counts())

    save_split_tables(df, out_dir)

    X_train = emb[df["split"] == "train"]   # Get the embeddings for the train set
    y_train = y[df["split"] == "train"]

    X_test = emb[df["split"] == "test"]
    y_test = y[df["split"] == "test"]

    # Logistic Regression
    clf = LogisticRegression(
        max_iter=5000,
        solver="lbfgs"              # Use the L-BFGS algorithm
    )

    clf.fit(X_train, y_train)       # Fit the model

    preds = clf.predict(X_test)     # Predict

    report = classification_report(
        y_test,
        preds,
        target_names=le.classes_,
        zero_division=0
    )
    (out_dir / "test_report.txt").write_text(report)

    test_macro_f1 = f1_score(y_test, preds, average="macro", zero_division=0)
    test_weighted_f1 = f1_score(y_test, preds, average="weighted", zero_division=0)
    test_acc = accuracy_score(y_test, preds)

    metrics = {
        "dataset_name": args.dataset_name,
        "model_name": "linear_probe",
        "seed": args.seed,
        "n_total": int(len(df)),
        "n_train": int((df["split"] == "train").sum()),
        "n_val": int((df["split"] == "val").sum()),
        "n_test": int((df["split"] == "test").sum()),
        "accuracy": float(test_acc),
        "macro_f1": float(test_macro_f1),
        "weighted_f1": float(test_weighted_f1),
        "classes": le.classes_.tolist(),
        "emb_path": args.emb,
        "idx_path": args.idx,
    }
    (out_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))

    predictions_df = pd.DataFrame(
        {
            "true_label": [le.classes_[i] for i in y_test],
            "pred_label": [le.classes_[i] for i in preds],
        }
    )
    predictions_df.to_csv(out_dir / "test_predictions.csv", index=False)

    print(report)
    print(f"✅ Saved: {out_dir / 'test_report.txt'}")
    print(f"✅ Saved: {out_dir / 'metrics.json'}")
    print(f"✅ Saved: {out_dir / 'split_counts_viruses.csv'}")
    print(f"✅ Saved: {out_dir / 'split_counts_proteins.csv'}")
    print(f"✅ Saved: {out_dir / 'test_predictions.csv'}")


if __name__ == "__main__":
    main()