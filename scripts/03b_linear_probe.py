from __future__ import annotations

import argparse
import torch
import pandas as pd
import numpy as np

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.preprocessing import LabelEncoder


def split_by_virus(df, train_frac=0.8, val_frac=0.1, seed=42):
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

    return split


def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("--emb", type=str)
    ap.add_argument("--idx", type=str)
    args = ap.parse_args()

    emb = torch.load(args.emb, map_location="cpu").numpy()
    df = pd.read_csv(args.idx)

    le = LabelEncoder()
    y = le.fit_transform(df["host_genus"])

    df["split"] = split_by_virus(df)

    X_train = emb[df["split"] == "train"]
    y_train = y[df["split"] == "train"]

    X_test = emb[df["split"] == "test"]
    y_test = y[df["split"] == "test"]

    clf = LogisticRegression(
        max_iter=5000,
        multi_class="multinomial",
        solver="lbfgs"
    )

    clf.fit(X_train, y_train)

    preds = clf.predict(X_test)

    print(classification_report(
        y_test,
        preds,
        target_names=le.classes_,
        zero_division=0
    ))


if __name__ == "__main__":
    main()