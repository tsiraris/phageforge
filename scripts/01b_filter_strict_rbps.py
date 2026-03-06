from __future__ import annotations
import re
from pathlib import Path
import pandas as pd

KEEP_PATTERNS = [
    r"tail\s*fib(e|r)e",          # tail fiber / fibre
    r"tail\s*spike",              # tail spike
    r"tailspike",
    r"receptor\s*binding",
    r"receptor-binding",
    r"host\s*recognition",
    r"adsorption",
]

DROP_PATTERNS = [
    r"chaperone",
    r"assembly",
    r"wedge",
    r"hub",
    r"baseplate",
    r"lysozyme",
    r"hypothetical",
]

def main():
    in_csv = Path("data/processed/rbp_dataset_eskapee_stage1.csv")
    out_csv = Path("data/processed/rbp_dataset_eskapee_strict.csv")

    df = pd.read_csv(in_csv)
    prod = df["product"].fillna("").str.lower()

    keep = False
    for p in KEEP_PATTERNS:
        keep = keep | prod.str.contains(p, regex=True)

    drop = False
    for p in DROP_PATTERNS:
        drop = drop | prod.str.contains(p, regex=True)

    df2 = df[keep & ~drop].copy()

    # remove very short sequences just in case
    df2 = df2[df2["aa_sequence"].astype(str).str.len() >= 60]

    df2.to_csv(out_csv, index=False)
    print("Wrote:", out_csv, "rows:", len(df2))
    print(df2["host_genus"].value_counts())
    print(df2["product"].value_counts().head(20))

if __name__ == "__main__":
    main()