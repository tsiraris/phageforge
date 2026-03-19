from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd


# High-confidence keep patterns: proteins likely to be directly involved in host recognition / adsorption specificity
KEEP_PATTERNS = [
    r"tail\s*fib(?:e|r)e",              # tail fiber / fibre
    r"tail\s*spike",                    # tail spike
    r"tailspike",
    r"receptor\s*binding",
    r"receptor-binding",
    r"host\s*recognition",
    r"adsorption",
    r"fiber\s*protein",
    r"fibre\s*protein",
]

# Exclude proteins that are often not themselves the main host-specific determinant
DROP_PATTERNS = [
    r"chaperone",
    r"assembly",
    r"baseplate\s*hub",
    r"baseplate\s*wedge",
    r"\bhub\b",
    r"\bwedge\b",
    r"lysozyme",
    r"muramidase",
    r"hypothetical",
    r"tail\s*lysozyme",
    r"peptidase",
]

# Optional “review” patterns (when running low in data): proteins that may be relevant but are annotation-ambiguous
REVIEW_PATTERNS = [
    r"depolymerase",
    r"tail\s*tip",
    r"baseplate",
]


def build_mask(series: pd.Series, patterns: list[str]) -> pd.Series:
    """ 
    Takes a Pandas column and a list of regex patterns and returns a boolean mask indicating whether each row matched any of the patterns in the list.
    True for rows that matched at least one pattern, False otherwise. 
    """
    # Create a pandas column of False values with the same index as series (= the product column)
    mask = pd.Series(False, index=series.index)
    
    # Loop over every regex pattern in the list 
    for p in patterns:
        # If the row was already marked as True or the regex pattern is found in the series column, set the corresponding index in the mask to True
        mask = mask | series.str.contains(p, regex=True, na=False)
    return mask


def main():
    # Setup arguments for inputs, outputs, and length limits
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", type=str, default="data/processed/rbp_dataset_eskapee_strict.csv")
    ap.add_argument("--out_csv", type=str, default="data/processed/rbp_dataset_eskapee_structural.csv")
    ap.add_argument("--review_csv", type=str, default="data/processed/rbp_dataset_eskapee_structural_review.csv")
    ap.add_argument("--min_len", type=int, default=80)
    ap.add_argument("--max_len", type=int, default=1400)
    args = ap.parse_args()

    # Input, output, and review CSV paths
    in_csv = Path(args.in_csv)          # "Strict" dataset we generated in Stage 1b 
    out_csv = Path(args.out_csv)
    review_csv = Path(args.review_csv)

    # Raise error if input CSV does not exist
    if not in_csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {in_csv}")

    # Convert input CSV to DataFrame
    df = pd.read_csv(in_csv)

    # Raise error if input CSV does not contain 'product' and 'aa_sequence' columns
    if "product" not in df.columns or "aa_sequence" not in df.columns:
        raise ValueError("Expected columns 'product' and 'aa_sequence' in input CSV.")

    # Fill NaN values in 'product' column with an empty string and convert to lowercase (for case-insensitive matching)
    prod = df["product"].fillna("").astype(str).str.lower() 
    
    # Convert the aa_sequence column to string and then count the number of characters in the string
    aa_len = df["aa_sequence"].astype(str).str.len()    # = a Series whose each row is the length of the corresponding aa_sequence 

    # Build True/False masks for keep, drop, and review
    keep_mask = build_mask(prod, KEEP_PATTERNS)
    drop_mask = build_mask(prod, DROP_PATTERNS)
    review_mask = build_mask(prod, REVIEW_PATTERNS)

    # Length mask: True only if the length is between min_len and max_len amino acids
    length_mask = (aa_len >= args.min_len) & (aa_len <= args.max_len)

    # Final structural set: Keep only rows (proteins) that match keep patterns, don't match drop patterns, and pass the length filter
    structural_df = df[keep_mask & ~drop_mask & length_mask].copy()

    # Review set: ambiguous but potentially useful proteins (match review patterns) for later inspection
    review_df = df[review_mask & ~keep_mask & ~drop_mask & length_mask].copy()

    # Write the structural and review datasets to CSV
    structural_df["aa_len"] = structural_df["aa_sequence"].astype(str).str.len()
    review_df["aa_len"] = review_df["aa_sequence"].astype(str).str.len()

    # Make sure the output directories exist
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    structural_df.to_csv(out_csv, index=False)
    review_df.to_csv(review_csv, index=False)

    # Print a summary
    print(f"✅ Saved structural dataset: {out_csv}")
    print(f"✅ Saved review dataset:     {review_csv}")

    print("\nStructural dataset size:")
    print(structural_df.shape)

    print("\nStructural dataset host distribution:")
    print(structural_df["host_genus"].value_counts())

    print("\nTop structural product strings:")
    print(structural_df["product"].value_counts().head(20))

    print("\nReview dataset size:")
    print(review_df.shape)

    if len(review_df) > 0:
        print("\nTop review product strings:")
        print(review_df["product"].value_counts().head(20))


if __name__ == "__main__":
    main()