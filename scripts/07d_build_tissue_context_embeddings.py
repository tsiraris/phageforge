"""
===============================================================================
07d: Build optional tissue-context embeddings for Stage 07 multimodal reranking.
===============================================================================

This script builds the optional tissue branch.
It accepts either:
- precomputed embedding columns, or
- raw numeric/categorical metadata that can be converted into one embedding table

The result is a tensor plus a metadata CSV that later scripts can merge into the
final multimodal ranking stage.
"""

from __future__ import annotations                                                    # Enable postponed annotation evaluation for cleaner typing.
import argparse                                                                       # Parse command-line arguments.
import json                                                                           # Save a small metadata JSON summary.
from pathlib import Path                                                              # Build output paths robustly.
import numpy as np                                                                    # Assemble numeric embedding blocks.
import pandas as pd                                                                   # Read the tissue-context CSV and write aligned metadata.
import torch                                                                          # Save the final embedding matrix as a tensor.
from sklearn.preprocessing import OneHotEncoder, StandardScaler                       # Encode categorical variables and normalize numeric ones.


def parse_args() -> argparse.Namespace:
    """Return command-line arguments for building optional tissue-context embeddings."""
    ap = argparse.ArgumentParser(description="Build optional tissue context embeddings for Stage 07.")                                              # Create the parser for this script.
    ap.add_argument("--input_csv", type=str, required=True, help="CSV with histopathology / omics / metadata context rows.")                        # Point to the source tissue table.
    ap.add_argument("--embedding_cols_prefix", type=str, default="", help="Optional prefix identifying precomputed embedding columns.")             # Enable direct use of precomputed embeddings.
    ap.add_argument("--context_id_col", type=str, default="tissue_context_id", help="Column that uniquely identifies each tissue context row.")     # Choose or create the context ID column.
    ap.add_argument("--out_pt", type=str, required=True, help="Where to store the output embedding tensor.")                                        # Point to the tensor output path.
    ap.add_argument("--out_csv", type=str, required=True, help="Where to store metadata aligned to the output tensor.")                             # Point to the aligned metadata CSV path.
    return ap.parse_args()                                                                                                                          # Parse the CLI and return the namespace.


# Main script: Loads tissue data, builds embeddings, and saves the tissue embedding matrix as a tensor, and a small metadata CSV.
def main() -> None:
    args = parse_args() 
    # Read the tissue-context input table csv, and ensure that every row has a stable tissue-context column ID.
    df = pd.read_csv(args.input_csv)                                                  # Read the tissue-context input table from disk.
    if args.context_id_col not in df.columns:                                         # Ensure that every row has a stable tissue-context identifier.
        df[args.context_id_col] = [f"ctx_{i:04d}" for i in range(len(df))]            # Create synthetic context IDs when the column is absent.

    # Either reuse precomputed embeddings (from a pathological FM for example) directly when a prefix was supplied.
    if args.embedding_cols_prefix:  
        # Scan the CSV for any columns starting with the used prefix, group them, and converts them into a float matrix.
        emb_cols = [c for c in df.columns if c.startswith(args.embedding_cols_prefix)]  # Collect all columns that match the requested embedding prefix.
        if not emb_cols:                                                                # Ensure the requested prefix actually matches some columns.
            raise ValueError(f"No columns found with prefix {args.embedding_cols_prefix!r}.")
        emb = df[emb_cols].to_numpy(dtype=np.float32)                                   # Materialize the precomputed embedding block as a float matrix.
    # Otherwise, build embeddings from numeric/categorical columns (of raw clinical metadata for example).
    else:         
        # Separate CSV columns into numeric and categorical blocks                                                                    
        numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c]) and c != args.context_id_col]  # Collect numeric columns excluding the ID.
        cat_cols = [c for c in df.columns if c != args.context_id_col and c not in numeric_cols]  # Collect non-numeric columns excluding the ID.
        blocks = []                                                                   # Collect numeric and categorical feature blocks before concatenation.

        # Standardize numeric and one-hot encode categorical columns
        if numeric_cols:                                                              # Standardize numeric metadata when numeric columns exist.
            # Substract the mean and divide by the standard deviation for each column.
            scaler = StandardScaler()                                                 # Create a standard scaler for numeric features.
            blocks.append(scaler.fit_transform(df[numeric_cols].fillna(0.0)))         # Fit-transform numeric columns after filling missing values with zeros.
        if cat_cols:                                                                  
            # Encode categorical columns as one-hot vectors, and fill missing data with a special "missing" category.
            enc = OneHotEncoder(handle_unknown="ignore", sparse_output=False)         
            blocks.append(enc.fit_transform(df[cat_cols].fillna("missing")))          

        # If no blocks were collected, create a one-dimensional zero embedding as a safe fallback (for example, for a CSV with only the context ID column).
        if not blocks:                                                                
            emb = np.zeros((len(df), 1), dtype=np.float32)                            
        # Otherwise, concatenate the numeric and categorical blocks to build the final tissue embedding matrix.
        else:                                                                         
            emb = np.concatenate(blocks, axis=1).astype(np.float32)                   

    # Save the tissue embedding matrix as a tensor, and the aligned metadata as a CSV.
    out_pt = Path(args.out_pt)                                                        # Convert the tensor output path into a Path object.
    out_csv = Path(args.out_csv)                                                      # Convert the metadata output path into a Path object.
    out_pt.parent.mkdir(parents=True, exist_ok=True)                                  # Create the tensor output directory if needed.
    out_csv.parent.mkdir(parents=True, exist_ok=True)                                 # Create the metadata output directory if needed.

    torch.save(torch.tensor(emb, dtype=torch.float32), out_pt)                        # Save the tissue embedding matrix as a float tensor.
    df.to_csv(out_csv, index=False)                                                   # Save the aligned tissue metadata CSV.

    meta = {                                                                          # Build a small metadata JSON summary for traceability.
        "rows": int(len(df)),                                                         # Record how many tissue-context rows were embedded.
        "embedding_dim": int(emb.shape[1]),                                           # Record the final tissue embedding dimensionality.
        "context_id_col": args.context_id_col,                                        # Record which column stores the stable tissue-context IDs.
    }
    (out_pt.parent / "tissue_context_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")  # Save the metadata summary as JSON.
    print(f"Wrote: {out_pt}")                                                         # Print the tensor output path.
    print(f"Wrote: {out_csv}")                                                        # Print the metadata CSV output path.


if __name__ == "__main__":                                                            # Standard Python entrypoint guard.
    main()                                                                            # Execute the tissue-embedding CLI.
