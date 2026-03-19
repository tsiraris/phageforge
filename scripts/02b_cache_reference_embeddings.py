from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmForMaskedLM


def embed_sequences(
    sequences: list[str],
    tokenizer,
    model,
    batch_size: int,
    max_aa: int,
    device: str,
) -> np.ndarray:
    """
    Embed amino-acid sequences using the encoder part of an MLM (Masked Language Model) checkpoint.
    Returns a numpy array of shape [N, H].
    """
    embs = []

    # Embed sequences in batches
    with torch.inference_mode():
        for start in range(0, len(sequences), batch_size):
            batch = sequences[start : start + batch_size]

            # Tokenize the input sequences batch
            toks = tokenizer(
                batch,
                return_tensors="pt",    # Return PyTorch tensors
                padding=True,
                truncation=True,
                max_length=max_aa,
            ).to(device)

            # Pass the tokens through the model
            out = model.esm(**toks)

            # Extract the last_hidden_state from the output to obtain the embeddings
            last_hidden_state = out.last_hidden_state  # [B, L, H]

            # Compute the mean of the last_hidden_state along axis 1
            mask = toks["attention_mask"].unsqueeze(-1).to(last_hidden_state.dtype)      # [B, L] -> [B, L, 1]
            summed = (last_hidden_state * mask).sum(dim=1)                               # [B, L, H] -> [B, H]: Sum along the axis 1 real tokens collapsing the aa sequence length of each batch protein into a single vector
            denom = mask.sum(dim=1).clamp(min=1.0)                                       # Number of real tokens in each batch: Sum the MASK along the sequence length and secure a non-zero denominator with clamp
            pooled = summed / denom                                                      # Mean along axis 1: Summed tokens divided by the number of real tokens

            embs.append(pooled.detach().cpu())                                           # Add the embeddings to the list and move them to the CPU to avoid memory issues

    emb = torch.cat(embs, dim=0).numpy()    # Concatenate the embeddings rowwise and convert them to a numpy array
    return emb


def normalize_rows(x: np.ndarray) -> np.ndarray:
    """ L2-normalize each embedding row for cosine-similarity computations. """
    return x / np.clip(np.linalg.norm(x, axis=1, keepdims=True), 1e-12, None)   # Clip the norm to avoid division by zero


def main():
    # Parse command-line arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference_csv", type=str, required=True)
    ap.add_argument("--esm_model", type=str, default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--batch_size", type=int, default=8)
    ap.add_argument("--max_aa", type=int, default=1022)
    ap.add_argument("--out_embeddings_npy", type=str, required=True)
    ap.add_argument("--out_index_csv", type=str, required=True)
    args = ap.parse_args()

    reference_csv = Path(args.reference_csv)
    out_embeddings_npy = Path(args.out_embeddings_npy)
    out_index_csv = Path(args.out_index_csv)

    if not reference_csv.exists():
        raise FileNotFoundError(f"Missing file: {reference_csv}")

    # Read the novelty-reference dataset that will be reused across many summary runs
    ref_df = pd.read_csv(reference_csv)                                                  
    if "aa_sequence" not in ref_df.columns:
        raise ValueError("Expected an 'aa_sequence' column in the reference CSV.")

    # Drop rows without sequences (cannot be embedded) and drop duplicate sequences (so cached novelty embeddings are not wasted on repeated proteins)
    ref_df = ref_df.dropna(subset=["aa_sequence"]).copy().reset_index(drop=True)         
    ref_df = ref_df.drop_duplicates(subset=["aa_sequence"]).reset_index(drop=True)       
    
    device = "cuda" if torch.cuda.is_available() else "cpu"

    # Load ESM MLM once and reuse it for the whole reference panel
    tokenizer = AutoTokenizer.from_pretrained(args.esm_model, do_lower_case=False)
    esm_mlm = EsmForMaskedLM.from_pretrained(args.esm_model)
    esm_mlm.to(device)
    esm_mlm.eval()

    # Compute one ESM embedding for each unique reference sequence
    embeddings = embed_sequences(
        sequences=ref_df["aa_sequence"].astype(str).tolist(),
        tokenizer=tokenizer,
        model=esm_mlm,
        batch_size=args.batch_size,
        max_aa=args.max_aa,
        device=device,
    )                                                                                   

    # Normalize and cast the cached reference embeddings once so downstream novelty comparisons only need matrix multiplication
    embeddings_normalized = normalize_rows(embeddings).astype(np.float32)     # cast to float32 to save space          

    # Create output directories
    out_embeddings_npy.parent.mkdir(parents=True, exist_ok=True)
    out_index_csv.parent.mkdir(parents=True, exist_ok=True)

    # Save the normalized cached reference embeddings to disk and the aligned reference index so the cached embeddings can be interpreted later
    np.save(out_embeddings_npy, embeddings_normalized)                                  
    ref_df.to_csv(out_index_csv, index=False)                                           

    # Print some stats
    print(f"✅ Saved cached reference embeddings: {out_embeddings_npy}")
    print(f"✅ Saved cached reference index:      {out_index_csv}")
    print("")
    print("Reference dataset size after deduplication:")
    print(ref_df.shape)


if __name__ == "__main__":
    main()
