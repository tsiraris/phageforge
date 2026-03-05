from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import torch
from tqdm import tqdm
from transformers import AutoTokenizer, AutoModel
import json

def mean_pool(last_hidden_state: torch.Tensor, attention_mask: torch.Tensor) -> torch.Tensor:
    """ Computes the mean of the last_hidden_state along axis 1. """
    # last_hidden_state size: [B, L, H] with B=batch_size, L=sequence_length, H=hidden_size (embedding dimension)
    # attention_mask to mark real tokens vs padding: [B, L] - 1 for real tokens and 0 for pad tokens
    mask = attention_mask.unsqueeze(-1).to(last_hidden_state.dtype)  # [B, L] -> [B, L, 1]
    summed = (last_hidden_state * mask).sum(dim=1)                   # [B, L, H] -> [B, H]; Sum along axis 1 the real tokens collapsing the aa sequence length of each batch protein into a single vector
    denom = mask.sum(dim=1).clamp(min=1.0)                           # Number of real tokens in each batch: Sum the MASK along the sequence length and secure a non-zero denominator with clamp
    return summed / denom                                            # Mean along axis 1: Summed tokens divided by the number of real tokens


def main():
    # Create an ArgumentParser object to run the same script with different parameters without editing code.
    ap = argparse.ArgumentParser()  
    ap.add_argument("--csv", type=str, default="data/processed/rbp_dataset_eskapee_stage1.csv")
    ap.add_argument("--out_dir", type=str, default="data/processed")
    ap.add_argument("--model", type=str, default="facebook/esm2_t12_35M_UR50D") # ESM-2 transformer encoder model
    ap.add_argument("--batch_size", type=int, default=8)
    ap.add_argument("--max_aa", type=int, default=1022)
    args = ap.parse_args()

    # Load data
    csv_path = Path(args.csv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Drop sequences with aa_sequence length longer than max_aa    
    df = pd.read_csv(csv_path)
    df["aa_len"] = df["aa_sequence"].astype(str).str.len()  # Convert the aa_sequence column to a string and then count the number of characters in the string
    keep = df["aa_len"] <= args.max_aa                      # Keep only the rows where the aa_len is less than or equal to max_aa   
    dropped = (~keep).sum()                                 # Count the number of rows that were dropped
    if dropped > 0:
        print(f"[INFO] Dropping {dropped} sequences longer than {args.max_aa} aa (ESM max length constraint).")
    df = df[keep].reset_index(drop=True)                    # Reset the index of the dataframe

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[INFO] device={device}")

    # Load the tokenizer and model, move the model to the device and set it to evaluation mode
    tokenizer = AutoTokenizer.from_pretrained(args.model, do_lower_case=False)  
    model = AutoModel.from_pretrained(args.model)                               
    model.to(device)                                                            
    model.eval()                                                                

    # To store the embeddings
    embs = []

    # Keep an index mapping so you can trace back the embeddings to the original dataframe
    index_df = df[["virus_accession", "host_genus", "protein_id", "product"]].copy()    # Copy the columns from df to index_df
    index_df["row_id"] = range(len(df))                                                 # Add a row_id column

    with torch.no_grad():
        for start in tqdm(range(0, len(df), args.batch_size)):
            batch = df.iloc[start : start + args.batch_size]  # Slice the dataframe into batches
            seqs = batch["aa_sequence"].tolist()              # Collect the aa_sequences column in a list

            # Tokenize the raw sequences (strings) into integers  
            toks = tokenizer(
                seqs,                       # List of sequences to tokenize - HF ESM tokenizer supports raw sequences.
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=args.max_aa,
            ).to(device)

            # Run the model forward pass to compute the embeddings
            out = model(**toks)    # shape: [B, L, H]
            
            # Mean-pool token embeddings into one vector per amino acid sequence
            pooled = mean_pool(out.last_hidden_state, toks["attention_mask"])  # [B, H]
            
            # Detach the tensor from the computation graph, move it to the CPU (so GPU memory doesn’t grow batch after batch) and append it to the embs list
            embs.append(pooled.detach().cpu())  

    # Concatenate the embeddings row-wise
    emb = torch.cat(embs, dim=0)  # [N, H]; N = number of different protein samples remaining in the dataframe
    print(f"[INFO] embeddings shape = {tuple(emb.shape)}")

    # Save the embeddings and index dataframe
    emb_path = out_dir / "esm2_embeddings.pt"
    idx_path = out_dir / "esm2_embeddings_index.csv"
    torch.save(emb, emb_path)
    index_df.to_csv(idx_path, index=False)

    # Save a small metadata JSON to be able to reproduce the results
    meta = {
        "model": args.model,
        "batch_size": args.batch_size,
        "max_aa": args.max_aa,
        "n_sequences": int(emb.shape[0]),
        "embedding_dim": int(emb.shape[1]),
    }
    (out_dir / "esm2_run_metadata.json").write_text(json.dumps(meta, indent=2))
    
    print(f"✅ Saved embeddings: {emb_path}")
    print(f"✅ Saved index:      {idx_path}")


if __name__ == "__main__":
    main()