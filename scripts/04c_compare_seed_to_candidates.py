from __future__ import annotations
import argparse
import json
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmForMaskedLM


def embed_sequences(
    sequences: list[str],
    model_name: str,
    batch_size: int,
    max_aa: int,
    device: str,
) -> np.ndarray:
    """
    Embed amino-acid sequences using the encoder part of an MLM checkpoint.
    Returns a numpy array of shape [N, H].
    """
    # Load the tokenizer and model, move the model to the device and set it to evaluation mode
    tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)
    model = EsmForMaskedLM.from_pretrained(model_name)
    model.to(device)
    model.eval()

    embs = []

    # Embed sequences in batches
    with torch.no_grad():
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

            # Pass the tokens through the model and extract the last hidden state to obtain the embeddings
            out = model.esm(**toks)
            last_hidden_state = out.last_hidden_state

            # Compute the mean of the last_hidden_state along axis 1
            mask = toks["attention_mask"].unsqueeze(-1).to(last_hidden_state.dtype)      # [B, L] -> [B, L, 1]
            summed = (last_hidden_state * mask).sum(dim=1)                               # [B, L, H] -> [B, H]: Sum along the axis 1 real tokens collapsing the aa sequence length of each batch protein into a single vector
            denom = mask.sum(dim=1).clamp(min=1.0)                                       # Number of real tokens in each batch: Sum the MASK along the sequence length and secure a non-zero denominator with clamp
            pooled = summed / denom                                                      # Mean along axis 1: Summed tokens divided by the number of real tokens

            embs.append(pooled.detach().cpu())                                           # Add the embeddings to the list and move them to the CPU to avoid memory issues

    emb = torch.cat(embs, dim=0).numpy()    # Concatenate the embeddings rowwise and convert them to a numpy array
    return emb


def load_seed_row(seed_csv: Path, seed_protein_id: str) -> pd.Series:
    """ Load the seed RBP sequence from a CSV file using protein_id. """
    # Read the seed CSV
    df = pd.read_csv(seed_csv)

    # Select all rows where the "protein_id" column is equal to the seed_protein_id
    row = df[df["protein_id"].astype(str) == str(seed_protein_id)]
    
    # If no rows were found, raise an error
    if len(row) == 0:
        raise ValueError(f"Protein ID '{seed_protein_id}' not found in {seed_csv}")
    
    return row.iloc[0] # Select and return only the first row


def hamming_like_distance(seq_a: str, seq_b: str) -> int:
    """
    Count the differing positions between two protein sequences.
    If protein lengths are equal the distance is the number of differing positions, 
    otherwise the distance is the number of differing positions plus the difference in lengths.
    """
    # Find the length of the shorter sequence
    n = min(len(seq_a), len(seq_b))
    # Count the number of differing positions between the sequences
    dist = sum(1 for i in range(n) if seq_a[i] != seq_b[i])
    # If the lengths are different, add the difference to the distance
    dist += abs(len(seq_a) - len(seq_b))
    return dist 


def main():
    # Parse command-line arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("--run_dir", type=str, required=True)
    ap.add_argument("--seed_csv", type=str, required=True)
    ap.add_argument("--model_path", type=str, required=True)
    ap.add_argument("--label_classes_path", type=str, required=True)
    ap.add_argument("--esm_model", type=str, default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--batch_size", type=int, default=4)
    ap.add_argument("--max_aa", type=int, default=1022)
    ap.add_argument("--top_k", type=int, default=10)
    args = ap.parse_args()

    # Paths
    run_dir = Path(args.run_dir)
    top_csv = run_dir / "top_candidates.csv"
    metadata_json = run_dir / "run_metadata.json"

    # Raise errors if top_csv or metadata_json don't exist
    if not top_csv.exists():
        raise FileNotFoundError(f"Missing file: {top_csv}")
    if not metadata_json.exists():
        raise FileNotFoundError(f"Missing file: {metadata_json}")

    # Load metadata and top candidates csv
    metadata = json.loads(metadata_json.read_text())
    top_df = pd.read_csv(top_csv)

    # Load seed RBP sequence, target host and convert to string
    seed_row = load_seed_row(Path(args.seed_csv), metadata["seed_protein_id"]) 
    seed_seq = str(seed_row["aa_sequence"]) 
    seed_host = str(seed_row["host_genus"])

    # Load model and label classes
    clf = joblib.load(args.model_path)
    label_classes = json.loads(Path(args.label_classes_path).read_text())

    # Check that the target host is in the label classes
    if metadata["target_host"] not in label_classes:
        raise ValueError(f"Target host '{metadata['target_host']}' not found in label classes.")

    device = "cuda" if torch.cuda.is_available() else "cpu"
    
    # Get the index of the target host in the label classes
    target_idx = label_classes.index(metadata["target_host"])

    # Score the seed sequence
    seed_emb = embed_sequences(
        sequences=[seed_seq],
        model_name=args.esm_model,
        batch_size=1,
        max_aa=args.max_aa,
        device=device,
    )
    
    seed_probs = clf.predict_proba(seed_emb)[0]                 # Get the predicted probabilities for the seed
    seed_pred_idx = int(np.argmax(seed_probs))                  # Get the index of the highest predicted label
    seed_pred_label = label_classes[seed_pred_idx]              # Get the predicted host label
    seed_target_score = float(seed_probs[target_idx])           # Get the predicted score for the target host

    # Compare top candidates against seed
    compare_df = top_df.head(args.top_k).copy()
    compare_df["seed_protein_id"] = metadata["seed_protein_id"]
    compare_df["seed_source_host"] = seed_host
    compare_df["seed_pred_label"] = seed_pred_label
    compare_df["seed_target_score"] = seed_target_score
    compare_df["delta_target_score"] = compare_df["target_score"] - seed_target_score 
    compare_df["distance_from_seed"] = compare_df["aa_sequence"].apply(lambda s: hamming_like_distance(seed_seq, str(s)))

    # Save results to CSV
    compare_df.to_csv(run_dir / "seed_vs_top_candidates.csv", index=False)

    # Create summary
    summary = {
        "seed_protein_id": metadata["seed_protein_id"],
        "seed_source_host": seed_host,
        "target_host": metadata["target_host"],
        "seed_pred_label": seed_pred_label,
        "seed_target_score": seed_target_score,
        "best_candidate_target_score": float(compare_df["target_score"].max()),
        "best_delta_target_score": float(compare_df["delta_target_score"].max()),
        "mean_delta_target_score_top_k": float(compare_df["delta_target_score"].mean()),
        "mean_distance_from_seed_top_k": float(compare_df["distance_from_seed"].mean()),
        "top_k": int(args.top_k),
    }

    # Save summary
    (run_dir / "seed_vs_top_candidates_summary.json").write_text(json.dumps(summary, indent=2))

    # Create Markdown report
    lines = []
    lines.append("# Seed vs Top Candidates")
    lines.append("")
    lines.append(f"- Seed protein ID: `{summary['seed_protein_id']}`")
    lines.append(f"- Seed source host: `{summary['seed_source_host']}`")
    lines.append(f"- Target host: `{summary['target_host']}`")
    lines.append(f"- Seed predicted label: `{summary['seed_pred_label']}`")
    lines.append(f"- Seed target score: `{summary['seed_target_score']:.6f}`")
    lines.append(f"- Best candidate target score: `{summary['best_candidate_target_score']:.6f}`")
    lines.append(f"- Best delta target score: `{summary['best_delta_target_score']:.6f}`")
    lines.append(f"- Mean delta target score (top-{summary['top_k']}): `{summary['mean_delta_target_score_top_k']:.6f}`")
    lines.append(f"- Mean distance from seed (top-{summary['top_k']}): `{summary['mean_distance_from_seed_top_k']:.2f}`")
    lines.append("")
    lines.append("## Top candidates vs seed")
    lines.append("")
    lines.append(
        compare_df[
            [
                "candidate_id",
                "pred_label",
                "target_score",
                "seed_target_score",
                "delta_target_score",
                "distance_from_seed",
                "mutations",
                "n_mutations",
            ]
        ].to_markdown(index=False)
    )
    lines.append("")

    # Save report
    report_path = run_dir / "seed_vs_top_candidates_report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")

    # Print success message
    print(f"✅ Saved: {run_dir / 'seed_vs_top_candidates.csv'}")
    print(f"✅ Saved: {run_dir / 'seed_vs_top_candidates_summary.json'}")
    print(f"✅ Saved: {report_path}")


if __name__ == "__main__":
    main()