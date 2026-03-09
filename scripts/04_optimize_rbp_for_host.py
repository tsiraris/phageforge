from __future__ import annotations
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import argparse
import json
import random
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmModel

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

def set_seed(seed: int):
    """ Set Python / NumPy / PyTorch random seeds for reproducibility. """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def mean_pool(last_hidden_state: torch.Tensor, attention_mask: torch.Tensor) -> torch.Tensor:
    """ Mean-pool token embeddings using the attention mask. """
    mask = attention_mask.unsqueeze(-1).to(last_hidden_state.dtype)  # [B, L, 1]
    summed = (last_hidden_state * mask).sum(dim=1)                   # [B, H]
    denom = mask.sum(dim=1).clamp(min=1.0)                           # [B, 1]
    return summed / denom


def load_seed_sequence(seed_csv: Path, seed_protein_id: str) -> tuple[str, str, str, str]:
    """
    Load the seed RBP sequence from a CSV file using protein_id.
    Returns:
        virus_accession, host_genus, protein_id, aa_sequence
    """
    df = pd.read_csv(seed_csv)                                                          # Read CSV
    row = df[df["protein_id"].astype(str) == str(seed_protein_id)]                      # Select all rows where the "protein_id" column is equal to the seed_protein_id
    if len(row) == 0:
        raise ValueError(f"Protein ID '{seed_protein_id}' not found in {seed_csv}")
    row = row.iloc[0]                                                                   # Select only the first row 
    return (
        str(row["virus_accession"]),
        str(row["host_genus"]),
        str(row["protein_id"]),
        str(row["aa_sequence"]),
    )


def mutate_sequence(seq: str, n_mutations: int, rng: random.Random) -> tuple[str, str]:
    """
    Create a mutated version of a sequence by applying n_mutations random amino-acid substitutions.
    Returns:
        mutated_sequence = the resulting mutated sequence, mutation_string = a string of all the mutations that were applied (aa changes) and in which positions of the sequence
    mutation_string example:
        "A15G;L42P" = mutation at position 15 from A to G and mutation at position 42 from L to P
    """
    # If n_mutations is zero, return the original sequence
    if n_mutations <= 0:
        return seq, ""

    seq_list = list(seq)                            # Convert the sequence to a list
    mutable_positions = list(range(len(seq_list)))  # List of positions that can be mutated

    # If there are fewer mutable positions than n_mutations, set n_mutations to the number of mutable positions
    if len(mutable_positions) < n_mutations:
        n_mutations = len(mutable_positions)

    # Randomly select n_mutations positions
    chosen_positions = rng.sample(mutable_positions, n_mutations)
    mutations = []

    # Apply the mutations
    for pos in chosen_positions:                                # For each chosen position
        old_aa = seq_list[pos]                                  # Get the old amino acid at the position
        choices = [aa for aa in AMINO_ACIDS if aa != old_aa]    # List of all amino acids except the old one
        new_aa = rng.choice(choices)                            # Randomly select a new amino acid from the choices
        seq_list[pos] = new_aa                                  # Replace the old amino acid with the new one
        mutations.append(f"{old_aa}{pos+1}{new_aa}")            # Add the mutation to the mutations track list

    mutated_seq = "".join(seq_list)                             # Join the mutated sequence to a string
    mutation_str = ";".join(mutations)                          # Join all the mutations track strings with a semicolon  
    return mutated_seq, mutation_str                                            


def embed_sequences( sequences: list[str], model_name: str, batch_size: int, max_aa: int, device: str) -> np.ndarray:
    """
    Embed amino-acid sequences with ESM-2 and return a numpy array of shape [N, H].
    """
    # Load the tokenizer and model, move the model to the device and set it to evaluation mode
    tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)
    model = EsmModel.from_pretrained(model_name)
    model.to(device)
    model.eval()

    embs = []

    # Embed sequences in batches
    with torch.no_grad():
        for start in range(0, len(sequences), batch_size):
            batch = sequences[start : start + batch_size]

            toks = tokenizer(
                batch,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=max_aa,
            ).to(device)

            out = model(**toks)
            pooled = mean_pool(out.last_hidden_state, toks["attention_mask"])
            embs.append(pooled.detach().cpu())

    emb = torch.cat(embs, dim=0).numpy()    # Concatenate the embeddings rowwise and convert them to a numpy array
    return emb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed_csv", type=str, default="data/processed/rbp_dataset_eskapee_strict.csv")
    ap.add_argument("--seed_protein_id", type=str, required=True)
    ap.add_argument("--target_host", type=str, required=True)                               # The host genus to optimize for
    ap.add_argument("--model_path", type=str, required=True)
    ap.add_argument("--label_classes_path", type=str, required=True)
    ap.add_argument("--esm_model", type=str, default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--batch_size", type=int, default=4)
    ap.add_argument("--max_aa", type=int, default=1022)
    ap.add_argument("--rounds", type=int, default=5)                                        # Number of optimization rounds
    ap.add_argument("--candidates_per_round", type=int, default=64)                         # Number of candidate sequences per round
    ap.add_argument("--mutations_per_candidate", type=int, default=2)                       # Number of mutations per candidate
    ap.add_argument("--keep_top_k", type=int, default=10)                                   # Number of top candidates to keep
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out_dir", type=str, default="results/design_runs")
    args = ap.parse_args()
 
    set_seed(args.seed)                   # Set the random seed for reproducibility (Numpy, PyTorch, random)
    rng = random.Random(args.seed)        # Create a random number generator

    # Create output directory
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load seed sequence, model, and label classes
    seed_csv = Path(args.seed_csv)
    model_path = Path(args.model_path)
    label_classes_path = Path(args.label_classes_path)

    # Load seed sequence
    virus_accession, source_host, protein_id, seed_sequence = load_seed_sequence(seed_csv=seed_csv, seed_protein_id=args.seed_protein_id)

    # Check if seed sequence length exceeds max_aa
    if len(seed_sequence) > args.max_aa:
        raise ValueError(
            f"Seed sequence length {len(seed_sequence)} exceeds max_aa={args.max_aa}. "
            f"Choose a shorter seed or add chunking later."
        )

    # Load model and label classes
    clf = joblib.load(model_path)
    label_classes = json.loads(label_classes_path.read_text())  # list of strings

    # Check if target host is in label classes
    if args.target_host not in label_classes:
        raise ValueError(
            f"Target host '{args.target_host}' not found in label classes: {label_classes}"
        )

    device = "cuda" if torch.cuda.is_available() else "cpu"

    # Track all generated candidates across rounds
    all_rows = []

    # A list of dictionaries, each representing the best so far protein sequences which are in the current pool.
    current_pool = [
        {
            "round": 0,
            "parent_rank": 0,
            "virus_accession": virus_accession,
            "source_host": source_host,
            "protein_id": protein_id,
            "candidate_id": "seed",
            "mutations": "",
            "aa_sequence": seed_sequence,
        }
    ]

    # Get the index of the target host
    target_idx = label_classes.index(args.target_host)              

    # Beam search optimization
    # In each round, we mutate the parents (= sequences currently in the pool), score the children, and keep only the top k children to become the parents for the next round.
    for round_idx in range(1, args.rounds + 1):                     # Iterate over the optimization rounds (= number of generations we run this process for)
        proposed_rows = []
        candidate_counter = 0

        # Generate mutations from the current pool
        for parent_rank, parent in enumerate(current_pool):         # Iterate over the parents in the current pool
            for _ in range(args.candidates_per_round):              # Iterate over the number of candidates per round (= how many random mutations for each parent)
                
                # Mutate one of the parent sequences
                mutated_seq, mutation_str = mutate_sequence(parent["aa_sequence"], n_mutations=args.mutations_per_candidate, rng=rng)   

                # Add the mutated sequence to the proposed rows
                proposed_rows.append(
                    {
                        "round": round_idx,
                        "parent_rank": parent_rank,
                        "virus_accession": virus_accession,
                        "source_host": source_host,
                        "protein_id": protein_id,
                        "candidate_id": f"round{round_idx}_cand{candidate_counter}",
                        "mutations": mutation_str,
                        "aa_sequence": mutated_seq,
                    }
                )
                candidate_counter += 1

        # Convert proposed rows to dataframe and drop duplicates
        proposed_df = pd.DataFrame(proposed_rows).drop_duplicates(subset=["aa_sequence"]).reset_index(drop=True)

        # Embed all proposed candidates
        embeddings = embed_sequences(
            sequences=proposed_df["aa_sequence"].tolist(),
            model_name=args.esm_model,
            batch_size=args.batch_size,
            max_aa=args.max_aa,
            device=device,
        )

        # Score each candidate with the trained linear probe
        probs = clf.predict_proba(embeddings)               # Predicted probabilities for all the candidates
        pred_idx = probs.argmax(axis=1)                     # Index of the predicted label (host)
        pred_label = [label_classes[i] for i in pred_idx]   # Predicted label (host)
        target_scores = probs[:, target_idx]                # Predicted score for the target host

        proposed_df["pred_label"] = pred_label              # Add predicted labels column to the dataframe
        proposed_df["target_host"] = args.target_host       # Add target host column
        proposed_df["target_score"] = target_scores         # Add target score column

        # Count number of substitutions we introduced this round (through the number of semicolons in the mutation string) and add them as a new column to the dataframe
        proposed_df["n_mutations"] = proposed_df["mutations"].apply(lambda s: 0 if s == "" else len(s.split(";")))

        # Sort proposed candidates by target score descending and number of mutations ascending
        proposed_df = proposed_df.sort_values(by=["target_score", "n_mutations"], ascending=[False, True]).reset_index(drop=True)

        # Assign a rank to each proposed candidate (= the order in which they were proposeds)
        proposed_df["rank_in_round"] = range(1, len(proposed_df) + 1)

        # Save all proposed candidates for each round to csv
        proposed_df.to_csv(out_dir / f"round_{round_idx}_candidates.csv", index=False)

        # Keep track of all proposed candidates across rounds
        all_rows.append(proposed_df)

        # Keep only the top k candidates of this round to seed the next round
        current_pool = proposed_df.head(args.keep_top_k).to_dict(orient="records")

    # Save all rounds together in a single csv
    all_candidates_df = pd.concat(all_rows, axis=0).reset_index(drop=True)
    all_candidates_df.to_csv(out_dir / "all_candidates.csv", index=False)

    # Sort the proposed candidates across all rounds by target score descending, number of mutations ascending, and save final top head (50) candidates
    final_top_df = (
        all_candidates_df.sort_values(by=["target_score", "n_mutations"], ascending=[False, True])
        .drop_duplicates(subset=["aa_sequence"])
        .head(50)                                   
        .reset_index(drop=True)
    )
    final_top_df.to_csv(out_dir / "top_candidates.csv", index=False)

    # Save a compact run metadata file
    metadata = {
        "seed_protein_id": args.seed_protein_id,
        "virus_accession": virus_accession,
        "source_host": source_host,
        "target_host": args.target_host,
        "esm_model": args.esm_model,
        "model_path": str(model_path),
        "label_classes_path": str(label_classes_path),
        "rounds": args.rounds,
        "candidates_per_round": args.candidates_per_round,
        "mutations_per_candidate": args.mutations_per_candidate,
        "keep_top_k": args.keep_top_k,
        "seed": args.seed,
        "n_total_candidates": int(len(all_candidates_df)),
        "n_top_candidates_saved": int(len(final_top_df)),
    }
    (out_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2))

    print(f"✅ Saved: {out_dir / 'all_candidates.csv'}")
    print(f"✅ Saved: {out_dir / 'top_candidates.csv'}")
    print(f"✅ Saved: {out_dir / 'run_metadata.json'}")


if __name__ == "__main__":
    main()