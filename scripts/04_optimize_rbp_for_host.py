from __future__ import annotations
import argparse
import json
import random
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmForMaskedLM

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


def set_seed(seed: int):
    """ Set Python / NumPy / PyTorch random seeds for reproducibility. """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


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


def choose_mutation_positions_randomly(seq: str, n_mutations: int, rng: random.Random) -> list[int]:
    """
    Randomly choose which sequence positions to mutate.
    Returns a list of positions in the protein sequence.
    """
    # If n_mutations is zero, return the original sequence (Adapted for this function: return empty list)
    if n_mutations <= 0:
        return []

    positions = list(range(len(seq)))               # List of positions that can be mutated

    # If there are fewer mutable positions than n_mutations, set n_mutations to the number of mutable positions
    if len(positions) < n_mutations:
        n_mutations = len(positions)

    # Randomly select n_mutations positions from the positions list and sort them
    return sorted(rng.sample(positions, n_mutations))


def get_amino_acid_token_ids(tokenizer) -> list[int]:
    """
    Collect tokenizer ids corresponding to one-letter amino-acid tokens.
    """
    aa_token_ids = []                                                                               # List of tokenizer ids
    for aa in AMINO_ACIDS:                                                                          # For each one-letter amino acid
        tok_id = tokenizer.convert_tokens_to_ids(aa)                                                # Get the tokenizer id (= the integer corresponding to the one-letter amino acid)
        if tok_id is not None and tok_id != tokenizer.unk_token_id:                                 # If the tokenizer id is not None and not the unknown token
            aa_token_ids.append(int(tok_id))                                                        # Add the tokenizer id to the list
    if len(aa_token_ids) == 0:                                                                      # If the list is empty
        raise ValueError("Could not resolve one-letter amino-acid token ids from tokenizer.")       # Raise an error
    return aa_token_ids                                                                             # Return the list


def compute_position_entropies(
    seq: str,
    tokenizer,
    model,
    max_aa: int,
    device: str,
    aa_token_ids: list[int],
) -> np.ndarray:
    """
    Compute an uncertainty score (Shannon entropy) for each residue position in the input sequence.
    The entropy is computed from the ESM logits restricted to valid amino-acid tokens (higher entropy --> multiple amino acids plausible at that position).
    Returns a numpy array of shape [L] with uncertainty scores.
    """
    # If the sequence is longer than max_aa, raise an error
    if len(seq) > max_aa:
        raise ValueError(f"Sequence length {len(seq)} exceeds max_aa={max_aa}")

    # Tokenize one sequence: Add special tokens, look up the integers in the vocabulary, and return a PyTorch tensor dictionary of tokens
    toks = tokenizer(
        [seq],
        return_tensors="pt",
        padding=True,
        truncation=True,
        max_length=max_aa,
    ).to(device)                # [1, L] dictionary of PyTorch tensors; L = tokenized sequence length

    with torch.inference_mode():
        # Feed the unmasked sequence into the model and inspect its token distribution at each residue position
        outputs = model(input_ids=toks["input_ids"], attention_mask=toks["attention_mask"])
        # Extract the logits (raw, un-normalized prediction scores) for every single token in the vocabulary
        logits = outputs.logits  # [1, L, vocab]

    aa_logits = logits[0, 1 : len(seq) + 1, aa_token_ids]           # Restrict the logits to valid one-letter amino-acid tokens only [seq_len, 20]
    aa_probs = torch.softmax(aa_logits, dim=-1)                     # Get the softmax probabilities in the vertical dimension [seq_len, 20]
    aa_log_probs = torch.log(aa_probs.clamp(min=1e-12))             # Get the log probabilities [seq_len, 20]
    entropies = -(aa_probs * aa_log_probs).sum(dim=-1)              # Compute the Shannon entropy [seq_len]
    
    return entropies.detach().cpu().numpy()                         # Convert the PyTorch tensor to a numpy array


def choose_mutation_positions_entropy_guided(
    seq: str,
    n_mutations: int,
    position_scores: np.ndarray,
    rng: random.Random,
    candidate_pool_size: int,
) -> list[int]:
    """
    Choose mutation positions using uncertainty / entropy scores.

    Strategy:
    1. rank positions by entropy
    2. keep only the top candidate_pool_size uncertain positions
    3. sample without replacement from that pool using entropy-weighted probabilities

    Returns a list of positions in the protein sequence.
    """
    # If n_mutations is zero, return the original sequence (Adapted for this function: return empty list)
    if n_mutations <= 0:
        return []

    positions = list(range(len(seq)))               # List of positions that can be mutated

    # If there are fewer mutable positions than n_mutations, set n_mutations to the number of mutable positions
    if len(positions) < n_mutations:
        n_mutations = len(positions)

    # If no scores are available, fall back to random choice
    if len(position_scores) != len(seq):
        return choose_mutation_positions_randomly(seq=seq, n_mutations=n_mutations, rng=rng)

    # Rank positions by entropy descending and keep a candidate pool of the most uncertain positions
    ranked_positions = np.argsort(-position_scores).tolist()
    pool_size = min(max(candidate_pool_size, n_mutations), len(ranked_positions))
    candidate_positions = ranked_positions[:pool_size]

    # Convert candidate scores to positive sampling weights
    candidate_scores = np.array([max(float(position_scores[p]), 1e-8) for p in candidate_positions], dtype=float)   # Make sure all scores are positive (min = 1e-8) 
    candidate_scores = candidate_scores / candidate_scores.sum()                                                    # Normalize

    # Create a list of selected positions and a list of remaining positions and scores to sample from
    selected_positions = []
    remaining_positions = candidate_positions.copy()
    remaining_scores = candidate_scores.copy()

    # Sample without replacement using entropy-weighted probabilities
    while len(selected_positions) < n_mutations and len(remaining_positions) > 0:           # Keep sampling until n_mutations positions have been selected
        chosen_idx = rng.choices(                                                           # Randomly sample
            population=list(range(len(remaining_positions))),                               # From the list of remaining positions
            weights=remaining_scores.tolist(),                                              # Using the remaining scores as sampling weights
            k=1,                                                                            # Sample one position
        )[0]                                                                                # Get the index of the selected position
        selected_positions.append(remaining_positions.pop(chosen_idx))                      # Add the selected position to the list
        remaining_scores = np.delete(remaining_scores, chosen_idx)                          # Remove the selected position from the list

        if len(remaining_scores) > 0:                                                       # If there are remaining positions
            remaining_scores = remaining_scores / remaining_scores.sum()                    # Normalize the remaining scores

    return sorted(selected_positions)        # Return the list of selected positions


def propose_esm_guided_mutation(
    seq: str,
    positions: list[int],
    tokenizer,
    model,
    max_aa: int,
    device: str,
    top_k: int,     # Number of top-k amino acid predictions from ESM to consider
    rng: random.Random,
) -> tuple[str, str]:
    """
    Propose a mutated sequence by masking selected positions and sampling plausible residues
    from the top-k ESM predictions at each position.

    Returns:
        mutated_sequence = the resulting mutated sequence, mutation_string = a string of all the mutations that were applied (aa changes) and in which positions of the sequence
    mutation_string example:
        "A15G;L42P" = mutation at position 15 from A to G and mutation at position 42 from L to P
    """
    # If the sequence is longer than max_aa, raise an error
    if len(seq) > max_aa:
        raise ValueError(f"Sequence length {len(seq)} exceeds max_aa={max_aa}")

    # Tokenize one sequence: Add special tokens, look up the integers in the vocabulary, and return a PyTorch tensor dictionary of tokens
    toks = tokenizer(
        [seq],
        return_tensors="pt",
        padding=True,
        truncation=True,
        max_length=max_aa,
    ).to(device)                # [1, L] dictionary of PyTorch tensors; L = tokenized sequence length

    # Clone the tokenized sequence to avoid modifying it
    input_ids = toks["input_ids"].clone()

    # Get the specific integer (ID) of the <mask> token (blank space) from the tokenizer vocabulary
    mask_token_id = tokenizer.mask_token_id
    if mask_token_id is None:
        raise ValueError("Tokenizer does not expose a mask token id.")

    # Convert chosen for mutation positions to token positions by shifting them by 1 (because ESM tokenizers add special tokens at the beginning)
    token_positions = [p + 1 for p in positions]

    # Apply the <mask> token at the chosen positions
    for tp in token_positions:

        if tp >= input_ids.shape[1] - 1:        # Skip if the token position is out of bounds
            continue

        input_ids[0, tp] = mask_token_id        # Masked tokenized sequence tensor

    with torch.inference_mode():
        # Feed the mask tensor into the model, which will predict the most likely amino acid at each position based on the surrounding context
        outputs = model(input_ids=input_ids, attention_mask=toks["attention_mask"])
        # Extract the logits (raw, un-normalized prediction scores) for every single token in the vocabulary
        logits = outputs.logits  # [1, L, vocab]

    seq_list = list(seq)                                        # Convert the protein sequence to a list
    mutations = []

    # Apply the mutations
    for pos, tp in zip(positions, token_positions):             # For each position and token position
        # If the token position is out of bounds, skip
        if tp >= logits.shape[1] - 1:
            continue

        # Get the logits only for masked (/to be mutated) positions and select the top-k residues
        residue_logits = logits[0, tp]  # [vocab]
        top_ids = torch.topk(residue_logits, k=top_k).indices.tolist() # A list of the top-k token ids (integers) with the highest logit score

        # Convert ids to tokens and keep only valid one-letter amino acids
        candidate_aas = []
        for tok_id in top_ids:
            tok = tokenizer.convert_ids_to_tokens(tok_id)       # Convert the token id to a token (string)
            if tok in AMINO_ACIDS and tok != seq_list[pos]:     # If the token is a valid one-letter amino acid and not the same as the original amino acid at the position
                candidate_aas.append(tok)                       # Add the token to the candidates

        # If there are no valid candidates, skip
        if not candidate_aas:
            continue

        # Randomly select a new amino acid from the candidates
        old_aa = seq_list[pos]                                  # Get the old amino acid at the position
        new_aa = rng.choice(candidate_aas)                      # Randomly select a new amino acid from the candidates list
        seq_list[pos] = new_aa                                  # Replace the old amino acid with the new one
        mutations.append(f"{old_aa}{pos+1}{new_aa}")            # Add the mutation to the mutations track list

    mutated_seq = "".join(seq_list)                             # Join the final mutated sequence to a string
    mutation_str = ";".join(mutations)                          # Join all the mutations track strings with a semicolon
    return mutated_seq, mutation_str


def sequence_distance(seq_a: str, seq_b: str) -> int:
    """ Count the differing positions between two protein sequences. """
    # If the lengths are different, add the difference to the distance
    if len(seq_a) != len(seq_b):
        return abs(len(seq_a) - len(seq_b)) + sum(
            a != b for a, b in zip(seq_a, seq_b)
        )
    # If the lengths are equal, just return the number of differing positions
    return sum(a != b for a, b in zip(seq_a, seq_b))


def select_diverse_top_candidates(
    df: pd.DataFrame,
    keep_top_k: int,
    min_distance: int,
) -> pd.DataFrame:
    """ Select top k candidates that are diverse enough from each other.
    Input args:
        df: dataframe of candidates sorted by target score
        keep_top_k: number of candidates to keep
        min_distance: minimum distance between candidates to keep them diverse enough
    Output:
        dataframe of top k candidates that are diverse enough from each other
    """

    selected_rows = []
    selected_sequences = []

    # Iterate over the dataframe and select top k candidates that are diverse enough
    for _, row in df.iterrows():
        seq = row["aa_sequence"]

        # If the sequence is diverse enough, select it
        if all(sequence_distance(seq, prev) >= min_distance for prev in selected_sequences):
            selected_rows.append(row)
            selected_sequences.append(seq)
        # When there are enough diverse candidates, stop
        if len(selected_rows) >= keep_top_k:
            break

    # If there are not enough diverse candidates, just select the first k that aren't already selected (dataframe is already sorted by target score)
    if len(selected_rows) < keep_top_k:

        selected_ids = {row["candidate_id"] for row in selected_rows}   # Set of selected candidate ids
        for _, row in df.iterrows():                                    # Iterate over the dataframe
            if row["candidate_id"] in selected_ids:                     # If the candidate id is already selected
                continue                                                # Skip it
            selected_rows.append(row)                                   # Add the row to the selected rows
            if len(selected_rows) >= keep_top_k:                        # If there are enough selected rows
                break                                                   # Stop
    # Return the selected rows as a dataframe and reset the index
    return pd.DataFrame(selected_rows).reset_index(drop=True)


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
    ap.add_argument("--min_mutations", type=int, default=1)                                 # Number of mutations per candidate (min)
    ap.add_argument("--max_mutations", type=int, default=3)                                 # Number of mutations per candidate (max)
    ap.add_argument("--keep_top_k", type=int, default=10)                                   # Number of top candidates to keep from each optimization round
    ap.add_argument("--proposal_top_k", type=int, default=8)                                # Number of top proposed by the ESM candidates to keep
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--diversity_on", action="store_true")
    ap.add_argument("--diversity_min_distance", type=int, default=8)
    ap.add_argument("--position_selection_strategy", type=str, default="random", choices=["random", "entropy"])     # Position selection strategy: random or entropy
    ap.add_argument("--entropy_guided_position_pool_size", type=int, default=32)                           # Number of high-entropy positions eligible for weighted sampling when using entropy-guided targeting
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
    virus_accession, source_host, protein_id, seed_sequence = load_seed_sequence(
        seed_csv=seed_csv,
        seed_protein_id=args.seed_protein_id
    )

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

    # Load ESM MLM once and reuse it across all rounds
    tokenizer = AutoTokenizer.from_pretrained(args.esm_model, do_lower_case=False)
    esm_mlm = EsmForMaskedLM.from_pretrained(args.esm_model)
    esm_mlm.to(device)
    esm_mlm.eval()

    # Resolve valid amino-acid token ids once and reuse them for entropy-guided mutation targeting
    aa_token_ids = get_amino_acid_token_ids(tokenizer)

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

        # Cache per-parent entropy scores so they are computed only once per parent sequence in each round
        parent_position_scores_cache = {}

        # Generate mutations from the current pool
        for parent_rank, parent in enumerate(current_pool):         # For each parent in the pool
            parent_seq = parent["aa_sequence"]

            # Compute parent-specific position scores once and reuse them across all children of that parent
            if args.position_selection_strategy == "entropy":                                       # If we're using entropy-guided targeting
                if parent_seq not in parent_position_scores_cache:                                  # If the parent sequence is not in the cache
                    parent_position_scores_cache[parent_seq] = compute_position_entropies(          # Compute the position entropies
                        seq=parent_seq,
                        tokenizer=tokenizer,
                        model=esm_mlm,
                        max_aa=args.max_aa,
                        device=device,
                        aa_token_ids=aa_token_ids,
                    )
                parent_position_scores = parent_position_scores_cache[parent_seq]                   # Get the (pre-computed) position entropies
            else:                                                                                   # If we're not using entropy-guided targeting
                parent_position_scores = None                                                       # Set the position entropies to None

            for _ in range(args.candidates_per_round):              # Generate multiple mutated candidates (# candidates_per_round)

                # Generate a random number of mutations to happen 
                n_mut = rng.randint(args.min_mutations, args.max_mutations)
                
                # Choose mutation positions (entropy-guided or random)
                if args.position_selection_strategy == "entropy":
                    positions = choose_mutation_positions_entropy_guided(
                        seq=parent_seq,
                        n_mutations=n_mut,
                        position_scores=parent_position_scores,
                        rng=rng,
                        candidate_pool_size=args.entropy_guided_position_pool_size,
                    )
                else:
                    positions = choose_mutation_positions_randomly(
                        parent_seq,
                        n_mutations=n_mut,
                        rng=rng
                    )

                # Propose a mutated sequence
                mutated_seq, mutation_str = propose_esm_guided_mutation(
                    seq=parent_seq,
                    positions=positions,
                    tokenizer=tokenizer,
                    model=esm_mlm,
                    max_aa=args.max_aa,
                    device=device,
                    top_k=args.proposal_top_k,
                    rng=rng,
                )

                # Skip if the mutation string is empty
                if mutation_str == "":
                    continue

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

        # Check if there are any proposed rows
        if len(proposed_rows) == 0:
            print(f"[WARN] No candidates proposed in round {round_idx}. Stopping early.")
            break

        # Convert proposed rows to dataframe and drop duplicates
        proposed_df = pd.DataFrame(proposed_rows).drop_duplicates(subset=["aa_sequence"]).reset_index(drop=True)

        # Embed all proposed candidates
        embeddings = embed_sequences(
            sequences=proposed_df["aa_sequence"].tolist(),
            tokenizer=tokenizer,
            model=esm_mlm,
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
        proposed_df["n_mutations"] = proposed_df["mutations"].apply(lambda s: len(s.split(";")) if s else 0)

        # Sort proposed candidates by target score descending and number of mutations ascending
        proposed_df = proposed_df.sort_values(
            by=["target_score", "n_mutations"],
            ascending=[False, True]
        ).reset_index(drop=True)

        # Assign rank within the round after sorting by target score and mutation count
        proposed_df["rank_in_round"] = range(1, len(proposed_df) + 1)

        # Save all proposed candidates for each round to csv
        proposed_df.to_csv(out_dir / f"round_{round_idx}_candidates.csv", index=False)

        # Keep track of all proposed candidates across rounds
        all_rows.append(proposed_df)

        # If diversity is enabled, select the most diverse top k candidates
        if args.diversity_on:
            selected_df = select_diverse_top_candidates(
                proposed_df,
                keep_top_k=args.keep_top_k,
                min_distance=args.diversity_min_distance,
            )
        # Else keep only the top k candidates based on target score of this round to seed the next round
        else:
            selected_df = proposed_df.head(args.keep_top_k).copy()

        # Update the current pool by adding the selected candidates
        current_pool = selected_df.to_dict(orient="records")    # Convert the dataframe to a list of dictionaries

    # Raise error if no candidates were generated
    if len(all_rows) == 0:
        raise RuntimeError("No candidate rounds were successfully generated.")

    # Save all rounds together in a single csv
    all_candidates_df = pd.concat(all_rows, axis=0).reset_index(drop=True)
    all_candidates_df.to_csv(out_dir / "all_candidates.csv", index=False)

    # Rank all candidates across all rounds by score and mutation count.
    all_candidates_ranked_df = (
        all_candidates_df
        .sort_values(by=["target_score", "n_mutations"], ascending=[False, True])
        .drop_duplicates(subset=["aa_sequence"])
        .reset_index(drop=True)
    )

    # If diversity is enabled, select a diverse final top-50 set.
    if args.diversity_on:
        final_top_df = select_diverse_top_candidates(
            all_candidates_ranked_df,
            keep_top_k=50,
            min_distance=args.diversity_min_distance,
        )
    else:
        final_top_df = all_candidates_ranked_df.head(50).reset_index(drop=True)

    # Save final top 50 candidates (diverse or not) to csv
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
        "min_mutations": args.min_mutations,
        "max_mutations": args.max_mutations,
        "keep_top_k": args.keep_top_k,
        "proposal_top_k": args.proposal_top_k,
        "seed": args.seed,
        "n_total_candidates": int(len(all_candidates_df)),
        "n_top_candidates_saved": int(len(final_top_df)),
        "diversity_on": args.diversity_on,
        "diversity_min_distance": args.diversity_min_distance,
        "position_selection_strategy": args.position_selection_strategy,
        "entropy_guided_position_pool_size": args.entropy_guided_position_pool_size,
    }
    (out_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2))

    print(f"✅ Saved: {out_dir / 'all_candidates.csv'}")
    print(f"✅ Saved: {out_dir / 'top_candidates.csv'}")
    print(f"✅ Saved: {out_dir / 'run_metadata.json'}")


if __name__ == "__main__":
    main()