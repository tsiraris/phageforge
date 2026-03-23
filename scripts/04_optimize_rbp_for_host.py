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


def choose_mutation_positions(seq: str, n_mutations: int, rng: random.Random) -> list[int]:
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
    aa_token_ids = []
    for aa in AMINO_ACIDS:
        tok_id = tokenizer.convert_tokens_to_ids(aa)                                     # Convert each amino-acid token string to its tokenizer id
        if tok_id is not None and tok_id != tokenizer.unk_token_id:                      # Keep only valid amino-acid ids and ignore unknown-token collisions
            aa_token_ids.append(int(tok_id))
    if len(aa_token_ids) == 0:
        raise ValueError("Could not resolve one-letter amino-acid token ids from tokenizer.")
    return aa_token_ids


def compute_position_entropies_batch(
    sequences: list[str],
    tokenizer,
    model,
    max_aa: int,
    device: str,
    aa_token_ids: list[int],
    batch_size: int,
) -> list[np.ndarray]:
    """
    Compute an uncertainty score (Shannon entropy) for each residue position in a batch of input sequences.

    The entropy is computed from the ESM logits restricted to valid amino-acid tokens.
    Higher entropy means the model considers multiple amino acids plausible at that position.
    """
    all_entropies = []

    # Compute entropies in batches so the current parent pool can be processed with fewer forward passes
    with torch.inference_mode():
        for start in range(0, len(sequences), batch_size):
            batch = sequences[start : start + batch_size]

            # Check that every sequence in the batch fits inside the ESM maximum sequence length
            for seq in batch:
                if len(seq) > max_aa:
                    raise ValueError(f"Sequence length {len(seq)} exceeds max_aa={max_aa}")

            # Tokenize the input sequences batch
            toks = tokenizer(
                batch,
                return_tensors="pt",    # Return PyTorch tensors
                padding=True,
                truncation=True,
                max_length=max_aa,
            ).to(device)

            # Feed the unmasked sequences into the model and inspect their token distributions at each residue position
            outputs = model(input_ids=toks["input_ids"], attention_mask=toks["attention_mask"])
            logits = outputs.logits                                                     # [B, L, vocab]

            # Compute one entropy vector per sequence while ignoring padding and special tokens
            for batch_idx, seq in enumerate(batch):
                seq_len = len(seq)                                                      # Real amino-acid sequence length
                aa_logits = logits[batch_idx, 1 : seq_len + 1, aa_token_ids]            # [seq_len, 20]: Restrict the logits to valid one-letter amino-acid tokens only
                aa_probs = torch.softmax(aa_logits, dim=-1)                             # [seq_len, 20]: Convert the logits to probabilities
                aa_log_probs = torch.log(aa_probs.clamp(min=1e-12))                     # [seq_len, 20]: Stable log-probabilities for the entropy formula
                entropies = -(aa_probs * aa_log_probs).sum(dim=-1)                      # [seq_len]: Shannon entropy per residue position
                all_entropies.append(entropies.detach().cpu().numpy())                  # Save one entropy vector per sequence in CPU memory

    return all_entropies


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
    1. rank sequence positions by entropy
    2. keep only the top uncertain positions equal to candidate_pool_size
    3. sample without replacement from that pool using entropy-weighted probabilities
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
        return choose_mutation_positions(seq=seq, n_mutations=n_mutations, rng=rng)

    # Rank positions by entropy descending and keep a pool of the most uncertain positions equal to candidate_pool_size
    ranked_positions = np.argsort(-position_scores).tolist()
    pool_size = min(max(candidate_pool_size, n_mutations), len(ranked_positions))
    candidate_positions = ranked_positions[:pool_size]

    # Convert candidate scores to positive sampling weights
    candidate_scores = np.array([max(float(position_scores[p]), 1e-8) for p in candidate_positions], dtype=float)   # Ensure non-zero scores
    candidate_scores = candidate_scores / candidate_scores.sum()                                                    # Normalize

    selected_positions = []
    remaining_positions = candidate_positions.copy()
    remaining_scores = candidate_scores.copy()

    # Sample from the candidate pool without replacement using entropy-weighted probabilities
    while len(selected_positions) < n_mutations and len(remaining_positions) > 0:       # As long as there are positions left in the pool
        chosen_idx = rng.choices(                                                       # Sample without replacement
            population=list(range(len(remaining_positions))),
            weights=remaining_scores.tolist(),
            k=1,
        )[0]                                                                            # Get the index of the chosen position
        selected_positions.append(remaining_positions.pop(chosen_idx))                  # Add the chosen position to the selected positions
        remaining_scores = np.delete(remaining_scores, chosen_idx)                      # Remove the chosen position from the pool

        if len(remaining_scores) > 0:                                                   # If there are still positions left in the pool
            remaining_scores = remaining_scores / remaining_scores.sum()                # Normalize the remaining scores

    return sorted(selected_positions)


def build_masked_candidate_batch(
    parent_sequence: str,
    positions_list: list[list[int]],
    tokenizer,
    max_aa: int,
    device: str,
) -> tuple[dict[str, torch.Tensor], list[list[int]]]:
    """
    Build a batch of masked tokenized sequences for a single parent sequence and many mutation-position sets.

    Returns:
        token_batch = tokenized batch of sequences with the requested positions masked
        valid_positions_list = mutation positions aligned with the returned batch rows
    """
    # Ensure the parent sequence is not too long
    if len(parent_sequence) > max_aa:
        raise ValueError(f"Sequence length {len(parent_sequence)} exceeds max_aa={max_aa}")

    # Tokenize one repeated parent sequence for each candidate so the masked proposals can be scored in one batched forward pass
    toks = tokenizer(
        [parent_sequence] * len(positions_list),        # Repeat the parent sequence as many times as there are candidates
        return_tensors="pt",                            # Return PyTorch tensors
        padding=True,
        truncation=True,
        max_length=max_aa,
    ).to(device)

    # Clone the tokenized batch to avoid modifying the original tokens
    input_ids = toks["input_ids"].clone()       
    
    # Get the specific integer (ID) of the <mask> token from the tokenizer vocabulary        
    mask_token_id = tokenizer.mask_token_id             
    if mask_token_id is None:
        raise ValueError("Tokenizer does not expose a mask token id.")

    valid_positions_list = []

    # Apply the <mask> token at the chosen positions for every candidate in the batch
    for row_idx, positions in enumerate(positions_list):                    # Loop over the candidates
        valid_positions = []

        for pos in positions:                                               # Loop over the mutation positions
            token_position = pos + 1                                        # Convert sequence positions to token positions by shifting them by 1 because ESM tokenizers add special tokens at the beginning
            if token_position >= input_ids.shape[1] - 1:                    # Skip if the token position is out of bounds
                continue
            input_ids[row_idx, token_position] = mask_token_id              # Mask the requested residue position for this candidate row
            valid_positions.append(pos)                                     # Add the mutation position  to the list

        valid_positions_list.append(valid_positions)                        # Add the list of valid mutation positions to the batch

    toks["input_ids"] = input_ids                                           # Replace the input ids in the token dictionary with the masked version
    return toks, valid_positions_list                                       # Return the tokenized batch and the list of valid mutation positions


def sample_mutations_from_batched_logits(
    parent_sequence: str,
    valid_positions_list: list[list[int]],
    logits: torch.Tensor,
    tokenizer,
    top_k: int,
    rng: random.Random,
) -> tuple[list[str], list[str]]:
    """
    Convert batched masked-token logits into mutated sequences and mutation strings.

    Returns:
        mutated_sequences = list of resulting mutated sequences
        mutation_strings = list of mutation-string annotations
    """

    mutated_sequences = []
    mutation_strings = []

    # Convert each batch row of logits to a mutated sequence and a mutation string
    for row_idx, positions in enumerate(valid_positions_list):
        seq_list = list(parent_sequence)                                                # Convert the protein sequence to a list so residues can be replaced in place
        mutations = []

        for pos in positions:                                                           # For each chosen mutation position
            token_position = pos + 1                                                    # Convert sequence position to token position by shifting it by 1
            if token_position >= logits.shape[1] - 1:                                   # Skip if the token position is out of bounds
                continue

            residue_logits = logits[row_idx, token_position]                            # [vocab]: Get the logits only for the masked (/to be mutated) position
            top_ids = torch.topk(residue_logits, k=top_k).indices.tolist()              # A list of the top-k token ids (integers) with the highest logit score

            candidate_aas = []
            for tok_id in top_ids:
                tok = tokenizer.convert_ids_to_tokens(tok_id)                           # Convert the token id to a token (string)
                if tok in AMINO_ACIDS and tok != seq_list[pos]:                         # If the token is a valid one-letter amino acid and not the same as the original amino acid at the position
                    candidate_aas.append(tok)                                           # Add the token to the candidate list

            if not candidate_aas:                                                       # If there are no valid candidates, skip this position
                continue

            old_aa = seq_list[pos]                                                      # Get the old amino acid at the position
            new_aa = rng.choice(candidate_aas)                                          # Randomly select a new amino acid from the candidate list
            seq_list[pos] = new_aa                                                      # Replace the old amino acid with the new one
            mutations.append(f"{old_aa}{pos+1}{new_aa}")                                # Add the mutation to the mutation-track list

        mutated_sequences.append("".join(seq_list))                                     # Join the final mutated sequence to a string
        mutation_strings.append(";".join(mutations))                                    # Join all the mutation strings with a semicolon

    return mutated_sequences, mutation_strings


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
        df: dataframe of candidates sorted by selection score
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

    # If there are not enough diverse candidates, just select the first k that aren't already selected (dataframe is already sorted by selection score)
    if len(selected_rows) < keep_top_k:

        selected_ids = {row["candidate_id"] for row in selected_rows}                   # Set of selected candidate ids
        for _, row in df.iterrows():                                                    # Iterate over the dataframe
            if row["candidate_id"] in selected_ids:                                     # If the candidate id is already selected
                continue                                                                # Skip it
            selected_rows.append(row)                                                   # Add the row to the selected rows
            if len(selected_rows) >= keep_top_k:                                        # If there are enough selected rows
                break                                                                   # Stop

    # Return the selected rows as a dataframe and reset the index
    return pd.DataFrame(selected_rows).reset_index(drop=True)


def normalize_rows(x: np.ndarray) -> np.ndarray:
    """ L2-normalize each embedding row for cosine-similarity computations. """
    return x / np.clip(np.linalg.norm(x, axis=1, keepdims=True), 1e-12, None)       # Clip the norm to avoid division by zero


def compute_seed_similarity_columns(
    candidate_embeddings: np.ndarray,
    seed_embedding_normalized: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute seed cosine similarity and novelty distance for a batch of candidate embeddings.

    Returns:
        seed_cosine_similarity = cosine similarity between each candidate and the seed embedding
        seed_novelty_distance = 1 - cosine similarity
    """
    candidate_embeddings_normalized = normalize_rows(candidate_embeddings)                  # Normalize candidate embeddings once so cosine similarity becomes a dot product
    seed_cosine_similarity = candidate_embeddings_normalized @ seed_embedding_normalized    # [N]: Cosine similarity between every candidate and the seed embedding
    seed_novelty_distance = 1.0 - seed_cosine_similarity                                    # [N]: Novelty distance away from the seed embedding
    return seed_cosine_similarity, seed_novelty_distance


def compute_novelty_penalty_weight(
    round_idx: int,
    total_rounds: int,
    novelty_penalty_lambda: float,
    novelty_penalty_schedule: str,
) -> float:
    """
    Compute the novelty-penalty weight used in selection_score for the current round.
    """
    if novelty_penalty_schedule == "constant":                                           # Keep the same novelty penalty in every round for easier interpretability
        return novelty_penalty_lambda

    if total_rounds <= 1:                                                                # Avoid division by zero if the run has only one round
        return novelty_penalty_lambda

    round_fraction = (round_idx - 1) / max(total_rounds - 1, 1)                          # Convert the round number to a fraction between 0 and 1
    return novelty_penalty_lambda * round_fraction                                       # Increase the novelty penalty gradually across rounds with a simple linear schedule


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed_csv", type=str, default="data/processed/rbp_dataset_eskapee_strict.csv")
    ap.add_argument("--seed_protein_id", type=str, required=True)
    ap.add_argument("--target_host", type=str, required=True)                               # The host genus to optimize for
    ap.add_argument("--model_path", type=str, required=True)
    ap.add_argument("--label_classes_path", type=str, required=True)
    ap.add_argument("--esm_model", type=str, default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--batch_size", type=int, default=8)                                    # Batch size for sequence embedding and parent-entropy computation
    ap.add_argument("--proposal_batch_size", type=int, default=16)                          # Batch size for masked mutation-proposal forward passes from the same parent sequence
    ap.add_argument("--max_aa", type=int, default=1022)
    ap.add_argument("--rounds", type=int, default=5)                                        # Number of optimization rounds
    ap.add_argument("--candidates_per_round", type=int, default=64)                         # Number of candidate sequences per round
    ap.add_argument("--min_mutations", type=int, default=1)                                 # Number of mutations per candidate (min)
    ap.add_argument("--max_mutations", type=int, default=3)                                 # Number of mutations per candidate (max)
    ap.add_argument("--keep_top_k", type=int, default=10)                                   # Number of top candidates to keep from each optimization round
    ap.add_argument("--proposal_top_k", type=int, default=8)                                # Number of top proposed by the ESM candidates to keep
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--diversity_on", action="store_true")
    ap.add_argument("--diversity_min_distance", type=int, default=20)                       # Minimum sequence distance between retained candidates when diversity-aware selection is enabled
    ap.add_argument("--position_selection_strategy", type=str, default="random", choices=["random", "entropy"])
    ap.add_argument("--position_pool_size", type=int, default=32)                           # Number of high-entropy positions eligible for weighted sampling when using entropy-guided targeting
    ap.add_argument("--novelty_penalty_lambda", type=float, default=0.0)                    # Strength of the penalty for remaining too close to the seed embedding during candidate selection
    ap.add_argument("--novelty_penalty_schedule", type=str, default="constant", choices=["constant", "linear"])
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

    # Embed the seed sequence once and normalize it so repeated seed-similarity calculations are cheap
    seed_embedding = embed_sequences(
        sequences=[seed_sequence],
        tokenizer=tokenizer,
        model=esm_mlm,
        batch_size=1,
        max_aa=args.max_aa,
        device=device,
    )
    seed_embedding_normalized = normalize_rows(seed_embedding)[0]         # [H]: One normalized seed embedding reused throughout the run

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
    # In each round, we mutate the parents (= sequences currently in the pool), score the children, and keep only the top k (most diverse or not) children to become the parents for the next round.
    for round_idx in range(1, args.rounds + 1):                     # Iterate over the optimization rounds (= number of generations we run this process for)
        proposed_rows = []
        candidate_counter = 0

        # Pre-compute parent entropy scores in batch once per round so all children of the same parent reuse the same position scores
        if args.position_selection_strategy == "entropy":
            parent_sequences = [parent["aa_sequence"] for parent in current_pool]         # List the current parent sequences
            parent_position_scores_list = compute_position_entropies_batch(               # Compute position entropies for each amino acid position for all parent sequences in parallel
                sequences=parent_sequences,
                tokenizer=tokenizer,
                model=esm_mlm,
                max_aa=args.max_aa,
                device=device,
                aa_token_ids=aa_token_ids,
                batch_size=args.batch_size,
            )
            parent_position_scores_cache = {                                              # Cache is a dict of aa parent position entropies; shape: {parent_seq: {pos: entropy}}
                seq: scores for seq, scores in zip(parent_sequences, parent_position_scores_list)
            }
        else:
            parent_position_scores_cache = {}                                             # If position selection strategy is not entropy, cache dict is empty 

        # Generate mutations for each parent in the current pool
        for parent_rank, parent in enumerate(current_pool):         # For each parent in the pool
            parent_seq = parent["aa_sequence"]
            positions_requests = []
            n_mut_list = []

            # Build the requested mutation-position sets (= which positions to mutate) for this parent before running any proposal-model forward passes
            for _ in range(args.candidates_per_round):                                      # Loop to generate multiple mutated candidates (# candidates_per_round)
                n_mut = rng.randint(args.min_mutations, args.max_mutations)                 # Select randomly the number of mutations between min_mutations and max_mutations
                n_mut_list.append(n_mut)                                                    # Save the number of mutations for this candidate

                if args.position_selection_strategy == "entropy":
                    parent_position_scores = parent_position_scores_cache[parent_seq]       # Reuse the entropy scores already computed for this parent sequence (from cache)
                    positions = choose_mutation_positions_entropy_guided(                   # Choose which sequence positions to mutate using entropy
                        seq=parent_seq,
                        n_mutations=n_mut,
                        position_scores=parent_position_scores,
                        rng=rng,
                        candidate_pool_size=args.position_pool_size,
                    )
                else:                                                                       # If position selection strategy is not entropy, choose positions randomly                        
                    positions = choose_mutation_positions(
                        parent_seq,
                        n_mutations=n_mut,
                        rng=rng
                    )

                positions_requests.append(positions)                                        # Save the requested mutation positions in a list so they can be proposed in batched forward passes

            # Run the masked mutation proposals for this parent in batches
            for start in range(0, len(positions_requests), args.proposal_batch_size):                           # For each batch of mutation-position sets
                batch_positions_requests = positions_requests[start : start + args.proposal_batch_size]         # Get the mutation-position sets for this batch

                toks, valid_positions_list = build_masked_candidate_batch(                                      # Build a batch of masked tokenized sequences
                    parent_sequence=parent_seq,                                                                 # Each row in the batch is a single parent sequence with one of the requested mutation-positions masked
                    positions_list=batch_positions_requests,
                    tokenizer=tokenizer,
                    max_aa=args.max_aa,
                    device=device,
                )                                                                                               # toks: [B=proposal_batch_size, max_aa] tensor matrix of differently masked tokenized sequences, valid_positions_list: mutation positions tracking list aligned with the returned batch rows

                # Run the proposal-model forward pass for the whole parent-specific batch and get the masked-token proposal logits
                with torch.inference_mode():
                    outputs = esm_mlm(input_ids=toks["input_ids"], attention_mask=toks["attention_mask"])       
                    logits = outputs.logits                                                                     # [B, L, vocab]: Masked-token proposal logits for the whole parent-specific batch

                mutated_sequences, mutation_strings = sample_mutations_from_batched_logits(                     # Convert each of the batch rows of masked-token logits into mutated sequences and mutation strings
                    parent_sequence=parent_seq,
                    valid_positions_list=valid_positions_list,
                    logits=logits,
                    tokenizer=tokenizer,
                    top_k=args.proposal_top_k,
                    rng=rng,
                )

                for mutated_seq, mutation_str in zip(mutated_sequences, mutation_strings):
                    if mutation_str == "":                                              # Skip if the mutation string is empty
                        continue

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
        probs = clf.predict_proba(embeddings)                                           # Predicted probabilities for all the candidates
        pred_idx = probs.argmax(axis=1)                                                 # Index of the predicted label (host)
        pred_label = [label_classes[i] for i in pred_idx]                               # Predicted label (host)
        target_scores = probs[:, target_idx]                                            # Predicted score for the target host

        # Compute seed similarity columns once for all candidates in the round so novelty-aware selection is cheap
        seed_cosine_similarity, seed_novelty_distance = compute_seed_similarity_columns(
            candidate_embeddings=embeddings,
            seed_embedding_normalized=seed_embedding_normalized,
        )
        novelty_penalty_weight = compute_novelty_penalty_weight(
            round_idx=round_idx,
            total_rounds=args.rounds,
            novelty_penalty_lambda=args.novelty_penalty_lambda,
            novelty_penalty_schedule=args.novelty_penalty_schedule,
        )

        proposed_df["pred_label"] = pred_label                                          # Add predicted labels column to the dataframe
        proposed_df["target_host"] = args.target_host                                   # Add target host column
        proposed_df["target_score"] = target_scores                                     # Add target score column
        proposed_df["seed_cosine_similarity"] = seed_cosine_similarity                  # Add seed cosine similarity column so each candidate records how close it remains to the original seed embedding
        proposed_df["seed_novelty_distance"] = seed_novelty_distance                    # Add seed novelty distance column so each candidate records how far it moved away from the original seed embedding
        proposed_df["selection_score"] = proposed_df["target_score"] - (
            novelty_penalty_weight * proposed_df["seed_cosine_similarity"]
        )                                                                               # Add the final selection score column that optionally penalizes candidates for staying too close to the seed embedding

        # Count number of substitutions we introduced this round (through the number of semicolons in the mutation string) and add them as a new column to the dataframe
        proposed_df["n_mutations"] = proposed_df["mutations"].apply(lambda s: len(s.split(";")) if s else 0)

        # Sort proposed candidates by selection score descending, target score descending, and number of mutations ascending
        proposed_df = proposed_df.sort_values(
            by=["selection_score", "target_score", "n_mutations"],
            ascending=[False, False, True]
        ).reset_index(drop=True)

        # Assign rank within the round after sorting by selection score, target score, and mutation count
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
        # Else keep only the top k candidates based on selection score of this round to seed the next round
        else:
            selected_df = proposed_df.head(args.keep_top_k).copy()

        # Update the current pool by adding the selected candidates
        current_pool = selected_df.to_dict(orient="records")

    # Raise error if no candidates were generated
    if len(all_rows) == 0:
        raise RuntimeError("No candidate rounds were successfully generated.")

    # Save all rounds together in a single csv
    all_candidates_df = pd.concat(all_rows, axis=0).reset_index(drop=True)
    all_candidates_df.to_csv(out_dir / "all_candidates.csv", index=False)

    # Rank all candidates across all rounds by selection score, target score, and mutation count
    all_candidates_ranked_df = (
        all_candidates_df
        .sort_values(by=["selection_score", "target_score", "n_mutations"], ascending=[False, False, True])
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
        "position_pool_size": args.position_pool_size,
        "proposal_batch_size": args.proposal_batch_size,
        "novelty_penalty_lambda": args.novelty_penalty_lambda,
        "novelty_penalty_schedule": args.novelty_penalty_schedule,
    }
    (out_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2))

    print(f"✅ Saved: {out_dir / 'all_candidates.csv'}")
    print(f"✅ Saved: {out_dir / 'top_candidates.csv'}")
    print(f"✅ Saved: {out_dir / 'run_metadata.json'}")


if __name__ == "__main__":
    main()
