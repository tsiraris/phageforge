"""
Phase A / Step 06b: Guided protein design under multiple constraints
===================

This script performs family-constrained, retrieval-guided host retargeting for one Phase A step.

Compared with the original Stage 04 optimizer, this script adds four key ideas:
1. it mutates only inside a host-ladder window derived from corrected hotspot evidence,
2. it samples mutation positions using target-specific priors,
3. it scores candidates with both host-probe score and family/target manifold terms,
4. it supports ladder mode where the seed can be either the canonical seed protein or a FASTA
   sequence produced by the previous ladder step.

"""
from __future__ import annotations
import argparse                                                                        # Parse command-line arguments.
import importlib.util                                                                  # Dynamically load helper functions from the existing Stage 04 script.
import json                                                                            # Read the Phase A plan and write run metadata.
import random                                                                          # Control reproducible random sampling.
import tempfile                                                                        # Create temporary files when the ladder seed comes from a FASTA input.
from pathlib import Path                                                               # Work with filesystem paths robustly.
from typing import Dict, List, Tuple                                                   # Improve readability of function signatures.
import joblib                                                                          # Load the trained host probe.
import numpy as np                                                                     # Numerical computations for embeddings and score terms.
import pandas as pd                                                                    # Candidate tables and CSV outputs.
import torch                                                                           # Load and run the ESM model.
from transformers import AutoTokenizer, EsmForMaskedLM                                 # Tokenizer and masked language model for mutation proposals.

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")                                             # Standard 20 amino acids accepted in generated proteins.


def load_stage04_helpers(stage04_script: Path):
    """Dynamically load the original Stage 04 optimizer module (source code) to reuse its utilities."""
    spec = importlib.util.spec_from_file_location("stage04_optimize_module", stage04_script)  # Create a module spec from the Stage 04 file path.
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load Stage 04 script from: {stage04_script}")
    module = importlib.util.module_from_spec(spec)                                      # Instantiate an importable module object from the spec.
    spec.loader.exec_module(module)                                                     # Execute the Stage 04 source code into the module namespace.
    return module


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for one constrained Phase A optimization step."""
    ap = argparse.ArgumentParser(description="Run one family-constrained Phase A host-retargeting step.")
    ap.add_argument("--phaseA_plan_json", type=str, required=True, help="Plan file produced by 06a_select_phaseA_family.py.")
    ap.add_argument("--stage04_script", type=str, default="scripts/04_optimize_rbp_for_host.py", help="Path to the existing Stage 04 optimizer script to reuse helper functions from.")
    ap.add_argument("--seed_csv", type=str, default="data/processed/rbp_dataset_eskapee_strict.csv", help="Strict RBP dataset used to recover the canonical seed when no FASTA seed is provided.")
    ap.add_argument("--seed_protein_id", type=str, default="", help="Seed protein ID for the first ladder step. Leave empty to use the canonical seed from the plan.")
    ap.add_argument("--seed_fasta", type=str, default="", help="Optional FASTA file containing the sequence that should seed the current ladder step.")
    ap.add_argument("--seed_label", type=str, default="", help="Optional human-readable label for FASTA-based ladder seeds.")
    ap.add_argument("--target_host", type=str, required=True, help="Target host genus for this ladder step.")
    ap.add_argument("--model_path", type=str, required=True, help="Path to the trained host probe joblib file.")
    ap.add_argument("--label_classes_path", type=str, required=True, help="JSON file containing the probe class label order.")
    ap.add_argument("--esm_model", type=str, default="facebook/esm2_t33_650M_UR50D", help="ESM masked language model checkpoint used for proposals and embeddings.")
    ap.add_argument("--embedding_pt", type=str, required=True, help="Cached embedding tensor used to recover family and target anchors.")
    ap.add_argument("--embedding_index_csv", type=str, required=True, help="Embedding index CSV used to align protein IDs to embedding rows.")
    ap.add_argument("--rounds", type=int, default=5, help="Number of optimization rounds.")
    ap.add_argument("--candidates_per_round", type=int, default=64, help="Total number of candidates to propose per round.")
    ap.add_argument("--min_mutations", type=int, default=2, help="Minimum number of mutations per candidate.")
    ap.add_argument("--max_mutations", type=int, default=6, help="Maximum number of mutations per candidate.")
    ap.add_argument("--keep_top_k", type=int, default=10, help="How many candidates seed the next round.")
    ap.add_argument("--proposal_top_k", type=int, default=8, help="Per-position amino-acid proposal top-k sampled from the MLM logits.")
    ap.add_argument("--batch_size", type=int, default=8, help="Batch size for candidate embedding forward passes.")
    ap.add_argument("--proposal_batch_size", type=int, default=16, help="Batch size for masked mutation-proposal passes.")
    ap.add_argument("--max_aa", type=int, default=1022, help="Maximum amino-acid sequence length supported by the ESM checkpoint.")
    ap.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    ap.add_argument("--diversity_min_distance", type=int, default=8, help="Minimum sequence distance between retained parents when diversity selection is applied.")
    ap.add_argument("--lambda_family", type=float, default=0.20, help="Weight for staying close to the chosen family centroid.")
    ap.add_argument("--lambda_target_anchor", type=float, default=0.15, help="Weight for moving toward the target-host centroid inside the family when available.")
    ap.add_argument("--lambda_seed_novelty", type=float, default=0.05, help="Penalty weight for staying too close to the immediate seed embedding.")
    ap.add_argument("--lambda_mutation_penalty", type=float, default=0.03, help="Penalty weight per mutation to discourage gratuitous drift.")
    ap.add_argument("--use_entropy", action="store_true", help="If set, combine target priors with Stage 04 entropy-guided position selection.")
    ap.add_argument("--out_dir", type=str, required=True, help="Directory where this ladder-step run will be written.")
    return ap.parse_args()


def set_seed(seed: int):
    """Set Python / NumPy / PyTorch seeds for reproducibility."""
    random.seed(seed)                                                                   # Seed Python's random module.
    np.random.seed(seed)                                                                # Seed NumPy's random number generator.
    torch.manual_seed(seed)                                                             # Seed PyTorch on CPU.
    if torch.cuda.is_available():                                                       # Seed CUDA as well when a GPU is available.
        torch.cuda.manual_seed_all(seed)


def read_fasta_sequence(seed_fasta: Path) -> Tuple[str, str]:
    """Read the first FASTA record and return `(label, sequence)`."""
    header = None                                                                       # Will store the first FASTA header line without the leading '>' (e.g., ">Candidate_r4_c17_Enterobacter_Winner").
    seq_chunks: List[str] = []                                                          # Empty list to collect all the sequence lines in case the sequence is broken across multiple lines.
    # Open the FASTA file, recording only the first record's header and sequence.
    with open(seed_fasta, "r", encoding="utf-8") as handle:                             # Open the FASTA file using explicit UTF-8 encoding.
        for line in handle:
            line = line.strip()                                                         # Remove whitespace and trailing newlines.
            if not line:
                continue                                                                # Skip blank lines.
            if line.startswith(">"):
                if header is not None:
                    break                                                               # Stop after the first FASTA record because ladder seeds are single-sequence files.
                header = line[1:].strip() or seed_fasta.stem                            # Use the FASTA header, or fall back to the filename stem.
            else:
                seq_chunks.append(line.upper())                                         # Accumulate sequence lines in uppercase.
    if header is None or len(seq_chunks) == 0:
        raise ValueError(f"Could not read a valid FASTA record from: {seed_fasta}")
    return header, "".join(seq_chunks)                                                  # Return the header and sequence for the first FASTA record only.


def l2_normalize(x: np.ndarray) -> np.ndarray:
    """Normalize rows to unit norm so cosine similarity becomes a dot product afterward."""
    norms = np.linalg.norm(x, axis=1, keepdims=True)                                    # Compute one L2 norm per row.
    norms = np.clip(norms, 1e-12, None)                                                 # Protect against division by zero.
    return x / norms


def load_embedding_bank(embedding_pt: Path, embedding_index_csv: Path) -> pd.DataFrame:
    """Load cached embeddings, attach the embedding matrix and return the entire index table."""
    emb = torch.load(embedding_pt, map_location="cpu")                                  # Load the cached embeddings onto CPU.
    if isinstance(emb, torch.Tensor):                                                   # Convert tensors into NumPy arrays for fast cosine math.
        emb_array = emb.detach().cpu().numpy()
    else:
        emb_array = np.asarray(emb)
    index_df = pd.read_csv(embedding_index_csv)                                         # Load the embedding metadata index.
    row_col = "row_id" if "row_id" in index_df.columns else "idx"                       # Support both index naming conventions.
    index_df["embedding_row"] = index_df[row_col].astype(int)                           # Normalize the row pointer into one canonical column.
    index_df["protein_id"] = index_df["protein_id"].astype(str)                         # Normalize protein IDs for joins.
    index_df.attrs["embedding_array"] = emb_array                                       # Attach the raw embedding array to the DataFrame metadata.
    return index_df


def build_family_runtime_objects(phaseA_plan: dict, embedding_index_df: pd.DataFrame) -> dict:
    """
    Recover runtime family/target centroid objects from the saved Phase A plan.
    Input: Phase A plan JSON and embedding index DataFrame.
    Output: Dictionary of family metadata, embeddings, and family/target centroids. 
    """
    # Recover family embeddings and centroid from the Phase A plan.
    emb_array = embedding_index_df.attrs["embedding_array"]                                 # Recover the raw embedding matrix.
    family_ids = set(str(x) for x in phaseA_plan["family_member_ids"])                      # Collect the family protein IDs defined in the Phase A plan.
    family_meta = embedding_index_df[embedding_index_df["protein_id"].isin(family_ids)].copy()  # Keep only cached embeddings belonging to the chosen family.
    family_emb = emb_array[family_meta["embedding_row"].to_numpy()]                         # Gather family embeddings aligned to the family metadata rows.
    family_emb = l2_normalize(np.asarray(family_emb, dtype=np.float32))                     # Normalize family embeddings for cosine operations.
    family_centroid = np.asarray(phaseA_plan["family_centroid"], dtype=np.float32)          # Load the family centroid saved by Step 06a.
    family_centroid = family_centroid / max(np.linalg.norm(family_centroid), 1e-12)         # Re-normalize defensively in case of float serialization drift.

    # Recover target-host centroids from the Phase A plan.
    target_centroids = {}                                                                   # Store one normalized centroid per target host when available.
    for target, centroid in phaseA_plan.get("target_reference_centroids", {}).items():      # Iterate over the saved target-host centroids from the plan.
        if centroid is None:
            target_centroids[target] = None                                                 # Some targets may have no in-family anchor and should be handled gracefully.
            continue
        arr = np.asarray(centroid, dtype=np.float32)                                        # Recover the centroid as a NumPy vector.
        target_centroids[target] = arr / max(np.linalg.norm(arr), 1e-12)                    # Normalize it again defensively.

    return {
        "family_meta": family_meta,
        "family_emb": family_emb,
        "family_centroid": family_centroid,
        "target_centroids": target_centroids,
    }                                                                                       # Return a compact dictionary of runtime family objects.


def choose_positions_with_target_prior(
    seq: str,
    n_mutations: int,
    allowed_positions: List[int],
    target_prior: np.ndarray,
    rng: random.Random,
    entropy_scores: np.ndarray | None = None,
) -> List[int]:
    """
    Choose mutation positions inside the allowed window using target-specific position priors (probabilities).

    If entropy scores are provided, they are multiplied by the target prior so the selection favors
    positions that are both biologically recurrent and currently uncertain under the language model.
    
    Returns: Sorted list of mutation positions that were chosen.
    """
    # If n_mutations is not zero, and the length of allowed_positions is not zero, store the allowed positions and number of mutations. 
    if n_mutations <= 0:
        return []                                                                           # No mutations requested means no positions to choose.
    allowed_positions = [p for p in allowed_positions if 0 <= p < len(seq)]                 # Clamp the allowed set to the actual sequence length.
    if len(allowed_positions) == 0:
        return []                                                                           # No allowed positions means no proposals can be made.
    n_mutations = min(n_mutations, len(allowed_positions))                                  # Never ask for more mutations than available allowed positions.

    # Construct a sampling distribution from the target-specific positional prior and optionally the entropy scores.
    weights = np.asarray([float(target_prior[p]) for p in allowed_positions], dtype=float)  # Start from the target-specific positional prior.
    weights = np.clip(weights, 1e-12, None)                                                 # Ensure numerical stability and non-zero sampling mass.
    if entropy_scores is not None and len(entropy_scores) == len(seq):                      # Optionally combine the target prior with entropy-guided uncertainty.
        entropy_part = np.asarray([float(max(entropy_scores[p], 1e-12)) for p in allowed_positions], dtype=float)
        weights = weights * entropy_part                                                    # Favor positions that are both target-recurrent and uncertain.
    weights = weights / weights.sum()                                                       # Normalize into a sampling distribution.

    # Sample n_mutations positions from the sampling distribution.
    selected: List[int] = []                                                                # Store the final chosen positions here.
    available_positions = allowed_positions.copy()                                          # Work on a mutable copy for sampling without replacement.
    available_weights = weights.copy()                                                      # Keep weights aligned to the mutable position list.
    while len(selected) < n_mutations and len(available_positions) > 0:                     # While we have not yet selected the requested number of positions...
        chosen_idx = rng.choices(range(len(available_positions)), weights=available_weights.tolist(), k=1)[0]  # Sample one weighted index.
        selected.append(available_positions.pop(chosen_idx))                                # Move the chosen position into the selected set.
        available_weights = np.delete(available_weights, chosen_idx)                        # Remove the sampled weight so we sample without replacement.
        if len(available_weights) > 0:
            available_weights = available_weights / available_weights.sum()                 # Re-normalize remaining weights.
    return sorted(selected)                                                                 # Return the selected positions in sorted order.


def compute_family_similarity_terms(candidate_embeddings: np.ndarray, runtime: dict) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the cosine similarity between a batch of candidate embeddings and the family centroid.
    Returns: L2-normalized candidate embeddings and candidate-family cosine similarities.
    """
    cand = l2_normalize(np.asarray(candidate_embeddings, dtype=np.float32))              # Normalize candidate embeddings once for cosine computations.
    family_centroid = runtime["family_centroid"].reshape(-1)                             # Recover the normalized family centroid.
    family_cos = cand @ family_centroid                                                  # Cosine similarity to the local family manifold center.
    return cand, family_cos                                                              # Return normalized embeddings because later terms reuse them.


def main():
    # Parse command-line arguments, seed all relevant RNGs, and create the output directory.
    args = parse_args()                                                                  # Parse command-line arguments.
    set_seed(args.seed)                                                                  # Seed all relevant RNGs for reproducibility.
    rng = random.Random(args.seed)                                                       # Create an explicit local RNG for sampling operations.

    out_dir = Path(args.out_dir)                                                         # Convert the run output path into a Path object.
    out_dir.mkdir(parents=True, exist_ok=True)                                           # Ensure the run output directory exists.

    # Load the Phase A plan produced by Step 06a, load the Stage 04 helper functions, 
    # and recover the seed metadata (from seed_fasta or phaseA_JSON), label classes and target index. 
    with open(args.phaseA_plan_json, "r", encoding="utf-8") as handle:                   # Read the Phase A plan produced by Step 06a.
        phaseA_plan = json.load(handle)

    stage04 = load_stage04_helpers(Path(args.stage04_script))                            # Load the Stage 04 helper functions dynamically from the repo.

    if args.seed_fasta:                                                                  # Ladder follow-up steps can provide a FASTA seed instead of a protein ID.
        fasta_label, seed_sequence = read_fasta_sequence(Path(args.seed_fasta))          # Recover the seed sequence from the FASTA file.
        seed_protein_id = args.seed_label or fasta_label                                 # Use an explicit seed label if one was provided.
        virus_accession = f"LADDER::{seed_protein_id}"                                   # FASTA-based seeds do not correspond to a canonical accession, so we create a readable synthetic label.
        source_host = "LADDER_STEP"                                                      # Mark the source host generically because the sequence may already be retargeted.
    else:
        canonical_seed = phaseA_plan["canonical_seed"]                                   # Recover the canonical seed metadata from the plan.
        effective_seed_protein_id = args.seed_protein_id or canonical_seed["seed_protein_id"]  # Default to the canonical seed unless the user overrides it.
        virus_accession, source_host, seed_protein_id, seed_sequence = stage04.load_seed_sequence(  # Reuse the original Stage 04 helper to load the seed sequence from the strict CSV.
            seed_csv=Path(args.seed_csv),
            seed_protein_id=effective_seed_protein_id,
        )

    if len(seed_sequence) > args.max_aa:
        raise ValueError(f"Seed sequence length {len(seed_sequence)} exceeds max_aa={args.max_aa}")  # Fail early if the seed is too long for the ESM checkpoint.

    with open(args.label_classes_path, "r", encoding="utf-8") as handle:                 # Load the host-class label order matching the trained probe.
        label_classes = json.load(handle)
    if args.target_host not in label_classes:
        raise ValueError(f"Target host '{args.target_host}' not found in label classes: {label_classes}")  # Ensure the target exists in the trained probe.
    target_idx = label_classes.index(args.target_host)                                   # Recover the classifier column corresponding to the target host.

    # 
    device = "cuda" if torch.cuda.is_available() else "cpu"                              # Use GPU automatically when available because Step 06b is the heavy part.
    print(f"Using device: {device}")                                                     # Surface the actual runtime device so the user can sanity-check the environment.

    tokenizer = AutoTokenizer.from_pretrained(args.esm_model)                            # Load the tokenizer that matches the chosen ESM checkpoint.
    esm_mlm = EsmForMaskedLM.from_pretrained(args.esm_model).to(device).eval()           # Load the masked language model used for both proposals and embeddings.
    clf = joblib.load(args.model_path)                                                   # Load the trained host probe.
    aa_token_ids = stage04.get_amino_acid_token_ids(tokenizer)                           # Reuse the Stage 04 helper that maps amino-acid tokens to tokenizer IDs.

    # Load the reference embeddings and the index table, family and target centroids, 
    # the target-specific positional prior vector, and the allowed mutation window.
    embedding_index_df = load_embedding_bank(Path(args.embedding_pt), Path(args.embedding_index_csv))   # Load cached reference embeddings and the index table.
    runtime = build_family_runtime_objects(phaseA_plan, embedding_index_df)                             # Recover family and target centroids defined in the Phase A plan.
    target_prior = np.asarray(phaseA_plan["target_position_priors"][args.target_host], dtype=float)     # Load the target-specific positional prior vector.
    allowed_positions = [int(x) for x in phaseA_plan["mutation_window_positions_0based"]]               # Recover the allowed mutation window from the plan.
    target_centroid = runtime["target_centroids"].get(args.target_host, None)                           # Recover the optional target-host centroid inside the family.
    
    # Embed the seed sequence and normalize it for cosine similarity.
    seed_emb = stage04.embed_sequences(                                                  # Embed the immediate seed so novelty penalties are defined relative to the current ladder step.
        sequences=[seed_sequence],
        tokenizer=tokenizer,
        model=esm_mlm,
        batch_size=1,
        max_aa=args.max_aa,
        device=device,
    )
    seed_emb = l2_normalize(np.asarray(seed_emb, dtype=np.float32))                      # Normalize the seed embedding for cosine computations.
    
    # Create the evolving parent pool as a list of simple records.
    current_pool = [                                                                    
        {
            "candidate_id": f"seed::{seed_protein_id}",
            "aa_sequence": seed_sequence,
            "selection_score": 0.0,
            "target_score": 0.0,
            "n_mutations": 0,
            "mutations": "",
        }
    ]
    all_rounds: List[pd.DataFrame] = []                                                  # Collect one DataFrame per round and concatenate at the end.
    candidate_counter = 0                                                                # Keep a running counter so every candidate ID is unique across rounds.

    # For each round: 
    for round_idx in range(1, args.rounds + 1):                                          
        proposals: List[dict] = []                                                       # Store all candidate proposal rows generated in this round.
        per_parent_budget = int(np.ceil(args.candidates_per_round / max(len(current_pool), 1)))  # Split the round budget (= number of candidates per parent) roughly evenly across current parents.

        entropy_by_parent: Dict[str, np.ndarray] = {}                                    # Cache optional entropy vectors for each parent sequence.
        # If the user requested entropy-aware targeting: Compute and cache one entropy vector per parent sequence.
        if args.use_entropy:
            parent_sequences = [row["aa_sequence"] for row in current_pool]              # Collect the current parent sequences once.
            entropies = stage04.compute_position_entropies_batch(                        # Reuse the Stage 04 entropy computation if the user requested entropy-aware targeting.
                sequences=parent_sequences,
                tokenizer=tokenizer,
                model=esm_mlm,
                max_aa=args.max_aa,
                device=device,
                aa_token_ids=aa_token_ids,
                batch_size=args.batch_size,
            )
            for row, ent in zip(current_pool, entropies):
                entropy_by_parent[row["candidate_id"]] = ent                             # Cache one entropy vector per parent candidate ID.

        # For each parent:
        for parent_rank, parent_row in enumerate(current_pool, start=1):                 # Expand proposals from each parent kept from the previous round.
            parent_seq = str(parent_row["aa_sequence"])                                  # Recover the parent amino-acid sequence.
            entropy_scores = entropy_by_parent.get(parent_row["candidate_id"], None)     # Recover optional entropy scores for this parent.
            positions_requests: List[List[int]] = []                                     # Collect the mutation-position sets requested for this parent.
            # For each of the candidates of each parent: 
            for _ in range(per_parent_budget):
                # Sample how many mutations to make and choose the positions inside the allowed window.
                n_mut = rng.randint(args.min_mutations, args.max_mutations)              
                positions = choose_positions_with_target_prior(                          
                    seq=parent_seq,
                    n_mutations=n_mut,
                    allowed_positions=allowed_positions,
                    target_prior=target_prior,
                    rng=rng,
                    entropy_scores=entropy_scores,
                )
                if len(positions) > 0:
                    # Create a list of all the valid position sets that were sampled for each candidate from this parent.
                    positions_requests.append(positions)                                  

            if len(positions_requests) == 0:
                continue                                                                  # If no valid positions were sampled for this parent, skip it.
            
            # For each batch of position sets:
            for batch_start in range(0, len(positions_requests), args.proposal_batch_size):     # Process parent-specific masked proposals in batches to control GPU memory.
                # Build the batch of masked candidates, and score them in one batched forward pass.
                batch_positions = positions_requests[batch_start : batch_start + args.proposal_batch_size]  # Slice one batch of position sets.
                toks, valid_positions_list = stage04.build_masked_candidate_batch(              # Reuse the Stage 04 batch builder that masks the chosen positions.
                    parent_sequence=parent_seq,
                    positions_list=batch_positions,
                    tokenizer=tokenizer,
                    max_aa=args.max_aa,
                    device=device,
                )
                with torch.inference_mode():                                              # Disable gradients because we are only generating and scoring proposals.
                    outputs = esm_mlm(input_ids=toks["input_ids"], attention_mask=toks["attention_mask"])
                    logits = outputs.logits                                               # Obtain MLM logits for every masked candidate row.
                # Sample mutations from the batched logits.
                mutated_sequences, mutation_strings = stage04.sample_mutations_from_batched_logits(  # Reuse the original Stage 04 mutation sampler from batched logits.
                    parent_sequence=parent_seq,
                    valid_positions_list=valid_positions_list,
                    logits=logits,
                    tokenizer=tokenizer,
                    top_k=args.proposal_top_k,
                    rng=rng,
                )

                # For each of the proposed mutated sequences from this batch, record the metadata needed downstream.
                for mutated_seq, mutation_str in zip(mutated_sequences, mutation_strings):
                    if mutation_str == "":
                        continue                                                          # Skip degenerate proposals that produced no real mutation.
                    proposals.append({
                        "round": round_idx,
                        "parent_rank": parent_rank,
                        "parent_candidate_id": parent_row["candidate_id"],
                        "virus_accession": virus_accession,
                        "source_host": source_host,
                        "seed_protein_id": seed_protein_id,
                        "protein_id": seed_protein_id,
                        "candidate_id": f"round{round_idx}_cand{candidate_counter}",
                        "target_host": args.target_host,
                        "mutations": mutation_str,
                        "aa_sequence": mutated_seq,
                    })                                                                    # Record the mutated candidate row with the metadata needed downstream.
                    candidate_counter += 1

        if len(proposals) == 0:
            print(f"[WARN] No proposals generated in round {round_idx}; stopping early.")
            break                                                                         # Stop early if the search produced no valid candidates.
        
        # Gather all the unique proposed candidates and get their embeddings and target-host score all at once.
        proposed_df = pd.DataFrame(proposals).drop_duplicates(subset=["aa_sequence"]).reset_index(drop=True)  # Deduplicate identical amino-acid sequences before expensive scoring.
        embeddings = stage04.embed_sequences(                                             # Reuse the Stage 04 embedding helper for all proposed candidates.
            sequences=proposed_df["aa_sequence"].tolist(),
            tokenizer=tokenizer,
            model=esm_mlm,
            batch_size=args.batch_size,
            max_aa=args.max_aa,
            device=device,
        )
        probs = clf.predict_proba(embeddings)                                             # Score all candidates with the trained host probe.
        target_scores = probs[:, target_idx]                                              # Extract the target-host probability column.
        pred_idx = probs.argmax(axis=1)                                                   # Recover the predicted host label per candidate.
        pred_label = [label_classes[i] for i in pred_idx]                                 # Convert predicted indices back into host-label strings.

        # Compute similarity between the candidate and the centroid of the seed's family, the centroid of the family of protein's-that-infect-the-target (target anchors), and the seed protein itself.
        normalized_cand, family_cos = compute_family_similarity_terms(embeddings, runtime)  # Compute family-preservation similarity for every candidate ("Stay within the cloud of stable, working tail fibers" - Structural stability "consensus").
        if target_centroid is None:
            target_anchor_cos = np.zeros(len(proposed_df), dtype=np.float32)              # If no target anchor exists in the family, set the attraction term to zero.
        else:
            target_anchor_cos = normalized_cand @ np.asarray(target_centroid, dtype=np.float32).reshape(-1)  # Compute attraction toward the in-family target centroid.
        seed_cos = (normalized_cand @ seed_emb.reshape(-1))                              # Measure how close each proposal candidate still is to the immediate seed of this ladder step.
        # Count the number of mutation proposals from the mutation string, and use that as the penalty input.
        n_mut = proposed_df["mutations"].fillna("").apply(lambda s: 0 if str(s).strip() == "" else len(str(s).split(";"))).astype(int).to_numpy()  # Count actual proposed mutations from the mutation string.
        mutation_penalty = n_mut.astype(np.float32)                                      # Use the raw count as the penalty input.

        # Create a dataframe with all the computed metadata of the proposed candidates, and compute the selection score.
        proposed_df["pred_label"] = pred_label                                          # Store the predicted host label on the candidate table.
        proposed_df["target_score"] = target_scores                                     # Store the target-host score.
        proposed_df["family_cosine"] = family_cos                                       # Store family-preservation similarity.
        proposed_df["target_anchor_cosine"] = target_anchor_cos                         # Store target-family attraction similarity.
        proposed_df["seed_cosine_similarity"] = seed_cos                                # Store similarity to the immediate seed sequence embedding.
        proposed_df["n_mutations"] = n_mut                                              # Store the mutation count explicitly for ranking and reporting.
        # Selection score = the sum of the target score and a weighted sum of the family, target and seed similarity terms and the mutation penalty.
        proposed_df["selection_score"] = (
            proposed_df["target_score"] +                                               # Main term: high host-probe score for the desired target.
            (args.lambda_family * proposed_df["family_cosine"]) +                       # Keep proposals inside the chosen scaffold family manifold.
            (args.lambda_target_anchor * proposed_df["target_anchor_cosine"]) -         # Pull proposals toward the target-host centroid when it exists.
            (args.lambda_seed_novelty * proposed_df["seed_cosine_similarity"]) -        # Lightly discourage proposals from remaining too close to the current seed.
            (args.lambda_mutation_penalty * proposed_df["n_mutations"])                 # Penalize gratuitous extra mutations.
        )                                                                               # Combine all score terms into one practical selection objective.
        # Sort the proposed candidates by their selection score, and rank them within the round.
        proposed_df = proposed_df.sort_values(
            by=["selection_score", "target_score", "family_cosine", "target_anchor_cosine", "n_mutations"],
            ascending=[False, False, False, False, True],
        ).reset_index(drop=True)                                                         # Sort candidates so the strongest family-consistent proposals rise to the top.
        proposed_df["rank_in_round"] = np.arange(1, len(proposed_df) + 1)                # Store the within-round rank for later reporting.
        proposed_df.to_csv(out_dir / f"round_{round_idx}_candidates.csv", index=False)   # Save the full scored proposal table for this round.
        all_rounds.append(proposed_df)                                                   # Append the round table to the run history.

        # Select the top k diverse candidates from the round for the next round of proposals.
        current_pool = stage04.select_diverse_top_candidates(                            # Reuse the Stage 04 diversity selector to seed the next round.
            df=proposed_df,
            keep_top_k=args.keep_top_k,
            min_distance=args.diversity_min_distance,
        ).to_dict(orient="records")                                                      # Convert the retained diverse top candidates back into parent records.

    if len(all_rounds) == 0:
        raise RuntimeError("No candidate rounds were successfully generated in Phase A optimization.")  # Fail loudly if the run produced nothing.

    # Combine the candidate tables from all rounds into a single master run table, and keep a final set of top diverse candidates.
    all_candidates_df = pd.concat(all_rounds, axis=0).reset_index(drop=True)            # Concatenate all round candidate tables into one master run table.
    all_candidates_df.to_csv(out_dir / "all_candidates.csv", index=False)               # Save the full run history.
    top_df = (
        all_candidates_df
        .sort_values(by=["selection_score", "target_score", "family_cosine", "target_anchor_cosine", "n_mutations"], ascending=[False, False, False, False, True])
        .drop_duplicates(subset=["aa_sequence"])                                        # Collapse duplicate sequences across rounds before final ranking.
        .reset_index(drop=True)
    )
    top_df = stage04.select_diverse_top_candidates(                                     # Keep a final diverse top set, not just the top rows by score.
        df=top_df,
        keep_top_k=50,
        min_distance=args.diversity_min_distance,
    )
    top_df.to_csv(out_dir / "top_candidates.csv", index=False)                          # Save the final top candidate set for this ladder step.

    # Write out the metadata for this Phase A optimization run.
    metadata = {
        "phase": "A",
        "seed_protein_id": seed_protein_id,
        "seed_label": args.seed_label or seed_protein_id,
        "source_host": source_host,
        "target_host": args.target_host,
        "seed_fasta": args.seed_fasta,
        "phaseA_plan_json": str(args.phaseA_plan_json),
        "stage04_script": str(args.stage04_script),
        "esm_model": args.esm_model,
        "model_path": str(args.model_path),
        "label_classes_path": str(args.label_classes_path),
        "embedding_pt": str(args.embedding_pt),
        "embedding_index_csv": str(args.embedding_index_csv),
        "rounds": int(args.rounds),
        "candidates_per_round": int(args.candidates_per_round),
        "min_mutations": int(args.min_mutations),
        "max_mutations": int(args.max_mutations),
        "keep_top_k": int(args.keep_top_k),
        "proposal_top_k": int(args.proposal_top_k),
        "proposal_batch_size": int(args.proposal_batch_size),
        "batch_size": int(args.batch_size),
        "max_aa": int(args.max_aa),
        "seed": int(args.seed),
        "diversity_min_distance": int(args.diversity_min_distance),
        "lambda_family": float(args.lambda_family),
        "lambda_target_anchor": float(args.lambda_target_anchor),
        "lambda_seed_novelty": float(args.lambda_seed_novelty),
        "lambda_mutation_penalty": float(args.lambda_mutation_penalty),
        "use_entropy": bool(args.use_entropy),
        "mutation_window_span_1based": phaseA_plan["mutation_window_span_1based"],
        "n_total_candidates": int(len(all_candidates_df)),
        "n_top_candidates_saved": int(len(top_df)),
    }                                                                                    # Collect the run configuration and summary in one metadata object.
    with open(out_dir / "run_metadata.json", "w", encoding="utf-8") as handle:           # Write the metadata JSON for this ladder step.
        json.dump(metadata, handle, indent=2)

    print(f"Saved: {out_dir / 'all_candidates.csv'}")                                   # Confirm the full candidate history path.
    print(f"Saved: {out_dir / 'top_candidates.csv'}")                                   # Confirm the final ranked candidate path.
    print(f"Saved: {out_dir / 'run_metadata.json'}")                                    # Confirm the metadata path.


if __name__ == "__main__":
    main()
