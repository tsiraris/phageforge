from __future__ import annotations
import argparse
import json
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


def sequence_distance(seq_a: str, seq_b: str) -> int:
    """ Count the differing positions between two protein sequences. """
    # If the lengths are different, add the difference to the distance
    if len(seq_a) != len(seq_b):
        return abs(len(seq_a) - len(seq_b)) + sum(
            a != b for a, b in zip(seq_a, seq_b)
        )
    # If the lengths are equal, just return the number of differing positions
    return sum(a != b for a, b in zip(seq_a, seq_b))


def average_pairwise_distance(sequences: list[str]) -> float:
    """
    Compute the mean pairwise sequence distance among a list of sequences.
    """
    if len(sequences) < 2:
        return 0.0

    distances = []
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            distances.append(sequence_distance(sequences[i], sequences[j]))
    return float(np.mean(distances)) if len(distances) > 0 else 0.0


def count_unique_mutation_positions(mutation_strings: list[str]) -> int:
    """
    Count the number of unique mutated sequence positions represented across the provided mutation strings.
    """
    positions = set()
    for mutation_str in mutation_strings:
        if pd.isna(mutation_str) or mutation_str == "":
            continue
        for token in str(mutation_str).split(";"):
            token = token.strip()
            if token == "":
                continue
            # Example token: S456A
            numeric_part = "".join(ch for ch in token if ch.isdigit())
            if numeric_part != "":
                positions.add(int(numeric_part))
    return len(positions)


def normalize_rows(x: np.ndarray) -> np.ndarray:
    """ L2-normalize each embedding row for cosine-similarity computations. """
    return x / np.clip(np.linalg.norm(x, axis=1, keepdims=True), 1e-12, None)


def compute_chunked_nearest_neighbors(
    query_embeddings_normalized: np.ndarray,
    reference_embeddings_normalized: np.ndarray,
    chunk_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute nearest-neighbor cosine similarity and reference index in chunks to reduce peak memory usage.
    """
    nn_idx_all = []
    nn_sim_all = []

    # Process query embeddings in chunks so cosine-similarity calculation stays memory-efficient
    for start in range(0, len(query_embeddings_normalized), chunk_size):
        query_chunk = query_embeddings_normalized[start : start + chunk_size]            # [chunk, H]: Slice the next chunk of normalized query embeddings
        sim_chunk = query_chunk @ reference_embeddings_normalized.T                       # [chunk, N_ref]: Cosine similarity becomes a matrix multiplication because both matrices are already normalized
        nn_idx_chunk = sim_chunk.argmax(axis=1)                                           # [chunk]: Index of the nearest reference sequence for every query sequence in this chunk
        nn_sim_chunk = sim_chunk[np.arange(len(query_chunk)), nn_idx_chunk]               # [chunk]: Nearest-neighbor cosine similarity for every query sequence in this chunk
        nn_idx_all.append(nn_idx_chunk)
        nn_sim_all.append(nn_sim_chunk)

    return np.concatenate(nn_idx_all), np.concatenate(nn_sim_all)


def compute_nearest_neighbor_novelty(
    top_df: pd.DataFrame,
    reference_df: pd.DataFrame | None,
    esm_model: str,
    batch_size: int,
    max_aa: int,
    cached_reference_embeddings_npy: str | None,
    cached_reference_index_csv: str | None,
    similarity_chunk_size: int,
) -> pd.DataFrame:
    """
    Compute nearest-neighbor novelty in ESM embedding space against a reference dataset.

    Returns a dataframe with one row per top candidate.
    """
    # Keep only rows with valid sequences
    top_valid = top_df.dropna(subset=["aa_sequence"]).copy().reset_index(drop=True)

    if len(top_valid) == 0:
        raise ValueError("No valid candidate sequences found for novelty calculation.")

    # Load cached normalized reference embeddings if available, otherwise build them from the provided reference dataframe
    if cached_reference_embeddings_npy is not None and cached_reference_index_csv is not None:
        cached_reference_embeddings_path = Path(cached_reference_embeddings_npy)
        cached_reference_index_path = Path(cached_reference_index_csv)

        if not cached_reference_embeddings_path.exists():
            raise FileNotFoundError(f"Missing file: {cached_reference_embeddings_path}")
        if not cached_reference_index_path.exists():
            raise FileNotFoundError(f"Missing file: {cached_reference_index_path}")

        ref_unique = pd.read_csv(cached_reference_index_path).dropna(subset=["aa_sequence"]).reset_index(drop=True)          # Load the cached reference index aligned to the cached embeddings
        ref_emb_normalized = np.load(cached_reference_embeddings_path).astype(np.float32)                                     # Load the cached normalized reference embeddings from disk
    else:
        if reference_df is None:
            raise ValueError("Either a novelty reference CSV or cached reference embeddings must be provided.")

        ref_valid = reference_df.dropna(subset=["aa_sequence"]).copy().reset_index(drop=True)
        if len(ref_valid) == 0:
            raise ValueError("No valid reference sequences found for novelty calculation.")

        ref_unique = ref_valid.drop_duplicates(subset=["aa_sequence"]).reset_index(drop=True)                                 # Deduplicate reference sequences so novelty computation is not slowed by exact duplicates

        device = "cuda" if torch.cuda.is_available() else "cpu"

        # Load ESM MLM once and reuse it for novelty calculation
        tokenizer = AutoTokenizer.from_pretrained(esm_model, do_lower_case=False)
        esm_mlm = EsmForMaskedLM.from_pretrained(esm_model)
        esm_mlm.to(device)
        esm_mlm.eval()

        ref_emb = embed_sequences(
            sequences=ref_unique["aa_sequence"].tolist(),
            tokenizer=tokenizer,
            model=esm_mlm,
            batch_size=batch_size,
            max_aa=max_aa,
            device=device,
        )
        ref_emb_normalized = normalize_rows(ref_emb).astype(np.float32)                                                       # Normalize the reference embeddings once so later cosine similarities are just dot products

    device = "cuda" if torch.cuda.is_available() else "cpu"

    # Load ESM MLM once and reuse it for candidate novelty calculation
    tokenizer = AutoTokenizer.from_pretrained(esm_model, do_lower_case=False)
    esm_mlm = EsmForMaskedLM.from_pretrained(esm_model)
    esm_mlm.to(device)
    esm_mlm.eval()

    # Embed top candidates
    top_emb = embed_sequences(
        sequences=top_valid["aa_sequence"].tolist(),
        tokenizer=tokenizer,
        model=esm_mlm,
        batch_size=batch_size,
        max_aa=max_aa,
        device=device,
    )
    top_emb_normalized = normalize_rows(top_emb).astype(np.float32)                                                            # Normalize candidate embeddings once before chunked cosine similarity calculation

    # Compute nearest neighbors with chunked cosine similarity to reduce peak memory usage
    nn_idx, nn_sim = compute_chunked_nearest_neighbors(
        query_embeddings_normalized=top_emb_normalized,
        reference_embeddings_normalized=ref_emb_normalized,
        chunk_size=similarity_chunk_size,
    )
    novelty_distance = 1.0 - nn_sim

    nn_df = pd.DataFrame(
        {
            "candidate_id": top_valid["candidate_id"].tolist(),
            "nearest_neighbor_reference_index": nn_idx.astype(int),
            "nearest_neighbor_cosine_similarity": nn_sim.astype(float),
            "novelty_distance": novelty_distance.astype(float),
            "nearest_neighbor_protein_id": ref_unique.iloc[nn_idx]["protein_id"].astype(str).tolist() if "protein_id" in ref_unique.columns else ["N/A"] * len(nn_idx),
            "nearest_neighbor_host_genus": ref_unique.iloc[nn_idx]["host_genus"].astype(str).tolist() if "host_genus" in ref_unique.columns else ["N/A"] * len(nn_idx),
            "nearest_neighbor_virus_accession": ref_unique.iloc[nn_idx]["virus_accession"].astype(str).tolist() if "virus_accession" in ref_unique.columns else ["N/A"] * len(nn_idx),
        }
    )

    return nn_df


def main():
    # Create an ArgumentParser object to run the same script with different parameters without editing code.
    ap = argparse.ArgumentParser()
    ap.add_argument("--run_dir", type=str, required=True)
    ap.add_argument("--top_k", type=int, default=10)
    ap.add_argument("--novelty_reference_csv", type=str, default=None)
    ap.add_argument("--cached_reference_embeddings_npy", type=str, default=None)
    ap.add_argument("--cached_reference_index_csv", type=str, default=None)
    ap.add_argument("--esm_model", type=str, default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--batch_size", type=int, default=8)
    ap.add_argument("--max_aa", type=int, default=1022)
    ap.add_argument("--similarity_chunk_size", type=int, default=64)
    args = ap.parse_args()

    # Load data
    run_dir = Path(args.run_dir)
    top_csv = run_dir / "top_candidates.csv"
    all_csv = run_dir / "all_candidates.csv"
    metadata_json = run_dir / "run_metadata.json"

    # Check that the files exist
    if not top_csv.exists():
        raise FileNotFoundError(f"Missing file: {top_csv}")
    if not all_csv.exists():
        raise FileNotFoundError(f"Missing file: {all_csv}")
    if not metadata_json.exists():
        raise FileNotFoundError(f"Missing file: {metadata_json}")

    # Load data and metadata
    top_df = pd.read_csv(top_csv)
    all_df = pd.read_csv(all_csv)
    metadata = json.loads(metadata_json.read_text())

    # Create summary files and preview
    pred_label_counts = (
        top_df["pred_label"]
        .value_counts()
        .rename_axis("pred_label")
        .reset_index(name="count")
    )

    # Count the distribution of mutation counts across the top candidates
    mutation_count_counts = (
        top_df["n_mutations"]
        .value_counts()
        .sort_index()
        .rename_axis("n_mutations")
        .reset_index(name="count")
    )

    # Compute summary statistics
    score_summary = pd.DataFrame(
        {
            "metric": ["min", "mean", "median", "max"],
            "target_score": [
                float(top_df["target_score"].min()),
                float(top_df["target_score"].mean()),
                float(top_df["target_score"].median()),
                float(top_df["target_score"].max()),
            ],
            "selection_score": [
                float(top_df["selection_score"].min()) if "selection_score" in top_df.columns else float("nan"),
                float(top_df["selection_score"].mean()) if "selection_score" in top_df.columns else float("nan"),
                float(top_df["selection_score"].median()) if "selection_score" in top_df.columns else float("nan"),
                float(top_df["selection_score"].max()) if "selection_score" in top_df.columns else float("nan"),
            ],
        }
    )

    # Compute additional diversity-focused metrics
    top_sequences = top_df["aa_sequence"].dropna().astype(str).tolist() if "aa_sequence" in top_df.columns else []
    avg_pairwise_dist = average_pairwise_distance(top_sequences)
    unique_mutation_site_count = count_unique_mutation_positions(top_df["mutations"].fillna("").astype(str).tolist()) if "mutations" in top_df.columns else 0

    diversity_summary = pd.DataFrame(
        {
            "metric": ["avg_pairwise_distance_top_candidates", "unique_mutation_site_count_top_candidates"],
            "value": [float(avg_pairwise_dist), int(unique_mutation_site_count)],
        }
    )

    seed_similarity_summary = None
    if "seed_cosine_similarity" in top_df.columns and "seed_novelty_distance" in top_df.columns:
        seed_similarity_summary = pd.DataFrame(
            {
                "metric": [
                    "min_seed_cosine_similarity",
                    "mean_seed_cosine_similarity",
                    "median_seed_cosine_similarity",
                    "max_seed_cosine_similarity",
                    "min_seed_novelty_distance",
                    "mean_seed_novelty_distance",
                    "median_seed_novelty_distance",
                    "max_seed_novelty_distance",
                ],
                "value": [
                    float(top_df["seed_cosine_similarity"].min()),
                    float(top_df["seed_cosine_similarity"].mean()),
                    float(top_df["seed_cosine_similarity"].median()),
                    float(top_df["seed_cosine_similarity"].max()),
                    float(top_df["seed_novelty_distance"].min()),
                    float(top_df["seed_novelty_distance"].mean()),
                    float(top_df["seed_novelty_distance"].median()),
                    float(top_df["seed_novelty_distance"].max()),
                ],
            }
        )                                                                                 # Summarize how close the final top candidates remain to the original seed embedding if the optimizer wrote those columns

    # Save summary files
    pred_label_counts.to_csv(run_dir / "pred_label_counts.csv", index=False)
    mutation_count_counts.to_csv(run_dir / "mutation_count_counts.csv", index=False)
    score_summary.to_csv(run_dir / "target_score_summary.csv", index=False)
    diversity_summary.to_csv(run_dir / "diversity_summary.csv", index=False)
    if seed_similarity_summary is not None:
        seed_similarity_summary.to_csv(run_dir / "seed_similarity_summary.csv", index=False)

    novelty_df = None
    novelty_summary = None

    # Optionally compute nearest-neighbor novelty against a reference dataset or cached novelty-reference embeddings
    if (
        args.novelty_reference_csv is not None
        or (args.cached_reference_embeddings_npy is not None and args.cached_reference_index_csv is not None)
    ):
        reference_df = None
        if args.novelty_reference_csv is not None:
            novelty_reference_csv = Path(args.novelty_reference_csv)
            if not novelty_reference_csv.exists():
                raise FileNotFoundError(f"Missing file: {novelty_reference_csv}")
            reference_df = pd.read_csv(novelty_reference_csv)

        novelty_df = compute_nearest_neighbor_novelty(
            top_df=top_df,
            reference_df=reference_df,
            esm_model=args.esm_model,
            batch_size=args.batch_size,
            max_aa=args.max_aa,
            cached_reference_embeddings_npy=args.cached_reference_embeddings_npy,
            cached_reference_index_csv=args.cached_reference_index_csv,
            similarity_chunk_size=args.similarity_chunk_size,
        )

        novelty_df.to_csv(run_dir / "nearest_neighbor_novelty.csv", index=False)

        novelty_summary = pd.DataFrame(
            {
                "metric": [
                    "min_novelty_distance",
                    "mean_novelty_distance",
                    "median_novelty_distance",
                    "max_novelty_distance",
                    "min_nearest_neighbor_cosine_similarity",
                    "mean_nearest_neighbor_cosine_similarity",
                    "median_nearest_neighbor_cosine_similarity",
                    "max_nearest_neighbor_cosine_similarity",
                ],
                "value": [
                    float(novelty_df["novelty_distance"].min()),
                    float(novelty_df["novelty_distance"].mean()),
                    float(novelty_df["novelty_distance"].median()),
                    float(novelty_df["novelty_distance"].max()),
                    float(novelty_df["nearest_neighbor_cosine_similarity"].min()),
                    float(novelty_df["nearest_neighbor_cosine_similarity"].mean()),
                    float(novelty_df["nearest_neighbor_cosine_similarity"].median()),
                    float(novelty_df["nearest_neighbor_cosine_similarity"].max()),
                ],
            }
        )

        novelty_summary.to_csv(run_dir / "novelty_summary.csv", index=False)

    # Save preview
    preview_cols = [
        c for c in [
            "candidate_id",
            "source_host",
            "target_host",
            "pred_label",
            "target_score",
            "selection_score",
            "seed_cosine_similarity",
            "seed_novelty_distance",
            "mutations",
            "n_mutations",
        ] if c in top_df.columns
    ]
    top_preview = top_df[preview_cols].head(args.top_k)

    # Create the summary file
    lines = []
    lines.append("# Design Run Summary")
    lines.append("")
    lines.append("## Run metadata")
    lines.append("")
    lines.append(f"- Seed protein ID: `{metadata.get('seed_protein_id')}`")
    lines.append(f"- Virus accession: `{metadata.get('virus_accession')}`")
    lines.append(f"- Source host: `{metadata.get('source_host')}`")
    lines.append(f"- Target host: `{metadata.get('target_host')}`")
    lines.append(f"- ESM model: `{metadata.get('esm_model')}`")
    lines.append(f"- Rounds: `{metadata.get('rounds')}`")
    lines.append(f"- Candidates per round: `{metadata.get('candidates_per_round')}`")
    lines.append(f"- Min mutations: `{metadata.get('min_mutations', metadata.get('mutations_per_candidate'))}`")
    lines.append(f"- Max mutations: `{metadata.get('max_mutations', metadata.get('mutations_per_candidate'))}`")
    lines.append(f"- Top-K kept per round: `{metadata.get('keep_top_k')}`")
    lines.append(f"- Proposal top-K: `{metadata.get('proposal_top_k', 'N/A')}`")
    lines.append(f"- Seed: `{metadata.get('seed')}`")
    lines.append(f"- Total candidates evaluated: `{metadata.get('n_total_candidates')}`")
    lines.append(f"- Top candidates saved: `{metadata.get('n_top_candidates_saved')}`")
    lines.append(f"- Diversity enabled: `{metadata.get('diversity_on', False)}`")
    lines.append(f"- Diversity min distance: `{metadata.get('diversity_min_distance', 'N/A')}`")
    lines.append(f"- Position selection strategy: `{metadata.get('position_selection_strategy', 'random')}`")
    lines.append(f"- Position pool size: `{metadata.get('position_pool_size', 'N/A')}`")
    lines.append(f"- Proposal batch size: `{metadata.get('proposal_batch_size', 'N/A')}`")
    lines.append(f"- Novelty penalty lambda: `{metadata.get('novelty_penalty_lambda', 0.0)}`")
    lines.append(f"- Novelty penalty schedule: `{metadata.get('novelty_penalty_schedule', 'constant')}`")
    lines.append("")

    lines.append("## Target score summary")
    lines.append("")
    lines.append(score_summary.to_markdown(index=False))
    lines.append("")

    lines.append("## Predicted label distribution in top candidates")
    lines.append("")
    lines.append(pred_label_counts.to_markdown(index=False))
    lines.append("")

    lines.append("## Mutation count distribution in top candidates")
    lines.append("")
    lines.append(mutation_count_counts.to_markdown(index=False))
    lines.append("")

    lines.append("## Diversity summary")
    lines.append("")
    lines.append(diversity_summary.to_markdown(index=False))
    lines.append("")

    if seed_similarity_summary is not None:
        lines.append("## Seed similarity summary")
        lines.append("")
        lines.append(seed_similarity_summary.to_markdown(index=False))
        lines.append("")

    if novelty_summary is not None:
        lines.append("## Nearest-neighbor novelty summary")
        lines.append("")
        lines.append(novelty_summary.to_markdown(index=False))
        lines.append("")

        novelty_preview_cols = [
            c for c in [
                "candidate_id",
                "nearest_neighbor_protein_id",
                "nearest_neighbor_host_genus",
                "nearest_neighbor_cosine_similarity",
                "novelty_distance",
            ] if c in novelty_df.columns
        ]
        novelty_preview = novelty_df[novelty_preview_cols].head(args.top_k)

        lines.append(f"## Nearest-neighbor novelty preview (top {args.top_k})")
        lines.append("")
        lines.append(novelty_preview.to_markdown(index=False))
        lines.append("")

    lines.append(f"## Top {args.top_k} candidates")
    lines.append("")
    lines.append(top_preview.to_markdown(index=False))
    lines.append("")

    # Save report
    report_path = run_dir / "design_run_report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"✅ Saved: {run_dir / 'pred_label_counts.csv'}")
    print(f"✅ Saved: {run_dir / 'mutation_count_counts.csv'}")
    print(f"✅ Saved: {run_dir / 'target_score_summary.csv'}")
    print(f"✅ Saved: {run_dir / 'diversity_summary.csv'}")
    if seed_similarity_summary is not None:
        print(f"✅ Saved: {run_dir / 'seed_similarity_summary.csv'}")
    if novelty_df is not None:
        print(f"✅ Saved: {run_dir / 'nearest_neighbor_novelty.csv'}")
        print(f"✅ Saved: {run_dir / 'novelty_summary.csv'}")
    print(f"✅ Saved: {report_path}")


if __name__ == "__main__":
    main()
