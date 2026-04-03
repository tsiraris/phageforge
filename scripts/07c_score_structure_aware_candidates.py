"""Stage 07c: Score generated candidates against the seed, family manifold, and target centroid."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from phageforge.stage07_utils import cosine_similarity, embed_sequences, read_json


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for structure-aware candidate scoring."""
    ap = argparse.ArgumentParser(description="Score Stage 07 candidates with ESM embeddings.")
    ap.add_argument("--context_json", type=str, required=True, help="Stage 07 context JSON produced by 07a_prepare_stage07_design_context.py.")
    ap.add_argument("--generated_csv", type=str, required=True, help="Generated candidate CSV produced by 07b_generate_rbps_with_esm3.py.")
    ap.add_argument("--scored_csv", type=str, required=True, help="Where to write the structure-aware scored CSV.")
    ap.add_argument("--esm_model", type=str, default="facebook/esm2_t33_650M_UR50D", help="ESM embedding model used for manifold scoring.")
    ap.add_argument("--batch_size", type=int, default=2, help="Batch size for embedding candidate sequences.")
    ap.add_argument("--max_aa", type=int, default=2048, help="Maximum sequence length passed into the embedding model.")
    return ap.parse_args()


def main() -> None:
    # Read the Stage 07 context and the generated candidates to prepare embedding-based scoring.
    args = parse_args()
    context = read_json(args.context_json)
    generated_df = pd.read_csv(args.generated_csv)
    if generated_df.empty:
        raise ValueError(f"Generated CSV is empty: {args.generated_csv}")

    valid_df = generated_df.loc[generated_df["generation_status"].astype(str) == "ok"].copy().reset_index(drop=True)
    if valid_df.empty:
        raise ValueError("No successful ESM3 candidates were found to score.")

    # Embed the selected seed and the candidate sequences using the requested ESM representation model.
    seed_sequence = str(context["selected_seed"]["seed_sequence"])
    sequences = [seed_sequence] + valid_df["candidate_sequence"].astype(str).tolist()
    embeddings = embed_sequences(sequences, model_name=args.esm_model, batch_size=args.batch_size, max_length=args.max_aa)
    seed_embedding = embeddings[0]
    candidate_embeddings = embeddings[1:]

    # Compare each candidate embedding to the seed, the family centroid, and the target-host centroid.
    family_centroid = np.asarray(context["family_context"].get("family_centroid", []), dtype=np.float32)
    target_centroid = np.asarray(context["target_context"].get("target_centroid", []), dtype=np.float32)
    valid_df["seed_cosine"] = [cosine_similarity(emb, seed_embedding) for emb in candidate_embeddings]
    valid_df["esm_embedding_model"] = args.esm_model
    valid_df["family_cosine"] = [cosine_similarity(emb, family_centroid) if family_centroid.size else 0.0 for emb in candidate_embeddings]
    valid_df["target_anchor_cosine"] = [cosine_similarity(emb, target_centroid) if target_centroid.size else 0.0 for emb in candidate_embeddings]
    valid_df["structure_score"] = (
        0.25 * valid_df["seed_cosine"]
        + 0.30 * valid_df["family_cosine"]
        + 0.45 * valid_df["target_anchor_cosine"]
    )

    # Merge the scores back onto the original generated table so provenance rows are preserved as well.
    scored_df = generated_df.merge(
        valid_df[["sample_id", "seed_cosine", "esm_embedding_model", "family_cosine", "target_anchor_cosine", "structure_score"]],
        on="sample_id",
        how="left",
    )
    out_path = Path(args.scored_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    scored_df.to_csv(out_path, index=False)
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
