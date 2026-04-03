"""Stage 07e: Merge Stage 07 scores, deduplicate candidates, and build a diversity-aware validation panel."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from phageforge.stage07_utils import embed_sequences, greedy_diverse_order


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for multimodal Stage 07 ranking."""
    ap = argparse.ArgumentParser(description="Rank Stage 07 candidates with multimodal and diversity-aware criteria.")
    ap.add_argument("--generated_csv", type=str, required=True, help="Candidate CSV from 07b_generate_rbps_with_esm3.py.")
    ap.add_argument("--structure_scored_csv", type=str, required=True, help="Scored CSV from 07c_score_structure_aware_candidates.py.")
    ap.add_argument("--host_validity_csv", type=str, default=None, help="Optional host-validity CSV to merge if available.")
    ap.add_argument("--tissue_embeddings_pt", type=str, default=None, help="Reserved optional argument for backward compatibility.")
    ap.add_argument("--tissue_metadata_csv", type=str, default=None, help="Reserved optional argument for backward compatibility.")
    ap.add_argument("--out_csv", type=str, required=True, help="Where to write the final multimodal ranked CSV.")
    ap.add_argument("--diverse_top_k", type=int, default=5, help="Number of candidates to mark for the diversity-aware validation panel.")
    ap.add_argument("--diversity_penalty_weight", type=float, default=0.25, help="How strongly to penalize similarity to already selected candidates.")
    ap.add_argument("--per_regime_pool", type=int, default=3, help="How many top candidates per generation regime to prioritize in the diversity pass.")
    ap.add_argument("--embedding_model", type=str, default="facebook/esm2_t12_35M_UR50D", help="Compact ESM model used for diversity selection embeddings.")
    return ap.parse_args()


def main() -> None:
    # Read the Stage 07 tables, deduplicate exact candidate sequences, and initialize optional placeholder modalities.
    args = parse_args()
    scored_df = pd.read_csv(args.structure_scored_csv)
    ranked_df = scored_df.loc[scored_df["generation_status"].astype(str) == "ok"].copy().reset_index(drop=True)
    ranked_df = ranked_df.drop_duplicates(subset=["candidate_sequence"]).reset_index(drop=True)
    ranked_df["tissue_score"] = 0.0
    ranked_df["target_score"] = ranked_df.get("target_anchor_cosine", 0.0)
    ranked_df["strict_manifold_score"] = 0.5 * ranked_df.get("seed_cosine", 0.0) + 0.5 * ranked_df.get("family_cosine", 0.0)
    ranked_df["used_esm3_api"] = ranked_df.get("used_esm3_api", False).astype(bool)
    ranked_df["used_local_generator"] = ranked_df.get("used_local_generator", ranked_df["generator_mode"].astype(str).str.contains("local", case=False, na=False)).astype(bool)
    ranked_df["used_local_esm3"] = ranked_df.get("used_local_esm3", ranked_df["generator_mode"].astype(str).str.contains("esm3_local|local:esm3", case=False, na=False)).astype(bool)
    ranked_df["used_esm2_fallback"] = ranked_df.get("used_esm2_fallback", False).astype(bool)
    ranked_df["provenance_complete"] = (~ranked_df["used_esm2_fallback"]).astype(bool)

    # Combine the available scores into a raw ranking score before diversity-aware reranking.
    ranked_df["base_rank_score"] = (
        0.42 * ranked_df["target_score"]
        + 0.28 * ranked_df["strict_manifold_score"]
        + 0.20 * ranked_df["structure_score"]
        + 0.05 * ranked_df["tissue_score"]
        + 0.05 * ranked_df.get("guided_mutation_score", 0.0).fillna(0.0)
        - 0.01 * ranked_df["mutation_penalty"].fillna(0.0)
    )
    ranked_df = ranked_df.sort_values(["base_rank_score", "target_score", "strict_manifold_score"], ascending=[False, False, False]).reset_index(drop=True)
    ranked_df["rank_raw"] = np.arange(1, len(ranked_df) + 1)
    ranked_df["regime_pool_rank"] = ranked_df.groupby("generation_regime", dropna=False)["base_rank_score"].rank(method="first", ascending=False)
    preferred_mask = ranked_df["regime_pool_rank"] <= max(1, int(args.per_regime_pool))

    # Build a diversity-aware panel by greedily penalizing candidates that are too close to already selected sequences.
    diversity_embeddings = embed_sequences(ranked_df["candidate_sequence"].astype(str).tolist(), model_name=args.embedding_model, batch_size=4)
    greedy_order, penalties = greedy_diverse_order(
        diversity_embeddings,
        ranked_df["base_rank_score"].to_numpy(dtype=np.float32),
        penalty_weight=args.diversity_penalty_weight,
        preferred_mask=preferred_mask.to_numpy(dtype=bool),
    )
    rank_diverse = np.empty(len(ranked_df), dtype=int)
    rank_diverse[greedy_order] = np.arange(1, len(ranked_df) + 1)
    ranked_df["nearest_selected_cosine"] = penalties
    ranked_df["diversity_penalty"] = penalties
    ranked_df["rank_diverse"] = rank_diverse
    ranked_df["selected_for_panel"] = ranked_df["rank_diverse"] <= int(args.diverse_top_k)
    ranked_df["final_multimodal_rank_score"] = ranked_df["base_rank_score"] - args.diversity_penalty_weight * ranked_df["diversity_penalty"]
    ranked_df = ranked_df.sort_values(["rank_diverse", "rank_raw"]).reset_index(drop=True)

    out_path = Path(args.out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ranked_df.drop(columns=["base_rank_score"]).to_csv(out_path, index=False)

    diverse_path = out_path.with_name("top_validation_panel.csv")
    ranked_df.loc[ranked_df["selected_for_panel"]].sort_values("rank_diverse").to_csv(diverse_path, index=False)
    print(f"Wrote: {out_path}")
    print(f"Wrote: {diverse_path}")


if __name__ == "__main__":
    main()
