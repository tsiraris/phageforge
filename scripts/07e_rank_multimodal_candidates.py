"""Stage 07e: merge Stage 07 scoring tables, deduplicate candidates, and build a diversity-aware panel.

This script takes the structure-aware Stage 07 scoring output, computes a final multimodal ranking
score, and then applies a diversity-aware greedy selection pass so the validation shortlist does not
collapse into near-duplicate candidates from one dominant generation regime.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from phageforge.stage07_utils import embed_sequences, greedy_diverse_order


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 07 multimodal ranking step."""
    parser = argparse.ArgumentParser(description="Rank Stage 07 candidates with multimodal and diversity-aware criteria.")
    parser.add_argument("--generated_csv", type=str, required=True, help="Candidate CSV from 07b_generate_rbps_with_esm3.py.")
    parser.add_argument("--structure_scored_csv", type=str, required=True, help="Scored CSV from 07c_score_structure_aware_candidates.py.")
    parser.add_argument("--host_validity_csv", type=str, default=None, help="Optional host-validity CSV to merge if available.")
    parser.add_argument("--tissue_embeddings_pt", type=str, default=None, help="Reserved optional argument for backward compatibility.")
    parser.add_argument("--tissue_metadata_csv", type=str, default=None, help="Reserved optional argument for backward compatibility.")
    parser.add_argument("--out_csv", type=str, required=True, help="Where to write the final multimodal ranked CSV.")
    parser.add_argument("--diverse_top_k", type=int, default=5, help="Number of candidates to mark for the diversity-aware validation panel.")
    parser.add_argument("--diversity_penalty_weight", type=float, default=0.25, help="How strongly to penalize similarity to already selected candidates.")
    parser.add_argument("--per_regime_pool", type=int, default=3, help="How many top candidates per generation regime to prioritize in the diversity pass.")
    parser.add_argument("--embedding_model", type=str, default="facebook/esm2_t12_35M_UR50D", help="Compact ESM model used for diversity selection embeddings.")
    return parser.parse_args()          # Return the parsed CLI namespace.


def main() -> None:
    # Parse the CLI once at startup so all downstream settings are centralized.
    args = parse_args()     # Read the command-line arguments for this ranking run.

    # Load the structure-aware scoring table and keep only candidates that completed successfully.
    scored_df = pd.read_csv(args.structure_scored_csv)  # Read the Stage 07 structure-scored candidate table.
    ranked_df = scored_df.loc[scored_df["generation_status"].astype(str) == "ok"].copy().reset_index(drop=True)     # Keep only successful generations.
    ranked_df = ranked_df.drop_duplicates(subset=["candidate_sequence"]).reset_index(drop=True)                     # Remove exact duplicate sequences before ranking.

    # Initialize optional modalities and provenance flags so later ranking columns are always present.
    ranked_df["tissue_score"] = 0.0                                                                                             # Fill the reserved tissue modality with zeros for compatibility.
    ranked_df["target_score"] = ranked_df.get("target_anchor_cosine", 0.0)                                                      # Reuse the target-anchor cosine as the target-side score.
    ranked_df["strict_manifold_score"] = 0.5 * ranked_df.get("seed_cosine", 0.0) + 0.5 * ranked_df.get("family_cosine", 0.0)    # Average seed and family preservation scores.
    ranked_df["used_esm3_api"] = ranked_df.get("used_esm3_api", False).astype(bool)                                             # Normalize the Forge/API provenance flag.
    ranked_df["used_local_generator"] = ranked_df.get("used_local_generator", ranked_df["generator_mode"].astype(str).str.contains("local", case=False, na=False)).astype(bool)  # Recover whether a local generator path was used.
    ranked_df["used_local_esm3"] = ranked_df.get("used_local_esm3", ranked_df["generator_mode"].astype(str).str.contains("esm3_local|local:esm3", case=False, na=False)).astype(bool)  # Recover whether local ESM3 was used.
    ranked_df["used_esm2_fallback"] = ranked_df.get("used_esm2_fallback", False).astype(bool)                                   # Normalize the ESM2 fallback provenance flag.
    ranked_df["provenance_complete"] = (~ranked_df["used_esm2_fallback"]).astype(bool)                                          # Mark candidates that stayed inside the preferred provenance path.

    # Combine the main Stage 07 ranking signals into one base score before diversity reranking.
    ranked_df["base_rank_score"] = (  # Build the raw ranking score used before diversity penalties.
        0.42 * ranked_df["target_score"]
        + 0.28 * ranked_df["strict_manifold_score"]
        + 0.20 * ranked_df["structure_score"]
        + 0.05 * ranked_df["tissue_score"]
        + 0.05 * ranked_df.get("guided_mutation_score", 0.0).fillna(0.0)
        - 0.01 * ranked_df["mutation_penalty"].fillna(0.0)
    )
    ranked_df = ranked_df.sort_values(["base_rank_score", "target_score", "strict_manifold_score"], ascending=[False, False, False]).reset_index(drop=True)  # Sort by the base score and tie-breakers.
    ranked_df["rank_raw"] = np.arange(1, len(ranked_df) + 1)  # Record the raw pre-diversity rank.
    ranked_df["regime_pool_rank"] = ranked_df.groupby("generation_regime", dropna=False)["base_rank_score"].rank(method="first", ascending=False)  # Rank candidates within each generation regime.
    preferred_mask = ranked_df["regime_pool_rank"] <= max(1, int(args.per_regime_pool))  # Mark the top candidates from each regime as preferred during diversity selection.

    # Compute compact sequence embeddings and then greedily build a diverse validation ordering.
    diversity_embeddings = embed_sequences(ranked_df["candidate_sequence"].astype(str).tolist(), model_name=args.embedding_model, batch_size=4)  # Embed candidates for diversity-aware ordering.
    greedy_order, penalties = greedy_diverse_order(  # Run the greedy diversity pass over the ranking table.
        diversity_embeddings,
        ranked_df["base_rank_score"].to_numpy(dtype=np.float32),
        penalty_weight=args.diversity_penalty_weight,
        preferred_mask=preferred_mask.to_numpy(dtype=bool),
    )
    rank_diverse = np.empty(len(ranked_df), dtype=int)              # Allocate the final diversity-aware rank array.
    rank_diverse[greedy_order] = np.arange(1, len(ranked_df) + 1)   # Convert the greedy order into rank positions.
    ranked_df["nearest_selected_cosine"] = penalties                # Store the nearest-selected similarity for inspection.
    ranked_df["diversity_penalty"] = penalties                      # Reuse the same similarity term as the diversity penalty.
    ranked_df["rank_diverse"] = rank_diverse                        # Store the final diversity-aware rank.
    ranked_df["selected_for_panel"] = ranked_df["rank_diverse"] <= int(args.diverse_top_k)  # Mark the intended validation panel members.
    ranked_df["final_multimodal_rank_score"] = ranked_df["base_rank_score"] - args.diversity_penalty_weight * ranked_df["diversity_penalty"]  # Compute the final score after diversity penalization.
    ranked_df = ranked_df.sort_values(["rank_diverse", "rank_raw"]).reset_index(drop=True)  # Present the output in the final diversity-aware order.

    # Write both the full ranked table and the compact validation panel export.
    out_path = Path(args.out_csv)                                               # Convert the output CSV argument into a path object.
    out_path.parent.mkdir(parents=True, exist_ok=True)                          # Create the output directory if it does not already exist.
    ranked_df.drop(columns=["base_rank_score"]).to_csv(out_path, index=False)   # Persist the full ranked table without the intermediate base score.

    diverse_path = out_path.with_name("top_validation_panel.csv")               # Define the companion validation-panel CSV path.
    ranked_df.loc[ranked_df["selected_for_panel"]].sort_values("rank_diverse").to_csv(diverse_path, index=False)  # Write the compact validation shortlist.
    print(f"Wrote: {out_path}")                                                 # Report the main ranked CSV path.
    print(f"Wrote: {diverse_path}")                                             # Report the validation-panel CSV path.


if __name__ == "__main__":
    main()  # Execute the CLI entrypoint when the script is run directly.
