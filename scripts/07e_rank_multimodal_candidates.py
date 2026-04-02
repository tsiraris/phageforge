
"""
=====================================================================================
07e: Merge host / validity / structure / tissue terms into one Stage 07 final ranking.
=====================================================================================

This script is the final Stage 07 reranker.
It merges:
- generated candidates
- structure-aware scores
- optional host/validity scores from the Stage 05 stack
- optional tissue-context scores (priors)

Then it computes one final multimodal ranking score.
"""

from __future__ import annotations                                                    # Enable postponed annotation evaluation for cleaner typing.
import argparse                                                                       # Parse command-line arguments.
from pathlib import Path                                                              # Build output paths robustly.
import numpy as np                                                                    # Aggregate tissue priors numerically.
import pandas as pd                                                                   # Merge and rank candidate tables.
import torch                                                                          # Load optional tissue embedding tensors.
from phageforge.eval.stage07_metrics import rank_score_dataframe                      # Reuse the shared final-ranking helper.


def parse_args() -> argparse.Namespace:
    """Return command-line arguments for final Stage 07 multimodal ranking."""
    ap = argparse.ArgumentParser(description="Rank Stage 07 candidates with optional tissue context.")  # Create the parser for the final reranking stage.
    ap.add_argument("--generated_csv", type=str, required=True, help="Generated candidate CSV from 07b.")  # Point to the raw generated-candidate table.
    ap.add_argument("--structure_scored_csv", type=str, required=True, help="Structure-aware scored CSV from 07c.")  # Point to the structure-scored candidate table.
    ap.add_argument("--host_validity_csv", type=str, default="", help="Optional Stage 05-style scoring table with target_score / strict_manifold_score columns.")  # Point to optional host/validity scores.
    ap.add_argument("--tissue_embeddings_pt", type=str, default="", help="Optional tissue embedding tensor from 07d.")  # Point to optional tissue embedding tensor.
    ap.add_argument("--tissue_metadata_csv", type=str, default="", help="Optional tissue metadata CSV from 07d.")  # Point to optional tissue metadata CSV.
    ap.add_argument("--out_csv", type=str, required=True, help="Where to write the final ranked candidate table.")  # Point to the final ranked output CSV.
    return ap.parse_args()                                                            # Parse the CLI and return the namespace.


def merge_optional_scores(base: pd.DataFrame, host_validity_csv: str) -> pd.DataFrame:
    """Attach optional Stage 05 host/validity scores to the current candidate table."""
    df = base.copy()                                                                  # Work on a copy so the caller keeps the base table unchanged.
    if host_validity_csv:                                                             # Only run the merge when an optional host/validity table was provided.
        hv = pd.read_csv(host_validity_csv).copy()                                    # Read the optional host/validity table from disk.
        merge_key = "candidate_sequence" if "candidate_sequence" in hv.columns else ("aa_sequence" if "aa_sequence" in hv.columns else None)  # Detect the sequence key.
        if merge_key is None:
            raise ValueError("host_validity_csv must contain candidate_sequence or aa_sequence column.")
        hv = hv.rename(columns={merge_key: "candidate_sequence"})                     # Normalize the sequence key name for merging.
        wanted = [c for c in ["candidate_sequence", "target_score", "strict_manifold_score", "family_cosine", "target_anchor_cosine"] if c in hv.columns]  # Keep only the ranking columns we care about.
        df = df.merge(
            hv[wanted].drop_duplicates(subset=["candidate_sequence"]),                # Deduplicate repeated sequence rows before merging.
            on="candidate_sequence",
            how="left",
            suffixes=("", "_hv"),
        )
        for col in ["target_score", "strict_manifold_score", "family_cosine", "target_anchor_cosine"]:
            alt = f"{col}_hv"
            if alt in df.columns:
                df[col] = df[col].fillna(df[alt]) if col in df.columns else df[alt]   # Prefer current values and backfill from the optional table.
                df = df.drop(columns=[alt])
    return df


def attach_tissue_scores(df: pd.DataFrame, tissue_embeddings_pt: str, tissue_metadata_csv: str) -> pd.DataFrame:
    """Attach a simple optional tissue prior to every candidate."""
    out = df.copy()
    if not tissue_embeddings_pt or not tissue_metadata_csv:
        out["tissue_score"] = 0.0                                                     # Fill the tissue score column with neutral zeros.
        return out

    tissue_emb = torch.load(tissue_embeddings_pt, map_location="cpu")                 # Load the optional tissue embedding tensor from disk.
    if isinstance(tissue_emb, torch.Tensor):
        tissue_emb = tissue_emb.numpy()
    tissue_emb = np.asarray(tissue_emb, dtype=np.float32)
    tissue_meta = pd.read_csv(tissue_metadata_csv).copy()                             # Read the tissue metadata CSV aligned to the tissue tensor.
    tissue_meta["tissue_score_prior"] = (
        tissue_meta["tissue_score_prior"] if "tissue_score_prior" in tissue_meta.columns else 0.0
    )
    priors = tissue_meta["tissue_score_prior"].to_numpy(dtype=float) if len(tissue_meta) else np.zeros(1, dtype=float)
    mean_prior = float(priors.mean()) if priors.size else 0.0
    out["tissue_score"] = mean_prior                                                 # Attach the same mean tissue prior to every candidate for now.
    return out


def ensure_provenance_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure ESM3-specific provenance columns always exist and derive report-friendly helpers."""
    out = df.copy()
    defaults = {
        "generator_mode": "unknown_generator",
        "esm3_model": pd.NA,
        "esm3_prompt_sequence": pd.NA,
        "esm3_hotspot_positions": pd.NA,
        "esm3_temperature": pd.NA,
        "esm3_top_k": pd.NA,
        "esm3_num_steps": pd.NA,
        "esm3_sampling_seed": pd.NA,
        "mutation_positions": pd.NA,
        "editable_hotspot_count": 0,
    }
    for col, default in defaults.items():
        if col not in out.columns:
            out[col] = default

    out["used_esm3_api"] = out["generator_mode"].astype(str).str.startswith("esm3_api:")     # Flag candidates produced by the real Forge ESM3 backend.
    out["used_local_generator"] = out["generator_mode"].astype(str).eq("local_conditional_generator")  # Flag candidates produced by the local conditioned generator.
    out["used_local_esm3"] = out["generator_mode"].astype(str).str.startswith("esm3_local:")  # Flag candidates produced by the local open-weights ESM3 backend.
    out["used_esm2_fallback"] = out["generator_mode"].astype(str).str.startswith("esm2_masked_lm:")  # Flag candidates produced by the ESM2 fallback backend.
    out["provenance_complete"] = (
        out["generator_mode"].notna() &
        out["esm3_prompt_sequence"].notna() &
        out["mutation_positions"].notna()
    )
    return out


def fill_stage07_score_fallbacks(df: pd.DataFrame) -> pd.DataFrame:
    """Populate host/validity score columns from Stage 07-native proxies when external Stage 05 tables are absent."""
    out = df.copy()
    numeric_defaults = ["target_score", "strict_manifold_score", "family_cosine", "target_anchor_cosine", "seed_cosine", "structure_score", "tissue_score"]
    for col in numeric_defaults:
        if col not in out.columns:
            out[col] = 0.0
        out[col] = pd.to_numeric(out[col], errors="coerce")                           # Normalize numeric score columns before fallback logic.

    # Build Stage 07-native host / manifold surrogates only when the corresponding columns are absent or effectively zero everywhere.
    if out["target_score"].fillna(0.0).abs().sum() == 0.0:                           # If no explicit host-transfer score exists, derive one from the target-anchor signal.
        out["target_score"] = out["target_anchor_cosine"].fillna(out["structure_score"]).fillna(0.0)
    else:
        out["target_score"] = out["target_score"].fillna(out["target_anchor_cosine"]).fillna(0.0)

    if out["strict_manifold_score"].fillna(0.0).abs().sum() == 0.0:                  # If no explicit strict-manifold score exists, derive one from seed/family preservation.
        out["strict_manifold_score"] = (
            0.60 * out["family_cosine"].fillna(0.0) +
            0.40 * out["seed_cosine"].fillna(0.0)
        )
    else:
        out["strict_manifold_score"] = out["strict_manifold_score"].fillna(
            0.60 * out["family_cosine"].fillna(0.0) + 0.40 * out["seed_cosine"].fillna(0.0)
        )

    out["family_cosine"] = out["family_cosine"].fillna(0.0)
    out["target_anchor_cosine"] = out["target_anchor_cosine"].fillna(0.0)
    out["tissue_score"] = out["tissue_score"].fillna(0.0)
    return out


# Main script: Reads Stage 07b generated candidates, merges inputs and writes the final ranked CSVs.
def main() -> None:
    args = parse_args()
    generated = pd.read_csv(args.generated_csv)                                       # Read the raw generated-candidate table.
    struct = pd.read_csv(args.structure_scored_csv)                                   # Read the structure-aware scored candidate table.

    df = generated.merge(
        struct.drop(columns=[c for c in generated.columns if c in struct.columns and c != "candidate_sequence"]),
        on="candidate_sequence",
        how="left",
    )
    df = merge_optional_scores(df, args.host_validity_csv)                            # Attach Stage 05-style host-transfer and validity terms when available.
    df = attach_tissue_scores(df, args.tissue_embeddings_pt, args.tissue_metadata_csv)  # Attach an optional global tissue prior when available.
    df = ensure_provenance_columns(df)                                                # Normalize ESM3-specific provenance columns before final ranking and reporting.
    df = fill_stage07_score_fallbacks(df)                                             # Repair missing ranking terms using Stage 07-native structure signals.

    if "mutation_penalty" not in df.columns:
        df["mutation_penalty"] = 0.0

    ranked_no_tissue = rank_score_dataframe(df, tissue_enabled=False)                 # Build a no-tissue comparison table for ablation reporting.
    ranked_with_tissue = rank_score_dataframe(
        df,
        tissue_enabled=bool(args.tissue_embeddings_pt and args.tissue_metadata_csv),
    )

    out_path = Path(args.out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ranked_with_tissue.to_csv(out_path, index=False)
    ranked_no_tissue.to_csv(out_path.with_name("top_candidates_no_tissue.csv"), index=False)
    ranked_with_tissue.to_csv(out_path.with_name("top_candidates_with_tissue.csv"), index=False)
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
