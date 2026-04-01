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
    # Only run the merge when an optional host/validity table was provided.
    if host_validity_csv:                                                             
        hv = pd.read_csv(host_validity_csv).copy()                                    # Read the optional host/validity table from disk.
        # Detect the "candidate_sequence" or "aa_sequence" column names if present in the host/validity table.
        merge_key = "candidate_sequence" if "candidate_sequence" in hv.columns else ("aa_sequence" if "aa_sequence" in hv.columns else None)  
        if merge_key is None:                                                         # Ensure the host/validity table has a sequence key.
            raise ValueError("host_validity_csv must contain candidate_sequence or aa_sequence column.")
        hv = hv.rename(columns={merge_key: "candidate_sequence"})                     # Normalize the sequence key name for merging.
        # Keep only the columns we care about from the host/validity table.
        wanted = [c for c in ["candidate_sequence", "target_score", "strict_manifold_score", "family_cosine", "target_anchor_cosine"] if c in hv.columns]  # Keep only the ranking columns we care about.
        # Merge the optional host/validity columns onto the candidate table.
        df = df.merge(                                                                # Merge the optional host/validity columns onto the candidate table.
            hv[wanted].drop_duplicates(subset=["candidate_sequence"]),                # Deduplicate repeated sequence rows before merging.
            on="candidate_sequence",                                                  # Merge by the normalized sequence column.
            how="left",                                                               # Keep every candidate even if some optional scores are missing.
            suffixes=("", "_hv"),                                                     # Protect existing columns during the merge.
        )
        # Fill missing primary columns from merged fallback columns.
        for col in ["target_score", "strict_manifold_score", "family_cosine", "target_anchor_cosine"]:  
            alt = f"{col}_hv"                                                         # Build the fallback merged column name.
            if alt in df.columns:                                                     # Only repair columns that actually appeared after the merge.
                df[col] = df[col].fillna(df[alt]) if col in df.columns else df[alt]   # Prefer existing values, else backfill from the host/validity table.
                df = df.drop(columns=[alt])                                           # Drop the temporary fallback column after reconciliation.
    return df                                                                         # Return the candidate table with optional host/validity scores attached.


def attach_tissue_scores(df: pd.DataFrame, tissue_embeddings_pt: str, tissue_metadata_csv: str) -> pd.DataFrame:
    """Attach a simple optional tissue prior to every candidate.

    When no explicit candidate-to-context mapping exists, it uses a global mean prior.
    That keeps the tissue branch usable now while leaving room for finer candidate-specific context scoring later.
    """
    out = df.copy()                                                                   # Work on a copy so the caller keeps the original table unchanged.
    # If no tissue inputs exist, default the tissue score to zero.
    if not tissue_embeddings_pt or not tissue_metadata_csv:                           
        out["tissue_score"] = 0.0                                                     # Fill the tissue score column with neutral zeros.
        return out                                                                    # Return early because there is nothing else to attach.

    # Read the optional tissue embedding tensor and metadata CSV.
    tissue_emb = torch.load(tissue_embeddings_pt, map_location="cpu")                 # Load the optional tissue embedding tensor from disk.
    if isinstance(tissue_emb, torch.Tensor):                                          # Convert tensors to NumPy arrays when necessary.
        tissue_emb = tissue_emb.numpy()                                               # Materialize the tensor on CPU as a NumPy array.
    tissue_emb = np.asarray(tissue_emb, dtype=np.float32)                             # Normalize the tissue matrix into a float NumPy array.
    tissue_meta = pd.read_csv(tissue_metadata_csv).copy()                             # Read the tissue metadata CSV aligned to the tissue tensor.
    # Look if a prior column exists in the tissue metadata (even when missing from the file), and default to zero if not.
    tissue_meta["tissue_score_prior"] = (                                             
        tissue_meta["tissue_score_prior"] if "tissue_score_prior" in tissue_meta.columns else 0.0
    )
    # Extract the tissue prior from the metadata into a NumPy array and calculate the mean, default to zero if not present.
    priors = tissue_meta["tissue_score_prior"].to_numpy(dtype=float) if len(tissue_meta) else np.zeros(1, dtype=float)  # Gather optional scalar tissue priors.
    mean_prior = float(priors.mean()) if priors.size else 0.0                         # Collapse the priors to one practical global score.
    # Attach the same mean tissue prior to every candidate for now.
    out["tissue_score"] = mean_prior                                                  
    return out                                                                        # Return the candidate table with an added tissue score column.


def ensure_provenance_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure ESM3-specific provenance columns always exist in the candidate table and derive a few report-friendly helpers.
    Returns the candidate table with added provenance columns.
    """
    out = df.copy()                                                                   # Work on a copy so the original merged table stays unchanged.

    # Add missing provenance columns with neutral placeholders to the candidate table so later reporting code can rely on a stable schema.
    defaults = {                                                                    
        "generator_mode": "unknown_generator",                                      # Record which backend produced the candidate.
        "esm3_model": pd.NA,                                                        # Record the Forge model when used; otherwise keep it missing.
        "esm3_prompt_sequence": pd.NA,                                              # Preserve the masked prompt passed to the generator.
        "esm3_hotspot_positions": pd.NA,                                            # Preserve hotspot positions under an ESM3-specific label.
        "esm3_temperature": pd.NA,                                                  # Preserve generation temperature for auditability.
        "esm3_top_k": pd.NA,                                                        # Preserve local / ESM2 top-k when relevant.
        "esm3_num_steps": pd.NA,                                                    # Preserve iterative Forge unmasking steps when relevant.
        "esm3_sampling_seed": pd.NA,                                                # Preserve the RNG seed used during generation.
        "mutation_positions": pd.NA,                                                # Preserve the compact mutation summary if available.
        "editable_hotspot_count": 0,                                                # Preserve how many sites were editable.
    }
    for col, default in defaults.items():                                             # Visit each provenance field that downstream reporting expects.
        if col not in out.columns:                                                    # Create only the columns that are still missing.
            out[col] = default                                                        # Fill missing provenance fields with a neutral placeholder.

    # Derive small helper fields (mostly flags) that make report summaries and filtering easier without changing the ranking formula itself.
    out["used_esm3_api"] = out["generator_mode"].astype(str).str.startswith("esm3_api:")  # Flag candidates produced by the real Forge ESM3 backend.
    out["used_local_generator"] = out["generator_mode"].astype(str).eq("local_conditional_generator")  # Flag candidates produced by the local conditioned generator.
    out["used_esm2_fallback"] = out["generator_mode"].astype(str).str.startswith("esm_masked_lm:")  # Flag candidates produced by the ESM2 fallback backend.
    out["provenance_complete"] = (                                                   # Mark rows that contain the minimum useful ESM-style provenance set.
        out["generator_mode"].notna() &
        out["esm3_prompt_sequence"].notna() &
        out["mutation_positions"].notna()
    )
    return out                                                                        # Return the candidate table with stable provenance columns and helper flags.


# Main script: Reads Stage 07b generated candidates, merges inputs and writes the final ranked CSVs.
def main() -> None:
    # Read the raw generated-candidate table and the structure-aware scored candidate table.
    args = parse_args()                                                               # Parse command-line arguments.
    generated = pd.read_csv(args.generated_csv)                                       # Read the raw generated-candidate table.
    struct = pd.read_csv(args.structure_scored_csv)                                   # Read the structure-aware scored candidate table.

    # Merge the generated candidates with the structure-aware scores.
    df = generated.merge(                                                             # Merge generated candidates with structure-aware scores by sequence.
        struct.drop(columns=[c for c in generated.columns if c in struct.columns and c != "candidate_sequence"]),  # Drop duplicate non-key columns before merging.
        on="candidate_sequence",                                                      # Merge on the generated candidate sequence.
        how="left",                                                                   # Keep all generated rows even if some structure scores are missing.
    )
    # Optionally attach Stage 05 host/validity terms and tissue context priors.
    df = merge_optional_scores(df, args.host_validity_csv)                            # Attach Stage 05-style host-transfer and validity terms when available.
    df = attach_tissue_scores(df, args.tissue_embeddings_pt, args.tissue_metadata_csv)  # Attach an optional global tissue prior when available.
    df = ensure_provenance_columns(df)                                                # Normalize ESM3-specific provenance columns before final ranking and reporting.

    # Ensure the mutation penalty column exists always, filling with a neutral zero if missing.
    if "mutation_penalty" not in df.columns:                                         
        df["mutation_penalty"] = 0.0                                                  

    # Compute the final rank with and without tissue disabled, and enabled when inputs exist.
    ranked_no_tissue = rank_score_dataframe(df, tissue_enabled=False)                 # Build a no-tissue comparison table for ablation reporting.
    ranked_with_tissue = rank_score_dataframe(                                        # Build the main ranked table using tissue only when supplied.
        df,
        tissue_enabled=bool(args.tissue_embeddings_pt and args.tissue_metadata_csv),
    )

    # Save the final ranked candidate tables to CSV.
    out_path = Path(args.out_csv)                                                     # Convert the main output path into a Path object.
    out_path.parent.mkdir(parents=True, exist_ok=True)                                # Create the output directory if needed.
    ranked_with_tissue.to_csv(out_path, index=False)                                  # Save the main final ranked candidate table.
    ranked_no_tissue.to_csv(out_path.with_name("top_candidates_no_tissue.csv"), index=False)   # Save the no-tissue comparison table.
    ranked_with_tissue.to_csv(out_path.with_name("top_candidates_with_tissue.csv"), index=False)  # Save the tissue-enabled comparison table.
    print(f"Wrote: {out_path}")                                                       # Print the main output path for quick confirmation.


if __name__ == "__main__":                                                            # Standard Python entrypoint guard.
    main()                                                                            # Execute the final ranking CLI.
