"""Ranking helpers used by the final Stage 07 multimodal reranker."""

from __future__ import annotations                                                   # Delay type annotation evaluation for compact typing.
import numpy as np                                                                   # Normalize score columns numerically with NumPy.
import pandas as pd                                                                  # Work with ranking tables stored as DataFrames.


# A helper function that standardizes score columns safely.
def safe_zscore(x: pd.Series | np.ndarray) -> np.ndarray:
    """Return a z-scored numeric array, falling back to zeros for constant inputs."""
    arr = np.asarray(x, dtype=float)                                                  # Convert the input series or array into a float NumPy array.
    if arr.size == 0:                                                                 # Handle empty arrays gracefully.
        return arr                                                                    # Return the empty array unchanged.
    std = arr.std()                                                                   # Compute the standard deviation used for z-scoring.
    if std < 1e-12:                                                                   # Detect effectively constant inputs that would divide by zero.
        return np.zeros_like(arr)                                                     # Return neutral zeros when standardization is not meaningful.
    return (arr - arr.mean()) / std                                                   # Return the standard z-score for each element.


# Final ranking helper: combines all score terms into the final Stage 07 rank.
def rank_score_dataframe(df: pd.DataFrame, tissue_enabled: bool = True) -> pd.DataFrame:
    """Compute the final multimodal Stage 07 score and return a sorted DataFrame."""
    out = df.copy()                                                                   # Work on a copy so the caller keeps the original table unchanged.
    defaults = {                                                                      # Define neutral fallback values for any missing score columns.
        "target_score": 0.0,
        "strict_manifold_score": 0.0,
        "family_cosine": 0.0,
        "target_anchor_cosine": 0.0,
        "structure_score": 0.0,
        "tissue_score": 0.0,
        "mutation_penalty": 0.0,
    }
    for col, default in defaults.items():                                             # Ensure every ranking component exists before z-scoring.
        if col not in out.columns:                                                    # Check whether the current score column is missing.
            out[col] = default                                                        # Fill missing columns with a neutral default value.

    z_target = safe_zscore(out["target_score"])                                       # Standardize the host-transfer score.
    z_validity = safe_zscore(out["strict_manifold_score"])                            # Standardize the manifold-validity score (= How close a candidate is to any known RBP in the entire strict dataset).
    z_family = safe_zscore(out["family_cosine"])                                      # Standardize the family-similarity score.
    z_anchor = safe_zscore(out["target_anchor_cosine"])                               # Standardize the target-anchor similarity score.
    z_structure = safe_zscore(out["structure_score"])                                 # Standardize the structure-aware plausibility score.
    z_tissue = safe_zscore(out["tissue_score"]) if tissue_enabled else np.zeros(len(out), dtype=float)  # Disable tissue when the branch is absent.
    z_penalty = safe_zscore(out["mutation_penalty"])                                  # Standardize the mutation penalty term.

    out["final_multimodal_rank_score"] = (                                            # Combine all normalized terms into one final ranking score.
        0.28 * z_target +                                                             # Give the largest weight to host-transfer performance.
        0.16 * z_validity +                                                           # Reward candidates that stay near the validated manifold.
        0.14 * z_family +                                                             # Reward candidates that remain close to the scaffold family.
        0.12 * z_anchor +                                                             # Reward proximity to target-host anchor proteins.
        0.20 * z_structure +                                                          # Reward structure-aware plausibility strongly.
        0.10 * z_tissue -                                                             # Add optional tissue-context support when enabled.
        0.08 * z_penalty                                                              # Penalize unnecessarily large mutation burdens.
    )

    return out.sort_values(                                                           # Return the table sorted from best to worst candidate.
        ["final_multimodal_rank_score", "target_score", "structure_score"],           # Break ties with host score then structure score.
        ascending=[False, False, False],                                              # Sort all ranking columns in descending order.
    ).reset_index(drop=True)                                                          # Reset the index so the ranked output is clean and contiguous.
