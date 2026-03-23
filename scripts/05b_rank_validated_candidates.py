# This script takes the massive table from the 05_compute_validity_metrics script, 
# throws away the non-validated_pass=True candidates, 
# and scientifically ranks the survivors to the most promising, diverse and physically plausible.

from __future__ import annotations                                                      
import argparse                                                                         
from pathlib import Path                                                                
import numpy as np                                                                      
import pandas as pd                                                                     


def method_family(run_name: str) -> str:                                                
    """
    This helper extracts a concise method label from the long run name, such as "diversity", "novelty", etc.
    It normalizes long run labels for later grouping.
    """
    text = str(run_name)                                                                # Convert the input to a string for robust handling.
    if "novelty" in text and "diversity" in text:                                       # Detect the novelty + diversity runs first because they are more specific.
        return "diversity_plus_novelty"                                                 # Return a compact label for those runs.
    if "diversity" in text:                                                             # Detect diversity-aware runs.
        return "diversity"                                                              # Return a compact label for diversity-aware runs.
    if "novelty" in text:                                                               # Detect novelty-aware runs.
        return "novelty"                                                                # Return a compact label for novelty-aware runs.
    if "entropy" in text:                                                               # Detect entropy-based position-selection runs.
        return "entropy"                                                                # Return a compact entropy label.
    if "random" in text:                                                                # Detect random baseline runs.
        return "random"                                                                 # Return a compact random label.
    return text                                                                         # Fall back to the original text when no pattern matches.


def safe_zscore(series: pd.Series) -> pd.Series:                                       
    """
    This helper z-scores a numeric series safely, avoiding divide-by-zero and returning zeros when the column is constant.
    """
    numeric = pd.to_numeric(series, errors="coerce")                                    # Force numeric conversion before standardization.
    std = numeric.std(ddof=0)                                                           # Compute population-style standard deviation.
    if pd.isna(std) or std == 0:                                                        # Avoid divide-by-zero and all-NaN cases.
        return pd.Series(np.zeros(len(numeric)), index=numeric.index)                   # Return neutral zeros when z-scoring is not meaningful.
    return (numeric - numeric.mean()) / std                                             # Return the standard z-score.


def build_rank_score(df: pd.DataFrame) -> pd.Series:                                    
    """
    This helper creates one scalar rank score after filtering has already happened.
    It ranks candidates by target score and number of mutations, and then based on strict manifold score, PLL, and density.
    """
    target_z = df.groupby("target_host")["target_score"].transform(safe_zscore)                 # Standardize target score within each target host instead of globally (judged only against other candidates that aim for the same host).
    manifold_z = df.groupby("target_host")["strict_manifold_score"].transform(safe_zscore)      # Standardize strict manifold score within each target host.
    pll_series = df["pll_delta_candidate_minus_seed"].fillna(df["pll_delta_candidate_minus_seed"].median())  # Build a PLL series with median imputation before within-host z-scoring.
    pll_z = pll_series.groupby(df["target_host"]).transform(safe_zscore)                        # Standardize PLL deltas within each target host.
    mutation_penalty = df["actual_n_mutations"].fillna(df["n_mutations"]).fillna(99).astype(float) * 0.10  # Penalize exact sequence-derived mutation burden first, falling back only when needed.
    density_col = next((c for c in df.columns if c.startswith("strict_top") and c.endswith("_mean_cosine")), None)  # Detect the strict density column dynamically.
    density_z = df.groupby("target_host")[density_col].transform(safe_zscore) if density_col is not None else 0.0  # Reward candidates that stay inside dense parts of the strict manifold.
    return target_z + 0.75 * manifold_z + 0.35 * pll_z + 0.35 * density_z - mutation_penalty    # Add strict-density support while still favoring host score first.


def take_balanced_top(df: pd.DataFrame, per_target: int, per_method_cap: int) -> pd.DataFrame:  
    """
    This helper builds a balanced shortlist per target host.
    It prevents one method from dominating the shortlist.
    """
    selected_rows = []                                                                  # Create a list that will collect selected rows.
    for target_host, target_df in df.groupby("target_host", dropna=False):              # Work independently within each target host.
        chosen_indices = []                                                             # Create a list of row indices selected for this target.
        method_counts = {}                                                              # Track how many rows have been taken from each method family.
        for idx, row in target_df.sort_values(["rank_score", "target_score", "strict_manifold_score"], ascending=False).iterrows():  # Walk through candidates from strongest to weakest.
            method_label = row["method_family"]                                         # Read the normalized method family label.
            current_count = method_counts.get(method_label, 0)                          # Look up how many candidates from this method are already selected.
            if current_count >= per_method_cap:                                         # Prevent one method from monopolizing the shortlist.
                continue                                                                # Skip this row and look at the next candidate.
            chosen_indices.append(idx)                                                  # Select this candidate row.
            method_counts[method_label] = current_count + 1                             # Update the per-method counter.
            if len(chosen_indices) >= per_target:                                       # Stop once the target-host quota is filled.
                break                                                                   # Exit the inner ranking loop.
        selected_rows.append(target_df.loc[chosen_indices])                             # Append the selected rows for this target host.
    if len(selected_rows) == 0:                                                         # Handle the edge case of no surviving candidates.
        return df.head(0).copy()                                                        # Return an empty DataFrame with the same schema.
    return pd.concat(selected_rows, axis=0).reset_index(drop=True)                      # Concatenate the per-target selections into one shortlist.


def main() -> None:                                                                     
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Rank and shortlist validated PhageForge candidates.")  
    parser.add_argument("--repo_root", type=str, default=".", help="Path to the phageforge repository root.")  
    parser.add_argument("--input_csv", type=str, default="results/analysis/validity/master_candidates_validated.csv", help="Master candidate table with validity metrics.")  
    parser.add_argument("--output_dir", type=str, default="results/analysis/validity", help="Directory where shortlist outputs will be written.")  
    parser.add_argument("--max_mutations", type=int, default=12, help="Maximum mutation count allowed in shortlist mode.")  # Default to the corrected actual-mutation closeout budget.
    parser.add_argument("--min_validated_strict_cosine", type=float, default=0.995, help="Optional absolute floor on strict nearest-neighbor cosine.")  # Use the corrected strict cosine floor by default.
    parser.add_argument("--max_negative_pll_delta", type=float, default=0.0, help="Most negative candidate-vs-seed PLL delta allowed.")  # Default to the corrected naturalness floor.
    parser.add_argument("--per_target", type=int, default=12, help="How many final shortlist candidates to keep per target host.")  
    parser.add_argument("--per_method_cap", type=int, default=4, help="Maximum shortlist count contributed by a single method family per target host.")  
    parser.add_argument("--drop_targets", nargs="*", default=["Pseudomonas"], help="Optional targets to remove from the final closeout shortlist.")  # Allow noisy targets to be excluded from the default closeout story.
    args = parser.parse_args()                                                          

    # Resolve paths
    repo_root = Path(args.repo_root).resolve()                                          # Resolve the repository root path.
    input_csv = (repo_root / args.input_csv).resolve()                                  # Resolve the input master table path.
    output_dir = (repo_root / args.output_dir).resolve()                                # Resolve the output directory path.
    output_dir.mkdir(parents=True, exist_ok=True)                                       # Create the output directory if needed.

    # Load the master candidate table, add a compact method-family label, and filter out invalid candidates.
    df = pd.read_csv(input_csv)                                                         # Load the master candidate table with validity metrics.
    df = df.copy()                                                                      # Work on a copy so the raw table remains untouched.
    df["method_family"] = df["run_name"].map(method_family)                             # Add a compact method-family label for balancing.
    df["strict_manifold_score"] = df["strict_manifold_score"].fillna(df["strict_nn_cosine"])  # Fall back to nearest strict cosine if the composite score is missing.
    df["pll_delta_candidate_minus_seed"] = df["pll_delta_candidate_minus_seed"].astype(float)  # Normalize the PLL delta column type.
    df["n_mutations"] = df["n_mutations"].fillna(999).astype(int)                      # Normalize the mutation-count column type.
    if "actual_n_mutations" not in df.columns:                                         # Be backward compatible with older 05 outputs that may not yet contain corrected mutation counts.
        df["actual_n_mutations"] = df["n_mutations"]                                   # Fall back to the logged mutation field when the corrected one is absent.
    validated_col = "validated_pass"                                                    # Keep the canonical validated column name in the merged 05 pipeline.

    filtered = df[                                                                      # Apply the fast hard filters before any weighted ranking.
        (df[validated_col] == True) &                                                   # Filter on the canonical corrected validated flag produced by the merged 05 script.
        (df["actual_n_mutations"] <= args.max_mutations) &                              # Apply the mutation budget to the exact sequence-derived mutation count.
        (df["strict_nn_cosine"] >= args.min_validated_strict_cosine) &                  # Optionally enforce an absolute manifold floor.
        (df["pll_delta_candidate_minus_seed"].fillna(0.0) >= args.max_negative_pll_delta)  # Reject severe local-naturalness collapses.
    ].copy()                                                                            # Materialize the filtered shortlist pool.

    if len(args.drop_targets) > 0:                                                      # Optionally remove noisy targets before deduplication and ranking.
        filtered = filtered[~filtered["target_host"].isin(args.drop_targets)].copy()     # Drop user-specified targets from the closeout shortlist pool.

    if len(filtered) == 0:                                                              # Fail early when no candidates survive filtering.
        raise ValueError("No candidates survived the shortlist filters. Relax the thresholds slightly and rerun.")  # Provide a useful action-oriented error.

    # Deduplicate by unique target-host sequence before ranking so repeated copies do not dominate the shortlist.
    filtered = filtered.sort_values(["target_host", "target_score", "strict_manifold_score", "pll_delta_candidate_minus_seed"], ascending=[True, False, False, False]).copy()  # Sort first so the strongest copy of each target/sequence pair is retained.
    filtered = filtered.drop_duplicates(subset=["target_host", "aa_sequence"], keep="first").copy()  # Keep only one row per unique sequence per target host.

    # Rank the remaining candidates by rank score, then target host, then target score, then strict manifold score
    filtered["rank_score"] = build_rank_score(filtered)                                  # Compute one scalar rank score for the surviving candidates.
    filtered = filtered.sort_values(["target_host", "rank_score", "target_score", "strict_manifold_score"], ascending=[True, False, False, False]).reset_index(drop=True)  # Sort candidates from strongest to weakest within each target host.

    # Shortlist the remaining candidates in a dataframe, keep two extra dataframes with one best candidate per run and per target, and save them all.
    final_shortlist = take_balanced_top(df=filtered, per_target=args.per_target, per_method_cap=args.per_method_cap)  # Build a method-balanced final shortlist.
    best_per_run = filtered.sort_values(["rank_score", "target_score"], ascending=False).groupby("run_name", dropna=False).head(1).reset_index(drop=True)  # Keep the strongest surviving candidate from each run.
    best_per_target = filtered.sort_values(["rank_score", "target_score"], ascending=False).groupby("target_host", dropna=False).head(args.per_target).reset_index(drop=True)  # Keep the top rows per target regardless of method balancing.

    filtered.to_csv(output_dir / "validated_pool_ranked.csv", index=False)              # Save the full filtered and ranked candidate pool.
    final_shortlist.to_csv(output_dir / "final_validated_shortlist.csv", index=False)   # Save the balanced final shortlist.
    best_per_run.to_csv(output_dir / "best_validated_candidate_per_run.csv", index=False)  # Save one best validated row per run.
    best_per_target.to_csv(output_dir / "best_validated_candidates_per_target.csv", index=False)  # Save the strongest rows per target host.

    # Build a concise summary of the final shortlist
    summary = final_shortlist.groupby(["target_host", "method_family"], dropna=False).agg(  # Build a concise shortlist summary for the writeup.
        n_selected=("aa_sequence", "size"),                                              # Count how many shortlist rows came from each method family.
        n_unique_sequences=("aa_sequence", pd.Series.nunique),                            # Report the number of unique sequences selected by each method family.
        mean_target_score=("target_score", "mean"),                                       # Average the target-host score within the selected rows.
        mean_strict_manifold=("strict_manifold_score", "mean"),                           # Average the strict manifold score within the selected rows.
        mean_pll_delta=("pll_delta_candidate_minus_seed", "mean"),                         # Average the local PLL delta within the selected rows.
        mean_actual_mutations=("actual_n_mutations", "mean"),                             # Report the exact sequence-derived mutation burden of the shortlist.
    ).reset_index()                                                                       # Flatten the grouped summary into a DataFrame.
    summary.to_csv(output_dir / "final_shortlist_summary.csv", index=False)              # Save the shortlist summary table.

    print("[INFO] Saved:")                                                               # Print a concise success summary.
    print(f"  - {output_dir / 'validated_pool_ranked.csv'}")                               # Show the ranked pool path.
    print(f"  - {output_dir / 'final_validated_shortlist.csv'}")                            # Show the final shortlist path.
    print(f"  - {output_dir / 'best_validated_candidate_per_run.csv'}")                     # Show the best-per-run path.
    print(f"  - {output_dir / 'best_validated_candidates_per_target.csv'}")                 # Show the best-per-target path.
    print(f"  - {output_dir / 'final_shortlist_summary.csv'}")                              # Show the summary path.


if __name__ == "__main__":                                                               # Run the CLI only when executed directly.
    main()                                                                                # Invoke the main function.
