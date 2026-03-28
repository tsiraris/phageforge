"""
Phase A / Step 06a: Building the design space
===================

This script selects the canonical seed and scaffold family for Phase A, then derives
host-ladder mutation priors from the corrected Stage 05 validity outputs.

Phase A in this project is intentionally narrow and fast:
    seed RBP -> Enterobacter -> Acinetobacter

The script does five things:
1. reads the corrected shortlist / master tables from Stage 05,
2. selects the best seed to carry into the host ladder,
3. constructs a seed-local scaffold family from the strict RBP bank,
4. extracts mutation hotspots from corrected shortlisted candidates,
5. writes a reusable `phaseA_plan.json` to be used by the constrained optimizer.
"""

from __future__ import annotations
import argparse                                                                        # Parse command-line arguments.
import json                                                                            # Read and write JSON planning artifacts.
import math                                                                            # Use a few small mathematical helpers.
import re                                                                              # Parse mutation strings safely.
from collections import Counter, defaultdict                                           # Count hosts, mutations, and hotspot frequencies.
from pathlib import Path                                                               # Work with filesystem paths cleanly.
from typing import Dict, List, Tuple                                                   # Make the code easier to read.
import numpy as np                                                                     # Numerical work for embeddings and priors.
import pandas as pd                                                                    # Table manipulation for shortlist and reference banks.
import torch                                                                           # Load the cached embedding tensor.


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for Phase A family selection."""
    ap = argparse.ArgumentParser(description="Select canonical Phase A family and mutation priors.")
    ap.add_argument("--corrected_shortlist_csv", type=str, required=True, help="Corrected final shortlist CSV from Stage 05.")
    ap.add_argument("--corrected_master_csv", type=str, required=True, help="Corrected master validated CSV from Stage 05.")
    ap.add_argument("--strict_csv", type=str, required=True, help="Strict RBP dataset CSV.")
    ap.add_argument("--embedding_pt", type=str, required=True, help="Torch tensor containing cached ESM embeddings.")
    ap.add_argument("--embedding_index_csv", type=str, required=True, help="CSV mapping proteins to embedding rows.")
    ap.add_argument("--preferred_targets", nargs="+", default=["Enterobacter", "Acinetobacter"], help="Targets to include in the Phase A ladder.")
    ap.add_argument("--family_cosine_floor", type=float, default=0.995, help="Minimum cosine similarity to include a strict RBP in the local scaffold family.")
    ap.add_argument("--family_top_k", type=int, default=128, help="Maximum number of family members to retain after cosine filtering.")
    ap.add_argument("--length_tolerance", type=float, default=0.15, help="Allowed relative seed-length deviation for family members.")
    ap.add_argument("--hotspot_min_count", type=int, default=2, help="Minimum occurrence count in shortlisted candidates for a position to become a hotspot.")
    ap.add_argument("--window_flank", type=int, default=16, help="Extra residues added on each side of the hotspot span to form a mutation window.")
    ap.add_argument("--target_anchor_top_k", type=int, default=16, help="How many nearest strict-bank references to retain per target host for target-anchor centroids.")
    ap.add_argument("--target_anchor_cosine_floor", type=float, default=0.94, help="Minimum seed cosine similarity for target-host anchor references when available.")
    ap.add_argument("--output_dir", type=str, required=True, help="Directory where Phase A planning outputs will be written.")
    return ap.parse_args()


def parse_mutation_positions(mutation_text: str) -> List[int]:
    """Convert mutation strings like 'W197R;K260P' into 0-based positions (1-based biological numbering to 0-based Python indexing)."""
    if pd.isna(mutation_text):                                                          # Missing mutation strings should simply produce no positions.
        return []
    positions: List[int] = []                                                           # Collect parsed 0-based residue indices here.
    for token in str(mutation_text).split(";"):                                        # One mutation token is expected per semicolon-separated segment.
        token = token.strip()                                                           # Remove whitespace around each token.
        if token == "":                                                                # Ignore accidental empty segments.
            continue
        match = re.search(r"(\d+)", token)                                             # Extract the numeric residue index from the mutation token.
        if match is None:                                                               # Skip malformed mutation tokens instead of crashing.
            continue
        positions.append(int(match.group(1)) - 1)                                       # Convert biological 1-based numbering into Python 0-based numbering.
    return sorted(set(p for p in positions if p >= 0))                                  # Deduplicate and keep only valid positions.


def l2_normalize(x: np.ndarray) -> np.ndarray:
    """L2-normalize input array rows so cosine similarity becomes a dot product."""
    norms = np.linalg.norm(x, axis=1, keepdims=True)                                    # Compute one vector norm per row.
    norms = np.clip(norms, 1e-12, None)                                                 # Prevent division by zero on degenerate rows.
    return x / norms                                                                    # Return normalized rows.


def load_embedding_bank(embedding_pt: Path, embedding_index_csv: Path) -> pd.DataFrame:
    """Load the cached embedding tensor of all (broad dataset) proteins and attach it along with row IDs to the returned index DataFrame, which already contains metadata such as host genus, protein ID, etc."""
    emb = torch.load(embedding_pt, map_location="cpu")                                  # Load the cached embeddings onto CPU so any machine can run Step 06a.
    if isinstance(emb, torch.Tensor):                                                   # The cache is expected to be a tensor in this repo.
        emb_array = emb.detach().cpu().numpy()                                          # Convert tensors into NumPy arrays for fast cosine computations.
    else:
        emb_array = np.asarray(emb)                                                     # Accept NumPy-like caches as well.
    index_df = pd.read_csv(embedding_index_csv)                                         # Load metadata for each embedding row.
    row_col = "row_id" if "row_id" in index_df.columns else "idx"                       # Support both row_id and idx column names.
    index_df["embedding_row"] = index_df[row_col].astype(int)                           # Normalize the row pointer into one canonical column.
    index_df["protein_id"] = index_df["protein_id"].astype(str)                         # Normalize protein IDs for joins.
    index_df["host_genus"] = index_df["host_genus"].astype(str)                         # Normalize host labels for grouping and reporting.
    index_df.attrs["embedding_array"] = emb_array                                       # Attach the embedding matrix so helper functions can recover it later.
    return index_df


def choose_canonical_seed(shortlist_df: pd.DataFrame, preferred_targets: List[str]) -> Tuple[str, dict]:
    """
    Read the shortlisted candidates and decide which protein seed that is the strongest across the preferred Phase A targets (= hosts).

    The selection score is intentionally simple and interpretable:
    - First by how many preferred targets a seed has,
    - Second by the mean corrected rank score across preferred targets,
    - Third by the mean target score across preferred targets.
    
    Returns a tuple of (seed_protein_id, seed_metadata dict). 
    """
    df = shortlist_df[shortlist_df["target_host"].isin(preferred_targets)].copy()           # Keep only the hosts that matter for Phase A.
    if df.empty:
        raise ValueError("No corrected shortlist rows found for the requested preferred targets.")

    # Groups shortlisted candidates by seed protein ID and summarize each seed's performance across the preferred targets.
    grouped = (
        df.groupby("seed_protein_id", dropna=False)
        .agg(
            n_rows=("aa_sequence", "size"),
            n_unique_targets=("target_host", "nunique"),                    # Count for how many different target hosts each seed generated valid candidates.
            mean_corrected_rank_score=("corrected_rank_score", "mean"),
            mean_target_score=("target_score", "mean"),
        )
        .reset_index()
    )  
    
    # Rank seeds by unique target coverage first, and then by the strongest overall score
    grouped = grouped.sort_values(
        by=["n_unique_targets", "n_rows", "mean_corrected_rank_score", "mean_target_score"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)                                                                # Rank seeds by target coverage first, then by strength.

    best_seed = str(grouped.iloc[0]["seed_protein_id"])                                     # Select the top-ranked seed protein ID.
    summary = grouped.iloc[0].to_dict()                                                     # Keep a JSON-serializable summary of why it was chosen.
    return best_seed, summary


def build_strict_family(
    strict_df: pd.DataFrame,
    embedding_index_df: pd.DataFrame,
    seed_protein_id: str,
    family_cosine_floor: float,
    family_top_k: int,
    length_tolerance: float,
) -> Tuple[pd.DataFrame, np.ndarray, dict]:
    """
    Construct a local scaffold family inside the strict RBP bank through a practical,
    local design manifold around the chosen seed that includes proteins that are: 
    high-similarity in embedding space, length-compatible, and still diverse enough to give target-host anchors.
    Returns a tuple of (family_df, family_emb, family_summary).
    """
    strict_df = strict_df.copy()                                                       # Work on a copy so the caller's table remains unchanged.
    strict_df["protein_id"] = strict_df["protein_id"].astype(str)                      # Normalize identifiers before joining to embeddings.
    strict_df["aa_sequence"] = strict_df["aa_sequence"].astype(str)                    # Normalize sequences before measuring lengths.
    strict_df["aa_len"] = strict_df["aa_sequence"].str.len()                           # Cache sequence lengths for the family-length filter.

    # Join the strict RBP bank dataframe to the embedding index dataframe
    merged = strict_df.merge(
        embedding_index_df[["protein_id", "embedding_row"]],
        on="protein_id",
        how="inner",
    )  
    # Keep only strict-bank rows that also have cached embeddings                                                                                  
    if merged.empty:
        raise ValueError("The strict RBP bank could not be aligned to the embedding index.")

    # Gather the embeddings and normalize them
    emb_array = embedding_index_df.attrs["embedding_array"]                             # Recover the full embedding matrix attached to the index DataFrame.
    merged_emb = emb_array[merged["embedding_row"].to_numpy()]                          # Gather embeddings aligned to the merged metadata rows.
    merged_emb = l2_normalize(np.asarray(merged_emb, dtype=np.float32))                 # Normalize once so later similarities are dot products.

    # Find the chosen seed inside the strict bank, cache its length, and slice its embedding as shape [1, H]
    seed_rows = merged[merged["protein_id"] == seed_protein_id]                         # Find the chosen seed inside the strict bank.
    if seed_rows.empty:
        raise ValueError(f"Chosen seed '{seed_protein_id}' was not found in the strict dataset with embeddings.")
    seed_row = seed_rows.iloc[0]                                                         # Use the first matching row as the canonical seed record.
    seed_idx = int(seed_row.name)                                                        # This row index matches the merged DataFrame indexing.
    seed_len = int(seed_row["aa_len"])                                                   # Cache seed length for the family length filter.
    seed_emb = merged_emb[seed_idx : seed_idx + 1]                                       # Slice the seed embedding as shape [1, H] for broadcasting.

    # Compute cosine similarity from every strict RBP to the seed, store it on the metadata table, and define the family manifold length bounds
    cosine = (merged_emb @ seed_emb.T).reshape(-1)                                       # Compute cosine similarity from every strict RBP to the seed.
    merged["seed_family_cosine"] = cosine                                                # Store the family similarity on the metadata table.
    min_len = math.floor(seed_len * (1.0 - length_tolerance))                            # Lower length bound of the proteins allowed in the family manifold.
    max_len = math.ceil(seed_len * (1.0 + length_tolerance))                             # Upper length bound allowed in the family manifold.
    
    # Find the family members: keep only close and length-compatible strict RBPs, sort by similarity score, and confine it to the family size
    family_df = merged[
        (merged["seed_family_cosine"] >= family_cosine_floor) &
        (merged["aa_len"] >= min_len) &
        (merged["aa_len"] <= max_len)
    ].copy()                                                                             
    if family_df.empty:
        family_df = merged.sort_values("seed_family_cosine", ascending=False).head(min(family_top_k, len(merged))).copy()  # Fall back to top-k nearest neighbors if the floor is too strict.
    else:
        family_df = family_df.sort_values("seed_family_cosine", ascending=False).head(family_top_k).copy()                 # Limit the number of rows to the top-k nearest neighbors (family size) to keep later scoring efficient.

    # Compute the family manifold embedding center and create a summary of the family
    family_emb = merged_emb[family_df.index.to_numpy()]                                  # Recover the family embedding rows aligned to family_df.
    family_centroid = l2_normalize(family_emb.mean(axis=0, keepdims=True)).reshape(-1)   # Compute one normalized centroid that defines the family manifold center.
    host_counts = family_df["host_genus"].value_counts().to_dict()                       # Summarize which host genera populate this local family.
    summary = {
        "seed_length": seed_len,
        "family_size": int(len(family_df)),
        "family_cosine_floor": float(family_cosine_floor),
        "family_top_k": int(family_top_k),
        "length_tolerance": float(length_tolerance),
        "host_counts": {str(k): int(v) for k, v in host_counts.items()},
    }                                                                                    # Prepare a compact JSON summary used in the plan file.
    return family_df.reset_index(drop=True), family_centroid, summary                    # Return the family DataFrame, its centroid, and the summary of the family.


def build_target_reference_summary(
    strict_df: pd.DataFrame,
    embedding_index_df: pd.DataFrame,
    seed_protein_id: str,
    preferred_targets: List[str],
    length_tolerance: float,
    target_anchor_top_k: int,
    target_anchor_cosine_floor: float,
) -> Tuple[dict, dict]:
    """
    Build per-target anchor information from the full strict bank, not only the seed-host family subset, searching for length-compatible, seed-near target-host RBPs.
    Returns a list of the centroids, and a summary of the metadata per target.
    """
    strict_df = strict_df.copy()                                                       # Work on a copy so the caller's table remains unchanged.
    strict_df["protein_id"] = strict_df["protein_id"].astype(str)                      # Normalize identifiers before joining to embeddings.
    strict_df["aa_sequence"] = strict_df["aa_sequence"].astype(str)                    # Normalize sequences before measuring lengths.
    strict_df["aa_len"] = strict_df["aa_sequence"].str.len()                           # Cache sequence lengths for filtering.
    
    # Merge strict-bank rows with cached embeddings, gather embeddings, and normalize
    merged = strict_df.merge(
        embedding_index_df[["protein_id", "embedding_row"]],
        on="protein_id",
        how="inner",
    )                                                                                   # Keep only strict-bank rows with cached embeddings.
    emb_array = embedding_index_df.attrs["embedding_array"]                             # Recover the raw embedding matrix.
    merged_emb = emb_array[merged["embedding_row"].to_numpy()]                          # Gather embeddings aligned to the merged strict-bank rows.
    merged_emb = l2_normalize(np.asarray(merged_emb, dtype=np.float32))                 # Normalize once for cosine similarity.

    # Find the canonical seed inside the strict bank, slice its embedding, cache its length, and define the target anchor length bounds
    seed_rows = merged[merged["protein_id"] == seed_protein_id]                         # Locate the canonical seed inside the strict bank.
    if seed_rows.empty:
        raise ValueError(f"Chosen seed '{seed_protein_id}' was not found in the strict bank when building target anchors.")
    seed_row = seed_rows.iloc[0]                                                         # Use the first seed row as the canonical anchor.
    seed_emb = merged_emb[int(seed_row.name) : int(seed_row.name) + 1]                   # Slice the seed embedding as shape [1, H].
    seed_len = int(seed_row["aa_len"])                                                   # Cache seed length for length-compatible target-anchor filtering.
    min_len = math.floor(seed_len * (1.0 - length_tolerance))                            # Lower length bound for plausible target anchors.
    max_len = math.ceil(seed_len * (1.0 + length_tolerance))                             # Upper length bound for plausible target anchors.

    # Measure similarity between the canonical seed and each strict-bank RBP, and keep only length-compatible target-anchor candidates
    merged["seed_anchor_cosine"] = (merged_emb @ seed_emb.T).reshape(-1)                 
    candidate_pool = merged[(merged["aa_len"] >= min_len) & (merged["aa_len"] <= max_len)].copy()  

    # For each of the preferred targets (i.e. ["Acinetobacter","Enterococcus"]): 
    # - Sort the target-host candidates by similarity to the canonical seed, and keep only the closest to the seed.
    # - Limit the number of candidates to the top-k nearest neighbors (target size) to keep later scoring efficient.
    # - Compute the target-host anchor embedding center (centroid) and create a summary of the target-host anchor candidates.
   
    target_centroids = {}                                                                # Store one normalized centroid vector per target host.
    target_refs = {}                                                                     # Store the reference rows that define each target centroid.
    for target in preferred_targets:                                                     # Build one target-anchor summary per ladder host.
        sub = candidate_pool[candidate_pool["host_genus"] == target].copy()              # Keep only rows belonging to the current target host.
        if sub.empty:
            target_centroids[target] = None                                              # Some targets may not exist in the strict bank subset.
            target_refs[target] = []
            continue
        sub = sub.sort_values("seed_anchor_cosine", ascending=False)                     # Rank the target-host references by similarity to the canonical seed.
        keep = sub[sub["seed_anchor_cosine"] >= target_anchor_cosine_floor].copy()       # Prefer only target-host references close enough to the seed.
        if keep.empty:
            keep = sub.head(target_anchor_top_k).copy()                                  # Fall back to the nearest target-host rows even if they miss the cosine floor.
        else:
            keep = keep.head(target_anchor_top_k).copy()                                 # Cap the number of target-host anchors to a manageable size.
        keep_emb = emb_array[keep["embedding_row"].astype(int).to_numpy()]               # Gather embeddings aligned to the chosen target-host anchor rows.
        centroid = l2_normalize(np.asarray(keep_emb, dtype=np.float32).mean(axis=0, keepdims=True)).reshape(-1)  # Compute a normalized target-host centroid.
        target_centroids[target] = centroid.tolist()                                     # Store as list for JSON serialization.
        target_refs[target] = (
            keep[["protein_id", "virus_accession", "host_genus", "product", "aa_len", "seed_anchor_cosine"]]
            .to_dict(orient="records")
        )                                                                                # Save the exact anchor proteins defining this target centroid.
    return target_centroids, target_refs


def derive_hotspot_priors(
    shortlist_df: pd.DataFrame,
    master_df: pd.DataFrame,
    seed_protein_id: str,
    preferred_targets: List[str],
    seed_length: int,
    hotspot_min_count: int,
    window_flank: int,
) -> Tuple[List[int], dict, pd.DataFrame]:
    """
    Derive mutation hotspots and per-target position priors (probability map) from corrected shortlisted candidates.

    Strategy:
    - use `actual_mutation_list` when present because it reflects sequence-vs-seed reality,
    - fall back to the legacy `mutations` column only if necessary,
    - summarize frequency per target host,
    - create a global mutation window spanning all recurrent positions plus a small flank,
    - create smoothed per-target position priors over the full seed length.
    """
    # Restrict the shortlist dataframe to the chosen seed and ladder targets.
    shortlist_seed = shortlist_df[
        (shortlist_df["seed_protein_id"].astype(str) == seed_protein_id) &
        (shortlist_df["target_host"].isin(preferred_targets))
    ].copy()   
    
    # Keep a broader master candidate view for context and frequency estimation.                                                                           
    master_seed = master_df[
        (master_df["seed_protein_id"].astype(str) == seed_protein_id) &
        (master_df["target_host"].isin(preferred_targets))
    ].copy()                                                                              
    if shortlist_seed.empty:
        raise ValueError("No corrected shortlist rows matched the chosen seed and preferred targets.")

    hotspot_rows: List[dict] = []                                                        # Collect a tidy hotspot table that will be written to CSV.
    target_priors: Dict[str, List[float]] = {}                                           # Store one full-length smoothed prior vector per target host.
    recurrent_union: set[int] = set()                                                    # Track the union of recurrent hotspot positions across targets.

    # For each of the preferred targets (i.e. ["Acinetobacter","Enterococcus"]):
    for target in preferred_targets:    
        # Keep only the shortlist and master rows for the current target host
        target_short = shortlist_seed[shortlist_seed["target_host"] == target].copy()    # Use shortlist rows as the primary evidence source.
        target_master = master_seed[master_seed["target_host"] == target].copy()         # Use master rows to stabilize priors with more examples.

        # Count recurrent positions in the shortlist and the master set
        short_counter: Counter = Counter()                                               # Count recurrent positions in the shortlist.
        master_counter: Counter = Counter()                                              # Count positions in the broader validated/master set.

        # Parse the shortlist and master mutation strings and count recurrent positions
        for _, row in target_short.iterrows():                                           # Parse the shortlist mutation strings first because they drive the host ladder.
            mutation_text = row.get("actual_mutation_list", row.get("mutations", ""))    # Prefer corrected mutations over legacy logged ones.
            for pos in parse_mutation_positions(mutation_text):                          # Counts how many times each position was mutated in the shortlist.
                short_counter[pos] += 1                                                  # Increment shortlist hotspot evidence.
        
        for _, row in target_master.iterrows():                                          # Parse the broader set to build smoother position priors.
            mutation_text = row.get("actual_mutation_list", row.get("mutations", ""))    # Again prefer corrected mutation lists.
            for pos in parse_mutation_positions(mutation_text):                          # Counts how many times each position was mutated in the master set.
                master_counter[pos] += 1                                                 # Increment broader evidence counts.

        # If an amino acid position was mutated more times than hotspot_min_count, it is a recurrent hotspot, and gets added to the global recurrent_union.
        recurrent = sorted([p for p, c in short_counter.items() if c >= hotspot_min_count])  
        recurrent_union.update(recurrent)                                                   

        # Build a smoothed per-target position prior (probability distribution map) over the full seed length
        prior = np.ones(seed_length, dtype=np.float64)                                   # Start with a small uniform prior so unseen positions never get zero probability.
        for pos, count in master_counter.items():                                        # For each position in the master set
            if 0 <= pos < seed_length:                                                   # If the position is within the seed length       
                prior[pos] += float(count)                                               # Add the master target-specific evidence counts onto the prior.
        if recurrent:                                                                    # Recurrent shortlist hotspots deserve an extra boost because they are strongest evidence.
            for pos in recurrent:
                prior[pos] += float(hotspot_min_count)                                   # Reinforce (=increase the value of) recurrent positions without hard-coding them as the only options.
        prior = prior / prior.sum()                                                      # Normalize into a proper probability distribution over positions.
        target_priors[target] = prior.tolist()                                           # Convert to a JSON-serializable list.

        # Write one row per target/position pair to the hotspot table
        for pos in sorted(set(list(short_counter.keys()) + list(master_counter.keys()))):
            hotspot_rows.append({
                "target_host": target,
                "position_0based": int(pos),
                "position_1based": int(pos + 1),
                "shortlist_count": int(short_counter.get(pos, 0)),
                "master_count": int(master_counter.get(pos, 0)),
                "is_recurrent_hotspot": bool(pos in recurrent),
            })                                                                            # Write one tidy row per target/position pair for later inspection.

    # If no recurrent hotspot positions were found from the shortlist, fall back to the strongest shortlist positions regardless the threshold.
    if not recurrent_union:                                                              # If no position met the recurrence cutoff, fall back to the strongest shortlist positions.
        fallback_counter: Counter = Counter()                                            # Aggregate shortlist evidence across the preferred targets.
        for _, row in shortlist_seed.iterrows():
            mutation_text = row.get("actual_mutation_list", row.get("mutations", ""))
            for pos in parse_mutation_positions(mutation_text):
                fallback_counter[pos] += 1
        top_fallback = [p for p, _ in fallback_counter.most_common(12)]                  # Keep a modest set of strongest fallback positions.
        recurrent_union.update(top_fallback)
    
    # If the shortlist hotspot union is still empty, fall back to the full seed
    if not recurrent_union:                                                              # If even the fallback is empty, widen all the way to the full seed.
        mutation_window_positions = list(range(seed_length))                             # This should be rare, but guarantees the optimizer still runs.
    else:
        # Define the mutation window based on the recurrent hotspot positions
        start = max(0, min(recurrent_union) - window_flank)                              # Expand the hotspot span to the left.
        end = min(seed_length - 1, max(recurrent_union) + window_flank)                  # Expand the hotspot span to the right.
        mutation_window_positions = list(range(start, end + 1))                          # Build the final allowed mutation window.

    # Write the hotspot table sorted by target, hotspot strength, and position
    hotspot_df = pd.DataFrame(hotspot_rows).sort_values(
        by=["target_host", "is_recurrent_hotspot", "shortlist_count", "master_count", "position_0based"],
        ascending=[True, False, False, False, True],
    ).reset_index(drop=True)                                                             # Sort the hotspot table for readable CSV output.
    return mutation_window_positions, target_priors, hotspot_df


def main():
    args = parse_args()                                                                  # Parse command-line arguments first.

    output_dir = Path(args.output_dir)                                                   # Convert the output directory into a Path object.
    output_dir.mkdir(parents=True, exist_ok=True)                                        # Ensure the directory exists before writing any files.

    shortlist_df = pd.read_csv(args.corrected_shortlist_csv)                             # Load the corrected final shortlist from Stage 05.
    master_df = pd.read_csv(args.corrected_master_csv)                                   # Load the corrected master validated table for broader evidence.
    strict_df = pd.read_csv(args.strict_csv)                                             # Load the strict RBP bank because Phase A stays inside this validated manifold.
    embedding_index_df = load_embedding_bank(Path(args.embedding_pt), Path(args.embedding_index_csv))  # Load cached embeddings and the matching index table.

    # Pick the seed that is strongest across the ladder targets.
    seed_protein_id, seed_selection_summary = choose_canonical_seed(shortlist_df, args.preferred_targets)   
    # Recover the seed metadata from the strict bank.
    strict_seed_rows = strict_df[strict_df["protein_id"].astype(str) == seed_protein_id].copy()             
    if strict_seed_rows.empty:
        raise ValueError(f"Chosen seed '{seed_protein_id}' was not found in the strict dataset.")
    seed_row = strict_seed_rows.iloc[0]                                                  # Use the first matching strict-bank seed row as the canonical record.

    # Construct the local family manifold around the chosen seed.
    family_df, family_centroid, family_summary = build_strict_family(
        strict_df=strict_df,
        embedding_index_df=embedding_index_df,
        seed_protein_id=seed_protein_id,
        family_cosine_floor=args.family_cosine_floor,
        family_top_k=args.family_top_k,
        length_tolerance=args.length_tolerance,
    )           
    # Build target-host anchors inside the chosen family whenever possible.                                                                         
    target_centroids, target_refs = build_target_reference_summary(
        strict_df=strict_df,
        embedding_index_df=embedding_index_df,
        seed_protein_id=seed_protein_id,
        preferred_targets=args.preferred_targets,
        length_tolerance=args.length_tolerance,
        target_anchor_top_k=args.target_anchor_top_k,
        target_anchor_cosine_floor=args.target_anchor_cosine_floor,
    )                                                                                    
    # Derive the mutation window and target-specific position priors from corrected shortlist evidence.
    mutation_window_positions, target_priors, hotspot_df = derive_hotspot_priors(
        shortlist_df=shortlist_df,
        master_df=master_df,
        seed_protein_id=seed_protein_id,
        preferred_targets=args.preferred_targets,
        seed_length=len(str(seed_row["aa_sequence"])),
        hotspot_min_count=args.hotspot_min_count,
        window_flank=args.window_flank,
    )                                                                                    
    # Save the local family members (which proteins define the manifold) and hotspot summary (for target-specific position evidence) for Phase A.
    family_df.to_csv(output_dir / "phaseA_family_members.csv", index=False)             
    hotspot_df.to_csv(output_dir / "phaseA_hotspots_by_target.csv", index=False)         

    # Save a compact target-anchor summary for quick inspection.
    pd.DataFrame([
        {
            "target_host": target,
            "n_target_refs": len(target_refs.get(target, [])),
            "top_reference_protein_id": target_refs.get(target, [{}])[0].get("protein_id", "") if len(target_refs.get(target, [])) > 0 else "",
            "top_reference_product": target_refs.get(target, [{}])[0].get("product", "") if len(target_refs.get(target, [])) > 0 else "",
        }
        for target in args.preferred_targets
    ]).to_csv(output_dir / "phaseA_target_reference_summary.csv", index=False)           

    # Build the Phase A plan by collecting everything the optimizer needs into one reusable JSON plan.
    phaseA_plan = {
        "phase_name": "Phase A",
        "goal": "Family-conditioned host ladder optimization from canonical seed to Enterobacter then Acinetobacter.",
        "preferred_targets": list(args.preferred_targets),
        "canonical_seed": {
            "seed_protein_id": str(seed_protein_id),
            "virus_accession": str(seed_row["virus_accession"]),
            "source_host": str(seed_row["host_genus"]),
            "product": str(seed_row["product"]),
            "aa_sequence": str(seed_row["aa_sequence"]),
            "selection_summary": {k: (float(v) if isinstance(v, (np.floating, float)) else int(v) if isinstance(v, (np.integer, int)) else v) for k, v in seed_selection_summary.items()},
        },
        "family_summary": family_summary,
        "family_member_ids": family_df["protein_id"].astype(str).tolist(),
        "family_centroid": family_centroid.tolist(),
        "target_reference_centroids": target_centroids,
        "target_reference_rows": target_refs,
        "mutation_window_positions_0based": [int(x) for x in mutation_window_positions],
        "mutation_window_span_1based": [int(min(mutation_window_positions) + 1), int(max(mutation_window_positions) + 1)] if len(mutation_window_positions) > 0 else [1, len(str(seed_row["aa_sequence"]))],
        "target_position_priors": target_priors,
        "hotspot_min_count": int(args.hotspot_min_count),
        "window_flank": int(args.window_flank),
    }                                                                                    

    # Write the canonical Phase A plan file.
    with open(output_dir / "phaseA_plan.json", "w", encoding="utf-8") as handle:      
        json.dump(phaseA_plan, handle, indent=2)

    # Print a summary of the Phase A plan.
    print(f"Phase A plan: {output_dir / 'phaseA_plan.json'}")
    print(f"Selected canonical seed: {seed_protein_id}")                                # Print the chosen seed so you can sanity-check the selection.
    print(f"Family size: {len(family_df)}")                                             # Print how many strict-bank proteins define the local family manifold.
    print(f"Mutation window (1-based): {phaseA_plan['mutation_window_span_1based']}")   # Print the mutation span that the constrained optimizer will use.
    print(f"Saved: {output_dir / 'phaseA_plan.json'}")                                  # Confirm where the reusable plan was written.


if __name__ == "__main__":
    main()
