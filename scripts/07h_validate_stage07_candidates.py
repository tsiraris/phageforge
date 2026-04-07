#!/usr/bin/env python3
"""
PhageForge Stage 07 Candidate Validation and Shortlisting Script
================================================================

What this script is used for:
This script acts as a rigorous post-generation filtering gateway for AI-designed
phage Receptor-Binding Proteins (RBPs). After Stage 07 generates and scores
mutated candidate sequences, this script validates their biological sanity, 
checks for sequence degeneracy, and enforces sequence diversity so that the 
wet lab receives a high-quality, varied shortlist of structural designs.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")                                              # Define a set of 20 standard valid amino acids


# Sequence metrics.
def normalized_shannon_entropy(sequence: str) -> float:
    """Return Shannon entropy normalized to [0, 1] for amino-acid composition."""
    if not sequence:                                                                # Check if the input sequence string is empty
        return 0.0                                                                  # Return an entropy of 0.0 for an empty sequence
    counts = pd.Series(list(sequence)).value_counts(normalize=True)                 # Compute the relative frequencies of each amino acid
    entropy = float(-(counts * np.log2(counts)).sum())                              # Calculate Shannon entropy using the formula: -sum(p * log2(p))
    return entropy / math.log2(len(VALID_AA))                                       # Normalize entropy to [0,1] dividing by max possible entropy


def max_residue_fraction(sequence: str) -> float:
    """Return the largest single-residue frequency fraction in the sequence."""
    if not sequence:                                                                # Check if the input sequence string is empty
        return 1.0                                                                  # Return 1.0 (max fraction) for an empty sequence safely
    counts = pd.Series(list(sequence)).value_counts(normalize=True)                 # Compute the relative frequencies of each amino acid
    return float(counts.max())                                                      # Return the frequency of the most abundant amino acid


def longest_homopolymer_run(sequence: str) -> int:
    """Return the longest run of the same residue."""
    if not sequence:                                                                # Check if the input sequence string is empty
        return 0                                                                    # Return 0 run length for an empty sequence
    best = run = 1                                                                  # Initialize the best and current run length variables to 1
    prev = sequence[0]                                                              # Set the previous character tracker to the first amino acid
    for aa in sequence[1:]:                                                         # Iterate through the rest of the sequence
        if aa == prev:                                                              # Check if the current amino acid matches the previous one
            run += 1                                                                # Increment the current contiguous run length counter
            best = max(best, run)                                                   # Update the best run length if the current run is longer
        else:                                                                       # If the current amino acid is different from the previous
            prev = aa                                                               # Update the previous amino acid tracker to the current one
            run = 1                                                                 # Reset the current contiguous run length counter back to 1
    return best                                                                     # Return the overall longest contiguous run length found


def low_complexity_fraction(sequence: str, k: int = 12, unique_threshold: int = 3) -> float:
    """Return fraction of windows with few unique residues."""
    if len(sequence) < k:                                                           # Check if sequence is shorter than the window size 'k'
        return float(len(set(sequence)) <= unique_threshold)                        # Return 1.0 if whole seq has few unique AAs, else 0.0
    flags = []                                                                      # Initialize a list to store boolean flags for each window
    for i in range(len(sequence) - k + 1):                                          # Iterate over all possible sliding windows of size 'k'
        flags.append(len(set(sequence[i : i + k])) <= unique_threshold)             # Append True if window has few unique AAs, else False
    return float(np.mean(flags)) if flags else 0.0                                  # Return the overall proportion of low-complexity windows


# Mutation parsing and geometry.
def parse_mutation_positions(mutation_positions: str) -> list[int]:
    """Extract integer positions from Stage 07 mutation annotation string."""
    if pd.isna(mutation_positions) or not str(mutation_positions).strip():          # Check if the input is NaN or an empty/whitespace string
        return []                                                                   # Return an empty list if there are no mutations
    positions: list[int] = []                                                       # Initialize an empty list to store parsed integer positions
    for token in str(mutation_positions).split(";"):                                # Split the annotation string by semicolon and iterate
        token = token.strip()                                                       # Remove leading and trailing whitespace from the token
        if not token or ":" not in token:                                           # Check if token is empty or lacks the expected ':' separator
            continue                                                                # Skip invalid or empty tokens
        pos_text = token.split(":", 1)[0].strip()                                   # Extract the position part before the colon
        try:                                                                        # Start a try block to handle potential integer parsing errors
            positions.append(int(pos_text))                                         # Convert text to int and append to the positions list
        except ValueError:                                                          # Catch the error if the extracted text is not a valid integer
            continue                                                                # Skip this token if integer conversion fails
    return positions                                                                # Return the final list of parsed integer positions


def parse_editable_positions(editable_hotspots: str) -> set[int]:
    """Extract editable hotspot positions from comma-separated annotation."""
    if pd.isna(editable_hotspots) or not str(editable_hotspots).strip():            # Check if input is NaN or an empty/whitespace string
        return set()                                                                # Return an empty set if there are no editable hotspots
    out: set[int] = set()                                                           # Initialize an empty set to store unique hotspot positions
    for token in str(editable_hotspots).split(","):                                 # Split the hotspot string by comma and iterate over tokens
        token = token.strip()                                                       # Remove leading and trailing whitespace from the token
        if not token:                                                               # Check if the token is empty after stripping whitespace
            continue                                                                # Skip empty tokens
        try:                                                                        # Start try block to safely parse the token as an integer
            out.add(int(token))                                                     # Convert token to integer and add it to the unique set
        except ValueError:                                                          # Catch the error if the token is not a valid integer
            continue                                                                # Skip this token if integer conversion fails
    return out                                                                      # Return the populated set of editable integer positions


def mutation_span(positions: Iterable[int]) -> int:
    """Return span covered by mutations."""
    pos = sorted(set(int(p) for p in positions))                                    # Deduplicate, convert to integers, and sort the positions
    if len(pos) < 2:                                                                # Check if there are fewer than 2 unique mutations
        return 0                                                                    # Return a span of 0 if a distance span cannot be formed
    return pos[-1] - pos[0]                                                         # Calculate and return difference between last and first pos


def mutation_clustering_score(positions: Iterable[int]) -> float:
    """Return a [0, 1] score where 1 is highly clustered and 0 is well spread."""
    pos = sorted(set(int(p) for p in positions))                                    # Deduplicate, convert to integers, and sort the positions
    if len(pos) < 3:                                                                # Check if there are fewer than 3 unique mutations
        return 0.0                                                                  # Return 0.0 as standard deviation isn't meaningful here
    gaps = np.diff(pos)                                                             # Calculate the adjacent differences (gaps) between positions
    if len(gaps) == 0:                                                              # Fallback check if there are no gaps
        return 1.0                                                                  # Return max clustering score if no gaps exist
    mean_gap = float(np.mean(gaps))                                                 # Calculate the mean of the gap sizes
    std_gap = float(np.std(gaps))                                                   # Calculate the standard deviation of the gap sizes
    if mean_gap == 0:                                                               # Check if the mean gap is zero
        return 1.0                                                                  # Return max clustering score
    return float(min(1.0, std_gap / mean_gap))                                      # Return coefficient of variation capped at 1.0 as the score


# Pairwise diversity metrics.
def sequence_identity(seq_a: str, seq_b: str) -> float:
    """Return position-wise identity fraction for same-length sequences."""
    if not seq_a or not seq_b or len(seq_a) != len(seq_b):                          # Check if either sequence is empty or lengths differ
        return 0.0                                                                  # Return 0.0 identity for invalid/mismatched pairs
    matches = sum(a == b for a, b in zip(seq_a, seq_b))                             # Count the number of identical characters at each position
    return matches / len(seq_a)                                                     # Calculate and return the fraction of matching positions


def mutation_overlap_fraction(positions_a: Iterable[int], positions_b: Iterable[int]) -> float:
    """Return overlap relative to the smaller mutation set."""
    a = set(int(x) for x in positions_a)                                            # Convert the first positions iterable to a set of integers
    b = set(int(x) for x in positions_b)                                            # Convert the second positions iterable to a set of integers
    if not a or not b:                                                              # Check if either of the generated sets is empty
        return 0.0                                                                  # Return 0.0 overlap if either set has no mutations
    return len(a & b) / min(len(a), len(b))                                         # Return intersection size divided by size of smaller set


@dataclass
class ValidationThresholds:
    min_mutations: int = 8                                                          # Minimum allowed number of mutations per sequence
    max_mutations: int = 40                                                         # Maximum allowed number of mutations per sequence
    min_normalized_entropy: float = 0.45                                            # Minimum allowed sequence Shannon entropy score
    max_single_residue_fraction: float = 0.22                                       # Maximum allowed frequency for any single amino acid
    max_homopolymer_run: int = 6                                                    # Maximum allowed length of identical repeating amino acids
    max_low_complexity_fraction: float = 0.20                                       # Maximum allowed fraction of low-complexity sequence windows
    min_mutation_span: int = 8                                                      # Minimum allowed distance between first and last mutation
    max_outside_editable_fraction: float = 0.0                                      # Maximum allowed fraction of mutations outside editable zones
    max_pairwise_identity: float = 0.985                                            # Maximum allowed sequence identity between any two candidates
    max_mutation_overlap: float = 0.80                                              # Maximum allowed mutation overlap between any two candidates


@dataclass
class CandidateDecision:
    sample_id: int | str                                                            # The unique identifier for the candidate sample
    keep: bool                                                                      # Boolean indicating whether the candidate passed validation
    stage: str                                                                      # The specific validation stage (e.g., hard_filters, diversity)
    reason: str                                                                     # Text description of the failure reason or "pass"


# Candidate-level validation.
def validate_candidate(row: pd.Series, seed_length: int, thresholds: ValidationThresholds) -> tuple[bool, list[str], dict[str, float | int | bool]]:
    """Apply hard candidate-level filters and return metrics."""
    seq = str(row.get("candidate_sequence", "") or "")                              # Extract candidate sequence from the row, defaulting to empty
    mutations = parse_mutation_positions(row.get("mutation_positions"))             # Parse the mutation positions string into a list of integers
    editable = parse_editable_positions(row.get("editable_hotspots"))               # Parse the editable hotspots string into a set of integers

    entropy = normalized_shannon_entropy(seq)                                       # Calculate the normalized Shannon entropy for the sequence
    max_frac = max_residue_fraction(seq)                                            # Calculate the highest single-residue fraction in the sequence
    longest_run = longest_homopolymer_run(seq)                                      # Find the length of the longest contiguous homopolymer run
    lc_frac = low_complexity_fraction(seq)                                          # Calculate the fraction of the sequence that has low complexity
    span = mutation_span(mutations)                                                 # Calculate the span distance between first and last mutation
    clustered = mutation_clustering_score(mutations)                                # Calculate the clustering score of the mutations
    outside_editable = [p for p in mutations if editable and p not in editable]     # Find mutations that fall outside the defined editable regions
    outside_frac = len(outside_editable) / len(mutations) if mutations else 0.0     # Calculate the fraction of mutations that are outside bounds

    failures: list[str] = []                                                        # Initialize an empty list to store failure reasons
    if not seq:                                                                     # Check if the sequence string is empty
        failures.append("empty_sequence")                                           # Append failure reason for empty sequence
    if any(aa not in VALID_AA for aa in seq):                                       # Check if the sequence contains any non-standard amino acids
        failures.append("invalid_amino_acids")                                      # Append failure reason for invalid amino acids
    if len(seq) != seed_length:                                                     # Check if the sequence length differs from the seed length
        failures.append("length_mismatch")                                          # Append failure reason for length mismatch
    if len(mutations) < thresholds.min_mutations:                                   # Check if the mutation count is below the minimum threshold
        failures.append("too_few_mutations")                                        # Append failure reason for insufficient mutations
    if len(mutations) > thresholds.max_mutations:                                   # Check if the mutation count exceeds the maximum threshold
        failures.append("too_many_mutations")                                       # Append failure reason for excessive mutations
    if entropy < thresholds.min_normalized_entropy:                                 # Check if the sequence entropy is below the minimum threshold
        failures.append("low_entropy")                                              # Append failure reason for low entropy
    if max_frac > thresholds.max_single_residue_fraction:                           # Check if maximum single-residue fraction exceeds limits
        failures.append("single_residue_enrichment")                                # Append failure reason for single residue over-enrichment
    if longest_run > thresholds.max_homopolymer_run:                                # Check if the longest homopolymer run exceeds the limit
        failures.append("long_homopolymer")                                         # Append failure reason for overly long homopolymer runs
    if lc_frac > thresholds.max_low_complexity_fraction:                            # Check if the low complexity fraction exceeds the limit
        failures.append("low_complexity")                                           # Append failure reason for too much low complexity regions
    if span < thresholds.min_mutation_span:                                         # Check if the overall mutation span is smaller than the limit
        failures.append("mutation_span_too_small")                                  # Append failure reason for an overly constrained mutation span
    if outside_frac > thresholds.max_outside_editable_fraction:                     # Check if the fraction of outside mutations exceeds the limit
        failures.append("outside_editable_region")                                  # Append failure reason for mutating non-editable regions

    metrics = {                                                                     # Initialize a dictionary to bundle all calculated metrics
        "sequence_length": len(seq),                                                # Store the length of the candidate sequence
        "mutation_count": len(mutations),                                           # Store the total number of valid mutations parsed
        "normalized_entropy": entropy,                                              # Store the calculated normalized Shannon entropy
        "max_single_residue_fraction": max_frac,                                    # Store the calculated maximum single-residue fraction
        "longest_homopolymer_run": longest_run,                                     # Store the calculated longest homopolymer run length
        "low_complexity_fraction": lc_frac,                                         # Store the calculated fraction of low-complexity regions
        "mutation_span": span,                                                      # Store the calculated mutation span
        "mutation_clustering_score": clustered,                                     # Store the calculated mutation clustering score
        "outside_editable_fraction": outside_frac,                                  # Store the calculated fraction of mutations outside hotpots
        "hard_pass": not failures,                                                  # Store boolean indicating if no hard filter failures occurred
    }                                                                               # Close the dictionary definition
    return not failures, failures, metrics                                          # Return the boolean pass status, failure list, and metrics dict


# Export helpers.
def write_fasta(df: pd.DataFrame, fasta_path: Path, score_col: str) -> None:
    """Write candidates to FASTA."""
    lines: list[str] = []                                                           # Initialize an empty list to accumulate lines for the FASTA file
    for _, row in df.iterrows():                                                    # Iterate over the rows of the provided pandas DataFrame
        score = float(row.get(score_col, float("nan")))                             # Extract the score from the specified column, defaulting to NaN
        lines.append(                                                               # Append the constructed FASTA header line to the list
            f">stage07_sample{row['sample_id']}|diverse_rank={row.get('validated_rank', '')}|score={score:.6f}|regime={row.get('generation_regime', 'na')}"
        )                                                                           # Close the string interpolation for the header
        lines.append(str(row["candidate_sequence"]))                                # Append the actual candidate sequence as the next line
    fasta_path.write_text("\n".join(lines) + ("\n" if lines else ""))               # Join lines with newlines and write them to the target file path


def build_parser() -> argparse.ArgumentParser:
    """Create CLI parser."""
    parser = argparse.ArgumentParser(description=__doc__)                           # Initialize the argument parser using the module's docstring
    parser.add_argument("--ranked_csv", required=True, help="Stage 07 final ranked candidate CSV.") # Add argument for input CSV path
    parser.add_argument("--out_dir", required=True, help="Output directory for validated exports.") # Add argument for output directory path
    parser.add_argument("--top_k", type=int, default=5, help="How many high-ranked rows to inspect initially.") # Add arg for initial sample cap
    parser.add_argument("--keep_top_k", type=int, default=5, help="How many validated rows to export as top set.") # Add arg for final export cap
    parser.add_argument("--keep_top3_k", type=int, default=3, help="How many validated rows to export as primary set.") # Add arg for top-3 export
    parser.add_argument("--min_mutations", type=int, default=8)                     # Add argument configuring minimum allowed mutations threshold
    parser.add_argument("--max_mutations", type=int, default=40)                    # Add argument configuring maximum allowed mutations threshold
    parser.add_argument("--min_normalized_entropy", type=float, default=0.45)       # Add argument configuring minimum normalized entropy threshold
    parser.add_argument("--max_single_residue_fraction", type=float, default=0.22)  # Add argument configuring max single residue fraction
    parser.add_argument("--max_homopolymer_run", type=int, default=6)               # Add argument configuring maximum allowed homopolymer run
    parser.add_argument("--max_low_complexity_fraction", type=float, default=0.20)  # Add argument configuring max low complexity fraction
    parser.add_argument("--min_mutation_span", type=int, default=8)                 # Add argument configuring minimum allowed mutation span
    parser.add_argument("--max_outside_editable_fraction", type=float, default=0.0) # Add argument configuring max mutations outside hotpots
    parser.add_argument("--max_pairwise_identity", type=float, default=0.985)       # Add argument configuring maximum pairwise sequence identity
    parser.add_argument("--max_mutation_overlap", type=float, default=0.80)         # Add argument configuring maximum pairwise mutation overlap
    return parser                                                                   # Return the fully configured argument parser object


def main() -> None:
    """Run validation workflow."""
    args = build_parser().parse_args()                                              # Parse command line arguments using the constructed parser
    out_dir = Path(args.out_dir)                                                    # Convert the output directory string path to a Path object
    out_dir.mkdir(parents=True, exist_ok=True)                                      # Create the output directory and missing parent dirs safely

    thresholds = ValidationThresholds(                                              # Instantiate the thresholds dataclass using parsed arguments
        min_mutations=args.min_mutations,                                           # Map parsed minimum mutations argument to the dataclass field
        max_mutations=args.max_mutations,                                           # Map parsed maximum mutations argument to the dataclass field
        min_normalized_entropy=args.min_normalized_entropy,                         # Map parsed minimum entropy argument to the dataclass field
        max_single_residue_fraction=args.max_single_residue_fraction,               # Map parsed max single residue fraction to the dataclass field
        max_homopolymer_run=args.max_homopolymer_run,                               # Map parsed max homopolymer run to the dataclass field
        max_low_complexity_fraction=args.max_low_complexity_fraction,               # Map parsed max low complexity fraction to the dataclass field
        min_mutation_span=args.min_mutation_span,                                   # Map parsed minimum mutation span to the dataclass field
        max_outside_editable_fraction=args.max_outside_editable_fraction,           # Map parsed max outside editable fraction to the dataclass field
        max_pairwise_identity=args.max_pairwise_identity,                           # Map parsed max pairwise identity to the dataclass field
        max_mutation_overlap=args.max_mutation_overlap,                             # Map parsed max mutation overlap to the dataclass field
    )                                                                               # Close the instantiation of ValidationThresholds

    ranked = pd.read_csv(args.ranked_csv).copy()                                    # Read the input CSV into a dataframe and make an explicit copy
    ranked = ranked.sort_values(["rank_diverse", "rank_raw"], kind="mergesort").reset_index(drop=True) # Sort stably by rankings and reset index
    inspect_df = ranked.head(args.top_k).copy()                                     # Take the top 'k' candidates for inspection and make a copy

    if inspect_df.empty:                                                            # Check if the extracted inspection dataframe is empty
        raise ValueError("No candidates found in ranked CSV.")                      # Raise an error if no candidate rows were found

    seed_length = len(str(inspect_df.iloc[0]["seed_sequence"]))                     # Determine reference length from first candidate's seed sequence
    decisions: list[CandidateDecision] = []                                         # Initialize a list to track keep/drop decisions for candidates
    metrics_rows: list[dict[str, object]] = []                                      # Initialize a list of dictionaries to store metrics for each row

    # Candidate-level hard filters.
    hard_pass_rows: list[pd.Series] = []                                            # Initialize list to hold row data for candidates passing hard filters
    for _, row in inspect_df.iterrows():                                            # Iterate over each candidate row in the inspection dataframe
        keep, failures, metrics = validate_candidate(row, seed_length, thresholds)  # Validate the current candidate against hard filters
        metrics_rows.append({**row.to_dict(), **metrics, "failure_reasons": ";".join(failures)}) # Append original data, metrics, and failures
        if keep:                                                                    # Check if the candidate passed all hard filter checks
            hard_pass_rows.append(row)                                              # Append the passing row to the hard_pass_rows list
        decisions.append(                                                           # Record the decision outcome for the current candidate
            CandidateDecision(                                                      # Instantiate a CandidateDecision object
                sample_id=row["sample_id"],                                         # Record the sample ID for the current candidate
                keep=keep,                                                          # Record whether the candidate was kept or dropped
                stage="hard_filters",                                               # Mark the current validation stage as "hard_filters"
                reason="pass" if keep else ";".join(failures),                      # Record "pass" or the semicolon-joined list of failure reasons
            )                                                                       # Close the instantiation of CandidateDecision
        )                                                                           # Close the append operation

    hard_pass_df = pd.DataFrame([r.to_dict() for r in hard_pass_rows])              # Convert the list of passing rows back into a Pandas DataFrame

    # Diversity pruning over surviving candidates.
    kept_rows: list[pd.Series] = []                                                 # Initialize list to track rows that pass the diversity filter
    diversity_notes: list[CandidateDecision] = []                                   # Initialize list to track diversity rejection/acceptance notes
    for _, row in hard_pass_df.sort_values(["rank_diverse", "rank_raw"], kind="mergesort").iterrows(): # Iterate over hard-passed candidates stably
        row_mut = parse_mutation_positions(row.get("mutation_positions"))           # Extract mutation positions for the current candidate
        row_seq = str(row["candidate_sequence"])                                    # Extract sequence string for the current candidate
        too_similar = False                                                         # Flag indicating if candidate is too similar to an already kept one
        reasons: list[str] = []                                                     # Initialize list to store specific similarity failure reasons
        for kept in kept_rows:                                                      # Iterate over candidates that have already been kept
            identity = sequence_identity(row_seq, str(kept["candidate_sequence"]))  # Calculate pairwise identity with the kept candidate
            overlap = mutation_overlap_fraction(row_mut, parse_mutation_positions(kept.get("mutation_positions"))) # Calculate mutation overlap
            if identity > thresholds.max_pairwise_identity:                         # Check if identity exceeds the maximum allowed threshold
                too_similar = True                                                  # Set the flag to true since it's too similar
                reasons.append(f"pairwise_identity>{thresholds.max_pairwise_identity:.3f}:sample{kept['sample_id']}") # Append exact identity reason
            if overlap > thresholds.max_mutation_overlap:                           # Check if mutation overlap exceeds the maximum allowed threshold
                too_similar = True                                                  # Set the flag to true since it's too similar
                reasons.append(f"mutation_overlap>{thresholds.max_mutation_overlap:.2f}:sample{kept['sample_id']}") # Append exact overlap reason
        if too_similar:                                                             # Check if the current candidate was flagged as too similar
            diversity_notes.append(                                                 # Record the rejection in the diversity notes list
                CandidateDecision(                                                  # Instantiate a CandidateDecision object for rejection
                    sample_id=row["sample_id"],                                     # Record the sample ID
                    keep=False,                                                     # Mark the candidate as dropped (not kept)
                    stage="diversity",                                              # Mark the current validation stage as "diversity"
                    reason=";".join(reasons),                                       # Record the semicolon-joined list of similarity failure reasons
                )                                                                   # Close the CandidateDecision instantiation
            )                                                                       # Close the append operation
            continue                                                                # Skip to the next candidate without adding this one to kept_rows
        kept_rows.append(row)                                                       # Add the unique candidate to the kept_rows list
        diversity_notes.append(                                                     # Record the success in the diversity notes list
            CandidateDecision(sample_id=row["sample_id"], keep=True, stage="diversity", reason="pass") # Instantiate and record a "pass" decision
        )                                                                           # Close the append operation

    validated = pd.DataFrame([r.to_dict() for r in kept_rows]).copy()               # Convert the uniquely kept rows into a Pandas DataFrame and copy it
    validated = validated.head(args.keep_top_k).reset_index(drop=True)              # Truncate to maximum requested number of top rows and reset index
    if not validated.empty:                                                         # Check if the final validated dataframe is not empty
        validated["validated_rank"] = np.arange(1, len(validated) + 1)              # Assign an integer ranking to each candidate from 1 to N

    validated_top3 = validated.head(args.keep_top3_k).copy()                        # Create a smaller subset dataframe containing just the top 3 items

    # Join validation metrics.
    metrics_df = pd.DataFrame(metrics_rows)                                         # Convert the raw accumulated metrics dictionaries into a DataFrame
    summary_df = metrics_df.merge(                                                  # Merge the metrics DataFrame with the decision records
        pd.DataFrame(asdict(x) for x in decisions + diversity_notes),               # Create a DataFrame from combined hard-filter and diversity decisions
        how="left",                                                                 # Perform a left join to keep all records from the metrics DataFrame
        on="sample_id",                                                             # Join the dataframes using the shared 'sample_id' column
        suffixes=("", "_decision"),                                                 # Append '_decision' to overlapping column names from right DataFrame
    )                                                                               # Close the merge operation

    # Write outputs.
    summary_csv = out_dir / "validation_summary.csv"                                # Construct the output filepath for the validation summary CSV
    validated_csv = out_dir / "validated_top5.csv"                                  # Construct the output filepath for the main validated candidates CSV
    validated_fasta = out_dir / "validated_top5.fasta"                              # Construct the output filepath for the main validated FASTA
    validated_top3_csv = out_dir / "validated_top3.csv"                             # Construct the output filepath for the top-3 validated CSV
    validated_top3_fasta = out_dir / "validated_top3.fasta"                         # Construct the output filepath for the top-3 validated FASTA
    thresholds_json = out_dir / "validation_thresholds.json"                        # Construct the output filepath for the thresholds JSON config

    summary_df.to_csv(summary_csv, index=False)                                     # Write the merged summary dataframe to a CSV file without row indices
    validated.to_csv(validated_csv, index=False)                                    # Write the top 5 validated dataframe to a CSV file without row indices
    write_fasta(validated, validated_fasta, score_col="final_multimodal_rank_score")# Delegate writing the top 5 candidates to a FASTA file
    validated_top3.to_csv(validated_top3_csv, index=False)                          # Write the top 3 validated dataframe to a CSV file without row indices
    write_fasta(validated_top3, validated_top3_fasta, score_col="final_multimodal_rank_score") # Delegate writing the top 3 candidates to a FASTA
    thresholds_json.write_text(json.dumps(asdict(thresholds), indent=2))            # Convert the thresholds dataclass to a formatted JSON string

    print(f"Wrote: {summary_csv}")                                                  # Print console confirmation for the summary CSV
    print(f"Wrote: {validated_csv}")                                                # Print console confirmation for the top 5 CSV
    print(f"Wrote: {validated_fasta}")                                              # Print console confirmation for the top 5 FASTA
    print(f"Wrote: {validated_top3_csv}")                                           # Print console confirmation for the top 3 CSV
    print(f"Wrote: {validated_top3_fasta}")                                         # Print console confirmation for the top 3 FASTA
    print(f"Wrote: {thresholds_json}")                                              # Print console confirmation for the thresholds JSON

    print("\nValidated top candidates:")                                            # Print a section header for the console summary output
    if validated.empty:                                                             # Check if no candidates survived the filtering process
        print("No candidates survived validation.")                                 # Print a warning message indicating total filtration
    else:                                                                           # If candidates did survive the validation process
        cols = [                                                                    # Initialize a list of column names to display in the console summary
            "validated_rank",                                                       # Add rank column to display list
            "sample_id",                                                            # Add sample ID column to display list
            "generation_regime",                                                    # Add generation regime column to display list
            "final_multimodal_rank_score",                                          # Add final score column to display list
            "target_score",                                                         # Add target score column to display list
            "strict_manifold_score",                                                # Add manifold score column to display list
            "structure_score",                                                      # Add structure score column to display list
            "mutation_positions",                                                   # Add mutation positions column to display list
        ]                                                                           # Close the column list definition
        existing = [c for c in cols if c in validated.columns]                      # Filter the display list to only include columns that actually exist
        print(validated[existing].to_string(index=False))                           # Print the selected columns of the validated dataframe as text


if __name__ == "__main__":
    main()                                                                          # Execute the main function, starting the script's workflow