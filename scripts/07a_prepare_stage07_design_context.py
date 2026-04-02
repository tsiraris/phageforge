"""
===============================================================================
07a: Prepare one self-contained Stage 07 context JSON from Stage 05/06 outputs.
===============================================================================

This script is the bridge from the existing validated pipeline to Stage 07.
It reads the exact outputs already produced by Stage 06 and converts them into
one reusable JSON context file that the later Stage 07 scripts can consume.

Supported upstream inputs:
- `phaseA_plan.json` from `06a_select_phaseA_family.py`
- either `top_candidates.csv` from `06b_optimize_family_constrained.py`
  or `phaseA_followup_seed_summary.json` from `06c_pick_phaseA_followup_seed.py`
- the strict RBP CSV for recovering the actual seed sequence by protein ID
- optional tissue metadata CSV
"""

from __future__ import annotations                                                    # Enable postponed evaluation of type annotations for cleaner typing.
import argparse                                                                       # Parse command-line arguments passed to the script.
import json                                                                           # Read upstream JSON files and write the final Stage 07 context JSON.
from pathlib import Path                                                              # Build filesystem paths safely across operating systems.
from typing import Any                                                                # Express helper return types that may vary across columns.
import pandas as pd                                                                   # Load CSV tables and work with them in a structured way.


def parse_args() -> argparse.Namespace:
    """Return command-line arguments for building the Stage 07 context JSON.

    The parser is intentionally flexible because the current repo can provide
    either a Phase 06b `top_candidates.csv` file or a Phase 06c follow-up JSON.
    That keeps Stage 07 compatible with the outputs you already have.
    """
    ap = argparse.ArgumentParser(                                                     # Create the main command-line parser for this script.
        description="Prepare Stage 07 generation context from Stage 05/06 outputs."   # Describe the script in `--help` output.
    )
    ap.add_argument(                                                                  # Add the path to the reusable Phase A plan JSON.
        "--phaseA_plan_json", type=str, required=True,
        help="Plan JSON from scripts/06a_select_phaseA_family.py.",
    )
    ap.add_argument(                                                                  # Add an optional direct path to the 06b ranked candidate table.
        "--phase06b_top_candidates_csv", type=str, default="",
        help="Optional top_candidates.csv from scripts/06b_optimize_family_constrained.py.",
    )
    ap.add_argument(                                                                  # Add an optional direct path to the 06c follow-up summary JSON.
        "--phase06c_followup_summary_json", type=str, default="",
        help="Optional phaseA_followup_seed_summary.json from scripts/06c_pick_phaseA_followup_seed.py.",
    )
    ap.add_argument(                                                                  # Add the strict sequence bank used to recover the actual amino-acid sequence.
        "--strict_csv", type=str, required=True,
        help="Strict RBP dataset CSV for recovering seed sequence when needed.",
    )
    ap.add_argument(                                                                  # Add the desired Stage 07 target host label.
        "--target_host", type=str, required=True,
        help="Target host for Stage 07, e.g. Acinetobacter.",
    )
    ap.add_argument(                                                                  # Allow the user to select a different seed rank from the 06b shortlist.
        "--seed_rank", type=int, default=0,
        help="Which ranked follow-up seed to use when reading a 06b top_candidates.csv file.",
    )
    ap.add_argument(                                                                  # Allow the editable region to be widened or narrowed manually.
        "--window_flank_override", type=int, default=-1,
        help="Optional override for hotspot flank size. Negative means use plan value.",
    )
    ap.add_argument(                                                                  # Accept optional pathology / tissue metadata for the multimodal branch.
        "--tissue_metadata_csv", type=str, default="",
        help="Optional CSV describing tissue / pathology / omics context.",
    )
    ap.add_argument(                                                                  # Accept the final JSON destination for Stage 07 downstream scripts.
        "--output_json", type=str, required=True,
        help="Where to write the Stage 07 context JSON.",
    )
    return ap.parse_args()                                                            # Parse the CLI arguments and return them as a namespace.


# ------------------------------ small normalization helpers ------------------------------ #
# The next helpers smooth over minor upstream naming differences so Stage 07
# can consume existing Stage 05/06 outputs without forcing you to rename files
# or rewrite older scripts.
def first_present(row: pd.Series, candidates: list[str], default: Any = None) -> Any:
    """Return the first non-null value from a candidate column list.

    Older and newer scripts sometimes use slightly different column names.
    This helper lets Stage 07 stay robust without forcing upstream rewrites.
    """
    for col in candidates:                                                            # Walk through the candidate column names in priority order.
        if col in row.index and pd.notna(row[col]):                                   # Keep the first column that exists and is not missing.
            return row[col]                                                           # Return the selected value immediately.
    return default                                                                    # Fall back to the provided default if nothing was found.


def nonempty_str(value: Any, default: str = "") -> str:
    """Return a stripped string unless it is missing / empty / None-like, then use default."""
    if value is None:                                                                 # Handle real Python None explicitly.
        return default                                                                # Return the requested default when no value exists.
    text = str(value).strip()                                                         # Normalize arbitrary input into a stripped string.
    if text == "" or text.lower() == "none" or text.lower() == "nan":               # Treat common missing-string spellings as missing values.
        return default                                                                # Fall back when the string is effectively empty.
    return text                                                                       # Keep the normalized non-empty string.


def load_json(path: Path) -> dict:
    """Load a JSON file and return it as a Python dictionary."""
    with open(path, "r", encoding="utf-8") as handle:                                 # Open the JSON file using UTF-8 for safe text handling.
        return json.load(handle)                                                      # Parse the JSON and return the resulting dictionary.


# ------------------------------ seed-sequence recovery helpers ------------------------------ #
# Stage 07 can start from three different kinds of upstream seed references:
# 1. a FASTA exported by Stage 06c
# 2. a generated candidate row stored in Stage 06b top_candidates.csv
# 3. a real protein_id that still exists in the strict RBP bank
#
# The helpers below resolve those cases explicitly instead of assuming that every
# Stage 06 seed identifier is always a strict-bank protein_id.
def load_fasta_sequence(path: Path) -> str:
    """Load a single FASTA sequence and return it as one amino-acid string."""
    lines = path.read_text(encoding="utf-8").splitlines()                             # Read the full FASTA file into memory line by line.
    seq_lines = [line.strip() for line in lines if line.strip() and not line.startswith(">") ]  # Keep only non-header sequence lines.
    if not seq_lines:                                                                 # Fail clearly when the FASTA has no sequence payload.
        raise ValueError(f"No sequence found in FASTA: {path}")
    return "".join(seq_lines)                                                         # Join multiline FASTA content into one contiguous sequence string.


def load_seed_sequence(strict_csv: Path, protein_id: str) -> tuple[str, str, str]:
    """Recover virus accession, source host, and amino-acid sequence for a protein ID.

    Stage 06 summary files sometimes keep only an upstream strict-bank protein ID.
    This helper looks the sequence back up in the strict bank so that generation
    can still start from a real validated seed.
    """
    df = pd.read_csv(strict_csv)                                                      # Read the strict RBP bank from disk.
    if "protein_id" not in df.columns:                                                # Validate that the lookup key exists.
        raise ValueError("Strict CSV must contain a protein_id column.")              # Fail early with a clear message if not.

    sub = df[df["protein_id"].astype(str) == str(protein_id)].copy()                  # Keep only rows that match the requested protein ID.
    if sub.empty:                                                                     # Handle the case where the seed is absent from the strict bank.
        raise ValueError(                                                             # Raise a useful error that includes both the ID and the file path.
            f"Could not find seed protein_id={protein_id} in strict_csv={strict_csv}."
        )

    row = sub.iloc[0]                                                                 # Take the first matching row as the seed record.
    seq = first_present(row, ["aa_sequence", "sequence"])                             # Accept either the newer or older sequence column name.
    if seq is None or str(seq) == "":                                                 # Ensure the recovered seed sequence is non-empty.
        raise ValueError(f"No sequence found for protein_id={protein_id}.")           # Stop if no valid sequence was recovered.

    return (                                                                          # Return the core seed metadata used downstream.
        str(row.get("virus_accession", "UNKNOWN")),                                   # Return the source phage accession if available.
        str(row.get("host_genus", "UNKNOWN")),                                        # Return the source host genus if available.
        str(seq),                                                                     # Return the amino-acid sequence as a string.
    )


def load_seed_row_from_top_candidates(top_csv: Path, seed_rank: int) -> dict:
    """Load one ranked candidate row from a Phase 06b `top_candidates.csv` file."""
    df = pd.read_csv(top_csv)                                                         # Read the top-candidate table created by the family-constrained optimizer.
    if len(df) == 0:                                                                  # Prevent later indexing into an empty file.
        raise ValueError(f"Phase 06b top candidate CSV is empty: {top_csv}")          # Raise a clear error if the file has no candidates.
    if seed_rank < 0 or seed_rank >= len(df):                                         # Guard against invalid rank selection.
        raise ValueError(f"--seed_rank={seed_rank} is outside the valid range [0, {len(df) - 1}].")
    return df.iloc[seed_rank].to_dict()                                               # Return the chosen row as a plain dictionary.


def load_seed_row_from_followup_summary(summary_json: Path) -> dict:
    """Construct a Stage 07 seed row from the 06c follow-up summary JSON."""
    summary = load_json(summary_json)                                                 # Read the 06c follow-up summary from disk.
    summary_series = pd.Series(summary)                                               # Wrap the raw JSON in a Series so `first_present` can be reused.
    return {                                                                          # Normalize the summary fields into a row-like dictionary.
        "candidate_id": first_present(                                                # Preserve whichever candidate-id key the summary actually uses.
            summary_series,
            ["chosen_candidate_id", "candidate_id", "selected_candidate_id"],
        ),
        "selected_seed_protein_id": first_present(                                    # Carry forward an explicit selected seed protein ID if one exists.
            summary_series,
            ["selected_seed_protein_id", "seed_protein_id"],
        ),
        "seed_protein_id": first_present(                                             # Keep the canonical strict-bank protein ID when present.
            summary_series,
            ["canonical_seed_protein_id", "seed_protein_id", "protein_id"],
        ),
        "protein_id": first_present(                                                  # Mirror the same ID under another common column name.
            summary_series,
            ["canonical_seed_protein_id", "protein_id", "seed_protein_id"],
        ),
        "chosen_target_score": summary.get("chosen_target_score"),                    # Carry forward the target score for traceability.
        "chosen_selection_score": summary.get("chosen_selection_score"),              # Carry forward the selection score for traceability.
        "chosen_family_cosine": summary.get("chosen_family_cosine"),                  # Carry forward the family similarity for traceability.
        "chosen_target_anchor_cosine": summary.get("chosen_target_anchor_cosine"),    # Carry forward the anchor similarity for traceability.
        "chosen_n_mutations": summary.get("chosen_n_mutations"),                      # Carry forward mutation burden for traceability.
        "followup_fasta_path": first_present(                                         # Record any FASTA path key used by the summary.
            summary_series,
            ["followup_fasta_path", "fasta_path", "seed_fasta"],
        ),
        "next_target": summary.get("next_target"),                                    # Preserve the next ladder target from the summary file.
        "source_host": summary.get("source_host"),                                    # Preserve the originating host if the summary recorded it.
    }


def infer_seed_metadata(phaseA: dict, seed_row: dict) -> tuple[str, str | None]:
    """Infer the chosen seed identifier and a human-readable source description.

    The returned identifier is only a *hint* at this point. It may refer to:
    - a strict-bank protein_id
    - a Stage 06 candidate_id
    - or a FASTA-exported candidate
    The concrete sequence is resolved later by `resolve_seed_sequence`.
    """
    seed_identifier = str(                                                            # Prefer explicit selected seed fields, then fall back to canonical plan data.
        first_present(
            pd.Series(seed_row),
            ["selected_seed_protein_id", "seed_protein_id", "protein_id", "candidate_id"],
            default=phaseA.get("canonical_seed", {}).get("seed_protein_id"),
        )
    )
    source_desc = None                                                                # Initialize the source description used for auditability.
    if "candidate_id" in seed_row and seed_row["candidate_id"] is not None:          # If the row came from 06b/06c, report that candidate ID.
        source_desc = f"candidate_id={seed_row['candidate_id']}"                      # Build a human-readable source descriptor.
    elif "followup_fasta_path" in seed_row and seed_row["followup_fasta_path"]:       # If the row came from 06c, report the FASTA path.
        source_desc = f"followup_fasta={seed_row['followup_fasta_path']}"             # Build a human-readable source descriptor.
    return seed_identifier, source_desc                                               # Return the normalized seed identifier and the source descriptor.


def resolve_relative_to_repo(raw_path: str | Path, repo_hint: Path) -> Path:
    """Resolve a possibly relative path using the current repo layout as a hint."""
    path = Path(raw_path)                                                             # Normalize the input into a Path object.
    if path.is_absolute():                                                            # Keep absolute paths untouched.
        return path                                                                   # Return immediately when the path is already absolute.
    candidates = [                                                                    # Try the most common relative anchors in your repo / Colab layout.
        repo_hint / path,                                                             # Resolve relative to the inferred repo root first.
        Path.cwd() / path,                                                            # Also try the current working directory.
        repo_hint.parent / path,                                                      # And try one level above the repo root if needed.
    ]
    for cand in candidates:                                                           # Return the first candidate path that actually exists.
        if cand.exists():
            return cand
    return repo_hint / path                                                           # Fall back to repo-relative resolution even if it does not exist yet.


def load_seed_row_from_top_candidates_by_candidate_id(top_csv: Path, candidate_id: str) -> dict:
    """Load one candidate row from a Phase 06b table using candidate_id instead of rank."""
    df = pd.read_csv(top_csv)                                                         # Read the optimizer output table from disk.
    if "candidate_id" not in df.columns:                                              # Ensure the required lookup key exists.
        raise ValueError(f"top_candidates.csv is missing a candidate_id column: {top_csv}")
    sub = df[df["candidate_id"].astype(str) == str(candidate_id)].copy()              # Keep only the candidate row matching the requested ID.
    if sub.empty:                                                                     # Fail clearly if the chosen ID is absent from the table.
        raise ValueError(f"Could not find candidate_id={candidate_id} in top_candidates_csv={top_csv}.")
    return sub.iloc[0].to_dict()                                                      # Return the matched row as a dictionary.


def resolve_seed_sequence(
    strict_csv: Path,
    seed_row: dict,
    repo_root: Path,
    top_candidates_csv: Path | None = None,
) -> tuple[str, str, str, str]:
    """Resolve the actual Stage 07 seed sequence from FASTA, candidate table, or strict CSV.

    Returns:
        virus_accession_like_id, source_host, seed_identifier, seed_sequence

    Resolution order:
    1. FASTA exported by Stage 06c (preferred, because it is the exact follow-up seed)
    2. candidate_id lookup in Stage 06b top_candidates.csv
    3. strict-bank protein_id lookup in the strict RBP CSV
    """
    # First, prefer the explicit FASTA exported by Stage 06c because it contains the
    # exact designed follow-up seed chosen for the next host jump.
    followup_fasta = seed_row.get("followup_fasta_path")                              # Read the follow-up FASTA path if the summary provided one.
    candidate_id = nonempty_str(seed_row.get("candidate_id"), default="")            # Keep the Stage 06 candidate ID handy for fallback lookups.
    candidate_row = None                                                              # Cache the resolved candidate table row when it can be found.
    if candidate_id and top_candidates_csv is not None and top_candidates_csv.exists():
        try:                                                                          # Best-effort lookup: useful for source_host and extra provenance even when FASTA exists.
            candidate_row = load_seed_row_from_top_candidates_by_candidate_id(top_candidates_csv, candidate_id)
        except Exception:
            candidate_row = None                                                      # Keep FASTA resolution alive even if the candidate table is unavailable.

    if followup_fasta:                                                                # Prefer FASTA when available because it is the concrete Stage 06 seed artifact.
        fasta_path = resolve_relative_to_repo(str(followup_fasta), repo_root)         # Resolve relative paths robustly in local and Colab runs.
        if fasta_path.exists():                                                       # Only use the FASTA branch when the file actually exists.
            seed_sequence = load_fasta_sequence(fasta_path)                           # Read the designed seed sequence directly from FASTA.
            source_host = nonempty_str(seed_row.get("source_host"), default="")       # First try the source_host stored in the summary itself.
            if source_host == "" and candidate_row is not None:                        # If summary metadata is incomplete, recover the source host from Stage 06b.
                source_host = nonempty_str(
                    first_present(pd.Series(candidate_row), ["source_host", "seed_host", "host_genus"], default=""),
                    default="",
                )
            if source_host == "":                                                    # If all else fails, keep a generic placeholder instead of a misleading missing string.
                source_host = "UNKNOWN"
            if candidate_row is not None:                                             # Recover a stable row identifier for provenance when the candidate table exists.
                virus_accession_like_id = nonempty_str(
                    first_present(pd.Series(candidate_row), ["candidate_id", "seed_protein_id", "protein_id"], default=candidate_id),
                    default=candidate_id or "stage06_followup_seed",
                )
            else:
                virus_accession_like_id = candidate_id or "stage06_followup_seed"     # Keep a candidate-like identifier even without the Stage 06b row.
            return virus_accession_like_id, source_host or "UNKNOWN", virus_accession_like_id, seed_sequence

    # Second, if the summary refers to a candidate_id, recover its sequence from the
    # Stage 06b optimizer output table rather than pretending it is a strict-bank protein.
    if candidate_id and candidate_row is not None:
        seed_sequence = str(first_present(pd.Series(candidate_row), ["aa_sequence", "candidate_sequence", "sequence"], default=""))  # Accept several possible sequence columns.
        if seed_sequence == "":                                                       # Stop if the candidate row did not contain a recoverable sequence.
            raise ValueError(
                f"Resolved candidate_id={candidate_id} in {top_candidates_csv} but found no sequence column."
            )
        virus_accession_like_id = nonempty_str(                                       # Prefer candidate-level identifiers for traceability.
            first_present(pd.Series(candidate_row), ["candidate_id", "seed_protein_id", "protein_id"], default=candidate_id),
            default=candidate_id,
        )
        source_host = nonempty_str(                                                   # Keep the optimizer's source-host annotation when present.
            first_present(pd.Series(candidate_row), ["source_host", "seed_host", "host_genus"], default=seed_row.get("source_host")),
            default="UNKNOWN",
        )
        return virus_accession_like_id, source_host, candidate_id, seed_sequence      # Return Stage 06 candidate-derived seed details.

    # Finally, fall back to the original strict-bank protein lookup when the seed is
    # genuinely an upstream strict RBP rather than a generated Stage 06 candidate.
    strict_protein_id = first_present(                                                # Recover whichever key points to the strict-bank protein ID.
        pd.Series(seed_row),
        ["selected_seed_protein_id", "seed_protein_id", "protein_id"],
    )
    if strict_protein_id is None:                                                     # Fail loudly if no usable seed identifier exists at all.
        raise ValueError("Could not resolve a Stage 07 seed from FASTA, top_candidates.csv, or strict protein_id.")
    virus_accession, source_host, seed_sequence = load_seed_sequence(
        strict_csv=strict_csv,
        protein_id=str(strict_protein_id),
    )
    return virus_accession, source_host, str(strict_protein_id), seed_sequence        # Return strict-bank seed details.


# ------------------------------ planning-field extraction helpers ------------------------------ #
# The Phase A plan evolved over the project. Some runs store hotspots and anchor
# metadata under newer names like `mutation_hotspots`, while your actual repo uses
# keys such as `mutation_window_positions_0based`, `family_size`, and
# `target_reference_rows`. The helpers below normalize those variations.
def extract_hotspot_positions(phaseA: dict) -> tuple[list[int], str]:
    """Return normalized hotspot positions plus a short source label."""
    if phaseA.get("mutation_hotspots"):                                              # Prefer an explicit hotspot list when the plan provides one.
        return [int(x) for x in phaseA.get("mutation_hotspots", [])], "mutation_hotspots"
    if phaseA.get("mutation_window_positions_0based"):                               # Otherwise use the Stage 06 position-prior window positions.
        return [int(x) for x in phaseA.get("mutation_window_positions_0based", [])], "mutation_window_positions_0based"
    return [], "none"                                                                 # Fall back to an empty hotspot list when no usable field exists.


def extract_edit_window(phaseA: dict, hotspot_positions: list[int], seed_sequence: str, flank: int) -> tuple[list[int], str]:
    """Return a 0-based half-open editable window plus a short source label."""
    base_window = phaseA.get("mutation_window", None)                                # Prefer any explicit legacy mutation window first.
    if isinstance(base_window, (list, tuple)) and len(base_window) == 2:
        return [int(base_window[0]), int(base_window[1])], "mutation_window"         # Normalize the explicit stored window to integers.
    if hotspot_positions:                                                             # Derive a compact editable window from hotspot positions when possible.
        return [
            max(0, min(hotspot_positions) - flank),                                   # Start slightly before the first hotspot while staying in bounds.
            min(len(seed_sequence), max(hotspot_positions) + flank + 1),              # End slightly after the last hotspot while staying in bounds.
        ], "derived_from_hotspots"
    return [0, len(seed_sequence)], "full_sequence_fallback"                         # Last resort: allow the full sequence to be editable.


def extract_family_context(phaseA: dict) -> dict:
    """Normalize family-level metadata from the Phase A plan into one stable dictionary."""
    family_summary = phaseA.get("family_summary", {})                                # Read the family summary block once for reuse below.
    family_member_count = first_present(                                              # Accept both the old and the current count field names.
        pd.Series(family_summary),
        ["n_family_members", "family_size"],
        default=0,
    )
    family_product_majority = first_present(                                          # Accept both the old and the current majority-product fields.
        pd.Series(family_summary),
        ["majority_product_label", "family_product_majority"],
        default=None,
    )
    if family_product_majority is None:                                               # If the family block has no product summary, inherit the canonical seed product.
        family_product_majority = first_present(
            pd.Series(phaseA.get("canonical_seed", {})),
            ["product", "product_label"],
            default=None,
        )
    return {
        "family_member_count": int(family_member_count),                            # Store the normalized family size as an integer.
        "family_cosine_floor": family_summary.get("family_cosine_floor", None),     # Preserve the cosine floor when present.
        "family_product_majority": family_product_majority,                          # Store the normalized majority product label.
        "family_member_ids": phaseA.get("family_member_ids", []),                   # Keep member IDs for downstream similarity reranking.
        "family_centroid": phaseA.get("family_centroid", None),                     # Keep the centroid if the plan already computed it.
    }


def extract_target_anchor_context(phaseA: dict, target_host: str) -> dict:
    """Normalize target-anchor references and target centroid fields for one host."""
    target_anchor_refs = phaseA.get("target_anchor_references", {}).get(target_host, None)  # Prefer the explicit newer anchor-reference field.
    if target_anchor_refs is None:
        target_anchor_refs = phaseA.get("target_reference_rows", {}).get(target_host, [])    # Fall back to the actual Stage 06 field name in your repo.
    target_reference_centroid = phaseA.get("target_reference_centroids", {}).get(target_host, None)  # Recover the target-host centroid if it exists.
    return {
        "target_anchor_references": target_anchor_refs,                             # Store normalized anchor-reference rows for downstream scoring.
        "target_reference_centroid": target_reference_centroid,                     # Store the target-host centroid when present.
    }


# Main script: Reads old outputs and writes one unified Stage 07 JSON.
def main() -> None:
    # Parse command-line arguments, load the Phase A plan, select the seed row, and infer the seed metadata.
    args = parse_args()                                                               # Parse all command-line arguments at the start of execution.

    phaseA = load_json(Path(args.phaseA_plan_json))                                   # Load the Phase A planning JSON produced by Stage 06a.

    if not args.phase06b_top_candidates_csv and not args.phase06c_followup_summary_json:  # Require at least one Stage 06 seed source.
        raise ValueError(                                                             # Raise an actionable error if the user forgot both inputs.
            "Provide either --phase06b_top_candidates_csv or --phase06c_followup_summary_json."
        )

    if args.phase06b_top_candidates_csv:                                              # Prefer the direct 06b shortlist if the user provides it.
        seed_row = load_seed_row_from_top_candidates(Path(args.phase06b_top_candidates_csv), args.seed_rank)
        source_kind = "06b_top_candidates_csv"                                        # Record which upstream artifact defined the Stage 07 seed.
    else:                                                                             # Otherwise fall back to the 06c follow-up summary JSON.
        seed_row = load_seed_row_from_followup_summary(Path(args.phase06c_followup_summary_json))
        source_kind = "06c_followup_summary_json"                                     # Record which upstream artifact defined the Stage 07 seed.

    seed_identifier_hint, source_desc = infer_seed_metadata(phaseA, seed_row)         # Infer a first-pass seed identifier and a source description string.

    # Infer the repo root from the Phase A plan path so relative FASTA / results paths
    # still resolve correctly inside Colab and in your local clone.
    repo_root = Path(args.phaseA_plan_json).resolve().parents[2]                      # results/phaseA/phaseA_plan.json -> repo root.

    # Recover the sequence for the chosen seed.
    # IMPORTANT:
    # - a Stage 06 follow-up seed may be a FASTA-exported generated candidate
    # - or a candidate_id in top_candidates.csv
    # - not necessarily a strict-bank protein_id
    # So we resolve those cases explicitly rather than always querying strict_csv.
    top_candidates_csv = None                                                         # Default to no explicit 06b top-candidate table.
    if args.phase06b_top_candidates_csv:                                              # Use the user-provided 06b top-candidate CSV when available.
        top_candidates_csv = Path(args.phase06b_top_candidates_csv)
    else:                                                                             # Otherwise try the canonical Phase A step-1 candidate table.
        guessed_top = repo_root / "results" / "phaseA" / "step1_enterobacter" / "top_candidates.csv"
        if guessed_top.exists():                                                      # Only use the guessed path if it actually exists.
            top_candidates_csv = guessed_top

    virus_accession, source_host, resolved_seed_identifier, seed_sequence = resolve_seed_sequence(
        strict_csv=Path(args.strict_csv),                                             # Use the strict bank for the final fallback case.
        seed_row=seed_row,                                                            # Pass the normalized Stage 06 seed row.
        repo_root=repo_root,                                                          # Pass the inferred repo root for robust relative-path resolution.
        top_candidates_csv=top_candidates_csv,                                        # Pass the optimizer table if present / guessable.
    )

    # If the FASTA / candidate summary did not carry a source host and the seed row also
    # lacked one, fall back to the canonical seed host so later reporting never sees "None".
    source_host = nonempty_str(                                                       # Clean up missing-string variants such as "None".
        source_host,
        default=nonempty_str(phaseA.get("canonical_seed", {}).get("source_host"), default="UNKNOWN"),
    )

    # Recover the hotspot positions and editable window from the actual field names in
    # your Phase A plan, then infer the final editable span used by Stage 07 generation.
    hotspot_positions, hotspot_source = extract_hotspot_positions(phaseA)             # Normalize hotspot positions from the plan.
    flank = (                                                                         # Decide which flank size to use around hotspots; Flank is the "breathing room" or "padding" we add around the core mutation hotspots to allow the protein to structurally compensate for those changes so that the new (mutated) protein doesn't fold incorrectly.
        int(phaseA.get("window_flank", 16))                                           # Prefer the original plan value if no override was supplied.
        if args.window_flank_override < 0                                             # Check whether the user requested a manual override.
        else int(args.window_flank_override)                                          # Otherwise use the manual override provided on the CLI.
    )
    window, window_source = extract_edit_window(phaseA, hotspot_positions, seed_sequence, flank)  # Normalize or derive the editable window.

    # Recover target-host centroid and anchor references from the Phase A plan, load optional tissue metadata, and assemble the final Stage 07 context.
    family_context = extract_family_context(phaseA)                                   # Normalize family-level metadata from the plan.
    target_anchor_context = extract_target_anchor_context(phaseA, args.target_host)   # Normalize target-anchor metadata for the selected host.

    tissue_payload = None                                                             # Default to no tissue metadata branch.
    if args.tissue_metadata_csv:                                                      # Load optional tissue metadata when the user provides it.
        tissue_df = pd.read_csv(args.tissue_metadata_csv)                             # Read the tissue CSV from disk.
        tissue_payload = tissue_df.to_dict(orient="records")                          # Convert it to JSON-serializable records.

    context = {                                                                       # Assemble the final Stage 07 context dictionary.
        "stage": "07",                                                                # Record the pipeline stage explicitly.
        "target_host": args.target_host,                                              # Store the final target host label for generation and ranking.
        "canonical_seed": phaseA.get("canonical_seed", {}),                           # Original canonical seed block from Stage 06a; The "Grandfather" seed of the ladder (i.e., the Klebsiella RBP).
        "selected_seed": {                                                            # Store the exact Stage 07 seed used for generation.
            "seed_rank": int(args.seed_rank),                                         # Record which ranked row was used when reading a shortlist.
            "seed_protein_id": resolved_seed_identifier,                              # Store the concrete resolved seed identifier used for generation.
            "seed_identifier_hint": seed_identifier_hint,                             # Preserve the original pre-resolution hint for auditability.
            "seed_source_kind": source_kind,                                          # Record whether the seed came from 06b or 06c outputs.
            "seed_source_desc": source_desc,                                          # Store a human-readable description of the upstream source.
            "virus_accession": virus_accession,                                       # Store the source phage accession or candidate-like ID for traceability.
            "source_host": source_host,                                               # Store the source host genus from the resolved seed source.
            "seed_sequence": seed_sequence,                                           # Store the amino-acid seed sequence for downstream generation.
            "sequence_length": len(seed_sequence),                                    # Store the seed length for convenience.
        },
        "family_context": family_context,                                             # Preserve normalized family-level metadata from Stage 06 planning.
        "editable_region": {                                                          # Store the editable window used by Stage 07 generation.
            "hotspot_positions": hotspot_positions,                                   # Keep the hotspot positions as discovered by Stage 06a.
            "hotspot_source": hotspot_source,                                         # Record which plan field provided the hotspots.
            "window_start": int(window[0]),                                           # Store the 0-based inclusive window start.
            "window_end": int(window[1]),                                             # Store the 0-based exclusive window end.
            "window_flank": int(flank),                                               # Store the flank size used to derive the window.
            "window_source": window_source,                                           # Record whether the window was explicit, hotspot-derived, or full-sequence fallback.
        },
        "target_anchor_context": target_anchor_context,                               # Store normalized target-anchor metadata for reranking.
        "stage06_seed_row": seed_row,                                                 # Preserve the selected upstream row for full traceability.
        "optional_tissue_context": tissue_payload,                                    # Include optional tissue metadata records when present.
    }

    out_path = Path(args.output_json)                                                 # Convert the requested output path to a Path object.
    out_path.parent.mkdir(parents=True, exist_ok=True)                                # Create the output directory if it does not yet exist.
    out_path.write_text(json.dumps(context, indent=2), encoding="utf-8")             # Serialize the Stage 07 context JSON to disk.
    print(f"Wrote: {out_path}")                                                       # Print the written path for quick confirmation.


if __name__ == "__main__":                                                            # Standard Python entrypoint guard.
    main()                                                                            # Run the CLI program.
