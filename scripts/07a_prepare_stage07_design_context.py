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


# Helpers that help normalize noisy upstream tables.
def first_present(row: pd.Series, candidates: list[str], default: Any = None) -> Any:
    """Return the first non-null value from a candidate column list.

    Older and newer scripts sometimes use slightly different column names.
    This helper lets Stage 07 stay robust without forcing upstream rewrites.
    """
    for col in candidates:                                                            # Walk through the candidate column names in priority order.
        if col in row.index and pd.notna(row[col]):                                   # Keep the first column that exists and is not missing.
            return row[col]                                                           # Return the selected value immediately.
    return default                                                                    # Fall back to the provided default if nothing was found.


def load_json(path: Path) -> dict:
    """Load a JSON file and return it as a Python dictionary."""
    with open(path, "r", encoding="utf-8") as handle:                                 # Open the JSON file using UTF-8 for safe text handling.
        return json.load(handle)                                                      # Parse the JSON and return the resulting dictionary.


def load_seed_sequence(strict_csv: Path, protein_id: str) -> tuple[str, str, str]:
    """Recover virus accession, source host, and amino-acid sequence for a protein ID.

    Stage 06 summary files usually keep the selected protein ID, but not always the
    full amino-acid sequence. This helper looks the sequence back up in the strict
    bank so that generation starts from a real validated seed.
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
    return {                                                                          # Normalize the summary fields into a row-like dictionary.
        "candidate_id": summary.get("chosen_candidate_id"),                           # Preserve the candidate ID selected by Stage 06c.
        "selected_seed_protein_id": summary.get("chosen_candidate_id"),               # Keep a placeholder key used elsewhere when available.
        "seed_protein_id": summary.get("canonical_seed_protein_id"),                  # Store the canonical seed protein ID chosen in the Phase A plan.
        "protein_id": summary.get("canonical_seed_protein_id"),                       # Mirror the same ID under another common column name.
        "chosen_target_score": summary.get("chosen_target_score"),                    # Carry forward the target score for traceability.
        "chosen_selection_score": summary.get("chosen_selection_score"),              # Carry forward the selection score for traceability.
        "chosen_family_cosine": summary.get("chosen_family_cosine"),                  # Carry forward the family similarity for traceability.
        "chosen_target_anchor_cosine": summary.get("chosen_target_anchor_cosine"),    # Carry forward the anchor similarity for traceability.
        "chosen_n_mutations": summary.get("chosen_n_mutations"),                      # Carry forward mutation burden for traceability.
        "followup_fasta_path": summary.get("fasta_path"),                             # Record the FASTA used by the next ladder step.
        "next_target": summary.get("next_target"),                                    # Preserve the next ladder target from the summary file.
    }


def infer_seed_metadata(phaseA: dict, seed_row: dict) -> tuple[str, str | None]:
    """Infer the chosen seed protein ID and the source filename description.

    This helper lets 07a accept either the Phase 06b candidate CSV or the Phase 06c
    summary JSON while keeping one common downstream representation.
    """
    seed_protein_id = str(                                                            # Prefer explicit selected seed fields, then fall back to canonical plan data.
        first_present(
            pd.Series(seed_row),
            ["selected_seed_protein_id", "seed_protein_id", "protein_id"],
            default=phaseA.get("canonical_seed", {}).get("seed_protein_id"),
        )
    )
    source_desc = None                                                                # Initialize the source description used for auditability.
    if "candidate_id" in seed_row:                                                    # If the row came from 06b, report that candidate ID.
        source_desc = f"candidate_id={seed_row['candidate_id']}"                      # Build a human-readable source descriptor.
    elif "followup_fasta_path" in seed_row:                                           # If the row came from 06c, report the FASTA path.
        source_desc = f"followup_fasta={seed_row['followup_fasta_path']}"             # Build a human-readable source descriptor.
    return seed_protein_id, source_desc                                               # Return the normalized seed protein ID and the source descriptor.


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

    seed_protein_id, source_desc = infer_seed_metadata(phaseA, seed_row)              # Infer the final seed protein ID and a source description string.
    
    # Recover the sequence for the chosen seed and the hotspot positions to use as the mutation window, then infer the window span.
    virus_accession, source_host, seed_sequence = load_seed_sequence(                 # Recover the concrete sequence for the selected protein ID.
        Path(args.strict_csv),
        seed_protein_id,
    )

    hotspot_positions = [int(x) for x in phaseA.get("mutation_hotspots", [])]         # Recover hotspot positions discovered during Stage 06 planning.
    base_window = phaseA.get("mutation_window", None)                                 # Recover any precomputed mutation window from the plan JSON.
    flank = (                                                                         # Decide which flank size to use around hotspots; Flank is the "breathing room" or "padding" we add around the core mutation hotspots to allow the protein to structurally compensate for those changes so that the new (mutated) protein doesn't fold incorrectly.
        int(phaseA.get("window_flank", 16))                                           # Prefer the original plan value if no override was supplied.
        if args.window_flank_override < 0                                             # Check whether the user requested a manual override.
        else int(args.window_flank_override)                                          # Otherwise use the manual override provided on the CLI.
    )

    # If the plan lacks a valid explicit window, derive one from hotspots.
    if base_window is None or len(base_window) != 2:                                  
        if len(hotspot_positions) > 0:                                                # Only derive a hotspot-centered window when hotspots exist.
            window = [                                                                # Build a 0-based half-open editable region.
                max(0, min(hotspot_positions) - flank),                               # Start slightly before the first hotspot while staying in bounds.
                min(len(seed_sequence), max(hotspot_positions) + flank + 1),          # End slightly after the last hotspot while staying in bounds.
            ]
        else:                                                                         # If no hotspots exist, allow the whole sequence to be editable.
            window = [0, len(seed_sequence)]                                          # Use the full seed sequence as the editable region.
    else:                                                                             # If a plan window exists, trust it.
        window = [int(base_window[0]), int(base_window[1])]                           # Normalize the stored values to integers.

    # Recover target-host centroid and anchor references from the Phase A plan, load optional tissue metadata, and assemble the final Stage 07 context.
    target_reference_centroid = (                                                     # Recover any target-host reference centroid stored in the plan.
        phaseA.get("target_reference_centroids", {}).get(args.target_host, None)
    )
    target_anchor_refs = phaseA.get("target_anchor_references", {}).get(args.target_host, [])  # Recover target-host anchor references for reranking.

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
            "seed_protein_id": seed_protein_id,                                       # Store the protein ID that defines the seed scaffold.
            "seed_source_kind": source_kind,                                          # Record whether the seed came from 06b or 06c outputs.
            "seed_source_desc": source_desc,                                          # Store a human-readable description of the upstream source.
            "virus_accession": virus_accession,                                       # Store the source phage accession for traceability.
            "source_host": source_host,                                               # Store the source host genus from the strict bank.
            "seed_sequence": seed_sequence,                                           # Store the amino-acid seed sequence for downstream generation.
            "sequence_length": len(seed_sequence),                                    # Store the seed length for convenience.
        },
        "family_context": {                                                          # Preserve family-level metadata from Stage 06 planning.
            "family_member_count": int(phaseA.get("family_summary", {}).get("n_family_members", 0)),
            "family_cosine_floor": phaseA.get("family_summary", {}).get("family_cosine_floor", None),
            "family_product_majority": phaseA.get("family_summary", {}).get("majority_product_label", None),
            "family_member_ids": phaseA.get("family_member_ids", []),                # Keep the family member IDs for family-similarity reranking later.
        },
        "editable_region": {                                                         # Store the editable window used by Stage 07 generation.
            "hotspot_positions": hotspot_positions,                                  # Keep the hotspot positions as discovered by Stage 06a.
            "window_start": int(window[0]),                                          # Store the 0-based inclusive window start.
            "window_end": int(window[1]),                                            # Store the 0-based exclusive window end.
            "window_flank": int(flank),                                              # Store the flank size used to derive the window.
        },
        "target_anchor_context": {                                                   # Store target-host anchor references for reranking.
            "target_anchor_references": target_anchor_refs,                          # Keep the anchor records as-is from the plan JSON.
            "target_reference_centroid": target_reference_centroid,                  # Keep the target reference centroid if one was computed.
        },
        "stage06_seed_row": seed_row,                                               # Preserve the selected upstream row for full traceability.
        "optional_tissue_context": tissue_payload,                                  # Include optional tissue metadata records when present.
    }

    out_path = Path(args.output_json)                                                # Convert the requested output path to a Path object.
    out_path.parent.mkdir(parents=True, exist_ok=True)                               # Create the output directory if it does not yet exist.
    out_path.write_text(json.dumps(context, indent=2), encoding="utf-8")             # Serialize the Stage 07 context JSON to disk.
    print(f"Wrote: {out_path}")                                                      # Print the written path for quick confirmation.


if __name__ == "__main__":                                                           # Standard Python entrypoint guard.
    main()                                                                            # Run the CLI program.
