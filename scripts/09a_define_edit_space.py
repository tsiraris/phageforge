#!/usr/bin/env python
"""Stage 09a: Define a stricter, structure-aware edit space around the selected Stage 07 seed.

This script narrows the Stage 07 editable region into a smaller Stage 09 edit space that is
more likely to preserve scaffold integrity. The output is a JSON artifact that marks:
- frozen positions,
- strongly editable positions,
- softly editable positions,
- mutation-budget recommendations,
- per-position substitution proposals derived from the family bank and target-host priors.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from phageforge.stage09_utils import (
    build_edit_proposals_from_context,
    choose_editable_positions,
    read_json,
    write_json,
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 09 edit-space builder."""
    ap = argparse.ArgumentParser(description="Build the Stage 09 structure-aware edit space JSON.")                                                                  # Initialize the argument parser object
    ap.add_argument("--context_json", type=str, required=True, help="Stage 07 context JSON produced by 07a_prepare_stage07_design_context.py.")                      # Define required argument for context JSON
    ap.add_argument("--strict_csv", type=str, default=None, help="Optional strict-bank CSV used to estimate target-host residue preferences at editable positions.") # Define optional argument for strict CSV
    ap.add_argument("--output_json", type=str, required=True, help="Where to write the Stage 09 edit-space JSON.")                                                   # Define required argument for output JSON
    ap.add_argument("--max_edit_positions", type=int, default=12, help="Maximum number of strongly editable positions to retain in the hard edit set.")            # Define argument for max hard positions
    ap.add_argument("--soft_buffer_positions", type=int, default=6, help="Additional lower-priority positions to retain as a soft edit set.")                      # Define argument for soft buffer size
    ap.add_argument("--min_mutations", type=int, default=3, help="Recommended minimum mutation count for Stage 09 localized search.")                              # Define argument for minimum mutations
    ap.add_argument("--max_mutations", type=int, default=8, help="Recommended maximum mutation count for Stage 09 localized search.")                              # Define argument for maximum mutations
    ap.add_argument("--seed", type=int, default=42, help="Random seed used only for deterministic tie breaking.")                                                  # Define argument for RNG seed
    return ap.parse_args()                                                                                                                                         # Parse and return the command-line arguments



def main() -> None:
    # Read the Stage 07 context and optional strict RBP bank that will be used to sharpen target-host residue proposals.
    args = parse_args()                                                                                                           # Execute argument parsing and store into args
    context = read_json(args.context_json)                                                                                        # Read the context dict from the JSON file
    strict_df = pd.read_csv(args.strict_csv) if args.strict_csv else None                                                         # Load CSV to dataframe if provided, else None

    # Build per-position substitution proposals from the Stage 07 position features, family rows, and optional target-host bank.
    proposals = build_edit_proposals_from_context(context=context, strict_df=strict_df)                                           # Build list of mutation proposals
    if not proposals:                                                                                                             # Check if the returned proposals list is empty
        raise ValueError("No Stage 09 edit proposals could be constructed from the provided Stage 07 context.")                   # Throw exception if no proposals were created

    # Choose a compact hard edit set, then retain a small soft buffer of secondary positions for optional search expansion.
    hard_positions = choose_editable_positions(proposals=proposals, max_positions=args.max_edit_positions, seed=args.seed)        # Determine the set of hard edit positions
    sorted_by_priority = sorted(proposals, key=lambda item: (item.functional_weight, -item.conservation_penalty, item.position), reverse=True) # Sort proposals by priority metrics
    soft_positions: list[int] = []                                                                                                # Initialize an empty list for soft positions
    for item in sorted_by_priority:                                                                                               # Iterate over sorted proposals sequentially
        if int(item.position) in hard_positions:                                                                                  # Check if position is already marked as hard
            continue                                                                                                              # Skip current loop iteration if it is hard
        soft_positions.append(int(item.position))                                                                                 # Add position integer to the soft edit list
        if len(soft_positions) >= args.soft_buffer_positions:                                                                     # Check if the soft list has reached its limit
            break                                                                                                                 # Terminate loop early if soft buffer is full

    # Freeze all other seed positions so the next stage searches a small scaffold-preserving neighborhood only.
    seed_sequence = str(context["selected_seed"]["seed_sequence"])                                                                # Extract original seed sequence as string
    all_positions = set(range(1, len(seed_sequence) + 1))                                                                         # Generate a set of all 1-indexed positions
    editable_positions = set(hard_positions) | set(soft_positions)                                                                # Union hard and soft sets into editable set
    frozen_positions = sorted(all_positions - editable_positions)                                                                 # Calculate and sort the remaining frozen positions

    # Serialize only plain JSON-friendly dictionaries so the edit-space artifact is easy to inspect and version.
    proposal_rows = []                                                                                                            # Initialize empty list for serialization payload
    for item in proposals:                                                                                                        # Iterate over each proposal object again
        proposal_rows.append(                                                                                                     # Append a new dictionary to the rows list
            {                                                                                                                     # Open dictionary definition for proposal
                "position": int(item.position),                                                                                   # Convert position attribute to standard int
                "seed_aa": item.seed_aa,                                                                                          # Extract the seed amino acid character
                "allowed_aas": list(item.allowed_aas),                                                                            # Cast the allowed amino acids set to a list
                "target_preference": dict(item.target_preference),                                                                # Convert target preference mapping to a dict
                "family_preference": dict(item.family_preference),                                                                # Convert family preference mapping to a dict
                "functional_weight": float(item.functional_weight),                                                               # Cast the functional weight value to a float
                "conservation_penalty": float(item.conservation_penalty),                                                         # Cast the conservation penalty value to float
                "region_name": item.region_name,                                                                                  # Extract region name identifier string
                "edit_tier": "hard" if int(item.position) in hard_positions else ("soft" if int(item.position) in soft_positions else "frozen"), # Determine specific edit tier logic
            }                                                                                                                     # Close dictionary definition for proposal
        )                                                                                                                         # Complete the append function call

    # Write the final edit-space JSON that the localized-search stage will use as its primary design constraint.
    out = {                                                                                                                       # Begin building final output payload dictionary
        "stage": "09a",                                                                                                           # Assign the pipeline stage identifier string
        "source_context_json": str(args.context_json),                                                                            # Record the path of the source context file
        "target_host": str(context["target_host"]),                                                                               # Extract and store the target host string
        "selected_seed": context["selected_seed"],                                                                                # Embed the selected seed metadata subset
        "mutation_budget": {                                                                                                      # Open the mutation budget sub-dictionary
            "recommended_min_mutations": int(args.min_mutations),                                                                 # Record integer recommended minimum mutations
            "recommended_max_mutations": int(args.max_mutations),                                                                 # Record integer recommended maximum mutations
            "hard_edit_position_count": int(len(hard_positions)),                                                                 # Compute integer count of hard edit positions
            "soft_edit_position_count": int(len(soft_positions)),                                                                 # Compute integer count of soft edit positions
        },                                                                                                                        # Close the mutation budget sub-dictionary
        "edit_space": {                                                                                                           # Open the edit space data sub-dictionary
            "hard_edit_positions": [int(x) for x in hard_positions],                                                              # Create an integer list of hard positions
            "soft_edit_positions": [int(x) for x in soft_positions],                                                              # Create an integer list of soft positions
            "editable_positions": [int(x) for x in sorted(editable_positions)],                                                   # Create sorted integer list of editable items
            "frozen_positions": [int(x) for x in frozen_positions],                                                               # Embed the list of frozen positions
            "proposal_rows": proposal_rows,                                                                                       # Attach the built proposal rows sequence
        },                                                                                                                        # Close the edit space data sub-dictionary
    }                                                                                                                             # Close final output payload dictionary
    write_json(out, args.output_json)                                                                                             # Serialize output dictionary to JSON path
    print(f"Wrote: {args.output_json}")                                                                                           # Log successfully written file path to terminal
    print(f"hard_edit_positions: {hard_positions}")                                                                               # Log the final hard edit positions subset
    print(f"soft_edit_positions: {soft_positions}")                                                                               # Log the final soft edit positions subset


if __name__ == "__main__":                                                                                                        # Check if module is being run as main script
    main()                                                                                                                        # Execute the main function of the application