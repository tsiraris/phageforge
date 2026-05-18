#!/usr/bin/env python
"""Stage 10a: Build the structure-conditioned redesign context.

This script takes the late Stage 07 design context and upgrades it into a much stricter
Stage 10 redesign specification. The main change is conceptual:
- Stage 09 still searched in sequence space and only approximated structure with proxies,
- Stage 10 explicitly anchors redesign to a fixed seed scaffold PDB.

The output JSON becomes the single source of truth for Stage 10. It records:
- the selected seed sequence,
- the seed scaffold PDB used for inverse folding,
- a compact set of editable positions,
- a strict mutation budget,
- and per-position allowed substitutions grounded in family and target-host priors.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from phageforge.stage10_utils import build_stage10_context, read_json, resolve_seed_pdb, write_json



def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 10 context builder."""
    ap = argparse.ArgumentParser(description="Build the Stage 10 structure-conditioned redesign context.")
    ap.add_argument("--context_json", type=str, required=True, help="Stage 07 context JSON produced by 07a_prepare_stage07_design_context.py.")
    ap.add_argument("--strict_csv", type=str, default=None, help="Optional strict RBP bank used to sharpen target-host residue priors.")
    ap.add_argument("--seed_pdb", type=str, default=None, help="Explicit seed scaffold PDB path. If omitted, a validation directory search will be used.")
    ap.add_argument("--validation_dir", type=str, default=None, help="Optional Stage 08/09 validation directory used to locate the seed scaffold PDB automatically.")
    ap.add_argument("--output_json", type=str, required=True, help="Where to write the Stage 10 context JSON.")
    ap.add_argument("--max_edit_positions", type=int, default=6, help="Maximum number of strongly editable positions retained in the hard edit set.")
    ap.add_argument("--soft_positions", type=int, default=3, help="Additional lower-priority positions retained as a soft edit buffer.")
    ap.add_argument("--min_mutations", type=int, default=1, help="Minimum recommended mutation count for Stage 10 redesign.")
    ap.add_argument("--max_mutations", type=int, default=4, help="Maximum recommended mutation count for Stage 10 redesign.")
    ap.add_argument("--seed", type=int, default=42, help="Random seed used only for deterministic tie-breaking when selecting edit positions.")
    return ap.parse_args()



def main() -> None:
    # Read the late Stage 07 design context and optionally the strict RBP bank used to estimate target-host residue preferences.
    args = parse_args()
    context = read_json(args.context_json)
    strict_df = pd.read_csv(args.strict_csv) if args.strict_csv else None

    # Resolve the physical seed scaffold PDB that Stage 10 will treat as fixed structure during redesign.
    seed_pdb_path = resolve_seed_pdb(seed_pdb=args.seed_pdb, validation_dir=args.validation_dir)

    # Build a smaller, scaffold-preserving edit space than the one used in Stage 09.
    stage10_context = build_stage10_context(
        context=context,
        seed_pdb_path=seed_pdb_path,
        strict_df=strict_df,
        max_edit_positions=args.max_edit_positions,
        soft_positions=args.soft_positions,
        min_mutations=args.min_mutations,
        max_mutations=args.max_mutations,
        seed=args.seed,
    )

    # Persist the finalized redesign context so all downstream Stage 10 steps use one consistent configuration.
    write_json(stage10_context, args.output_json)
    print(f"Wrote: {args.output_json}")
    print(f"Seed scaffold: {seed_pdb_path}")
    print(f"Editable positions: {stage10_context['editable_region']['editable_positions']}")


if __name__ == "__main__":
    main()
