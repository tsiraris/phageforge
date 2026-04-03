"""Stage 07a: Build the compact design context for the final Acinetobacter retargeting step."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from phageforge.stage07_utils import (
    build_position_feature_table,
    build_structured_windows,
    read_json,
    write_json,
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 07 context builder."""
    ap = argparse.ArgumentParser(description="Prepare the Stage 07 design context JSON.")
    ap.add_argument("--phaseA_plan_json", type=str, required=True, help="Plan JSON from 06a_select_phaseA_family.py.")
    ap.add_argument("--phase06c_followup_summary_json", type=str, required=True, help="Follow-up seed JSON from 06c_pick_phaseA_followup_seed.py.")
    ap.add_argument("--strict_csv", type=str, required=True, help="Strict processed RBP dataset used to recover family context rows.")
    ap.add_argument("--target_host", type=str, required=True, help="Target host genus for Stage 07, such as Acinetobacter.")
    ap.add_argument("--output_json", type=str, required=True, help="Where to write the Stage 07 context JSON.")
    return ap.parse_args()


def main() -> None:
    # Read the upstream planning artifacts and the strict family dataset.
    args = parse_args()
    plan = read_json(args.phaseA_plan_json)
    followup = read_json(args.phase06c_followup_summary_json)
    strict_df = pd.read_csv(args.strict_csv)

    # Recover the selected seed row from the previous ladder step so the exact sequence is propagated into Stage 07.
    source_top_candidates_csv = Path(followup["source_top_candidates_csv"])
    if not source_top_candidates_csv.exists():
        raise FileNotFoundError(f"Missing file: {source_top_candidates_csv}")
    source_df = pd.read_csv(source_top_candidates_csv)
    chosen_id = str(followup["chosen_candidate_id"])
    selected_rows = source_df.loc[source_df["candidate_id"].astype(str) == chosen_id].reset_index(drop=True)
    if selected_rows.empty:
        raise ValueError(f"Could not find chosen candidate_id={chosen_id} in {source_top_candidates_csv}")
    selected_row = selected_rows.iloc[0]
    seed_sequence = str(selected_row["aa_sequence"])

    # Recover family rows from the strict dataset so the family context remains grounded in real proteins.
    family_member_ids = list(plan.get("family_member_ids", []))
    family_rows = strict_df.loc[strict_df["protein_id"].astype(str).isin([str(x) for x in family_member_ids])].copy().reset_index(drop=True)
    family_sequences = family_rows.get("aa_sequence", pd.Series(dtype=str)).astype(str).tolist()

    # Build a functionally informed editable region using family variability and the target host position priors.
    hotspots_0based = [int(x) for x in plan.get("mutation_window_positions_0based", [])]
    hotspots_1based = [pos + 1 for pos in hotspots_0based]
    position_features = build_position_feature_table(
        seed_sequence=seed_sequence,
        family_sequences=family_sequences,
        hotspots_1based=hotspots_1based,
        target_position_priors=plan.get("target_position_priors", {}).get(args.target_host, []),
    )
    ranked_hotspots = [int(row["position"]) for row in position_features]
    priority_weights = {str(int(row["position"])): float(row["functional_weight"]) for row in position_features}
    default_mask_count = max(24, int(plan.get("hotspot_min_count", 24)))
    structured_windows = build_structured_windows(position_features, seed_sequence=seed_sequence, window_size=default_mask_count, top_k=3)
    active_positions = ranked_hotspots or hotspots_1based
    window_start = min(active_positions) if active_positions else 1
    window_end = max(active_positions) + 1 if active_positions else len(seed_sequence) + 1

    context = {
        "stage": "07",
        "target_host": args.target_host,
        "canonical_seed": plan.get("canonical_seed", {}),
        "selected_seed": {
            "seed_rank": 0,
            "seed_protein_id": str(selected_row.get("candidate_id", chosen_id)),
            "seed_identifier_hint": str(selected_row.get("seed_protein_id", plan.get("canonical_seed", {}).get("seed_protein_id", ""))),
            "seed_source_kind": "06c_followup_summary_json",
            "seed_source_desc": f"candidate_id={chosen_id}",
            "virus_accession": str(selected_row.get("candidate_id", selected_row.get("virus_accession", ""))),
            "source_host": str(selected_row.get("source_host", plan.get("canonical_seed", {}).get("source_host", ""))),
            "seed_sequence": seed_sequence,
            "sequence_length": int(len(seed_sequence)),
        },
        "family_context": {
            "family_member_count": int(len(family_rows)),
            "family_cosine_floor": float(plan.get("family_summary", {}).get("family_cosine_floor", 0.995)),
            "family_product_majority": str(family_rows["product"].mode().iat[0] if len(family_rows) else plan.get("canonical_seed", {}).get("product", "receptor-binding protein")),
            "family_member_ids": family_member_ids,
            "family_centroid": plan.get("family_centroid", []),
            "family_rows": family_rows.to_dict(orient="records"),
        },
        "target_context": {
            "target_centroid": plan.get("target_reference_centroids", {}).get(args.target_host, []),
            "target_reference_count": int(len(plan.get("target_reference_rows", {}).get(args.target_host, []))) if isinstance(plan.get("target_reference_rows", {}).get(args.target_host, []), list) else int(plan.get("target_reference_rows", {}).get(args.target_host, 0)),
        },
        "editable_region": {
            "window_start": int(window_start),
            "window_end": int(window_end),
            "hotspot_positions": ranked_hotspots,
            "hotspot_priority_weights": priority_weights,
            "position_features": position_features,
            "structured_windows": structured_windows,
            "target_position_priors": plan.get("target_position_priors", {}).get(args.target_host, []),
            "default_max_masked_positions": int(default_mask_count),
        },
        "upstream_artifacts": {
            "phaseA_plan_json": str(args.phaseA_plan_json),
            "phase06c_followup_summary_json": str(args.phase06c_followup_summary_json),
            "source_top_candidates_csv": str(source_top_candidates_csv),
            "strict_csv": str(args.strict_csv),
        },
    }

    write_json(context, args.output_json)
    print(f"Wrote: {args.output_json}")


if __name__ == "__main__":
    main()
