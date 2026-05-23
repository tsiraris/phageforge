#!/usr/bin/env python
"""Stage 11a: Build the self-contained Stage 11 redesign context.

Stage 11 starts from any wild-type RBP row in `rbp_dataset_eskapee_strict.csv`
and computes everything it needs to redesign that seed toward a different host
(family/target centroids, edit space, hotspot priors) **inside this script**.
It deliberately reads NO context, edit space, surrogate model, or candidate
sequence from Stages 06/07/09/10 — the Stage 10 post-mortem traced every
downstream structural failure to a corrupted Stage 06 chassis, and Stage 11's
core invariant is that the seed must be foldable in its wild-type form.

This script is therefore the gate stage. It:

1. Picks one seed protein_id from the strict CSV (with source-host validation).
2. ESMFolds the seed sequence in-process and writes `seed_wt.pdb`.
3. Runs the Baseline Qualification gate (default pLDDT >= 70.0). If the seed
   fails, the pipeline aborts with exit code 2 so the operator picks another
   seed — no compute is spent on the inverse-folding search.
4. Loads (or, with --no_reuse_cached_embeddings, recomputes) ESM-2 strict
   embeddings, finds the top-N family members + top-M target-host members,
   and computes their centroids.
5. Builds the from-scratch edit space (Shannon entropy + AA preferences),
   selects hard/soft/frozen positions, and packs everything into
   `stage11_context.json` (08a-compatible).
6. Writes a `run_metadata.json` capturing every CLI arg, GPU info, package
   versions, and a UTC ISO-8601 timestamp.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from phageforge.stage11_utils import (
    EXIT_BASELINE_QUALIFICATION_FAILED,
    EXIT_INFERENCE_ERROR,
    EXIT_INPUT_ERROR,
    EXIT_OK,
    baseline_qualification_gate,
    build_stage11_context,
    build_stage11_edit_proposals,
    compute_centroid,
    default_run_name,
    embed_strict_set_fresh,
    esmfold_single_sequence,
    find_family_members,
    find_target_members,
    load_embedding_backend,
    load_strict_dataset,
    load_strict_embeddings,
    seed_everything,
    select_seed_row,
    stage11_choose_editable_positions,
    write_json,
    write_run_metadata,
)


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for the Stage 11 context builder."""
    ap = argparse.ArgumentParser(description="Build the Stage 11 self-contained redesign context (gate stage).")        # CLI description.
    # ----- Required inputs -----
    ap.add_argument("--strict_csv", type=str, required=True,                                                            # Strict RBP CSV path.
                    help="Path to data/processed/rbp_dataset_eskapee_strict.csv.")
    ap.add_argument("--seed_protein_id", type=str, required=True,                                                       # Required: which strict row to use as the chassis.
                    help="protein_id of the strict RBP to use as the Stage 11 seed (no auto-pick).")
    ap.add_argument("--source_host", type=str, required=True,                                                           # Required: the seed's native host (validated against the CSV).
                    help="Native host genus of the seed (validated against the strict CSV).")
    ap.add_argument("--target_host", type=str, required=True,                                                           # Required: the host we are flipping toward.
                    help="Host genus the redesign is aiming at (must differ from --source_host).")
    ap.add_argument("--out_dir", type=str, required=True,                                                               # Required: where the run tree lives.
                    help="Root directory for the Stage 11 run (e.g. results/stage11/<run>/).")
    # ----- Cached embeddings (toggleable) -----
    ap.add_argument("--strict_embeddings", type=str,                                                                    # Optional: cached embeddings path.
                    default="data/processed/strict/esm2_embeddings.pt",
                    help="Cached ESM-2 strict embeddings (.pt). Used only when --reuse_cached_embeddings is set.")
    ap.add_argument("--strict_embeddings_index", type=str,                                                              # Optional: index CSV for the cached embeddings.
                    default="data/processed/strict/esm2_embeddings_index.csv",
                    help="Index CSV (row_id → protein_id) for the cached strict embeddings.")
    ap.add_argument("--reuse_cached_embeddings", action="store_true",                                                   # Opt-in flag for the ~5-minute speedup.
                    help="Reuse the cached ESM-2 strict embeddings instead of recomputing them inside Stage 11.")
    # ----- Trained host probe (treated as foundational input) -----
    ap.add_argument("--predictor_model", type=str,                                                                      # Trained logistic regression .joblib.
                    default="results/broad/linear_probe/seed_42/model.joblib",
                    help="Trained host-prediction probe used by Stage 11b (recorded in provenance).")
    ap.add_argument("--predictor_label_classes", type=str,                                                              # Label classes JSON.
                    default="results/broad/linear_probe/seed_42/label_classes.json",
                    help="Label classes JSON for the trained host probe.")
    # ----- Embedding model -----
    ap.add_argument("--embedding_model", type=str,                                                                      # ESM-2 model id.
                    default="facebook/esm2_t33_650M_UR50D",
                    help="HuggingFace model id of the ESM-2 backbone (must match the probe's training backbone).")
    # ----- ESMFold knobs -----
    ap.add_argument("--esmfold_device", type=str, default="auto",                                                       # cuda / cpu / auto.
                    help="Device for ESMFold (auto = cuda if available, else cpu).")
    ap.add_argument("--esmfold_chunk_size", type=int, default=128,                                                      # VRAM-vs-throughput knob.
                    help="ESMFold chunk size (reduce on OOM).")
    ap.add_argument("--esmfold_num_recycles", type=int, default=1,                                                      # Recycle depth.
                    help="ESMFold recycle iterations (1 is the Stage 08 default).")
    # ----- Baseline Qualification gate -----
    ap.add_argument("--min_seed_plddt", type=float, default=70.0,                                                       # The single most important guard.
                    help="Baseline Qualification threshold; seeds folding below this abort the run with exit 2.")
    # ----- Family / target context -----
    ap.add_argument("--family_top_n", type=int, default=32,                                                             # Family pool size.
                    help="Number of family members to keep (top-N by ESM-2 cosine to the seed).")
    ap.add_argument("--target_top_m", type=int, default=8,                                                              # Target pool size.
                    help="Number of target-host members to keep (top-M by ESM-2 cosine to the seed).")
    ap.add_argument("--family_cosine_floor", type=float, default=0.85,                                                  # Cosine quality gate.
                    help="Minimum cosine to seed for a row to qualify as a family member.")
    ap.add_argument("--length_tolerance", type=float, default=0.05,                                                     # Relative length gate.
                    help="Maximum relative length difference (|len_other-len_seed|/len_seed) for family/target inclusion.")
    # ----- Edit-space construction -----
    ap.add_argument("--max_edit_positions", type=int, default=6,                                                        # Hard cap.
                    help="Maximum number of hard editable positions.")
    ap.add_argument("--soft_positions", type=int, default=3,                                                            # Soft buffer.
                    help="Number of soft (lower-priority) editable positions.")
    ap.add_argument("--min_mutations", type=int, default=1,                                                             # Lower mutation budget.
                    help="Minimum number of mutations per Stage 11 candidate.")
    ap.add_argument("--max_mutations", type=int, default=4,                                                             # Upper mutation budget.
                    help="Maximum number of mutations per Stage 11 candidate.")
    ap.add_argument("--entropy_floor", type=float, default=0.20,                                                        # Conservation gate.
                    help="Per-position entropy floor for inclusion in the edit space.")
    ap.add_argument("--family_top_k", type=int, default=4,                                                              # Per-position AA preference depth (family).
                    help="Top-K family-side AAs to consider per position when building the allowed-AA list.")
    ap.add_argument("--target_top_k", type=int, default=4,                                                              # Per-position AA preference depth (target).
                    help="Top-K target-side AAs to consider per position when building the allowed-AA list.")
    ap.add_argument("--max_allowed_aas_per_pos", type=int, default=6,                                                   # Branching cap.
                    help="Max distinct AAs allowed per position (after combining family and target preferences).")
    ap.add_argument("--region_block", type=int, default=50,                                                             # Region bucket size.
                    help="Region bucket size in residues for spreading hard positions.")
    # ----- Misc -----
    ap.add_argument("--seed", type=int, default=42,                                                                     # Deterministic seed.
                    help="Random seed used for deterministic tie-breaking when selecting positions.")
    ap.add_argument("--run_name", type=str, default=None,                                                                # Optional override.
                    help="Override the auto-generated run name; default is <source>_to_<target>_<seedid>_seed<N>_<UTC>.")
    return ap.parse_args()                                                                                              # Parse and return.


def main() -> None:                                                                                                     # Main orchestration entry point.
    args = parse_args()                                                                                                 # CLI parsing.
    seed_everything(args.seed)                                                                                          # Lock all RNGs for reproducibility.

    # ----- Resolve the run directory (auto-suffix with UTC timestamp unless overridden) -----
    run_name = args.run_name or default_run_name(args.source_host, args.target_host, args.seed_protein_id, args.seed)   # Default run name uses the design-doc convention.
    out_root = Path(args.out_dir)                                                                                       # Root directory the operator passed in.
    if out_root.name != run_name and not any(part == run_name for part in out_root.parts):                              # Only nest inside the operator's path if not already there.
        out_root = out_root / run_name                                                                                  # Nest the run name under the operator's choice.
    context_dir = out_root / "context"                                                                                  # Subdirectory for the context artifacts.
    context_dir.mkdir(parents=True, exist_ok=True)                                                                      # Create the tree if missing.

    # ----- Validate basic CLI invariants before doing any compute -----
    if str(args.source_host) == str(args.target_host):                                                                  # Refuse a no-op redesign.
        print(f"[ERROR] --source_host and --target_host must differ (got '{args.source_host}').", file=sys.stderr)      # Loud error.
        sys.exit(EXIT_INPUT_ERROR)                                                                                      # Exit 1.

    # ----- Load and validate the strict CSV; pick the seed row -----
    try:                                                                                                                # Wrap CSV / seed selection in a try/except so we exit cleanly on bad input.
        strict_df = load_strict_dataset(args.strict_csv)                                                                # Load & validate the strict CSV.
        seed_row = select_seed_row(strict_df, args.seed_protein_id, args.source_host)                                   # Pick the seed row and validate host.
    except (FileNotFoundError, ValueError) as exc:                                                                      # CSV missing or seed invalid.
        print(f"[ERROR] {exc}", file=sys.stderr)                                                                        # Surface the message.
        sys.exit(EXIT_INPUT_ERROR)                                                                                      # Exit 1.

    seed_sequence = str(seed_row["aa_sequence"])                                                                        # Wild-type sequence.
    print(f"[INFO] Stage 11a — seed={args.seed_protein_id} ({args.source_host} → {args.target_host}), "                 # Loud progress line.
          f"length={len(seed_sequence)} AA, run_name={run_name}", flush=True)

    # ----- ESMFold the seed and run the Baseline Qualification gate -----
    seed_pdb_path = context_dir / "seed_wt.pdb"                                                                         # Destination for the wild-type PDB.
    print(f"[INFO] Folding seed with ESMFold (device={args.esmfold_device}, chunk={args.esmfold_chunk_size}, "          # Heads-up before the heavy step.
          f"recycles={args.esmfold_num_recycles})…", flush=True)
    try:                                                                                                                # Wrap heavy folding in a try/except to translate OOM into exit 3.
        seed_metrics = esmfold_single_sequence(                                                                         # Single-sequence ESMFold call.
            sequence=seed_sequence,                                                                                     # Seed sequence.
            out_pdb_path=seed_pdb_path,                                                                                 # Where to write the PDB.
            device=args.esmfold_device,                                                                                 # Device knob.
            chunk_size=args.esmfold_chunk_size,                                                                         # Chunk size knob.
            num_recycles=args.esmfold_num_recycles,                                                                     # Recycle knob.
        )
    except Exception as exc:                                                                                            # ESMFold failed for some reason.
        print(f"[ERROR] ESMFold failed on the seed: {exc}", file=sys.stderr)                                            # Surface the message.
        sys.exit(EXIT_INFERENCE_ERROR)                                                                                  # Exit 3.

    passed, reason = baseline_qualification_gate(seed_metrics, args.min_seed_plddt)                                     # Run the gate.
    bq = {                                                                                                              # Build the baseline_qualification dict.
        "passed": bool(passed),                                                                                         # Pass / fail flag.
        "reason": reason,                                                                                               # Human-readable explanation.
        "min_plddt_threshold": float(args.min_seed_plddt),                                                              # The threshold used.
        "seed_mean_plddt": float(seed_metrics["mean_plddt"]),                                                           # Actual seed pLDDT.
        "seed_per_residue_plddt_summary": dict(seed_metrics["per_residue_plddt_summary"]),                              # Summary stats.
        "n_residues": int(seed_metrics["n_residues"]),                                                                  # Folded residue count.
        "device_used": str(seed_metrics["device_used"]),                                                                # Where the work happened.
        "elapsed_seconds": float(seed_metrics["elapsed_seconds"]),                                                      # Wall-clock cost.
        "seed_pdb_path": str(seed_pdb_path.resolve()),                                                                  # Where the PDB lives.
        "seed_protein_id": str(seed_row["protein_id"]),                                                                 # Which seed was tested.
    }
    write_json(bq, context_dir / "baseline_qualification.json")                                                         # Persist the gate result immediately so even failed runs leave evidence.
    if not passed:                                                                                                      # Gate failed — abort before any further compute.
        print(f"[GATE] Baseline Qualification FAILED: {reason}", file=sys.stderr, flush=True)                           # Loud message to stderr so wrapping shells see it.
        print(f"[GATE] Wrote gate evidence to: {context_dir / 'baseline_qualification.json'}", file=sys.stderr)         # Tell the operator where to look.
        sys.exit(EXIT_BASELINE_QUALIFICATION_FAILED)                                                                    # Exit 2: separate code for wrapper dispatch.
    print(f"[GATE] Baseline Qualification PASSED: {reason}", flush=True)                                                # Loud pass message.

    # ----- Load the ESM-2 backbone (used for both embedding & centroid math) -----
    print(f"[INFO] Loading ESM-2 backbone '{args.embedding_model}' for centroid math…", flush=True)                     # Heads-up.
    try:                                                                                                                # Wrap backend load.
        torch_emb, tokenizer_emb, model_emb, emb_device = load_embedding_backend(args.embedding_model)                  # Load tokenizer + model + device.
    except Exception as exc:                                                                                            # Capture import / load failures.
        print(f"[ERROR] Failed to load ESM-2 backbone: {exc}", file=sys.stderr)                                         # Surface.
        sys.exit(EXIT_INFERENCE_ERROR)                                                                                  # Exit 3.

    # ----- Strict embeddings: cached (fast) or fresh (strict purity) -----
    if args.reuse_cached_embeddings:                                                                                    # Operator opted into caching.
        try:                                                                                                            # Try to load cached.
            strict_embeddings, strict_index = load_strict_embeddings(                                                   # Load the cached .pt + index.
                args.strict_embeddings, args.strict_embeddings_index)
            strict_emb_source = "cached"                                                                                # Provenance tag.
        except Exception as exc:                                                                                        # Fall through with a clear message.
            print(f"[ERROR] Failed to load cached strict embeddings: {exc}", file=sys.stderr)                           # Surface.
            sys.exit(EXIT_INPUT_ERROR)                                                                                  # Exit 1 — bad input config.
    else:                                                                                                               # Default: recompute fresh inside Stage 11.
        print("[INFO] --reuse_cached_embeddings not set; recomputing ESM-2 embeddings of the strict CSV from scratch…", # Loud heads-up.
              flush=True)
        strict_embeddings, strict_index = embed_strict_set_fresh(                                                       # Recompute fresh.
            strict_df=strict_df,                                                                                        # The validated strict CSV.
            embedding_model=args.embedding_model,                                                                       # Same backbone.
            batch_size=4,                                                                                               # Safe default for SageMaker g5.xlarge.
            torch=torch_emb, tokenizer=tokenizer_emb, model=model_emb, device=emb_device,                               # Reuse the persistent backend.
        )
        strict_emb_source = "fresh"                                                                                     # Provenance tag.

    # ----- Locate the seed's embedding row -----
    seed_idx_hits = strict_index.index[strict_index["protein_id"].astype(str) == str(args.seed_protein_id)].tolist()    # Find the seed's row in the index.
    if not seed_idx_hits:                                                                                               # Defensive: the seed must be in the embedding cache too.
        print(f"[ERROR] Seed '{args.seed_protein_id}' is not present in the strict embeddings index.", file=sys.stderr) # Surface.
        sys.exit(EXIT_INPUT_ERROR)                                                                                      # Exit 1.
    seed_row_in_index = int(seed_idx_hits[0])                                                                           # Row id in the embedding matrix.
    seed_embedding = np.asarray(strict_embeddings[seed_row_in_index], dtype=np.float32)                                 # Pull the seed vector.

    # ----- Family + target pools and their centroids -----
    print(f"[INFO] Finding family (top-{args.family_top_n}) and target (top-{args.target_top_m}) members…", flush=True) # Heads-up.
    family_df = find_family_members(                                                                                    # Family pool: cosine + length-matched, includes seed.
        seed_embedding=seed_embedding,
        seed_protein_id=args.seed_protein_id,
        strict_embeddings=strict_embeddings,
        strict_index=strict_index,
        strict_df=strict_df,
        top_n=args.family_top_n,
        cosine_floor=args.family_cosine_floor,
        length_tolerance=args.length_tolerance,
    )
    if len(family_df) < 3:                                                                                              # Refuse a degenerate family pool.
        print(f"[ERROR] Family pool has only {len(family_df)} members after length+cosine filtering; "                  # Surface.
              f"consider lowering --family_cosine_floor or widening --length_tolerance.", file=sys.stderr)
        sys.exit(EXIT_INPUT_ERROR)                                                                                      # Exit 1.

    target_df = find_target_members(                                                                                    # Target pool: cosine + length-matched, excludes seed.
        seed_embedding=seed_embedding,
        seed_protein_id=args.seed_protein_id,
        target_host=args.target_host,
        strict_embeddings=strict_embeddings,
        strict_index=strict_index,
        strict_df=strict_df,
        top_m=args.target_top_m,
        length_tolerance=args.length_tolerance,
    )
    if len(target_df) < 1:                                                                                              # Refuse an empty target pool.
        print(f"[ERROR] Target pool for host '{args.target_host}' is empty after length filtering; "                    # Surface.
              f"consider widening --length_tolerance.", file=sys.stderr)
        sys.exit(EXIT_INPUT_ERROR)                                                                                      # Exit 1.

    family_member_indices = [int(strict_index.index[strict_index["protein_id"].astype(str) == pid].tolist()[0])         # Map family protein_ids → embedding rows.
                              for pid in family_df["protein_id"].astype(str).tolist()]
    target_member_indices = [int(strict_index.index[strict_index["protein_id"].astype(str) == pid].tolist()[0])         # Map target protein_ids → embedding rows.
                              for pid in target_df["protein_id"].astype(str).tolist()]
    family_centroid = compute_centroid(strict_embeddings[family_member_indices])                                        # L2-normalized mean of family vectors.
    target_centroid = compute_centroid(strict_embeddings[target_member_indices])                                        # L2-normalized mean of target vectors.

    # ----- Build the edit space from family + target alignment columns -----
    family_aligned = family_df["aa_sequence"].astype(str).tolist()                                                      # Length-truncate happens inside the helper.
    target_aligned = target_df["aa_sequence"].astype(str).tolist()                                                      # Same handling for target sequences.
    print(f"[INFO] Building edit space (entropy_floor={args.entropy_floor}, "                                            # Heads-up.
          f"family_top_k={args.family_top_k}, target_top_k={args.target_top_k}, "
          f"max_allowed={args.max_allowed_aas_per_pos})…", flush=True)
    proposals = build_stage11_edit_proposals(                                                                           # The from-scratch edit-space builder.
        seed_sequence=seed_sequence,
        family_aligned_seqs=family_aligned,
        target_aligned_seqs=target_aligned,
        entropy_floor=args.entropy_floor,
        family_top_k=args.family_top_k,
        target_top_k=args.target_top_k,
        max_allowed_aas_per_pos=args.max_allowed_aas_per_pos,
        region_block=args.region_block,
    )
    if not proposals:                                                                                                   # Refuse a fully frozen seed.
        print("[ERROR] No editable proposals survived the entropy floor; "                                              # Surface.
              "consider lowering --entropy_floor.", file=sys.stderr)
        sys.exit(EXIT_INPUT_ERROR)                                                                                      # Exit 1.

    hard_positions, soft_positions, _ = stage11_choose_editable_positions(                                              # Pick hard + soft positions with regional spreading.
        proposals=proposals,
        max_hard=args.max_edit_positions,
        max_soft=args.soft_positions,
        seed=args.seed,
    )
    print(f"[INFO] Editable positions chosen: hard={hard_positions}, soft={soft_positions}", flush=True)                # Loud summary.

    # ----- Assemble the full Stage 11 context dict and persist it -----
    context = build_stage11_context(                                                                                    # Master context builder.
        seed_row=seed_row,
        target_host=args.target_host,
        source_host=args.source_host,
        seed_pdb_path=seed_pdb_path,
        baseline_qualification=bq,
        family_df=family_df,
        target_df=target_df,
        family_centroid=family_centroid,
        target_centroid=target_centroid,
        proposals=proposals,
        hard_positions=hard_positions,
        soft_positions=soft_positions,
        min_mutations=args.min_mutations,
        max_mutations=args.max_mutations,
        cli_args=vars(args),
        embedding_model=args.embedding_model,
        predictor_model_path=args.predictor_model,
        predictor_label_classes_path=args.predictor_label_classes,
        strict_csv_path=args.strict_csv,
        strict_embeddings_path=(args.strict_embeddings if args.reuse_cached_embeddings else None),                      # Only record cached path when actually used.
    )
    context_path = context_dir / "stage11_context.json"                                                                 # Final destination.
    write_json(context, context_path)                                                                                   # Persist.
    print(f"[OK] Wrote: {context_path}", flush=True)                                                                    # Confirmation.

    # ----- Write the run-level metadata file -----
    write_run_metadata(                                                                                                 # Capture every CLI knob + GPU + git + package versions.
        out_path=out_root / "run_metadata.json",                                                                        # Top-level location for cross-phase consumption.
        cli_args=vars(args),                                                                                            # Pass every parsed CLI flag.
        extras={                                                                                                        # Stage-11-specific extras.
            "run_name": run_name,                                                                                       # Resolved run name.
            "strict_embeddings_source": strict_emb_source,                                                              # "cached" or "fresh".
            "family_member_count": int(len(family_df)),                                                                 # For quick auditing.
            "target_member_count": int(len(target_df)),                                                                 # For quick auditing.
            "hard_positions": hard_positions,                                                                           # Final hard list.
            "soft_positions": soft_positions,                                                                           # Final soft list.
            "seed_mean_plddt": float(seed_metrics["mean_plddt"]),                                                       # Baseline Qualification result.
            "seed_pdb_path": str(seed_pdb_path.resolve()),                                                              # Where the seed PDB lives.
            "context_json_path": str(context_path.resolve()),                                                           # Where the context JSON lives.
        },
    )
    print(f"[OK] Wrote: {out_root / 'run_metadata.json'}", flush=True)                                                  # Confirmation.
    print(f"[OK] Stage 11a complete. Next: scripts/11b_run_inverse_folding_beam_search.py "                             # Tell the operator the next step.
          f"--stage11_context_json {context_path}", flush=True)
    sys.exit(EXIT_OK)                                                                                                   # Explicit success exit.


if __name__ == "__main__":                                                                                              # Standard CLI guard.
    main()                                                                                                              # Run.
