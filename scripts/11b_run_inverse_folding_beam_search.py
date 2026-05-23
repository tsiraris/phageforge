#!/usr/bin/env python
"""Stage 11b: Run inverse-folding beam search anchored to the Stage 11 seed scaffold.

This is the structural core of Stage 11. It mirrors Stage 10b's algorithm
(ESM-IF1 beam search with composite scoring) but reads its context exclusively
from `stage11_context.json` written by 11a and does not touch any artifact
produced by Stages 06/07/09/10.

Operationally it:
- loads the freshly-folded wild-type seed PDB (the 11a Baseline Qualification gate
  has already confirmed it is structurally healthy),
- expands a compact set of allowed local substitutions (from the from-scratch
  Stage 11 edit space),
- scores every proposal for target-host compatibility (trained linear probe),
- scores every proposal for backbone compatibility (ESM-IF1 log-likelihood),
- and keeps only the strongest, structurally diverse beam of candidates per round.

The output CSV deliberately includes both Stage 10's column names (`mutated_positions`,
`stage10_composite_score`) and 08a's legacy aliases (`mutation_positions`,
`final_multimodal_rank_score`, `generation_regime`) so the unmodified
08a_structural_fasttrack_validation.py can consume it via the Adapter pattern.
"""

from __future__ import annotations

# --- SageMaker / Python 3.12 compatibility shim -------------------------------------------------------------------------
# fair-esm's inverse-folding GVP layers import `torch_scatter` and `torch_sparse`, which fail to compile from source on
# Python 3.12 + recent CUDA toolchains. The only operation fair-esm invokes is `scatter_add`, so when the real compiled
# libraries are unavailable we register a tiny pure-PyTorch stand-in. This MUST run before any `import esm`. On a machine
# where the genuine torch_scatter is installed, the real library is used and this shim is a no-op.
import sys as _sys, types as _types                                                                                       # Minimal imports for module registration.
try:                                                                                                                      # Prefer the real, compiled extension when present.
    import torch_scatter as _real_ts                                                                                      # Attempt to load the genuine build.
    _ = _real_ts.scatter_add                                                                                              # Touch the symbol to confirm the binary actually loaded.
except Exception:                                                                                                         # Fall back only when torch_scatter is missing / fails to load.
    import torch as _torch                                                                                                # Local torch handle for the stand-in.
    def _mock_scatter_add(src, index, dim=-1, out=None, dim_size=None):                                                   # Pure-PyTorch replacement for scatter_add.
        if out is None:                                                                                                   # When the caller does not preallocate output.
            size = list(src.size())                                                                                       # Copy source shape.
            size[dim] = dim_size if dim_size is not None else (int(index.max()) + 1 if index.numel() else 0)              # Resolve scatter-dim length.
            out = _torch.zeros(size, dtype=src.dtype, device=src.device)                                                  # Allocate zeros on the same device/dtype.
        return out.scatter_add_(dim, index, src)                                                                          # Delegate to PyTorch's native in-place scatter-add.
    def _mock_scatter(src, index, dim=-1, out=None, dim_size=None, reduce="sum"):                                         # Replacement for scatter (sum/add only).
        if reduce in ("sum", "add"):                                                                                      # fair-esm only requests additive reduction.
            return _mock_scatter_add(src, index, dim, out, dim_size)                                                      # Reuse the scatter_add path.
        raise NotImplementedError(f"torch_scatter shim does not implement reduce={reduce!r}")                             # Fail loudly on anything else.
    _ts = _types.ModuleType("torch_scatter")                                                                              # Synthesize the module object.
    _ts.scatter_add = _mock_scatter_add; _ts.scatter = _mock_scatter                                                      # Expose the API fair-esm uses.
    _sys.modules["torch_scatter"] = _ts                                                                                   # Register the stand-in.
    _sys.modules.setdefault("torch_sparse", _types.ModuleType("torch_sparse"))                                           # fair-esm only imports the torch_sparse namespace; an empty stub suffices.
# ------------------------------------------------------------------------------------------------------------------------

import argparse
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd

from phageforge.stage11_utils import (
    EXIT_INFERENCE_ERROR,
    EXIT_INPUT_ERROR,
    EXIT_OK,
    Stage10Candidate,
    apply_mutation,
    composite_stage10_score,
    embed_sequences_with_backend,
    evaluate_candidate_table,
    greedy_diverse_subset,
    load_embedding_backend,
    load_inverse_folding_model,
    load_inverse_folding_structure,
    load_target_predictor,
    mutation_list,
    read_json,
    seed_everything,
    write_json,
)


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for the Stage 11 inverse-folding beam search."""
    ap = argparse.ArgumentParser(description="Run Stage 11 structure-conditioned inverse-folding search.")              # CLI description.
    # ----- Required inputs -----
    ap.add_argument("--stage11_context_json", type=str, required=True,                                                  # Required context path.
                    help="Path to stage11_context.json produced by 11a.")
    ap.add_argument("--predictor_model", type=str, required=True,                                                       # Trained host probe.
                    help="Trained logistic-regression host probe (.joblib).")
    ap.add_argument("--label_classes_json", type=str, required=True,                                                    # Label classes file.
                    help="Label classes JSON for the trained host probe.")
    ap.add_argument("--out_csv", type=str, required=True,                                                               # Output CSV.
                    help="Where to write the full Stage 11 search table.")
    ap.add_argument("--out_json", type=str, required=True,                                                              # Output summary JSON.
                    help="Where to write the compact Stage 11 search summary JSON.")
    # ----- Embedding + IF1 knobs -----
    ap.add_argument("--embedding_model", type=str, default="facebook/esm2_t33_650M_UR50D",                              # ESM-2 backbone.
                    help="Embedding backbone used by the trained predictor and family centroid scoring.")
    ap.add_argument("--if_device", type=str, default="cuda",                                                            # IF1 device.
                    help="Device for the ESM-IF1 model.")
    ap.add_argument("--if_chain_id", type=str, default="A",                                                             # Chain ID inside the PDB.
                    help="Chain identifier used when loading the seed scaffold.")
    # ----- Search knobs -----
    ap.add_argument("--rounds", type=int, default=4,                                                                    # Beam depth.
                    help="Number of sequential redesign rounds.")
    ap.add_argument("--beam_width", type=int, default=24,                                                               # Survivors per round.
                    help="Number of candidates retained after each redesign round.")
    ap.add_argument("--proposals_per_parent", type=int, default=8,                                                      # Branching factor per parent.
                    help="Maximum number of proposal branches created from each parent.")
    ap.add_argument("--substitutions_per_position", type=int, default=3,                                                # Per-position breadth.
                    help="Maximum number of amino-acid substitutions tried at each editable position.")
    ap.add_argument("--batch_size", type=int, default=4,                                                                # VRAM safety.
                    help="Batch size for embedding candidates with the predictor backbone.")
    # ----- Composite scoring weights (locked at Stage 10 defaults for v1) -----
    ap.add_argument("--w_target", type=float, default=0.30,                                                             # Target-host weight.
                    help="Composite weight on target_probability (default locked at Stage 10).")
    ap.add_argument("--w_if1", type=float, default=0.45,                                                                # IF1 weight.
                    help="Composite weight on if1_log_likelihood (default locked at Stage 10).")
    ap.add_argument("--w_family", type=float, default=0.15,                                                             # Family cosine weight.
                    help="Composite weight on family_cosine (default locked at Stage 10).")
    ap.add_argument("--w_identity", type=float, default=0.10,                                                           # Identity weight.
                    help="Composite weight on seed_identity (default locked at Stage 10).")
    ap.add_argument("--w_mut_penalty", type=float, default=0.10,                                                        # Mutation count penalty (subtracted).
                    help="Composite penalty weight on mutation_count (subtracted; default locked at Stage 10).")
    # ----- Misc -----
    ap.add_argument("--seed", type=int, default=42,                                                                     # Determinism.
                    help="Random seed used for deterministic search ordering.")
    return ap.parse_args()                                                                                              # Done.


def build_round_children(
    parents: list[Stage10Candidate],
    proposal_rows: list[dict],
    seed_sequence: str,
    max_mutations: int,
    proposals_per_parent: int,
    substitutions_per_position: int,
) -> list[Stage10Candidate]:
    """Spawn one generation of edits from the current beam.

    Sorts proposals by functional priority, refuses to double-edit the same
    position within a single parent sequence, scores each allowed AA via the
    same blended formula Stage 10 uses (0.55*target + 0.30*family + 0.15*functional),
    and emits up to `proposals_per_parent` children per parent.
    """
    children: list[Stage10Candidate] = []                                                                               # Accumulator for this round's children.
    proposal_rows = sorted(                                                                                             # Sort the rulebook by descending priority.
        proposal_rows,                                                                                                  # Pass the rulebook.
        key=lambda row: (                                                                                               # Sort hierarchically.
            float(row.get("functional_weight", 0.0)),                                                                   # Highest functional weight first.
            -float(row.get("conservation_penalty", 0.0)),                                                               # Lowest conservation penalty next.
            int(row["position"]),                                                                                       # Then by position for determinism.
        ),
        reverse=True,                                                                                                   # Descending.
    )

    for parent in parents:                                                                                              # Walk every elite parent.
        if len(parent.mutated_positions) >= max_mutations:                                                              # Parent already at mutation ceiling.
            continue                                                                                                    # Freeze it.
        parent_mutated = set(int(p) for p in parent.mutated_positions)                                                  # Fast lookup for double-edit prevention.
        emitted = 0                                                                                                     # Counter for branches spawned by this parent.
        for row in proposal_rows:                                                                                       # Walk the prioritized rulebook.
            pos = int(row["position"])                                                                                  # Position.
            if pos in parent_mutated:                                                                                   # Forbid double-editing.
                continue                                                                                                # Skip.

            allowed = list(row.get("allowed_aas", []))                                                                  # Allowed AAs at this position.
            target_preference = dict(row.get("target_preference", {}))                                                  # Target preferences.
            family_preference = dict(row.get("family_preference", {}))                                                  # Family preferences.
            functional_weight = float(row.get("functional_weight", 0.0))                                                # Position priority.

            aa_scores: list[tuple[str, float]] = []                                                                     # Rank allowed AAs by blended priority.
            for aa in allowed:                                                                                          # Walk every legal AA.
                score = (                                                                                               # Blended priority — matches Stage 10b.
                    0.55 * float(target_preference.get(aa, 0.0))                                                        # Heavy target weight.
                    + 0.30 * float(family_preference.get(aa, 0.0))                                                      # Moderate family weight.
                    + 0.15 * functional_weight                                                                          # Light functional weight.
                )
                aa_scores.append((aa, score))                                                                           # Save.
            aa_scores.sort(key=lambda item: (item[1], item[0]), reverse=True)                                           # Descending by score, then alphabetical tie-break.

            for aa, _ in aa_scores[: int(substitutions_per_position)]:                                                  # Take top-N AAs at this position.
                next_seq = apply_mutation(parent.candidate_sequence, pos, aa)                                           # Apply the substitution.
                if next_seq == parent.candidate_sequence:                                                               # No-op (seed AA was already there) — skip.
                    continue                                                                                            # Skip.
                child = Stage10Candidate(                                                                               # Build the child object.
                    candidate_sequence=next_seq,                                                                        # New sequence.
                    parent_sequence=parent.parent_sequence,                                                             # Keep the original ancestor for audit.
                    mutations=mutation_list(seed_sequence, next_seq),                                                   # Full mutation list vs seed.
                    mutated_positions=sorted({*parent_mutated, pos}),                                                   # Union with parent's positions.
                    proposal_trace=list(parent.proposal_trace) + [f"{pos}:{row['seed_aa']}→{aa}"],                      # Audit trail.
                    round_index=int(parent.round_index) + 1,                                                            # Increment depth.
                )
                children.append(child)                                                                                  # Accept the child.
                emitted += 1                                                                                            # Count it.
                if emitted >= int(proposals_per_parent):                                                                # Per-parent cap reached.
                    break                                                                                               # Stop inner loop.
            if emitted >= int(proposals_per_parent):                                                                    # Per-parent cap reached (outer).
                break                                                                                                   # Stop position loop.
    return children                                                                                                     # Done.


def main() -> None:                                                                                                     # Orchestrator.
    import sys                                                                                                          # Used for exit codes.
    args = parse_args()                                                                                                 # CLI.
    seed_everything(args.seed)                                                                                          # Deterministic RNG.
    start_time = time.time()                                                                                            # Wall-clock start.

    # ----- Load the Stage 11 context -----
    context_path = Path(args.stage11_context_json)                                                                      # Path object.
    if not context_path.exists():                                                                                       # Defensive check.
        print(f"[ERROR] Missing Stage 11 context JSON: {context_path}", file=sys.stderr)                                # Surface.
        sys.exit(EXIT_INPUT_ERROR)                                                                                      # Exit 1.
    context = read_json(context_path)                                                                                   # Load.
    if context.get("stage") != "11a":                                                                                   # Defensive: must be Stage 11a output.
        print(f"[ERROR] Expected stage=='11a' in context, got stage={context.get('stage')!r}", file=sys.stderr)         # Surface.
        sys.exit(EXIT_INPUT_ERROR)                                                                                      # Exit 1.
    if not context.get("baseline_qualification", {}).get("passed", False):                                              # Defensive: 11a should never write a failed-gate context, but check anyway.
        print("[ERROR] Stage 11 context has baseline_qualification.passed=false; cannot proceed.", file=sys.stderr)     # Surface.
        sys.exit(EXIT_INPUT_ERROR)                                                                                      # Exit 1.

    # ----- Resolve scaffold + seed sequence -----
    seed_pdb_path = Path(context["seed_pdb_path"])                                                                      # Path to the freshly-folded seed PDB.
    if not seed_pdb_path.exists():                                                                                      # Defensive existence check.
        print(f"[ERROR] Seed PDB referenced by context is missing: {seed_pdb_path}", file=sys.stderr)                   # Surface.
        sys.exit(EXIT_INPUT_ERROR)                                                                                      # Exit 1.
    seed_sequence = str(context["selected_seed"]["seed_sequence"])                                                      # Seed sequence.
    target_host = str(context["target_host"])                                                                           # Target host.
    source_host = str(context.get("source_host", ""))                                                                   # Source host (for output provenance).
    seed_protein_id = str(context["selected_seed"].get("seed_protein_id", ""))                                          # Seed ID.

    # ----- Load IF1 + ESM-2 backbones (once, reused across rounds) -----
    print(f"[INFO] Loading inverse-folding model on device={args.if_device}…", flush=True)                              # Heads-up.
    try:                                                                                                                # Wrap heavy loads.
        coords, native_sequence = load_inverse_folding_structure(seed_pdb_path, chain_id=args.if_chain_id)              # Parse PDB.
        _, if_model, if_alphabet = load_inverse_folding_model(device=args.if_device)                                    # Load IF1.
        torch_emb, tokenizer_emb, model_emb, emb_device = load_embedding_backend(args.embedding_model)                  # Load ESM-2.
    except Exception as exc:                                                                                            # Surface.
        print(f"[ERROR] Failed to load inverse-folding or embedding backbone: {exc}", file=sys.stderr)                  # Loud error.
        sys.exit(EXIT_INFERENCE_ERROR)                                                                                  # Exit 3.

    if len(native_sequence) != len(seed_sequence):                                                                      # Sanity check the PDB-vs-context length.
        print(                                                                                                          # Loud warning but not fatal — Stage 10 has the same compensation.
            f"[WARN] IF1 native_sequence length {len(native_sequence)} != context seed length {len(seed_sequence)}; "
            f"continuing with the context seed sequence as the reference.",
            flush=True,
        )

    # ----- Load trained host probe -----
    try:                                                                                                                # Wrap predictor load.
        predictor, label_classes = load_target_predictor(args.predictor_model, args.label_classes_json)                 # Logistic regression + label order.
    except Exception as exc:                                                                                            # Surface.
        print(f"[ERROR] Failed to load target predictor: {exc}", file=sys.stderr)                                       # Loud.
        sys.exit(EXIT_INPUT_ERROR)                                                                                      # Exit 1.
    if target_host not in label_classes:                                                                                # Defensive: probe must know about this host.
        print(f"[ERROR] target_host='{target_host}' is not in predictor label_classes={label_classes}", file=sys.stderr)# Loud.
        sys.exit(EXIT_INPUT_ERROR)                                                                                      # Exit 1.

    # ----- Family centroid (recompute from the context's family rows on each backend pass) -----
    family_rows = context.get("family_context", {}).get("family_rows", [])                                              # The full family rows are embedded in the context.
    family_sequences = [str(r["aa_sequence"]) for r in family_rows] if family_rows else []                              # Pull just the sequences.
    family_centroid_list = context.get("family_context", {}).get("family_centroid", []) or []                           # Stage 11a already computed and persisted the centroid.
    family_centroid = np.asarray(family_centroid_list, dtype=np.float32) if family_centroid_list else None              # Convert to numpy if present.
    if family_centroid is None and family_sequences:                                                                    # Recompute on-the-fly if 11a somehow didn't persist it.
        print("[WARN] family_centroid missing from context; recomputing on the fly.", flush=True)                       # Loud.
        from phageforge.stage11_utils import compute_family_centroid as _compute_family_centroid                        # Local import to keep top-of-file lean.
        family_centroid = _compute_family_centroid(                                                                     # Centroid via persistent backend.
            family_sequences=family_sequences,
            embedding_model=args.embedding_model,
            batch_size=args.batch_size,
            torch=torch_emb, tokenizer=tokenizer_emb, model=model_emb, device=emb_device,
        )

    # ----- Extract edit-space constraints -----
    editable_region = dict(context["editable_region"])                                                                  # Edit-space block.
    proposal_rows = list(editable_region.get("proposal_rows", []))                                                      # Proposal rulebook.
    min_mutations = int(editable_region["min_mutations"])                                                               # Lower budget.
    max_mutations = int(editable_region["max_mutations"])                                                               # Upper budget.

    # ----- Initialize beam with the wild-type seed -----
    beam: list[Stage10Candidate] = [                                                                                    # Round-0 beam.
        Stage10Candidate(
            candidate_sequence=seed_sequence,                                                                           # Wild-type sequence as the only round-0 parent.
            parent_sequence=seed_sequence,                                                                              # Self as parent.
            mutations=[],                                                                                               # No mutations yet.
            mutated_positions=[],                                                                                       # No positions yet.
            proposal_trace=[],                                                                                          # Empty audit trail.
            round_index=0,                                                                                              # Round 0.
        )
    ]
    all_round_frames: list[pd.DataFrame] = []                                                                           # Per-round dataframes accumulated here.

    # ----- Beam-search loop -----
    for round_idx in range(1, int(args.rounds) + 1):                                                                    # Walk every round.
        print(f"[INFO] Round {round_idx}/{int(args.rounds)} — expanding from {len(beam)} parents…", flush=True)         # Heads-up.
        children = build_round_children(                                                                                # Spawn children.
            parents=beam,
            proposal_rows=proposal_rows,
            seed_sequence=seed_sequence,
            max_mutations=max_mutations,
            proposals_per_parent=args.proposals_per_parent,
            substitutions_per_position=args.substitutions_per_position,
        )
        if not children:                                                                                                # No children — beam exhausted.
            print(f"[INFO] Round {round_idx} produced no new children; terminating early.", flush=True)                 # Loud.
            break                                                                                                       # Exit the round loop.

        # Deduplicate by exact sequence (keep the first instance).
        dedup: dict[str, Stage10Candidate] = {}                                                                         # Dict to deduplicate.
        for item in children:                                                                                           # Walk children.
            dedup.setdefault(item.candidate_sequence, item)                                                              # First-wins dedup.
        children = list(dedup.values())                                                                                 # Replace list.

        # Score every unique child via the multi-modal evaluation gauntlet.
        try:                                                                                                            # Wrap heavy inference.
            score_frame = evaluate_candidate_table(                                                                     # Target prob + IF1 ll + family cos.
                sequences=[item.candidate_sequence for item in children],
                target_host=target_host,
                predictor_model_path=args.predictor_model,
                predictor_label_classes_path=args.label_classes_json,
                embedding_model=args.embedding_model,
                family_centroid=family_centroid,
                coords=coords,
                if_model=if_model,
                if_alphabet=if_alphabet,
                batch_size=args.batch_size,
                predictor=predictor,
                label_classes=label_classes,
                torch=torch_emb,
                tokenizer=tokenizer_emb,
                model=model_emb,
                device=emb_device,
            )
        except RuntimeError as exc:                                                                                     # OOM or other runtime.
            print(f"[ERROR] Round {round_idx} scoring failed: {exc}", file=sys.stderr)                                  # Surface.
            sys.exit(EXIT_INFERENCE_ERROR)                                                                              # Exit 3.

        # Attach trajectory metadata + seed identity.
        meta_frame = pd.DataFrame(                                                                                      # Metadata frame.
            {
                "candidate_sequence": [item.candidate_sequence for item in children],                                   # Primary key.
                "mutation_count": [len(item.mutated_positions) for item in children],                                   # Edit count.
                "mutated_positions": [";".join(map(str, item.mutated_positions)) for item in children],                 # Stage-10 column name.
                "mutation_positions": [";".join(item.mutations) for item in children],                                  # 08a-compatible column name (mutation strings).
                "mutation_text": [";".join(item.mutations) for item in children],                                       # Stage-10 column name.
                "proposal_trace": [";".join(item.proposal_trace) for item in children],                                 # Audit trail.
                "round_index": [int(item.round_index) for item in children],                                            # Round depth.
                "seed_identity": [                                                                                      # Fractional identity vs wild-type.
                    sum(a == b for a, b in zip(seed_sequence, item.candidate_sequence))
                    / max(len(seed_sequence), 1)
                    for item in children
                ],
            }
        )
        round_frame = score_frame.merge(meta_frame, on="candidate_sequence", how="inner")                               # Join scores with metadata.

        # ----- Composite score (locked weights, with CLI overrides for ablation) -----
        target_norm = round_frame["target_probability"].to_numpy(dtype=np.float32)                                      # Raw target prob.
        if1_norm = round_frame["if1_log_likelihood"].to_numpy(dtype=np.float32)                                         # Raw IF1.
        family_norm = round_frame["family_cosine"].to_numpy(dtype=np.float32)                                           # Raw family cosine.
        identity_norm = round_frame["seed_identity"].to_numpy(dtype=np.float32)                                         # Identity.
        mut_count = round_frame["mutation_count"].to_numpy(dtype=np.float32)                                            # Mutation count.

        # If the operator passed Stage-10 defaults, use the canonical helper exactly.
        if (math.isclose(args.w_target, 0.30) and math.isclose(args.w_if1, 0.45) and                                    # Detect default weights.
                math.isclose(args.w_family, 0.15) and math.isclose(args.w_identity, 0.10) and
                math.isclose(args.w_mut_penalty, 0.10)):
            round_frame["stage10_composite_score"] = composite_stage10_score(                                           # Use the canonical Stage-10 calculation.
                target_probability=target_norm, if1_log_likelihood=if1_norm,
                family_cosine=family_norm, seed_identity=identity_norm,
                mutation_count=mut_count,
            )
        else:                                                                                                           # Operator is ablating — apply weights manually.
            from phageforge.stage11_utils import robust_minmax                                                          # Local import.
            tn = robust_minmax(target_norm)                                                                             # Normalize each input.
            i1 = robust_minmax(if1_norm)                                                                                # Same.
            fc = robust_minmax(family_norm)                                                                             # Same.
            sid = robust_minmax(identity_norm)                                                                          # Same.
            mp = robust_minmax(mut_count)                                                                               # Same.
            round_frame["stage10_composite_score"] = (                                                                  # Custom weights.
                float(args.w_target) * tn
                + float(args.w_if1) * i1
                + float(args.w_family) * fc
                + float(args.w_identity) * sid
                - float(args.w_mut_penalty) * mp
            )

        # ----- 08a-compatible aliases on the same column space -----
        round_frame["stage11_composite_score"] = round_frame["stage10_composite_score"]                                 # Stage-11 alias of the composite.
        round_frame["final_multimodal_rank_score"] = round_frame["stage10_composite_score"]                             # 08a sorts by this column.
        round_frame["generation_regime"] = "stage11_inverse_folding"                                                    # 08a's expected regime column.
        round_frame["seed_pdb_path"] = str(seed_pdb_path.resolve())                                                     # Anchor traceability.
        round_frame["target_host"] = target_host                                                                        # Target host column.
        round_frame["source_host"] = source_host                                                                        # Source host column.
        round_frame["seed_protein_id"] = seed_protein_id                                                                # Seed ID column.

        # Append to the historical log (sorted by composite descending).
        all_round_frames.append(                                                                                        # Save this round.
            round_frame.sort_values("stage10_composite_score", ascending=False).reset_index(drop=True)
        )

        # Embed and pick diverse top-`beam_width` survivors for the next round.
        cand_emb = embed_sequences_with_backend(                                                                        # Embed every scored candidate.
            round_frame["candidate_sequence"].astype(str).tolist(),
            torch_emb, tokenizer_emb, model_emb,
            device=emb_device,
            batch_size=args.batch_size,
        )
        keep_idx = greedy_diverse_subset(                                                                               # Diversity reranking.
            cand_emb,
            round_frame["stage10_composite_score"].to_numpy(dtype=np.float32),
            top_k=args.beam_width,
        )
        keep_sequences = set(round_frame.iloc[keep_idx]["candidate_sequence"].astype(str).tolist())                     # Survivors as a fast lookup set.

        next_beam: list[Stage10Candidate] = []                                                                          # Build the next-round beam.
        for item in children:                                                                                           # Walk the children of this round.
            if item.candidate_sequence in keep_sequences:                                                               # Surviving children become parents.
                next_beam.append(item)
        beam = next_beam                                                                                                # Update the active beam.
        if not beam:                                                                                                    # No survivors — stop.
            print(f"[INFO] No survivors after round {round_idx}; terminating early.", flush=True)                       # Loud.
            break                                                                                                       # Exit round loop.

    # ----- Concatenate all round outputs and persist -----
    if all_round_frames:                                                                                                # We have at least one round of data.
        search_df = pd.concat(all_round_frames, ignore_index=True)                                                      # Stack vertically.
        search_df = search_df.drop_duplicates(subset=["candidate_sequence"], keep="first").reset_index(drop=True)       # Dedup by sequence (keep best-scored instance).
        # Inject a 1-based sample_id (08a wants this).
        search_df = search_df.sort_values(                                                                              # Sort by composite descending so sample_id 1 is the best.
            ["stage10_composite_score", "target_probability", "if1_log_likelihood"],
            ascending=False,
        ).reset_index(drop=True)
        search_df.insert(0, "sample_id", np.arange(1, len(search_df) + 1))                                              # Inject 1-based ID at front.
    else:                                                                                                               # Defensive: empty pipeline (shouldn't happen).
        search_df = pd.DataFrame(columns=[                                                                              # Empty but correctly-shaped frame.
            "sample_id", "candidate_sequence", "target_probability", "if1_log_likelihood",
            "family_cosine", "mutation_count", "mutated_positions", "mutation_positions",
            "mutation_text", "proposal_trace", "round_index", "seed_identity",
            "stage10_composite_score", "stage11_composite_score", "final_multimodal_rank_score",
            "generation_regime", "seed_pdb_path", "target_host", "source_host", "seed_protein_id",
        ])

    out_csv = Path(args.out_csv)                                                                                        # Destination CSV.
    out_csv.parent.mkdir(parents=True, exist_ok=True)                                                                   # Create directory.
    search_df.to_csv(out_csv, index=False)                                                                              # Write.
    print(f"[OK] Wrote: {out_csv}", flush=True)                                                                         # Confirm.

    # ----- Summary JSON -----
    elapsed = float(time.time() - start_time)                                                                           # Wall-clock cost.
    gpu_name = ""                                                                                                       # Best-effort GPU name.
    try:                                                                                                                # Wrap optional torch import.
        import torch as _torch                                                                                          # Local import.
        if _torch.cuda.is_available():                                                                                  # Only if GPU is present.
            gpu_name = str(_torch.cuda.get_device_name(0))                                                              # First device.
    except Exception:                                                                                                   # Silently skip if torch missing.
        pass                                                                                                            # Continue.

    summary = {                                                                                                         # Compact run manifest.
        "stage": "11b",                                                                                                 # Origin tag.
        "stage11_context_json": str(context_path.resolve()),                                                            # Provenance.
        "search_rows": int(len(search_df)),                                                                             # Total unique candidates evaluated.
        "rounds_completed": int(search_df["round_index"].max()) if len(search_df) else 0,                               # Deepest round reached.
        "top_candidate_sequence": (str(search_df.iloc[0]["candidate_sequence"]) if len(search_df) else ""),             # #1 candidate.
        "top_stage10_score": (float(search_df.iloc[0]["stage10_composite_score"]) if len(search_df) else float("nan")), # Best composite.
        "top_target_probability": (                                                                                      # Best target prob.
            float(search_df.iloc[0]["target_probability"]) if len(search_df) else float("nan")
        ),
        "top_if1_log_likelihood": (                                                                                     # Best IF1.
            float(search_df.iloc[0]["if1_log_likelihood"]) if len(search_df) else float("nan")
        ),
        "elapsed_seconds": elapsed,                                                                                     # Wall-clock cost.
        "gpu_name": gpu_name,                                                                                           # GPU info.
        "composite_weights": {                                                                                          # Record weights used.
            "w_target": float(args.w_target),
            "w_if1": float(args.w_if1),
            "w_family": float(args.w_family),
            "w_identity": float(args.w_identity),
            "w_mut_penalty": float(args.w_mut_penalty),
        },
        "out_csv": str(out_csv.resolve()),                                                                              # Output path.
    }
    write_json(summary, args.out_json)                                                                                  # Persist.
    print(f"[OK] Wrote: {args.out_json}", flush=True)                                                                   # Confirm.
    sys.exit(EXIT_OK)                                                                                                   # Explicit success.


if __name__ == "__main__":                                                                                              # CLI guard.
    main()                                                                                                              # Run.
