"""Stage 07b: Generate family-constrained RBP candidates with local ESM3 or Forge-backed ESM3."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from phageforge.stage07_utils import (
    Stage07Regime,
    candidate_guidance_score,
    choose_hotspots,
    make_masked_prompt,
    mutation_list,
    mutation_penalty,
    parse_regimes_json,
    read_json,
    seed_everything,
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for Stage 07 ESM3 generation."""
    ap = argparse.ArgumentParser(description="Generate Stage 07 RBP candidates with ESM3.")
    ap.add_argument("--context_json", type=str, required=True, help="Stage 07 context JSON produced by 07a_prepare_stage07_design_context.py.")
    ap.add_argument("--out_csv", type=str, required=True, help="Where to write the generated-candidate CSV.")
    ap.add_argument("--n_samples", type=int, default=12, help="Total number of candidate rows to request across all regimes.")
    ap.add_argument("--temperature", type=float, default=0.7, help="Default generation temperature when explicit regimes are not provided.")
    ap.add_argument("--top_k", type=int, default=5, help="Default top-k metadata value when explicit regimes are not provided.")
    ap.add_argument("--max_esm3_masked_positions", type=int, default=24, help="Maximum number of masked sequence positions per prompt.")
    ap.add_argument("--sampling_seed", type=int, default=42, help="Base random seed for hotspot choice and generation reproducibility.")
    ap.add_argument("--esm3_backend", type=str, default="local", choices=["local", "forge"], help="Run local esm3-open weights or call the Forge API.")
    ap.add_argument("--esm3_model", type=str, default="esm3-open", help="Model label for local or forge generation.")
    ap.add_argument("--esm3_num_steps", type=int, default=8, help="Number of ESM3 iterative unmasking steps.")
    ap.add_argument("--esm3_error_fallback", type=str, default="none", choices=["none", "esm2"], help="Optional fallback; 'none' keeps the run pure ESM3.")
    ap.add_argument("--regimes_json", type=str, default=None, help="Optional JSON list describing multiple generation regimes.")
    ap.add_argument("--hotspot_strategy", type=str, default="mixed", choices=["even", "priority", "mixed"], help="Fallback hotspot selection strategy when regimes_json is omitted.")
    ap.add_argument("--max_attempts_per_sample", type=int, default=3, help="How many guided ESM3 attempts to try before keeping the best candidate.")
    return ap.parse_args()


def load_esm3_client(backend: str, model_name: str):
    """Instantiate a local or Forge ESM3 client using the public ESM python interface."""
    import esm
    import torch
    from esm.models.esm3 import ESM3

    if backend == "local":
        client = ESM3.from_pretrained(model_name)
        if hasattr(client, "to"):
            client = client.to("cuda" if torch.cuda.is_available() else "cpu")
        return client

    token = os.environ.get("ESM_API_KEY") or os.environ.get("FORGE_API_TOKEN") or os.environ.get("ESM3_FORGE_TOKEN")
    if not token:
        raise RuntimeError("Forge generation was requested but no Forge token was found in ESM_API_KEY / FORGE_API_TOKEN / ESM3_FORGE_TOKEN.")
    return esm.sdk.client(model_name, token=token)


def generate_sequence(client, prompt_sequence: str, num_steps: int, temperature: float, top_k: int, sample_seed: int) -> str:
    """Generate a sequence completion from an ESM3 prompt string."""
    import torch
    from esm.sdk.api import ESMProtein, GenerationConfig

    seed_everything(sample_seed)
    torch.manual_seed(sample_seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(sample_seed)

    protein = ESMProtein(sequence=prompt_sequence)
    cfg = GenerationConfig(track="sequence", num_steps=num_steps, temperature=float(temperature))
    if hasattr(cfg, "top_k"):
        cfg.top_k = int(top_k)
    elif hasattr(cfg, "topk"):
        cfg.topk = int(top_k)
    out = client.generate(protein, cfg)
    sequence = getattr(out, "sequence", None)
    if not sequence:
        raise RuntimeError("ESM3 returned no sequence output.")
    return str(sequence)


def _attempt_generation(client, seed_sequence: str, prompt_sequence: str, masked_positions: list[int], regime: Stage07Regime, args: argparse.Namespace, position_features: list[dict], sample_seed: int) -> tuple[str, str | None, int, float]:
    """Run a few guided ESM3 attempts and keep the best successful candidate."""
    best_sequence = ""
    best_error = None
    best_score = float("-inf")
    best_seed = sample_seed

    for attempt_idx in range(max(1, int(args.max_attempts_per_sample))):
        attempt_seed = int(sample_seed + attempt_idx)
        try:
            candidate_sequence = generate_sequence(
                client=client,
                prompt_sequence=prompt_sequence,
                num_steps=regime.num_steps,
                temperature=regime.temperature,
                top_k=regime.top_k,
                sample_seed=attempt_seed,
            )
            mutated_positions = [pos for pos in masked_positions if pos <= len(candidate_sequence) and candidate_sequence[pos - 1] != seed_sequence[pos - 1]]
            guidance_score = candidate_guidance_score(seed_sequence, candidate_sequence, position_features, mutated_positions)
            novelty_bonus = 0.005 * len(mutated_positions)
            total_score = guidance_score + novelty_bonus
            if total_score > best_score:
                best_sequence = candidate_sequence
                best_error = None
                best_score = total_score
                best_seed = attempt_seed
        except Exception as exc:
            best_error = str(exc)

    return best_sequence, best_error, best_seed, best_score


def build_rows(seed_sequence: str, target_host: str, family_label: str, regimes: list[Stage07Regime], context: dict, client, args: argparse.Namespace) -> list[dict]:
    """Generate candidate rows across regimes, deduplicate sequences, and retain provenance."""
    hotspot_positions = context["editable_region"].get("hotspot_positions", [])
    hotspot_weights = context["editable_region"].get("hotspot_priority_weights", {})
    position_features = context["editable_region"].get("position_features", [])
    per_regime = max(1, (args.n_samples + len(regimes) - 1) // len(regimes))

    rows = []
    seen_sequences = set()
    sample_id = 0
    for regime_rank, regime in enumerate(regimes):
        for local_idx in tqdm(range(per_regime), total=per_regime):
            if len(rows) >= args.n_samples:
                break
            base_seed = int(args.sampling_seed + 1000 * regime_rank + 100 * local_idx)
            masked_positions = choose_hotspots(
                hotspots_1based=hotspot_positions,
                priority_weights=hotspot_weights,
                max_positions=regime.max_masked_positions,
                strategy=regime.hotspot_strategy,
                sample_seed=base_seed,
            )
            prompt_sequence = make_masked_prompt(seed_sequence, masked_positions)
            try:
                candidate_sequence, best_error, chosen_seed, guidance_score = _attempt_generation(
                    client=client,
                    seed_sequence=seed_sequence,
                    prompt_sequence=prompt_sequence,
                    masked_positions=masked_positions,
                    regime=regime,
                    args=args,
                    position_features=position_features,
                    sample_seed=base_seed,
                )
                if not candidate_sequence:
                    raise RuntimeError(best_error or "ESM3 generation did not return a candidate sequence.")
                status = "ok"
                error_text = ""
                generator_mode = f"esm3_{args.esm3_backend}:{args.esm3_model}"
            except Exception as exc:
                if args.esm3_error_fallback != "none":
                    raise RuntimeError("ESM2 fallback is intentionally not implemented in the upgraded pure-ESM3 Stage 07 path.") from exc
                candidate_sequence = ""
                chosen_seed = base_seed
                guidance_score = float("-inf")
                status = "error"
                error_text = str(exc)
                generator_mode = f"esm3_{args.esm3_backend}:{args.esm3_model}"

            if status == "ok" and candidate_sequence in seen_sequences:
                continue
            if status == "ok":
                seen_sequences.add(candidate_sequence)

            row = {
                "sample_id": sample_id,
                "generator_mode": generator_mode,
                "generation_regime": regime.name,
                "generation_regime_rank": regime_rank,
                "generation_status": status,
                "generation_error": error_text or None,
                "target_host": target_host,
                "family_id": family_label,
                "seed_sequence": seed_sequence,
                "candidate_sequence": candidate_sequence if candidate_sequence else seed_sequence,
                "editable_hotspots": ",".join(str(pos) for pos in masked_positions),
                "editable_hotspot_count": int(len(masked_positions)),
                "mutation_positions": ";".join(mutation_list(seed_sequence, candidate_sequence if candidate_sequence else seed_sequence)),
                "mutation_penalty": int(mutation_penalty(seed_sequence, candidate_sequence if candidate_sequence else seed_sequence)),
                "esm3_prompt_sequence": prompt_sequence,
                "esm3_hotspot_positions": ",".join(str(pos) for pos in masked_positions),
                "esm3_temperature": float(regime.temperature),
                "esm3_top_k": int(regime.top_k),
                "esm3_sampling_seed": int(chosen_seed),
                "esm3_model": args.esm3_model,
                "esm3_num_steps": int(regime.num_steps),
                "hotspot_strategy": regime.hotspot_strategy,
                "guided_mutation_score": float(guidance_score) if np.isfinite(guidance_score) else None,
                "used_esm3_api": bool(args.esm3_backend == "forge"),
                "used_local_generator": bool(args.esm3_backend == "local"),
                "used_local_esm3": bool(args.esm3_backend == "local"),
                "used_esm2_fallback": False,
            }
            rows.append(row)
            sample_id += 1
    return rows


def main() -> None:
    # Read the Stage 07 context, build generation regimes, and load the selected ESM3 client.
    args = parse_args()
    context = read_json(args.context_json)
    regimes = parse_regimes_json(args.regimes_json, args.temperature, args.top_k, args.max_esm3_masked_positions, args.esm3_num_steps)
    if args.regimes_json is None:
        for regime in regimes:
            regime.hotspot_strategy = args.hotspot_strategy if regime.name == "balanced" else regime.hotspot_strategy
    client = load_esm3_client(args.esm3_backend, args.esm3_model)

    seed_sequence = str(context["selected_seed"]["seed_sequence"])
    target_host = str(context["target_host"])
    family_label = str(context["family_context"].get("family_product_majority", "receptor-binding protein"))
    rows = build_rows(seed_sequence, target_host, family_label, regimes, context, client, args)
    if not rows:
        raise RuntimeError("ESM3 generation did not produce any rows.")

    # Write the generated table in a single place so later Stage 07 steps can read it directly.
    out_path = Path(args.out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).sort_values(["generation_status", "guided_mutation_score", "mutation_penalty", "sample_id"], ascending=[True, False, True, True]).to_csv(out_path, index=False)
    print(f"Wrote: {out_path}")
    print(f"used_backend: {args.esm3_backend}")
    print(f"rows_written: {len(rows)}")


if __name__ == "__main__":
    main()
