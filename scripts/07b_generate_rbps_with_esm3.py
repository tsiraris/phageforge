"""
==========================================================================================
Stage 07b: ESM3-guided masked editing over the validated Stage 07 window (clean version).
==========================================================================================

Hotfix:
- fixes the local/backend sampling_seed NameError by threading sampling_seed through
  generate_with_esm3_backend() and generate_with_esm2_masked_lm()
"""

from __future__ import annotations
import argparse
import json
import os
import random
from pathlib import Path
from typing import Any

import pandas as pd
import torch
from transformers import AutoTokenizer, EsmForMaskedLM

from phageforge.generation.dataset import VOCAB
from phageforge.generation.generator_model import ConditionalMaskedRBPGenerator
from phageforge.generation.guidance import (
    decode_ids,
    encode_seed,
    hotspot_masked_ids,
    mutation_burden,
    sample_hotspot_edits,
)


def parse_args() -> argparse.Namespace:
    """Return command-line arguments for Stage 07 candidate generation."""
    ap = argparse.ArgumentParser(description="Generate Stage 07 candidate RBP edits cleanly and reproducibly.")
    ap.add_argument("--context_json", type=str, required=True, help="Context JSON from 07a_prepare_stage07_design_context.py.")
    ap.add_argument("--out_csv", type=str, required=True, help="Where to write generated candidates.")
    ap.add_argument("--n_samples", type=int, default=8, help="How many candidates to generate in this run.")
    ap.add_argument("--temperature", type=float, default=0.7, help="Sampling temperature used at editable positions.")
    ap.add_argument("--top_k", type=int, default=5, help="Top-k token sampling at editable positions when supported.")
    ap.add_argument("--sampling_seed", type=int, default=42, help="Random seed for reproducibility.")

    ap.add_argument("--esm3_backend", type=str, default="forge", choices=["forge", "local", "auto"], help="How to run ESM3: Forge API, local open weights, or auto.")
    ap.add_argument("--use_esm3_api", action="store_true", help="Backward-compatible alias for --esm3_backend forge.")
    ap.add_argument("--esm3_model", type=str, default="esm3-open-2024-03", help="Forge ESM3 model name.")
    ap.add_argument("--esm3_num_steps", type=int, default=8, help="Number of ESM3 iterative unmasking steps.")
    ap.add_argument("--esm_api_key", type=str, default="", help="Optional Forge API key. Falls back to ESM_API_KEY env var if omitted.")
    ap.add_argument("--max_esm3_masked_positions", type=int, default=24, help="Maximum editable positions sent to ESM3 in one prompt.")
    ap.add_argument("--allow_full_window_if_no_hotspots", action="store_true", help="Allow full editable window if no hotspots exist.")

    ap.add_argument("--esm3_error_fallback", type=str, default="none", choices=["none", "esm2", "local_generator"], help="Optional explicit fallback mode.")
    ap.add_argument("--esm2_model", type=str, default="facebook/esm2_t33_650M_UR50D", help="Masked-LM fallback model.")
    ap.add_argument("--generator_checkpoint", type=str, default="", help="Local conditional generator checkpoint.")
    ap.add_argument("--stop_on_quota_exhausted", action="store_true", help="Stop early and write partial output when Forge credits are exhausted.")
    ap.add_argument("--resume_from_existing", action="store_true", help="Append to an existing CSV and continue from the next sample_id.")
    return ap.parse_args()


def set_generation_seed(seed: int) -> None:
    """Seed Python and PyTorch RNGs for reproducible generation behavior."""
    random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def load_context(path: Path) -> dict:
    """Load the Stage 07 context JSON from disk."""
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def resolve_hotspots(context: dict, seed_seq: str) -> list[int]:
    """Return explicit editable residue positions from hotspot list or editable window."""
    hotspot_positions = list(context.get("editable_region", {}).get("hotspot_positions", []))
    if hotspot_positions:
        return [int(x) for x in hotspot_positions]
    window_start = int(context.get("editable_region", {}).get("window_start", 0))
    window_end = int(context.get("editable_region", {}).get("window_end", len(seed_seq)))
    window_end = min(window_end, len(seed_seq))
    return list(range(max(0, window_start), max(0, window_end)))


def cap_positions_evenly(positions: list[int], max_positions: int) -> list[int]:
    """Evenly downsample editable positions to a maximum count while preserving span."""
    if max_positions <= 0 or len(positions) <= max_positions:
        return list(positions)
    if max_positions == 1:
        return [positions[len(positions) // 2]]
    step = (len(positions) - 1) / float(max_positions - 1)
    chosen = sorted({positions[round(i * step)] for i in range(max_positions)})
    return chosen


def build_esm3_masked_prompt(seed_seq: str, hotspot_positions: list[int]) -> str:
    """Return an ESM3 prompt where editable positions are replaced by underscores."""
    chars = list(seed_seq)
    for pos in hotspot_positions:
        if 0 <= pos < len(chars):
            chars[pos] = "_"
    return "".join(chars)


def describe_mutations(seed_seq: str, candidate_seq: str) -> str:
    """Return a concise semicolon-separated description of residue changes."""
    changes = []
    for i, (a, b) in enumerate(zip(seed_seq, candidate_seq)):
        if a != b:
            changes.append(f"{i}:{a}→{b}")
    return ";".join(changes)


def build_generation_record(
    *,
    sample_id: int,
    generator_mode: str,
    generation_status: str,
    generation_error: str | None,
    target_host: str,
    family_id: str,
    seed_sequence: str,
    candidate_sequence: str,
    hotspot_positions: list[int],
    masked_prompt: str,
    temperature: float,
    top_k: int,
    sampling_seed: int,
    esm3_model: str,
    esm3_num_steps: int,
) -> dict[str, Any]:
    """Build one standardized provenance-rich output row."""
    return {
        "sample_id": int(sample_id),
        "generator_mode": generator_mode,
        "generation_status": generation_status,
        "generation_error": generation_error,
        "target_host": target_host,
        "family_id": family_id,
        "seed_sequence": seed_sequence,
        "candidate_sequence": candidate_sequence,
        "editable_hotspots": ",".join(str(x) for x in hotspot_positions),
        "editable_hotspot_count": int(len(hotspot_positions)),
        "mutation_positions": describe_mutations(seed_sequence, candidate_sequence),
        "mutation_penalty": mutation_burden(seed_sequence, candidate_sequence),
        "esm3_prompt_sequence": masked_prompt,
        "esm3_hotspot_positions": ",".join(str(x) for x in hotspot_positions),
        "esm3_temperature": float(temperature),
        "esm3_top_k": int(top_k),
        "esm3_sampling_seed": int(sampling_seed),
        "esm3_model": esm3_model,
        "esm3_num_steps": int(esm3_num_steps),
    }


def parse_credit_or_access_error(message: str) -> str:
    """Classify common Forge error patterns for cleaner control flow."""
    text = str(message).lower()
    if "daily credit limit" in text or "429" in text or "quota" in text:
        return "quota_exhausted"
    if "does not have access to model" in text or "access to model" in text:
        return "model_access_denied"
    return "other"


def maybe_load_existing_rows(out_csv: Path, resume: bool) -> list[dict[str, Any]]:
    """Load existing output rows when resume mode is enabled."""
    if resume and out_csv.exists():
        return pd.read_csv(out_csv).to_dict(orient="records")
    return []


def try_import_local_esm3():
    """Import local ESM3-open classes only when needed."""
    from esm.models.esm3 import ESM3
    from esm.sdk.api import ESMProtein, GenerationConfig
    return ESM3, ESMProtein, GenerationConfig


def try_import_forge_esm3():
    """Import Forge client classes only when needed."""
    from esm.sdk import client as esm_client
    from esm.sdk.api import ESMProtein, GenerationConfig
    return esm_client, ESMProtein, GenerationConfig


def generate_with_esm3_backend(
    *,
    context: dict,
    backend: str,
    model_name: str,
    esm_api_key: str,
    n_samples: int,
    temperature: float,
    top_k: int,
    num_steps: int,
    max_masked_positions: int,
    allow_full_window_if_no_hotspots: bool,
    stop_on_quota_exhausted: bool,
    start_sample_id: int,
    sampling_seed: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Generate candidates with either Forge ESM3-open or local ESM3-open."""
    seed_seq = context["selected_seed"]["seed_sequence"]
    target_host = context["target_host"]
    family_id = context.get("family_context", {}).get("family_product_majority") or "unknown_family"
    hotspot_positions = resolve_hotspots(context, seed_seq)

    if len(hotspot_positions) == 0 and not allow_full_window_if_no_hotspots:
        raise ValueError("No hotspot positions were found in context, and full-window editing is disabled.")
    hotspot_positions = cap_positions_evenly(hotspot_positions, max_masked_positions)
    masked_prompt = build_esm3_masked_prompt(seed_seq, hotspot_positions)

    if backend == "local":
        ESM3, ESMProtein, GenerationConfig = try_import_local_esm3()
        device = "cuda" if torch.cuda.is_available() else "cpu"
        model = ESM3.from_pretrained("esm3-open").to(device)
        backend_label = "esm3_local:esm3-open"
        effective_model_name = "esm3-open"
    elif backend == "forge":
        esm_client, ESMProtein, GenerationConfig = try_import_forge_esm3()
        if not esm_api_key:
            raise ValueError("ESM API key is required for Forge backend. Set --esm_api_key or ESM_API_KEY.")
        model = esm_client(model_name, token=esm_api_key)
        backend_label = f"esm3_api:{model_name}"
        effective_model_name = model_name
    else:
        raise ValueError(f"Unsupported ESM3 backend: {backend}")

    rows = []
    stats = {"quota_exhausted": False, "backend": backend, "masked_positions": len(hotspot_positions)}

    for offset in range(n_samples):
        sample_id = start_sample_id + offset
        protein = ESMProtein(sequence=masked_prompt)
        try:
            generated = model.generate(
                protein,
                GenerationConfig(track="sequence", num_steps=num_steps, temperature=temperature),
            )
            candidate_seq = str(generated.sequence).replace(" ", "")
            rows.append(build_generation_record(
                sample_id=sample_id,
                generator_mode=backend_label,
                generation_status="ok",
                generation_error=None,
                target_host=target_host,
                family_id=family_id,
                seed_sequence=seed_seq,
                candidate_sequence=candidate_seq,
                hotspot_positions=hotspot_positions,
                masked_prompt=masked_prompt,
                temperature=temperature,
                top_k=top_k,
                sampling_seed=sampling_seed,
                esm3_model=effective_model_name,
                esm3_num_steps=num_steps,
            ))
        except Exception as exc:
            err = str(exc)
            kind = parse_credit_or_access_error(err)
            if kind == "quota_exhausted" and stop_on_quota_exhausted:
                stats["quota_exhausted"] = True
                stats["stop_reason"] = err
                break
            rows.append(build_generation_record(
                sample_id=sample_id,
                generator_mode=backend_label,
                generation_status="error",
                generation_error=err,
                target_host=target_host,
                family_id=family_id,
                seed_sequence=seed_seq,
                candidate_sequence=seed_seq,
                hotspot_positions=hotspot_positions,
                masked_prompt=masked_prompt,
                temperature=temperature,
                top_k=top_k,
                sampling_seed=sampling_seed,
                esm3_model=effective_model_name,
                esm3_num_steps=num_steps,
            ))
    return pd.DataFrame(rows), stats


def generate_with_esm2_masked_lm(
    context: dict,
    model_name: str,
    n_samples: int,
    temperature: float,
    top_k: int,
    start_sample_id: int,
    sampling_seed: int,
) -> pd.DataFrame:
    """Generate candidates with an ESM2 masked language model when explicitly requested."""
    device = "cuda" if torch.cuda.is_available() else "cpu"
    tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)
    model = EsmForMaskedLM.from_pretrained(model_name).to(device).eval()

    seed_seq = context["selected_seed"]["seed_sequence"]
    family_id = context.get("family_context", {}).get("family_product_majority") or "unknown_family"
    hotspot_positions = resolve_hotspots(context, seed_seq)
    chars = list(seed_seq)
    for pos in hotspot_positions:
        if 0 <= pos < len(chars):
            chars[pos] = tokenizer.mask_token
    masked_sequence = "".join(chars)
    rows = []

    for offset in range(n_samples):
        sample_id = start_sample_id + offset
        toks = tokenizer(masked_sequence, return_tensors="pt").to(device)
        with torch.no_grad():
            logits = model(**toks).logits[0]
        sampled_ids = toks["input_ids"][0].clone()
        mask_positions = (toks["input_ids"][0] == tokenizer.mask_token_id).nonzero(as_tuple=False).view(-1)
        for mask_pos in mask_positions.tolist():
            row = logits[mask_pos] / max(float(temperature), 1e-6)
            probs = torch.softmax(row, dim=-1)
            if 0 < top_k < probs.numel():
                top_vals, top_idx = torch.topk(probs, k=top_k)
                top_vals = top_vals / top_vals.sum()
                chosen = top_idx[torch.multinomial(top_vals, num_samples=1)[0]].item()
            else:
                chosen = torch.multinomial(probs, num_samples=1)[0].item()
            sampled_ids[mask_pos] = int(chosen)
        candidate_seq = tokenizer.decode(sampled_ids, skip_special_tokens=True).replace(" ", "")
        rows.append(build_generation_record(
            sample_id=sample_id,
            generator_mode=f"esm2_masked_lm:{model_name}",
            generation_status="ok",
            generation_error=None,
            target_host=context["target_host"],
            family_id=family_id,
            seed_sequence=seed_seq,
            candidate_sequence=candidate_seq,
            hotspot_positions=hotspot_positions,
            masked_prompt=masked_sequence,
            temperature=temperature,
            top_k=top_k,
            sampling_seed=sampling_seed,
            esm3_model=model_name,
            esm3_num_steps=0,
        ))
    return pd.DataFrame(rows)


def main() -> None:
    """Run Stage 07 generation with strict, explicit ESM3 behavior."""
    args = parse_args()
    if args.use_esm3_api:
        args.esm3_backend = "forge"
    set_generation_seed(args.sampling_seed)

    context = load_context(Path(args.context_json))
    out_path = Path(args.out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    existing_rows = maybe_load_existing_rows(out_path, args.resume_from_existing)
    start_sample_id = len(existing_rows)
    remaining = max(0, args.n_samples - start_sample_id)
    all_rows = list(existing_rows)

    esm_api_key = args.esm_api_key or os.environ.get("ESM_API_KEY", "")
    backend_order = [args.esm3_backend] if args.esm3_backend != "auto" else ["local", "forge"]

    generated_df = pd.DataFrame()
    used_backend = None
    last_error = None
    quota_exhausted = False

    for backend in backend_order:
        try:
            generated_df, stats = generate_with_esm3_backend(
                context=context,
                backend=backend,
                model_name=args.esm3_model,
                esm_api_key=esm_api_key,
                n_samples=remaining,
                temperature=args.temperature,
                top_k=args.top_k,
                num_steps=args.esm3_num_steps,
                max_masked_positions=args.max_esm3_masked_positions,
                allow_full_window_if_no_hotspots=args.allow_full_window_if_no_hotspots,
                stop_on_quota_exhausted=args.stop_on_quota_exhausted,
                start_sample_id=start_sample_id,
                sampling_seed=args.sampling_seed,
            )
            used_backend = backend
            quota_exhausted = bool(stats.get("quota_exhausted", False))
            break
        except Exception as exc:
            last_error = str(exc)
            continue

    if used_backend is None:
        if args.esm3_error_fallback == "none":
            raise RuntimeError(f"ESM3 generation did not start successfully. Last error: {last_error}")
        if args.esm3_error_fallback == "esm2":
            generated_df = generate_with_esm2_masked_lm(
                context=context,
                model_name=args.esm2_model,
                n_samples=remaining,
                temperature=args.temperature,
                top_k=args.top_k,
                start_sample_id=start_sample_id,
                sampling_seed=args.sampling_seed,
            )
        else:
            raise ValueError("This hotfix file only supports fallback mode 'none' or 'esm2'.")

    all_rows.extend(generated_df.to_dict(orient="records"))
    final_df = pd.DataFrame(all_rows)
    final_df.to_csv(out_path, index=False)

    print(f"Wrote: {out_path}")
    print(f"used_backend: {used_backend}")
    print(f"rows_written: {len(final_df)}")
    if quota_exhausted:
        print("quota_exhausted: True")
        print("tip: rerun later with --resume_from_existing to continue the same CSV.")


if __name__ == "__main__":
    main()
