"""Shared helpers for the Stage 07 local-ESM3 design workflow."""

from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import torch


AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class Stage07Regime:
    """Compact description of one generation regime for Stage 07."""

    name: str
    temperature: float
    top_k: int
    max_masked_positions: int
    hotspot_strategy: str = "mixed"
    num_steps: int = 8


def seed_everything(seed: int) -> None:
    """Set deterministic random seeds for python, numpy, and torch."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def read_json(path: str | Path) -> dict:
    """Read a JSON file and return the parsed object."""
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(obj: dict, path: str | Path) -> None:
    """Write a JSON object with stable indentation and UTF-8 encoding."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2)


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """Return cosine similarity between two vectors, guarding against zero norms."""
    a = np.asarray(a, dtype=np.float32)
    b = np.asarray(b, dtype=np.float32)
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom <= 0:
        return 0.0
    return float(np.dot(a, b) / denom)


def normalize_rows(x: np.ndarray) -> np.ndarray:
    """L2-normalize each row for cosine-based comparisons."""
    x = np.asarray(x, dtype=np.float32)
    denom = np.linalg.norm(x, axis=1, keepdims=True)
    return x / np.clip(denom, 1e-12, None)


def mutation_list(seed_sequence: str, candidate_sequence: str) -> list[str]:
    """Return human-readable residue substitutions using 1-based positions."""
    mutations = []
    for i, (seed_aa, cand_aa) in enumerate(zip(seed_sequence, candidate_sequence), start=1):
        if seed_aa != cand_aa:
            mutations.append(f"{i}:{seed_aa}→{cand_aa}")
    if len(candidate_sequence) > len(seed_sequence):
        for i, cand_aa in enumerate(candidate_sequence[len(seed_sequence) :], start=len(seed_sequence) + 1):
            mutations.append(f"{i}:∅→{cand_aa}")
    return mutations


def mutation_penalty(seed_sequence: str, candidate_sequence: str) -> int:
    """Count the number of changed positions between seed and candidate sequence."""
    shared = sum(a != b for a, b in zip(seed_sequence, candidate_sequence))
    return int(shared + abs(len(seed_sequence) - len(candidate_sequence)))


def parse_regimes_json(regimes_json: str | None, default_temperature: float, default_top_k: int, default_max_masked_positions: int, num_steps: int) -> list[Stage07Regime]:
    """Parse regime JSON or return a small default set that balances quality and diversity."""
    if regimes_json:
        raw = json.loads(regimes_json)
        return [
            Stage07Regime(
                name=str(item.get("name", f"regime_{idx}")),
                temperature=float(item.get("temperature", default_temperature)),
                top_k=int(item.get("top_k", default_top_k)),
                max_masked_positions=int(item.get("max_masked_positions", default_max_masked_positions)),
                hotspot_strategy=str(item.get("hotspot_strategy", "mixed")),
                num_steps=int(item.get("num_steps", num_steps)),
            )
            for idx, item in enumerate(raw)
        ]
    return [
        Stage07Regime("conservative", max(0.55, default_temperature - 0.1), max(4, default_top_k - 1), max(8, min(default_max_masked_positions, 16)), "priority", num_steps),
        Stage07Regime("balanced", default_temperature, default_top_k, default_max_masked_positions, "mixed", num_steps),
        Stage07Regime("exploratory", min(1.1, default_temperature + 0.15), max(default_top_k, 8), max(default_max_masked_positions, 32), "mixed", num_steps),
    ]


def _normalize_dict(values: dict[int, float]) -> dict[int, float]:
    """Min-max normalize a position keyed mapping into the [0, 1] interval."""
    if not values:
        return {}
    arr = np.asarray(list(values.values()), dtype=np.float32)
    lo = float(arr.min())
    hi = float(arr.max())
    if hi <= lo:
        return {int(k): 0.0 for k in values}
    return {int(k): float((v - lo) / (hi - lo)) for k, v in values.items()}


def build_position_feature_table(seed_sequence: str, family_sequences: list[str], hotspots_1based: list[int], target_position_priors: list[float] | dict[str, float] | None = None) -> list[dict]:
    """Summarize mutability and guidance features for each editable position."""
    hotspots = sorted({int(pos) for pos in hotspots_1based if 1 <= int(pos) <= len(seed_sequence)})
    if not hotspots:
        return []

    if isinstance(target_position_priors, list):
        target_lookup = {i + 1: float(v) for i, v in enumerate(target_position_priors)}
    else:
        target_lookup = {int(k): float(v) for k, v in (target_position_priors or {}).items()}

    midpoint = sum(hotspots) / max(len(hotspots), 1)
    span = max(max(hotspots) - min(hotspots), 1)
    family_weights: dict[int, float] = {}
    center_weights: dict[int, float] = {}
    target_weights: dict[int, float] = {pos: target_lookup.get(pos, 0.0) for pos in hotspots}
    rows = []

    for pos in hotspots:
        idx = pos - 1
        counts: dict[str, int] = {}
        total = 0
        for seq in family_sequences:
            if idx < len(seq):
                aa = str(seq[idx])
                counts[aa] = counts.get(aa, 0) + 1
                total += 1
        total = max(total, 1)
        probs = np.asarray([count / total for count in counts.values()], dtype=np.float32)
        entropy = float(-(probs * np.log(probs + 1e-12)).sum() / math.log(max(len(counts), 2))) if len(counts) > 1 else 0.0
        seed_aa = seed_sequence[idx]
        seed_freq = float(counts.get(seed_aa, 0) / total)
        consensus_aa, consensus_count = max(counts.items(), key=lambda item: item[1]) if counts else (seed_aa, 0)
        consensus_freq = float(consensus_count / total)
        family_mutability = max(0.0, min(1.0, 0.5 * entropy + 0.5 * (1.0 - seed_freq)))
        family_weights[pos] = family_mutability
        center_weights[pos] = float(1.0 - abs(pos - midpoint) / span)
        rows.append(
            {
                "position": int(pos),
                "seed_aa": seed_aa,
                "consensus_aa": consensus_aa,
                "seed_freq": seed_freq,
                "consensus_freq": consensus_freq,
                "family_entropy": entropy,
                "family_mutability": family_mutability,
                "target_prior": float(target_lookup.get(pos, 0.0)),
                "center_weight": center_weights[pos],
            }
        )

    family_norm = _normalize_dict(family_weights)
    target_norm = _normalize_dict(target_weights)
    center_norm = _normalize_dict(center_weights)
    for row in rows:
        pos = int(row["position"])
        row["functional_weight"] = float(0.50 * target_norm.get(pos, 0.0) + 0.35 * family_norm.get(pos, 0.0) + 0.15 * center_norm.get(pos, 0.0))
    return sorted(rows, key=lambda row: (-row["functional_weight"], row["position"]))


def build_structured_windows(position_features: list[dict], seed_sequence: str, window_size: int = 24, top_k: int = 3) -> list[dict]:
    """Create compact contiguous windows centered on the strongest functional regions."""
    if not position_features:
        return []
    pos_to_score = {int(row["position"]): float(row["functional_weight"]) for row in position_features}
    positions = sorted(pos_to_score)
    min_pos, max_pos = min(positions), max(positions)
    half = max(window_size // 2, 1)
    window_rows = []
    for center in positions:
        start = max(1, center - half)
        end = min(len(seed_sequence), start + window_size - 1)
        start = max(1, end - window_size + 1)
        window_positions = list(range(start, end + 1))
        covered_scores = [pos_to_score.get(pos, 0.0) for pos in window_positions]
        window_rows.append(
            {
                "name": f"window_{start}_{end}",
                "window_start": int(start),
                "window_end": int(end + 1),
                "positions": window_positions,
                "mean_functional_weight": float(np.mean(covered_scores)) if covered_scores else 0.0,
                "max_functional_weight": float(max(covered_scores, default=0.0)),
            }
        )
    window_rows = sorted(window_rows, key=lambda row: (-row["mean_functional_weight"], -row["max_functional_weight"], row["window_start"]))

    selected = []
    selected_positions: list[set[int]] = []
    for row in window_rows:
        pos_set = set(row["positions"])
        if any(len(pos_set & chosen) / max(len(pos_set), 1) > 0.6 for chosen in selected_positions):
            continue
        selected.append(row)
        selected_positions.append(pos_set)
        if len(selected) >= top_k:
            break
    return selected


def choose_hotspots(hotspots_1based: list[int], priority_weights: dict[str, float] | dict[int, float] | None, max_positions: int, strategy: str, sample_seed: int) -> list[int]:
    """Choose a compact subset of editable positions using even, priority, or mixed sampling."""
    unique_hotspots = sorted({int(pos) for pos in hotspots_1based})
    if len(unique_hotspots) <= max_positions:
        return unique_hotspots

    weights = {int(k): float(v) for k, v in (priority_weights or {}).items()}
    ranked = sorted(unique_hotspots, key=lambda pos: (weights.get(pos, 0.0), -pos), reverse=True)

    def spaced_pick(source: list[int], count: int) -> list[int]:
        if count <= 0:
            return []
        if len(source) <= count:
            return sorted(source)
        idxs = np.linspace(0, len(source) - 1, num=count, dtype=int)
        return sorted({source[i] for i in idxs})

    if strategy == "even":
        return spaced_pick(unique_hotspots, max_positions)
    if strategy == "priority":
        return spaced_pick(ranked, max_positions)

    rng = random.Random(sample_seed)
    head = ranked[: max_positions // 2]
    remaining = [pos for pos in unique_hotspots if pos not in head]
    rng.shuffle(remaining)
    mixed = sorted(set(head + remaining[: max_positions - len(head)]))
    return spaced_pick(mixed, max_positions)


def make_masked_prompt(sequence: str, masked_positions_1based: Iterable[int]) -> str:
    """Replace selected 1-based positions with underscores for ESM3 prompting."""
    chars = list(sequence)
    for pos in masked_positions_1based:
        if 1 <= pos <= len(chars):
            chars[pos - 1] = "_"
    return "".join(chars)


def candidate_guidance_score(seed_sequence: str, candidate_sequence: str, position_features: list[dict], mutated_positions: list[int]) -> float:
    """Score a candidate using target-aware mutability and family-compatible substitutions."""
    if not mutated_positions:
        return -1.0
    feature_lookup = {int(row["position"]): row for row in position_features}
    scores = []
    for pos in mutated_positions:
        feature = feature_lookup.get(int(pos))
        if feature is None or pos > len(candidate_sequence) or pos > len(seed_sequence):
            continue
        candidate_aa = candidate_sequence[pos - 1]
        seed_aa = seed_sequence[pos - 1]
        family_bonus = 1.0 if candidate_aa == feature["consensus_aa"] and candidate_aa != seed_aa else 0.0
        support_bonus = float(feature["consensus_freq"] * (candidate_aa == feature["consensus_aa"]))
        scores.append(float(0.55 * feature["functional_weight"] + 0.20 * feature["family_mutability"] + 0.15 * family_bonus + 0.10 * support_bonus))
    return float(np.mean(scores)) if scores else -1.0


def embed_sequences(sequences: list[str], model_name: str, batch_size: int = 4, max_length: int = 2048) -> np.ndarray:
    """Embed protein sequences with an ESM model using masked-mean pooling over token states."""
    from transformers import AutoModel, AutoTokenizer

    device = "cuda" if torch.cuda.is_available() else "cpu"
    tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)
    model = AutoModel.from_pretrained(model_name)
    model.to(device)
    model.eval()

    embeddings = []
    with torch.no_grad():
        for start in range(0, len(sequences), batch_size):
            batch = sequences[start : start + batch_size]
            toks = tokenizer(batch, padding=True, truncation=True, max_length=max_length, return_tensors="pt")
            toks = {key: value.to(device) for key, value in toks.items()}
            hidden = model(**toks).last_hidden_state
            mask = toks["attention_mask"].unsqueeze(-1)
            pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1)
            embeddings.append(pooled.cpu().numpy())
    return np.vstack(embeddings).astype(np.float32)


def greedy_diverse_order(embeddings: np.ndarray, base_scores: np.ndarray, penalty_weight: float = 0.25, preferred_mask: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Produce a diversity-aware greedy order and the associated penalties."""
    n = len(base_scores)
    if n == 0:
        return np.array([], dtype=int), np.array([], dtype=np.float32)
    normed = normalize_rows(embeddings)
    remaining = list(range(n))
    order = []
    penalties = np.zeros(n, dtype=np.float32)
    preferred_mask = np.ones(n, dtype=bool) if preferred_mask is None else np.asarray(preferred_mask, dtype=bool)

    first_pool = np.where(preferred_mask)[0]
    first = int(first_pool[np.argmax(base_scores[first_pool])]) if len(first_pool) else int(np.argmax(base_scores))
    order.append(first)
    remaining.remove(first)

    while remaining:
        best_idx = None
        best_value = -math.inf
        for idx in remaining:
            nearest = max(float(np.dot(normed[idx], normed[chosen])) for chosen in order)
            preference_bonus = 0.02 if preferred_mask[idx] else 0.0
            penalized = float(base_scores[idx] - penalty_weight * nearest + preference_bonus)
            if penalized > best_value:
                best_value = penalized
                best_idx = idx
                penalties[idx] = nearest
        order.append(int(best_idx))
        remaining.remove(int(best_idx))
    return np.asarray(order, dtype=int), penalties


def write_fasta(records: list[tuple[str, str]], path: str | Path) -> None:
    """Write a simple FASTA file from (header, sequence) tuples."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")
