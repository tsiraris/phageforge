"""Shared helpers for Stage 09 structure-aware localized redesign.

This module centralizes the reusable logic needed by the new Stage 09 scripts:
sequence-level diagnostics, edit-space parsing, substitution priors, target-model
scoring, surrogate loading, and diversity-aware filtering.
"""

from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from phageforge.stage07_utils import AMINO_ACIDS, cosine_similarity, embed_sequences, normalize_rows


VALID_AA = set(AMINO_ACIDS)


@dataclass
class EditProposal:
    """Compact container describing the allowed substitutions for one residue position."""

    position: int
    seed_aa: str
    allowed_aas: list[str]
    target_preference: dict[str, float]
    family_preference: dict[str, float]
    functional_weight: float
    conservation_penalty: float
    region_name: str


@dataclass
class SearchCandidate:
    """Container for one Stage 09 search state and its current metadata."""

    candidate_sequence: str
    mutations: list[str]
    mutated_positions: list[int]
    proposal_trace: list[str]
    round_index: int
    base_parent_id: str


# ----------------------------- JSON and path helpers ----------------------------- #


def read_json(path: str | Path) -> dict:
    """Read a JSON file and return the parsed Python object."""
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)



def write_json(obj: dict, path: str | Path) -> None:
    """Write a JSON object using stable indentation and UTF-8 encoding."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2)


# ----------------------------- Sequence metrics and parsing ----------------------------- #


def mutation_list(seed_sequence: str, candidate_sequence: str) -> list[str]:
    """Return human-readable 1-based mutation strings such as 197:W→R."""
    out: list[str] = []
    for idx, (seed_aa, cand_aa) in enumerate(zip(seed_sequence, candidate_sequence), start=1):
        if seed_aa != cand_aa:
            out.append(f"{idx}:{seed_aa}→{cand_aa}")
    return out



def parse_mutation_positions(mutation_text: str | Iterable[str] | None) -> list[int]:
    """Extract 1-based mutation positions from a string or list of mutation annotations."""
    if mutation_text is None:
        return []
    if isinstance(mutation_text, str):
        tokens = [tok.strip() for tok in mutation_text.split(";") if tok.strip()]
    else:
        tokens = [str(tok).strip() for tok in mutation_text if str(tok).strip()]
    positions: list[int] = []
    for token in tokens:
        if ":" not in token:
            continue
        try:
            positions.append(int(token.split(":", 1)[0]))
        except ValueError:
            continue
    return sorted(set(positions))



def normalized_shannon_entropy(sequence: str) -> float:
    """Return Shannon entropy normalized to [0, 1] for amino-acid composition."""
    if not sequence:
        return 0.0
    counts = pd.Series(list(sequence)).value_counts(normalize=True)
    entropy = float(-(counts * np.log2(counts)).sum())
    return entropy / math.log2(len(VALID_AA))



def max_residue_fraction(sequence: str) -> float:
    """Return the largest single-residue frequency fraction in the sequence."""
    if not sequence:
        return 1.0
    counts = pd.Series(list(sequence)).value_counts(normalize=True)
    return float(counts.max())



def longest_homopolymer_run(sequence: str) -> int:
    """Return the longest run of the same residue."""
    if not sequence:
        return 0
    best = run = 1
    prev = sequence[0]
    for aa in sequence[1:]:
        if aa == prev:
            run += 1
            best = max(best, run)
        else:
            prev = aa
            run = 1
    return best



def low_complexity_fraction(sequence: str, k: int = 12, unique_threshold: int = 3) -> float:
    """Return the fraction of windows with very few unique residues."""
    if len(sequence) < k:
        return float(len(set(sequence)) <= unique_threshold)
    flags = []
    for i in range(len(sequence) - k + 1):
        flags.append(len(set(sequence[i : i + k])) <= unique_threshold)
    return float(np.mean(flags)) if flags else 0.0



def mutation_span(positions: Iterable[int]) -> int:
    """Return the span covered by the selected mutations."""
    pos = sorted(set(int(p) for p in positions))
    if len(pos) < 2:
        return 0
    return pos[-1] - pos[0]



def outside_editable_fraction(positions: Iterable[int], editable_positions: set[int]) -> float:
    """Return the fraction of mutated positions that fall outside the allowed edit set."""
    pos = [int(p) for p in positions]
    if not pos:
        return 0.0
    return float(np.mean([p not in editable_positions for p in pos]))



def sequence_identity(seq_a: str, seq_b: str) -> float:
    """Return position-wise identity fraction for two same-length sequences."""
    if not seq_a or not seq_b or len(seq_a) != len(seq_b):
        return 0.0
    return float(sum(a == b for a, b in zip(seq_a, seq_b)) / len(seq_a))



def build_basic_sequence_features(seed_sequence: str, candidate_sequence: str, editable_positions: set[int]) -> dict[str, float | int]:
    """Compute lightweight candidate features used by Stage 09 search and surrogate scoring."""
    mutations = mutation_list(seed_sequence, candidate_sequence)
    positions = parse_mutation_positions(mutations)
    return {
        "mutation_count": len(positions),
        "mutation_span": mutation_span(positions),
        "normalized_entropy": normalized_shannon_entropy(candidate_sequence),
        "max_single_residue_fraction": max_residue_fraction(candidate_sequence),
        "longest_homopolymer_run": longest_homopolymer_run(candidate_sequence),
        "low_complexity_fraction": low_complexity_fraction(candidate_sequence),
        "outside_editable_fraction": outside_editable_fraction(positions, editable_positions),
        "seed_identity": sequence_identity(seed_sequence, candidate_sequence),
    }


# ----------------------------- Edit-space and substitution-prior helpers ----------------------------- #


def infer_target_residue_preferences(strict_df: pd.DataFrame, target_host: str, positions_1based: list[int]) -> dict[int, dict[str, float]]:
    """Estimate per-position residue preferences from target-host family rows in the strict bank."""
    target_rows = strict_df.loc[strict_df["host_genus"].astype(str) == str(target_host)].copy()
    if target_rows.empty:
        return {int(pos): {} for pos in positions_1based}

    preferences: dict[int, dict[str, float]] = {}
    for pos in positions_1based:
        residues = []
        idx = int(pos) - 1
        for seq in target_rows["aa_sequence"].astype(str):
            if idx < len(seq):
                aa = seq[idx]
                if aa in VALID_AA:
                    residues.append(aa)
        if not residues:
            preferences[int(pos)] = {}
            continue
        counts = pd.Series(residues).value_counts(normalize=True)
        preferences[int(pos)] = {str(aa): float(freq) for aa, freq in counts.items()}
    return preferences



def build_edit_proposals_from_context(context: dict, strict_df: pd.DataFrame | None = None) -> list[EditProposal]:
    """Build per-position substitution proposals from the Stage 07 context and optional strict bank."""
    seed_sequence = str(context["selected_seed"]["seed_sequence"])
    target_host = str(context["target_host"])
    position_features = list(context["editable_region"].get("position_features", []))
    structured_windows = list(context["editable_region"].get("structured_windows", []))

    # Map positions to their strongest structured window so proposal records can keep region provenance.
    position_to_region: dict[int, str] = {}
    for window in structured_windows:
        region_name = str(window.get("name", f"window_{window.get('window_start', 'na')}_{window.get('window_end', 'na')}"))
        for pos in window.get("positions", []):
            pos = int(pos)
            if pos not in position_to_region:
                position_to_region[pos] = region_name

    # Recover family residue frequencies directly from the family rows stored in the context.
    family_rows = pd.DataFrame(context["family_context"].get("family_rows", []))
    family_sequences = family_rows.get("aa_sequence", pd.Series(dtype=str)).astype(str).tolist() if not family_rows.empty else []
    positions_1based = [int(row["position"]) for row in position_features]
    target_preferences = infer_target_residue_preferences(strict_df, target_host, positions_1based) if strict_df is not None else {int(pos): {} for pos in positions_1based}

    proposals: list[EditProposal] = []
    for row in position_features:
        pos = int(row["position"])
        idx = pos - 1
        seed_aa = seed_sequence[idx]

        # Gather the family-supported residues observed at this seed-local position.
        residues = []
        for seq in family_sequences:
            if idx < len(seq):
                aa = seq[idx]
                if aa in VALID_AA:
                    residues.append(aa)
        family_counts = pd.Series(residues).value_counts(normalize=True) if residues else pd.Series(dtype=float)
        family_preference = {str(aa): float(freq) for aa, freq in family_counts.items()}

        # Build a compact allowed set anchored in both family support and target-host residue prevalence.
        allowed = set()
        allowed.update(list(family_preference.keys())[:4])
        allowed.update(list(target_preferences.get(pos, {}).keys())[:4])
        allowed.discard(seed_aa)
        allowed = {aa for aa in allowed if aa in VALID_AA}

        # Keep only positions that have at least one plausible alternative residue.
        if not allowed:
            continue

        proposals.append(
            EditProposal(
                position=pos,
                seed_aa=seed_aa,
                allowed_aas=sorted(allowed),
                target_preference=target_preferences.get(pos, {}),
                family_preference=family_preference,
                functional_weight=float(row.get("functional_weight", 0.0)),
                conservation_penalty=float(row.get("seed_freq", 0.0)),
                region_name=position_to_region.get(pos, "ungrouped"),
            )
        )
    return proposals



def choose_editable_positions(proposals: list[EditProposal], max_positions: int, seed: int) -> list[int]:
    """Choose a compact seed-local editable subset while preserving window diversity."""
    if len(proposals) <= max_positions:
        return [int(item.position) for item in proposals]

    rng = random.Random(seed)
    grouped: dict[str, list[EditProposal]] = {}
    for item in proposals:
        grouped.setdefault(item.region_name, []).append(item)
    for items in grouped.values():
        items.sort(key=lambda x: (x.functional_weight, -x.conservation_penalty, x.position), reverse=True)

    chosen: list[int] = []

    # First cover the strongest position from each region so the search does not over-concentrate immediately.
    for region_name in sorted(grouped, key=lambda name: max(item.functional_weight for item in grouped[name]), reverse=True):
        chosen.append(int(grouped[region_name][0].position))
        if len(chosen) >= max_positions:
            return sorted(chosen[:max_positions])

    # Then fill the remaining budget with the strongest global positions.
    remaining = [item for item in proposals if int(item.position) not in chosen]
    remaining.sort(key=lambda x: (x.functional_weight, -x.conservation_penalty, x.position), reverse=True)
    rng.shuffle(remaining[max_positions:])
    chosen.extend(int(item.position) for item in remaining[: max_positions - len(chosen)])
    return sorted(chosen[:max_positions])



def substitution_priority(item: EditProposal, aa: str) -> float:
    """Combine target preference, family support, and position weight into one substitution priority."""
    target = float(item.target_preference.get(aa, 0.0))
    family = float(item.family_preference.get(aa, 0.0))
    return float(0.50 * target + 0.30 * family + 0.20 * item.functional_weight)


# ----------------------------- Predictor and surrogate helpers ----------------------------- #


def _patch_logistic_regression_compat(model) -> object:
    """Restore LogisticRegression attributes that may be missing after sklearn version drift."""
    if not isinstance(model, LogisticRegression):
        return model
    if not hasattr(model, "multi_class"):
        model.multi_class = "auto"
    if not hasattr(model, "n_features_in_") and hasattr(model, "coef_"):
        model.n_features_in_ = int(model.coef_.shape[1])
    return model



def load_target_predictor(model_path: str | Path, label_classes_path: str | Path | None) -> tuple[object, list[str]]:
    """Load a host predictor and the associated label order used for class probabilities."""
    model = _patch_logistic_regression_compat(joblib.load(model_path))
    if label_classes_path is None:
        classes = [str(x) for x in getattr(model, "classes_", [])]
    else:
        with open(label_classes_path, "r", encoding="utf-8") as handle:
            classes = [str(x) for x in json.load(handle)]
    if not classes:
        raise ValueError("Could not recover label classes for the target predictor.")
    return model, classes



def predict_target_probability(model, label_classes: list[str], target_host: str, embeddings: np.ndarray) -> np.ndarray:
    """Return target-host probability for each embedding row using the loaded host predictor."""
    if target_host not in label_classes:
        raise ValueError(f"Target host '{target_host}' not found in predictor label order: {label_classes}")
    embeddings = np.asarray(embeddings, dtype=np.float32)
    if hasattr(model, "n_features_in_") and embeddings.shape[1] != int(model.n_features_in_):
        raise ValueError(f"Predictor expects {int(model.n_features_in_)} embedding features, but received {embeddings.shape[1]}. Use the embedding backbone that matches the trained predictor.")
    probs = model.predict_proba(embeddings)
    target_idx = label_classes.index(target_host)
    return np.asarray(probs[:, target_idx], dtype=np.float32)



def maybe_load_surrogate(path: str | Path | None) -> dict | None:
    """Load a saved Stage 09 structural surrogate bundle if a path was provided."""
    if path is None:
        return None
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing surrogate model bundle: {path}")
    return joblib.load(path)



def surrogate_structural_risk(bundle: dict | None, feature_frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return structural risk, predicted mean pLDDT, and predicted RMSD for each candidate.

    When a fitted surrogate bundle is unavailable, this function falls back to a rule-based score
    that penalizes large drift, weak entropy, and edit patterns that already failed in Stage 08.
    """
    df = feature_frame.copy()
    if bundle is None:
        risk = (
            0.40 * (1.0 - np.clip(df["seed_identity"].to_numpy(dtype=float), 0.0, 1.0))
            + 0.15 * np.clip(df["outside_editable_fraction"].to_numpy(dtype=float), 0.0, 1.0)
            + 0.15 * np.clip(df["low_complexity_fraction"].to_numpy(dtype=float), 0.0, 1.0)
            + 0.10 * np.clip((df["longest_homopolymer_run"].to_numpy(dtype=float) - 4.0) / 4.0, 0.0, 1.0)
            + 0.20 * np.clip(np.abs(df["mutation_count"].to_numpy(dtype=float) - 6.0) / 6.0, 0.0, 1.0)
        )
        pred_plddt = 90.0 - 35.0 * risk
        pred_rmsd = 1.0 + 6.0 * risk
        return np.asarray(risk, dtype=np.float32), np.asarray(pred_plddt, dtype=np.float32), np.asarray(pred_rmsd, dtype=np.float32)

    feature_cols = list(bundle["feature_columns"])
    X = df.reindex(columns=feature_cols, fill_value=0.0).to_numpy(dtype=np.float32)

    pass_model = bundle.get("pass_model")
    plddt_model = bundle.get("plddt_model")
    rmsd_model = bundle.get("rmsd_model")

    pass_prob = pass_model.predict_proba(X)[:, 1] if pass_model is not None else np.full(len(df), 0.5, dtype=np.float32)
    pred_plddt = plddt_model.predict(X) if plddt_model is not None else np.full(len(df), 65.0, dtype=np.float32)
    pred_rmsd = rmsd_model.predict(X) if rmsd_model is not None else np.full(len(df), 3.5, dtype=np.float32)
    risk = 1.0 - np.clip(pass_prob, 0.0, 1.0)
    return np.asarray(risk, dtype=np.float32), np.asarray(pred_plddt, dtype=np.float32), np.asarray(pred_rmsd, dtype=np.float32)


# ----------------------------- Diversity and selection helpers ----------------------------- #


def greedy_diverse_pick(embeddings: np.ndarray, scores: np.ndarray, top_k: int, penalty_weight: float = 0.25) -> tuple[list[int], np.ndarray]:
    """Select a diversity-aware top-k set by greedily penalizing near-duplicates."""
    if len(scores) == 0:
        return [], np.asarray([], dtype=np.float32)
    order: list[int] = []
    penalties = np.zeros(len(scores), dtype=np.float32)
    normed = normalize_rows(np.asarray(embeddings, dtype=np.float32))
    remaining = list(range(len(scores)))

    # Seed the diverse panel with the strongest-scoring candidate.
    first = int(np.argmax(scores))
    order.append(first)
    remaining.remove(first)

    # Add candidates greedily while penalizing those too close to the already selected ones.
    while remaining and len(order) < top_k:
        best_idx = None
        best_value = -math.inf
        for idx in remaining:
            nearest = max(float(np.dot(normed[idx], normed[chosen])) for chosen in order)
            penalized = float(scores[idx] - penalty_weight * nearest)
            if penalized > best_value:
                best_value = penalized
                best_idx = idx
                penalties[idx] = nearest
        order.append(int(best_idx))
        remaining.remove(int(best_idx))
    return order, penalties
