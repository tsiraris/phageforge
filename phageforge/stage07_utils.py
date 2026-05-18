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


AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")                                                                      # Creates a list of the 20 standard amino acid characters for sequence processing.


@dataclass                                                                                                      # Decorator that automatically generates special methods like __init__ for the class.
class Stage07Regime:
    """Compact description of one generation regime for Stage 07."""

    name: str                                                                                                   # Defines the string name identifier for this specific generation strategy.
    temperature: float                                                                                          # Defines the sampling temperature, controlling randomness in model predictions.
    top_k: int                                                                                                  # Defines the token restriction limit, restricting generation to the top K probable tokens.
    max_masked_positions: int                                                                                   # Defines the maximum number of residues that can be simultaneously mutated.
    hotspot_strategy: str = "mixed"                                                                             # Defines the method used to select mutation sites, defaulting to a mixed approach.
    num_steps: int = 8                                                                                          # Defines the number of iterative unmasking steps ESM3 should perform during generation.


def seed_everything(seed: int) -> None:
    """Set deterministic random seeds for python, numpy, and torch."""
    random.seed(seed)                                                                                           # Sets the seed for Python's built-in random number generator to ensure reproducibility.
    np.random.seed(seed)                                                                                        # Sets the seed for the NumPy library's random number operations.
    torch.manual_seed(seed)                                                                                     # Sets the seed for PyTorch's CPU-based random number generation.
    if torch.cuda.is_available():                                                                               # Checks if a CUDA-compatible GPU is currently available to PyTorch.
        torch.cuda.manual_seed_all(seed)                                                                        # Sets the deterministic seed across all available GPUs if they are present.


def read_json(path: str | Path) -> dict:
    """Read a JSON file and return the parsed object."""
    with open(path, "r", encoding="utf-8") as handle:                                                           # Opens the specified file in read mode using UTF-8 encoding safely via a context manager.
        return json.load(handle)                                                                                # Parses the JSON file content into a Python dictionary and returns it.


def write_json(obj: dict, path: str | Path) -> None:
    """Write a JSON object with stable indentation and UTF-8 encoding."""
    path = Path(path)                                                                                           # Converts the input path string into a flexible Pathlib object.
    path.parent.mkdir(parents=True, exist_ok=True)                                                              # Creates any missing parent directories for the file, ignoring errors if they exist.
    with open(path, "w", encoding="utf-8") as handle:                                                           # Opens the target file in write mode using UTF-8 encoding safely via a context manager.
        json.dump(obj, handle, indent=2)                                                                        # Serializes the Python dictionary to JSON with a 2-space indentation for readability.


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """Return cosine similarity between two vectors, guarding against zero norms."""
    a = np.asarray(a, dtype=np.float32)                                                                         # Converts the first input to a 32-bit float NumPy array for consistent math.
    b = np.asarray(b, dtype=np.float32)                                                                         # Converts the second input to a 32-bit float NumPy array for consistent math.
    denom = np.linalg.norm(a) * np.linalg.norm(b)                                                               # Calculates the denominator by multiplying the Euclidean lengths (L2 norms) of both vectors.
    if denom <= 0:                                                                                              # Checks if either vector has a magnitude of zero to prevent division-by-zero errors.
        return 0.0                                                                                              # Returns a similarity of 0.0 if a zero-magnitude vector is encountered.
    return float(np.dot(a, b) / denom)                                                                          # Computes the dot product, divides by the magnitude denominator, and casts to standard float.


def normalize_rows(x: np.ndarray) -> np.ndarray:
    """L2-normalize each row for cosine-based comparisons."""
    x = np.asarray(x, dtype=np.float32)                                                                         # Converts the input matrix to a 32-bit float NumPy array.
    denom = np.linalg.norm(x, axis=1, keepdims=True)                                                            # Calculates the L2 norm for each row independently, keeping the matrix dimensions intact.
    return x / np.clip(denom, 1e-12, None)                                                                      # Divides the matrix by the norms, clamping the minimum denominator to 1e-12 to prevent zero-division.


def mutation_list(seed_sequence: str, candidate_sequence: str) -> list[str]:
    """Return human-readable residue substitutions using 1-based positions."""
    mutations = []                                                                                              # Initializes an empty list to store the formatted mutation strings.
    for i, (seed_aa, cand_aa) in enumerate(zip(seed_sequence, candidate_sequence), start=1):                    # Iterates through both sequences simultaneously, tracking the 1-indexed position.
        if seed_aa != cand_aa:                                                                                  # Checks if the amino acid has changed at the current position.
            mutations.append(f"{i}:{seed_aa}→{cand_aa}")                                                        # Appends a formatted string (e.g., "15:A→T") to the mutations list.
    if len(candidate_sequence) > len(seed_sequence):                                                            # Checks if the candidate sequence has additional residues appended at the end.
        for i, cand_aa in enumerate(candidate_sequence[len(seed_sequence) :], start=len(seed_sequence) + 1):    # Iterates over the extra candidate residues, starting the index from the end of the seed.
            mutations.append(f"{i}:∅→{cand_aa}")                                                                # Appends a formatted string indicating an insertion (e.g., "100:∅→G").
    return mutations                                                                                            # Returns the complete list of formatted mutation strings.


def mutation_penalty(seed_sequence: str, candidate_sequence: str) -> int:
    """Count the number of changed positions between seed and candidate sequence."""
    shared = sum(a != b for a, b in zip(seed_sequence, candidate_sequence))                                     # Counts how many aligned residues differ between the two sequences.
    return int(shared + abs(len(seed_sequence) - len(candidate_sequence)))                                      # Adds the length difference (insertions/deletions) to the substitution count and returns the integer total.


def parse_regimes_json(regimes_json: str | None, default_temperature: float, default_top_k: int, default_max_masked_positions: int, num_steps: int) -> list[Stage07Regime]:
    """Parse regime JSON or return a small default set that balances quality and diversity."""
    if regimes_json:                                                                                                                                                                          # Checks if a custom JSON string of regimes was provided as an argument.
        raw = json.loads(regimes_json)                                                                                                                                                        # Parses the raw JSON string into a list of Python dictionaries.
        return [                                                                                                                                                                              # Begins a list comprehension to build Stage07Regime objects from the parsed data.
            Stage07Regime(                                                                                                                                                                    # Instantiates a new Stage07Regime data class.
                name=str(item.get("name", f"regime_{idx}")),                                                                                                                                  # Extracts the regime name, falling back to an indexed default if missing.
                temperature=float(item.get("temperature", default_temperature)),                                                                                                              # Extracts the temperature, falling back to the globally provided default.
                top_k=int(item.get("top_k", default_top_k)),                                                                                                                                  # Extracts the top-k parameter, falling back to the globally provided default.
                max_masked_positions=int(item.get("max_masked_positions", default_max_masked_positions)),                                                                                     # Extracts the masking limit, falling back to the globally provided default.
                hotspot_strategy=str(item.get("hotspot_strategy", "mixed")),                                                                                                                  # Extracts the selection strategy, defaulting to "mixed" if missing.
                num_steps=int(item.get("num_steps", num_steps)),                                                                                                                              # Extracts the ESM3 step count, falling back to the globally provided default.
            )                                                                                                                                                                                 # Closes the object instantiation.
            for idx, item in enumerate(raw)                                                                                                                                                   # Iterates over every dictionary in the parsed JSON list.
        ]                                                                                                                                                                                     # Closes the list comprehension and returns the custom regimes.
    return [                                                                                                                                                                                  # Begins the fallback return list containing three default archetypal regimes.
        Stage07Regime("conservative", max(0.55, default_temperature - 0.1), max(4, default_top_k - 1), max(8, min(default_max_masked_positions, 16)), "priority", num_steps),                 # Creates a cautious regime that strictly adheres to high-priority functional sites.
        Stage07Regime("balanced", default_temperature, default_top_k, default_max_masked_positions, "mixed", num_steps),                                                                      # Creates a moderate regime using standard defaults and a mixed mutation strategy.
        Stage07Regime("exploratory", min(1.1, default_temperature + 0.15), max(default_top_k, 8), max(default_max_masked_positions, 32), "mixed", num_steps),                                 # Creates an aggressive regime that generates riskier, highly mutated candidate variations.
    ]                                                                                                                                                                                         # Closes the fallback list.


def _normalize_dict(values: dict[int, float]) -> dict[int, float]:
    """Min-max normalize a position keyed mapping into the [0, 1] interval."""
    if not values:                                                                                              # Checks if the input dictionary is empty to prevent processing errors.
        return {}                                                                                               # Returns an empty dictionary immediately if no data is present.
    arr = np.asarray(list(values.values()), dtype=np.float32)                                                   # Extracts all the numeric values from the dictionary and converts them into a NumPy array.
    lo = float(arr.min())                                                                                       # Finds the absolute minimum value across the array.
    hi = float(arr.max())                                                                                       # Finds the absolute maximum value across the array.
    if hi <= lo:                                                                                                # Checks if all values are identical (max equals min), which would cause zero-division.
        return {int(k): 0.0 for k in values}                                                                    # Returns a dictionary where all keys are mapped to 0.0 to handle uniform data safely.
    return {int(k): float((v - lo) / (hi - lo)) for k, v in values.items()}                                     # Applies the min-max normalization formula to every value and returns the new dictionary.


def build_position_feature_table(seed_sequence: str, family_sequences: list[str], hotspots_1based: list[int], target_position_priors: list[float] | dict[str, float] | None = None) -> list[dict]:
    """Summarize mutability and guidance features for each editable position."""
    hotspots = sorted({int(pos) for pos in hotspots_1based if 1 <= int(pos) <= len(seed_sequence)})             # Deduplicates, validates against sequence length, casts to int, and sorts the input positions.
    if not hotspots:                                                                                            # Checks if the sanitization process resulted in an empty list of valid hotspots.
        return []                                                                                               # Escapes the function early by returning an empty list if no valid sites exist.

    if isinstance(target_position_priors, list):                                                                # Checks if the provided host-specific priors were passed in as a sequential list format.
        target_lookup = {i + 1: float(v) for i, v in enumerate(target_position_priors)}                         # Converts the list into a 1-indexed dictionary mapping position to prior weight.
    else:                                                                                                       # Executes if the priors were provided as a dictionary or None.
        target_lookup = {int(k): float(v) for k, v in (target_position_priors or {}).items()}                   # Casts keys to integers and values to floats, handling the None case safely with an empty dict.

    midpoint = sum(hotspots) / max(len(hotspots), 1)                                                            # Calculates the mathematical average location of all allowed hotspots on the sequence.
    span = max(max(hotspots) - min(hotspots), 1)                                                                # Calculates the total distance between the earliest and latest permitted hotspots.
    family_weights: dict[int, float] = {}                                                                       # Initializes dictionary to track evolutionary variability for each position.
    center_weights: dict[int, float] = {}                                                                       # Initializes dictionary to track spatial centrality biases for each position.
    target_weights: dict[int, float] = {pos: target_lookup.get(pos, 0.0) for pos in hotspots}                   # Populates target bias tracking using the provided lookup map, defaulting missing spots to 0.0.
    rows = []                                                                                                   # Initializes master list to hold the final generated feature rows.

    for pos in hotspots:                                                                                        # Iterates sequentially through every valid mutation site.
        idx = pos - 1                                                                                           # Converts the 1-based biological coordinate back to a 0-based computer index.
        counts: dict[str, int] = {}                                                                             # Initializes a tracker to count how often each amino acid appears at this spot in nature.
        total = 0                                                                                               # Initializes a counter for total valid evolutionary observations at this site.
        for seq in family_sequences:                                                                            # Iterates over the multiple sequence alignment of related proteins.
            if idx < len(seq):                                                                                  # Ensures the current family sequence is long enough to have an amino acid at this index.
                aa = str(seq[idx])                                                                              # Extracts the observed amino acid character from the relative.
                counts[aa] = counts.get(aa, 0) + 1                                                              # Increments the observation count for this specific amino acid in the tracker.
                total += 1                                                                                      # Increments the total valid observations counter.
        total = max(total, 1)                                                                                   # Clamps the total to a minimum of 1 to prevent division-by-zero on completely gapped columns.
        probs = np.asarray([count / total for count in counts.values()], dtype=np.float32)                      # Calculates the raw probability frequency of each observed amino acid character.
        entropy = float(-(probs * np.log(probs + 1e-12)).sum() / math.log(max(len(counts), 2))) if len(counts) > 1 else 0.0 # Calculates Shannon entropy (variability metric) scaled between 0 and 1.
        seed_aa = seed_sequence[idx]                                                                            # Extracts the amino acid character present in the original starting sequence.
        seed_freq = float(counts.get(seed_aa, 0) / total)                                                       # Calculates how common the starting amino acid is in natural evolution.
        consensus_aa, consensus_count = max(counts.items(), key=lambda item: item[1]) if counts else (seed_aa, 0) # Identifies the most common (consensus) amino acid naturally found at this position.
        consensus_freq = float(consensus_count / total)                                                         # Calculates the probability frequency of the consensus amino acid.
        family_mutability = max(0.0, min(1.0, 0.5 * entropy + 0.5 * (1.0 - seed_freq)))                         # Calculates a composite mutability score: high if site is variable and seed is naturally rare.
        family_weights[pos] = family_mutability                                                                 # Stores the computed evolutionary mutability in the tracking dictionary.
        center_weights[pos] = float(1.0 - abs(pos - midpoint) / span)                                           # Calculates a spatial bias score rewarding positions closer to the center of the hotspot cluster.
        rows.append(                                                                                            # Begins adding the compiled data for this position to the output matrix.
            {                                                                                                   # Opens dictionary structure for the row.
                "position": int(pos),                                                                           # Records the 1-based index location.
                "seed_aa": seed_aa,                                                                             # Records the original character.
                "consensus_aa": consensus_aa,                                                                   # Records the most common natural character.
                "seed_freq": seed_freq,                                                                         # Records the natural frequency of the original character.
                "consensus_freq": consensus_freq,                                                               # Records the natural frequency of the consensus character.
                "family_entropy": entropy,                                                                      # Records the calculated evolutionary variability metric.
                "family_mutability": family_mutability,                                                         # Records the finalized mutability score.
                "target_prior": float(target_lookup.get(pos, 0.0)),                                             # Records the input target-specific importance score.
                "center_weight": center_weights[pos],                                                           # Records the spatial centrality score.
            }                                                                                                   # Closes the dictionary structure.
        )                                                                                                       # Closes the append operation.

    family_norm = _normalize_dict(family_weights)                                                               # Rescales all evolutionary scores across the sequence to sit perfectly between 0.0 and 1.0.
    target_norm = _normalize_dict(target_weights)                                                               # Rescales all target importance scores to sit perfectly between 0.0 and 1.0.
    center_norm = _normalize_dict(center_weights)                                                               # Rescales all spatial scores to sit perfectly between 0.0 and 1.0.
    for row in rows:                                                                                            # Iterates through the compiled feature rows a second time to calculate the final aggregate score.
        pos = int(row["position"])                                                                              # Extracts the position coordinate from the current dictionary row.
        row["functional_weight"] = float(0.50 * target_norm.get(pos, 0.0) + 0.35 * family_norm.get(pos, 0.0) + 0.15 * center_norm.get(pos, 0.0)) # Computes the final weighted priority score combining all logic.
    return sorted(rows, key=lambda row: (-row["functional_weight"], row["position"]))                           # Returns the feature table, strongly sorted so the most critical functional positions appear first.


def build_structured_windows(position_features: list[dict], seed_sequence: str, window_size: int = 24, top_k: int = 3) -> list[dict]:
    """Create compact contiguous windows centered on the strongest functional regions."""
    if not position_features:                                                                                   # Checks if the input feature table is empty to prevent execution errors.
        return []                                                                                               # Returns an empty list immediately if no features are available to process.
    pos_to_score = {int(row["position"]): float(row["functional_weight"]) for row in position_features}         # Converts the list of dictionaries into a fast mapping of position integer to final importance score.
    positions = sorted(pos_to_score)                                                                            # Creates an ordered list of all available coordinates to establish spatial boundaries.
    min_pos, max_pos = min(positions), max(positions)                                                           # Finds the absolute extreme edges of the provided coordinate space.
    half = max(window_size // 2, 1)                                                                             # Calculates the radius of the window, ensuring it is at least 1.
    window_rows = []                                                                                            # Initializes list to store evaluated window boundary proposals.
    for center in positions:                                                                                    # Iterates, treating every available hotspot as the potential center of a new editing window.
        start = max(1, center - half)                                                                           # Calculates the left edge of the window, capping it at 1 so it doesn't fall off the sequence.
        end = min(len(seed_sequence), start + window_size - 1)                                                  # Calculates the right edge, capping it at the maximum sequence length to prevent out-of-bounds.
        start = max(1, end - window_size + 1)                                                                   # Recalculates left edge to ensure the window is as large as possible if pushed against the right boundary.
        window_positions = list(range(start, end + 1))                                                          # Generates the exhaustive list of sequential 1-based coordinates covered by this window.
        covered_scores = [pos_to_score.get(pos, 0.0) for pos in window_positions]                               # Looks up the functional importance scores for all residues captured inside this window.
        window_rows.append(                                                                                     # Begins appending the evaluated window metadata to the candidate list.
            {                                                                                                   # Opens window dictionary block.
                "name": f"window_{start}_{end}",                                                                # Assigns a unique naming label based on the boundary coordinates.
                "window_start": int(start),                                                                     # Records the left boundary index.
                "window_end": int(end + 1),                                                                     # Records the pythonic right boundary index (exclusive logic).
                "positions": window_positions,                                                                  # Records the explicit list of captured sequence positions.
                "mean_functional_weight": float(np.mean(covered_scores)) if covered_scores else 0.0,            # Computes the average biological importance of the residues in this specific window.
                "max_functional_weight": float(max(covered_scores, default=0.0)),                               # Identifies the single most important residue captured in this window.
            }                                                                                                   # Closes window dictionary block.
        )                                                                                                       # Closes the append operation.
    window_rows = sorted(window_rows, key=lambda row: (-row["mean_functional_weight"], -row["max_functional_weight"], row["window_start"])) # Sorts proposed windows heavily prioritizing highest average importance.

    selected = []                                                                                               # Initializes list to hold the finalized non-redundant window choices.
    selected_positions: list[set[int]] = []                                                                     # Initializes list of sets to quickly track which coordinates have already been claimed by selected windows.
    for row in window_rows:                                                                                     # Iterates over the sorted, prioritized window proposals.
        pos_set = set(row["positions"])                                                                         # Converts the current proposal's coordinate list into a set for fast intersection math.
        if any(len(pos_set & chosen) / max(len(pos_set), 1) > 0.6 for chosen in selected_positions):            # Checks if this proposed window overlaps heavily (more than 60%) with a window we already selected.
            continue                                                                                            # Discards the proposal if it is too redundant, skipping to the next loop iteration.
        selected.append(row)                                                                                    # Approves the non-redundant window and adds it to the final output list.
        selected_positions.append(pos_set)                                                                      # Records the newly claimed coordinates so future windows can be checked against them.
        if len(selected) >= top_k:                                                                              # Checks if we have successfully found the requested number of distinct windows.
            break                                                                                               # Escapes the loop early once the quota is met to save computation.
    return selected                                                                                             # Returns the final curated list of distinct, highly-functional editing windows.


def choose_hotspots(hotspots_1based: list[int], priority_weights: dict[str, float] | dict[int, float] | None, max_positions: int, strategy: str, sample_seed: int) -> list[int]:
    """Choose a compact subset of editable positions using even, priority, or mixed sampling."""
    unique_hotspots = sorted({int(pos) for pos in hotspots_1based})                                             # Removes duplicates, casts to integers, and spatially orders the available coordinates.
    if len(unique_hotspots) <= max_positions:                                                                   # Checks if the provided pool is already smaller than the model's masking limit.
        return unique_hotspots                                                                                  # Returns the entire pool immediately, as no subset selection is necessary.

    weights = {int(k): float(v) for k, v in (priority_weights or {}).items()}                                   # Sanitizes the weight dictionary, handling missing arguments and converting keys to integers.
    ranked = sorted(unique_hotspots, key=lambda pos: (weights.get(pos, 0.0), -pos), reverse=True)               # Re-orders the available coordinates strictly from most functionally important to least important.

    def spaced_pick(source: list[int], count: int) -> list[int]:                                                # Defines an internal helper to pluck widely distributed items from a list.
        if count <= 0:                                                                                          # Checks if the request asked for zero or negative items.
            return []                                                                                           # Returns an empty list safely without breaking math.
        if len(source) <= count:                                                                                # Checks if the source list is smaller than the requested extraction amount.
            return sorted(source)                                                                               # Returns the entire sorted source to prevent indexing errors.
        idxs = np.linspace(0, len(source) - 1, num=count, dtype=int)                                            # Calculates evenly spaced array indices across the entire length of the source list.
        return sorted({source[i] for i in idxs})                                                                # Extracts the items at those indices, deduplicates them just in case, and returns them sorted.

    if strategy == "even":                                                                                      # Checks if the generation regime specifically requested spatially distributed masking.
        return spaced_pick(unique_hotspots, max_positions)                                                      # Returns evenly distributed positions purely based on spatial coordinates, ignoring functional scores.
    if strategy == "priority":                                                                                  # Checks if the generation regime strictly requested the most critical functional sites.
        return spaced_pick(ranked, max_positions)                                                               # Returns evenly distributed positions from the score-ranked list, capturing top hits across the spectrum.

    rng = random.Random(sample_seed)                                                                            # Instantiates an isolated random number generator bound to the reproducible seed.
    head = ranked[: max_positions // 2]                                                                         # Identifies the absolute highest-scoring half of the required quota to guarantee their inclusion.
    remaining = [pos for pos in unique_hotspots if pos not in head]                                             # Isolates all leftover coordinates that were not captured in the elite "head" group.
    rng.shuffle(remaining)                                                                                      # Randomizes the order of the leftover, lower-tier coordinates.
    mixed = sorted(set(head + remaining[: max_positions - len(head)]))                                          # Combines the guaranteed elite sites with a random sampling of lower-tier sites to fill the quota.
    return spaced_pick(mixed, max_positions)                                                                    # Performs a final spacing pass over the mixed pool to ensure a clean sequence distribution.


def make_masked_prompt(sequence: str, masked_positions_1based: Iterable[int]) -> str:
    """Replace selected 1-based positions with underscores for ESM3 prompting."""
    chars = list(sequence)                                                                                      # Converts the raw sequence string into a mutable list of individual characters.
    for pos in masked_positions_1based:                                                                         # Iterates through the provided list of coordinates slated for mutation.
        if 1 <= pos <= len(chars):                                                                              # Verifies that the coordinate is within the valid biological bounds of the sequence.
            chars[pos - 1] = "_"                                                                                # Replaces the character at the converted 0-based index with the ESM3 mask token.
    return "".join(chars)                                                                                       # Reconstructs the character list back into a continuous string and returns it.


def candidate_guidance_score(seed_sequence: str, candidate_sequence: str, position_features: list[dict], mutated_positions: list[int]) -> float:
    """
    Score a generated candidate using target-aware mutability and family-compatible substitutions.
    
    It compares the candidate_sequence to the original seed_sequence to find exactly what changed. 
    Then, it checks those specific changes against the position_features (the evolutionary data). 
    If the model mutated the protein to an amino acid that naturally appears in the broader protein family (consensus_aa), it gets a massive bonus (family_bonus).
    
    Returns an average of all the scores, or -1.0 if no valid scores could be calculated. 
    """
    if not mutated_positions:                                                                                   # Checks if the candidate sequence is completely identical to the seed.
        return -1.0                                                                                             # Returns a heavily penalized score, as unmodified sequences fail the generation objective.
    feature_lookup = {int(row["position"]): row for row in position_features}                                   # Creates a quick dictionary mapping to access functional feature logic for specific sites.
    scores = []                                                                                                 # Initializes list to track the individual biological logic scores of each induced mutation.
    
    for pos in mutated_positions:                                                                               # Iterates only over the sites where the candidate actually changed an amino acid.
        feature = feature_lookup.get(int(pos))                                                                  # Retrieves the pre-calculated functional metadata for this specific mutation coordinate.
        if feature is None or pos > len(candidate_sequence) or pos > len(seed_sequence):                        # Validates that the metadata exists and the mutation doesn't exceed sequence bounds.
            continue                                                                                            # Skips the scoring for this residue if it is malformed or out of bounds.
        candidate_aa = candidate_sequence[pos - 1]                                                              # Identifies exactly which new amino acid was inserted at this location.
        seed_aa = seed_sequence[pos - 1]                                                                        # Identifies the original amino acid that was replaced.
        # Grants a boolean logic reward if the model smartly mutated this position to the family's natural consensus (broader protein family).
        family_bonus = 1.0 if candidate_aa == feature["consensus_aa"] and candidate_aa != seed_aa else 0.0
        # Scales an additional mathematical reward based on how dominant that consensus amino acid is.
        support_bonus = float(feature["consensus_freq"] * (candidate_aa == feature["consensus_aa"]))            
        # Calculates and records the heavily weighted composite score for this specific mutation.
        scores.append(float(0.55 * feature["functional_weight"] + 0.20 * feature["family_mutability"] + 0.15 * family_bonus + 0.10 * support_bonus)) 
    return float(np.mean(scores)) if scores else -1.0                                                           # Averages all recorded mutation scores, returning a failure flag if no valid scores were computed.


def embed_sequences(sequences: list[str], model_name: str, batch_size: int = 4, max_length: int = 2048) -> np.ndarray:
    """Embed protein sequences with an ESM model using masked-mean pooling over token states."""
    from transformers import AutoModel, AutoTokenizer

    device = "cuda" if torch.cuda.is_available() else "cpu"                                                     # Dynamically selects GPU acceleration if hardware supports it, falling back to CPU.
    tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)                                  # Loads the HuggingFace tokenizer associated with the specified ESM model without altering case.
    model = AutoModel.from_pretrained(model_name)                                                               # Loads the heavy HuggingFace neural network weights into memory.
    model.to(device)                                                                                            # Transfers the massive model weights to the selected hardware processing unit.
    model.eval()                                                                                                # Locks the model into inference mode, disabling dropout and gradient tracking to boost speed.

    embeddings = []                                                                                             # Initializes a master list to collect the high-dimensional mathematical representations.
    with torch.no_grad():                                                                                       # Temporarily turns off PyTorch's gradient calculation engine to save massive amounts of RAM.
        for start in range(0, len(sequences), batch_size):                                                      # Iterates over the massive sequence list in safe, manageable chunks defined by batch_size.
            batch = sequences[start : start + batch_size]                                                       # Slices out the current subset of sequences for processing.
            toks = tokenizer(batch, padding=True, truncation=True, max_length=max_length, return_tensors="pt")  # Tokenizes the strings, adding padding so matrices align, and formats as PyTorch tensors.
            toks = {key: value.to(device) for key, value in toks.items()}                                       # Transfers all the input tensor data dictionaries onto the compute hardware.
            hidden = model(**toks).last_hidden_state                                                            # Runs the neural network forward pass and extracts the raw, multi-dimensional feature states.
            mask = toks["attention_mask"].unsqueeze(-1)                                                         # Isolates the padding mask and adds a dimension so it perfectly aligns with the hidden features.
            pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1)                                  # Collapses the sequence length dimension by averaging the unpadded tokens to create one vector per sequence.
            embeddings.append(pooled.cpu().numpy())                                                             # Moves the resulting mathematical vector back to system RAM, converts to NumPy, and stores it.
    return np.vstack(embeddings).astype(np.float32)                                                             # Stacks all batch lists into one massive continuous 2D NumPy array matrix of floats.


def greedy_diverse_order(embeddings: np.ndarray, base_scores: np.ndarray, penalty_weight: float = 0.25, preferred_mask: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Produce a diversity-aware greedy order and the associated penalties."""
    n = len(base_scores)                                                                                        # Determines the total number of candidate items to be ranked.
    if n == 0:                                                                                                  # Checks if the input lists are completely empty.
        return np.array([], dtype=int), np.array([], dtype=np.float32)                                          # Bails out early, returning safely formatted empty arrays to prevent downstream math errors.
    normed = normalize_rows(embeddings)                                                                         # Scales all structural embedding vectors so cosine similarity calculations are fast and accurate.
    remaining = list(range(n))                                                                                  # Initializes an index tracking list of candidates that have not yet been placed in the final rank.
    order = []                                                                                                  # Initializes the master list that will hold the final sorted rankings.
    penalties = np.zeros(n, dtype=np.float32)                                                                   # Initializes a zeroed tracking array to record the structural redundancy penalty hit for each item.
    # Formats the boolean preference override mask, defaulting to all-true if none was provided.
    preferred_mask = np.ones(n, dtype=bool) if preferred_mask is None else np.asarray(preferred_mask, dtype=bool) 

    # Extracts the raw indices of all candidates explicitly marked as highly preferred.
    first_pool = np.where(preferred_mask)[0]                                                                    
    # Finds the absolute highest scoring candidate, restricted to the preferred pool if possible.
    first = int(first_pool[np.argmax(base_scores[first_pool])]) if len(first_pool) else int(np.argmax(base_scores)) 
    # Anchors the ranking list by placing the undisputed top candidate in the number one spot.
    order.append(first)  
    # Deletes the anchor candidate from the pending pool so it cannot be selected again.                                                                                       
    remaining.remove(first)                                                                                     

    while remaining:                                                                                            # Begins the iterative selection loop, running until every candidate has been assigned a rank.
        best_idx = None                                                                                         # Initializes tracker for the winner of the current iteration round.
        best_value = -math.inf                                                                                  # Initializes tracker for the highest penalized score observed in the current round.
        for idx in remaining:                                                                                   # Scans through every candidate still waiting in the pending pool.
            nearest = max(float(np.dot(normed[idx], normed[chosen])) for chosen in order)                       # Calculates how structurally similar the candidate is to the most similar item already placed in the rankings.
            preference_bonus = 0.02 if preferred_mask[idx] else 0.0                                             # Applies a tiny artificial score boost if the item was marked by the external preference mask.
            penalized = float(base_scores[idx] - penalty_weight * nearest + preference_bonus)                   # Computes the true ranking score: base quality minus a penalty for looking too much like higher-ranked items.
            if penalized > best_value:                                                                          # Checks if this candidate's adjusted score beats the current round leader.
                best_value = penalized                                                                          # Updates the highest penalized score tracker.
                best_idx = idx                                                                                  # Assigns the current candidate as the new round leader.
                penalties[idx] = nearest                                                                        # Records the specific redundancy penalty that this candidate absorbed.
        order.append(int(best_idx))                                                                             # Locks the round winner into the next available spot in the final rankings.
        remaining.remove(int(best_idx))                                                                         # Deletes the newly ranked winner from the pending pool.
    return np.asarray(order, dtype=int), penalties                                                              # Returns the final ordered indices and their associated diversity penalties as NumPy arrays.


def write_fasta(records: list[tuple[str, str]], path: str | Path) -> None:
    """Write a simple FASTA file from (header, sequence) tuples."""
    path = Path(path)                                                                                           # Casts the input target location to a standard Pathlib object for easy directory manipulation.
    path.parent.mkdir(parents=True, exist_ok=True)                                                              # Ensures the folder structure exists for the target file, creating missing folders silently.
    with open(path, "w", encoding="utf-8") as handle:                                                           # Safely opens the target file in write mode using standard UTF-8 encoding.
        for header, sequence in records:                                                                        # Iterates over all provided data pairs.
            handle.write(f">{header}\n{sequence}\n")                                                            # Writes the biological standard FASTA format: a line starting with '>' for the ID, followed by the raw sequence string.