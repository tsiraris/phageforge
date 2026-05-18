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


VALID_AA = set(AMINO_ACIDS)                                                                                                 # Creates a fast lookup set of valid amino acid characters for quick membership testing.


@dataclass                                                                                                                  # Decorator to automatically generate special methods like __init__ and __repr__ for the class.
class EditProposal:                                                                                                         # Defines a data structure to hold substitution rules for a specific position.
    """
    Compact container describing the allowed substitutions for one residue position.
    
    This dataclass acts as a strict biological rulebook for a single amino acid 
    during the generation search. It prevents the AI from hallucinating arbitrary 
    mutations by explicitly defining which amino acids are allowed at this index 
    based on family conservation and target-host probabilities.
    
    Example:
        proposal = EditProposal(
            position=42, 
            seed_aa='A', 
            allowed_aas=['S', 'T'], 
            functional_weight=0.85, ...
        )
    """

    position: int                                                                                                           # Stores the 1-indexed position of the amino acid in the sequence.
    seed_aa: str                                                                                                            # Stores the original amino acid character at this position.
    allowed_aas: list[str]                                                                                                  # Stores a list of valid amino acid characters that can be swapped in.
    target_preference: dict[str, float]                                                                                     # Stores probability weights of amino acids for the target host genus.
    family_preference: dict[str, float]                                                                                     # Stores probability weights of amino acids based on evolutionary family data.
    functional_weight: float                                                                                                # Stores the computed biological importance score of mutating this position.
    conservation_penalty: float                                                                                             # Stores a penalty score based on how conserved (unchanging) the original amino acid is.
    region_name: str                                                                                                        # Stores the name of the contiguous structural window this position belongs to.


@dataclass                                                                                                                  # Decorator to automatically generate class boilerplate methods.
class SearchCandidate:                                                                                                      # Defines a data structure tracking a specific mutant sequence during the generation search.
    """
    Container for one Stage 09 search state and its current metadata.
    
    This dataclass represents a "node" in the beam search tree. It stores the 
    current mutated sequence string, tracks the exact history of edits made to reach 
    this state, and records which algorithmic round produced it. This ensures perfect 
    provenance tracing for every generated protein.
    
    Example:
        candidate = SearchCandidate(
            candidate_sequence="MKAAG...", 
            mutations=["42:A->S"], 
            round_index=2, ...
        )
    """

    candidate_sequence: str                                                                                                 # Stores the full string of the mutated amino acid sequence.
    mutations: list[str]                                                                                                    # Stores a human-readable list of specific changes made (e.g., ['197:W→R']).
    mutated_positions: list[int]                                                                                            # Stores the integer indices where mutations have occurred.
    proposal_trace: list[str]                                                                                               # Stores the history of algorithmic decisions that led to this sequence.
    round_index: int                                                                                                        # Stores the iteration number of the search algorithm when this was generated.
    base_parent_id: str                                                                                                     # Stores a unique identifier mapping back to the original unmodified seed.


# ----------------------------- JSON and path helpers ----------------------------- #


def read_json(path: str | Path) -> dict:
    """
    Reads a JSON file from disk and parses it into a Python dictionary.
    
    It opens the file at the specified path using UTF-8 encoding to prevent 
    character corruption, parses the JSON string structure using the standard 
    library, and returns the resulting dictionary object.
    
    Example:
        config = read_json("stage07_context.json")
        -> {"target_host": "Acinetobacter", "stage": "07", ...}
    """
    with open(path, "r", encoding="utf-8") as handle:                                                                       # Opens the file at the given path in read mode with standard UTF-8 encoding.
        return json.load(handle)                                                                                            # Parses the file stream into a Python dictionary and returns it.



def write_json(obj: dict, path: str | Path) -> None:
    """
    Serializes a Python dictionary into a JSON file on disk.
    
    It converts the path to a robust pathlib object, aggressively creates any 
    missing parent directories to prevent "folder not found" crashes, and writes 
    the dictionary to disk with a human-readable 2-space indentation.
    
    Example:
        write_json({"status": "success"}, "results/output.json")
        -> (Creates the file with properly indented JSON text)
    """
    path = Path(path)                                                                                                       # Converts the string path into a robust pathlib.Path object.
    path.parent.mkdir(parents=True, exist_ok=True)                                                                          # Creates necessary parent directories if they don't exist, suppressing errors if they do.
    with open(path, "w", encoding="utf-8") as handle:                                                                       # Opens the destination file in write mode with standard UTF-8 encoding.
        json.dump(obj, handle, indent=2)                                                                                    # Serializes the Python dictionary to JSON format with 2-space indentation for readability.


# ----------------------------- Sequence metrics and parsing ----------------------------- #


def mutation_list(seed_sequence: str, candidate_sequence: str) -> list[str]:
    """
    Calculates a human-readable list of point mutations between a seed and a candidate.
    
    It zips the two sequence strings together side-by-side. It iterates through 
    the pairs, tracking the 1-based biological index. Whenever the characters differ, 
    it formats a standard biological mutation string and appends it to the output list.
    
    Example:
        mutation_list("MKAAG", "MKASG")
        -> ["4:A→S"]
    """
    out: list[str] = []                                                                                                     # Initializes an empty list to store the formatted mutation strings.
    for idx, (seed_aa, cand_aa) in enumerate(zip(seed_sequence, candidate_sequence), start=1):                              # Iterates through both sequences simultaneously, starting the index at 1.
        if seed_aa != cand_aa:                                                                                              # Checks if the amino acids differ at the current position.
            out.append(f"{idx}:{seed_aa}→{cand_aa}")                                                                        # Formats and appends the human-readable mutation string to the output list.
    return out                                                                                                              # Returns the populated list of mutation strings.



def parse_mutation_positions(mutation_text: str | Iterable[str] | None) -> list[int]:
    """
    Extracts raw integer positions from a formatted mutation annotation string/list.
    
    It normalizes the input into a list of strings, splitting by semicolons if 
    necessary. It iterates through the tokens, safely splits them at the colon, 
    casts the prefix to an integer, and returns a sorted, deduplicated list of 
    those integers.
    
    Example:
        parse_mutation_positions("42:A->C; 105:T->G")
        -> [42, 105]
    """
    if mutation_text is None:                                                                                               # Checks if the input is completely empty (None).
        return []                                                                                                           # Returns an empty list safely without crashing.
    if isinstance(mutation_text, str):                                                                                      # Checks if the input is a single concatenated string.
        tokens = [tok.strip() for tok in mutation_text.split(";") if tok.strip()]                                           # Splits by semicolon and removes whitespace, discarding any empty tokens.
    else:                                                                                                                   # Handles the case where input is already an iterable (like a list).
        tokens = [str(tok).strip() for tok in mutation_text if str(tok).strip()]                                            # Converts each item to string, strips whitespace, and discards empty tokens.
    positions: list[int] = []                                                                                               # Initializes a list to hold the parsed integer positions.
    for token in tokens:                                                                                                    # Iterates through the cleaned mutation string segments.
        if ":" not in token:                                                                                                # Skips tokens that lack the expected colon separator.
            continue                                                                                                        # Moves to the next token.
        try:                                                                                                                # Starts an error-handling block in case integer conversion fails.
            positions.append(int(token.split(":", 1)[0]))                                                                   # Splits by colon, takes the first part (the number), converts to int, and appends.
        except ValueError:                                                                                                  # Catches errors if the part before the colon wasn't a valid number.
            continue                                                                                                        # Fails silently and moves to the next token.
    return sorted(set(positions))                                                                                           # Removes duplicate positions via set(), sorts them numerically, and returns.



def normalized_shannon_entropy(sequence: str) -> float:
    """
    Calculates the Shannon entropy of an amino acid sequence, normalized to [0, 1].
    
    It converts the sequence into a pandas Series to compute the fractional 
    abundance of each character. It applies the Shannon entropy formula 
    (-sum(p * log2(p))) and divides by the theoretical maximum entropy (log2(20)) 
    to output a percentage score representing the sequence's amino acid diversity.
    
    Example:
        normalized_shannon_entropy("AAAAAAAAAA") 
        -> 0.0 (Zero diversity / Degenerate)
    """
    if not sequence:                                                                                                        # Returns 0 immediately if the sequence string is empty.
        return 0.0                                                                                                          # Returns baseline zero entropy.
    counts = pd.Series(list(sequence)).value_counts(normalize=True)                                                         # Converts string to list, then to pandas Series to calculate relative frequencies of each character.
    entropy = float(-(counts * np.log2(counts)).sum())                                                                      # Calculates standard Shannon entropy using the formula: -sum(p * log2(p)).
    return entropy / math.log2(len(VALID_AA))                                                                               # Normalizes the entropy by dividing by the maximum possible entropy (log2 of alphabet size).



def max_residue_fraction(sequence: str) -> float:
    """
    Identifies the single most abundant amino acid and returns its frequency fraction.
    
    It tallies all characters in the sequence, normalizes the counts into 
    percentages, and returns the absolute maximum percentage value. This is 
    a fast diagnostic to detect AI models over-enriching a sequence with one residue.
    
    Example:
        max_residue_fraction("AABBCCC") 
        -> 0.428 (Because 'C' appears 3 times out of 7)
    """
    if not sequence:                                                                                                        # Returns 1.0 safely if the sequence is empty to avoid math errors.
        return 1.0                                                                                                          # Returns maximum fraction fallback.
    counts = pd.Series(list(sequence)).value_counts(normalize=True)                                                         # Calculates the fractional abundance of every unique amino acid in the sequence.
    return float(counts.max())                                                                                              # Extracts and returns the fraction of the single most abundant residue.



def longest_homopolymer_run(sequence: str) -> int:
    """
    Finds the length of the longest contiguous homopolymer run in the sequence.
    
    It iterates through the string character by character. If the current character 
    matches the previous one, it increments a counter. If it differs, the counter 
    resets. It returns the highest recorded value. This detects pathological 
    "spaghetti" loops hallucinated by generative models.
    
    Example:
        longest_homopolymer_run("MKAAGGGGTAY") 
        -> 4 (Because 'G' repeats 4 times consecutively)
    """
    if not sequence:                                                                                                        # Returns 0 immediately if sequence is empty.
        return 0                                                                                                            # Returns 0 length.
    best = run = 1                                                                                                          # Initializes tracking variables for the overall longest run and current active run.
    prev = sequence[0]                                                                                                      # Sets the initial 'previous' character to the first character of the sequence.
    for aa in sequence[1:]:                                                                                                 # Iterates through the sequence starting from the second character.
        if aa == prev:                                                                                                      # Checks if the current character matches the immediately preceding one.
            run += 1                                                                                                        # Increments the current run counter.
            best = max(best, run)                                                                                           # Updates the best run counter if the current run exceeds it.
        else:                                                                                                               # Triggers if the character sequence breaks the repetition.
            prev = aa                                                                                                       # Updates the previous character tracker to the new character.
            run = 1                                                                                                         # Resets the current run counter back to 1.
    return best                                                                                                             # Returns the highest recorded consecutive identical amino acid count.



def low_complexity_fraction(sequence: str, k: int = 12, unique_threshold: int = 3) -> float:
    """
    Measures the percentage of the protein made up of low-complexity motifs.
    
    It slides a window of size `k` across the sequence. If a window contains 
    fewer unique amino acids than the `unique_threshold`, it flags the window 
    as "low complexity". It returns the fraction of total windows that were flagged.
    
    Example:
        low_complexity_fraction("GGSGGSGGSGGS") 
        -> 1.0 (100% of windows fail the threshold, highly repetitive)
    """
    if len(sequence) < k:                                                                                                   # Checks if the total sequence is shorter than the scanning window size.
        return float(len(set(sequence)) <= unique_threshold)                                                                # Evaluates the whole sequence at once, returning 1.0 if it lacks diversity, 0.0 otherwise.
    flags = []                                                                                                              # Initializes a list to hold boolean results for each sliding window.
    for i in range(len(sequence) - k + 1):                                                                                  # Iterates over every possible starting index for a window of size k.
        flags.append(len(set(sequence[i : i + k])) <= unique_threshold)                                                     # Slices the window, counts unique chars via set(), and appends True if below threshold.
    return float(np.mean(flags)) if flags else 0.0                                                                          # Calculates the ratio of True (low complexity) windows to total windows.



def mutation_span(positions: Iterable[int]) -> int:
    """
    Calculates the total linear sequence distance covered by a set of mutations.
    
    It deduplicates and sorts the integer positions from lowest to highest, 
    then subtracts the first position from the last position. This determines 
    how widely dispersed the AI's edits are across the chassis.
    
    Example:
        mutation_span([10, 12, 15, 25])
        -> 15 (Because 25 - 10 = 15)
    """
    pos = sorted(set(int(p) for p in positions))                                                                            # Casts all inputs to integers, removes duplicates, and sorts them sequentially.
    if len(pos) < 2:                                                                                                        # Checks if there are fewer than 2 mutations, meaning no span exists.
        return 0                                                                                                            # Returns 0 distance.
    return pos[-1] - pos[0]                                                                                                 # Calculates the sequence distance between the highest and lowest position indices.



def outside_editable_fraction(positions: Iterable[int], editable_positions: set[int]) -> float:
    """
    Calculates the fraction of mutations that occurred outside authorized hotspots.
    
    It iterates through the provided mutation positions and checks if they exist 
    within the pre-approved `editable_positions` set. It calculates the average 
    of these boolean checks to output the violation percentage.
    
    Example:
        outside_editable_fraction([10, 20], {10, 15})
        -> 0.5 (Because position 20 is illegal, representing 50% of the edits)
    """
    pos = [int(p) for p in positions]                                                                                       # Converts all mutated position inputs into integers.
    if not pos:                                                                                                             # Checks if the mutation list is empty.
        return 0.0                                                                                                          # Returns 0 fraction as there are no violations.
    return float(np.mean([p not in editable_positions for p in pos]))                                                       # Checks each mutation against the allowed set, returning the average of the booleans.



def sequence_identity(seq_a: str, seq_b: str) -> float:
    """
    Calculates the exact positional sequence identity fraction between two strings.
    
    It validates that both strings are the same length, pairs their characters 
    using `zip()`, counts the exact positional matches, and divides by the total 
    length to return a similarity percentage.
    
    Example:
        sequence_identity("ACTG", "ACTA")
        -> 0.75 (Because 3 out of 4 positions match exactly)
    """
    if not seq_a or not seq_b or len(seq_a) != len(seq_b):                                                                  # Validates that both sequences exist and are exactly the same length.
        return 0.0                                                                                                          # Returns 0 similarity if the validation fails.
    return float(sum(a == b for a, b in zip(seq_a, seq_b)) / len(seq_a))                                                    # Pairs up characters, counts exact matches, and divides by total length.



def build_basic_sequence_features(seed_sequence: str, candidate_sequence: str, editable_positions: set[int]) -> dict[str, float | int]:
    """
    Compiles a comprehensive dictionary of lightweight biological sequence metrics.
    
    It runs the candidate sequence through all the fast diagnostic functions 
    (entropy, homopolymers, mutation counts) and bundles the results into a 
    single dictionary. This is heavily used by the Structural Surrogate to assess risk.
    
    Example:
        build_basic_sequence_features("MK...", "MK...", {10,11})
        -> {"mutation_count": 2, "normalized_entropy": 0.85, ...}
    """
    mutations = mutation_list(seed_sequence, candidate_sequence)                                                            # Computes the list of specific substitution strings.
    positions = parse_mutation_positions(mutations)                                                                         # Extracts the integer locations of those substitutions.
    return {                                                                                                                # Begins constructing a dictionary of calculated metrics.
        "mutation_count": len(positions),                                                                                   # Stores the raw number of mutations made.
        "mutation_span": mutation_span(positions),                                                                          # Stores the sequence distance between the first and last mutation.
        "normalized_entropy": normalized_shannon_entropy(candidate_sequence),                                               # Stores the overall diversity score of the candidate string.
        "max_single_residue_fraction": max_residue_fraction(candidate_sequence),                                            # Stores the frequency of the most overrepresented amino acid.
        "longest_homopolymer_run": longest_homopolymer_run(candidate_sequence),                                             # Stores the length of the longest repeating amino acid stretch.
        "low_complexity_fraction": low_complexity_fraction(candidate_sequence),                                             # Stores the ratio of sequence segments lacking residue diversity.
        "outside_editable_fraction": outside_editable_fraction(positions, editable_positions),                              # Stores the ratio of illegal mutations to total mutations.
        "seed_identity": sequence_identity(seed_sequence, candidate_sequence),                                              # Stores the global similarity score against the parent sequence.
    }                                                                                                                       # Closes the dictionary of metrics.
            

def infer_target_residue_preferences(strict_df: pd.DataFrame, target_host: str, positions_1based: list[int]) -> dict[int, dict[str, float]]:
    """
    Calculates evolutionary substitution priors specific to the target host bacteria.
    
    It filters the strict protein database for rows matching the target host. For 
    every specified editable position, it counts which amino acids naturally occur 
    there across those filtered rows, returning a dictionary of probability weights.
    
    Example:
        infer_target_residue_preferences(df, "Acinetobacter", [42])
        -> {42: {'A': 0.80, 'S': 0.20}}
    """
    target_rows = strict_df.loc[strict_df["host_genus"].astype(str) == str(target_host)].copy()                             # Filters the dataframe to retain only rows where the bacteria matches the target host.
    if target_rows.empty:                                                                                                   # Checks if the filtered dataframe contains no data.
        return {int(pos): {} for pos in positions_1based}                                                                   # Returns an empty dictionary for every requested position.

    preferences: dict[int, dict[str, float]] = {}                                                                           # Initializes the master dictionary to store positional probabilities.
    for pos in positions_1based:                                                                                            # Iterates through each biologically requested index.
        residues = []                                                                                                       # Initializes a list to collect all amino acids found at this position.
        idx = int(pos) - 1                                                                                                  # Converts the 1-based biological position to a 0-based list index.
        for seq in target_rows["aa_sequence"].astype(str):                                                                  # Iterates through every sequence string in the filtered dataset.
            if idx < len(seq):                                                                                              # Ensures the sequence is actually long enough to have this position.
                aa = seq[idx]                                                                                               # Extracts the amino acid character at the target index.
                if aa in VALID_AA:                                                                                          # Verifies the character is a standard biological amino acid.
                    residues.append(aa)                                                                                     # Adds the valid character to our collection list.
        if not residues:                                                                                                    # Checks if no valid amino acids were found at this position across all rows.
            preferences[int(pos)] = {}                                                                                      # Assigns an empty dictionary since there is no data.
            continue                                                                                                        # Skips to the next position.
        counts = pd.Series(residues).value_counts(normalize=True)                                                           # Converts the list to a Series and calculates the fractional frequency of each residue.
        preferences[int(pos)] = {str(aa): float(freq) for aa, freq in counts.items()}                                       # Maps the fractional frequencies back to standard string/float pairs in the dictionary.
    return preferences                                                                                                      # Returns the fully populated preferences mapping.



def build_edit_proposals_from_context(context: dict, strict_df: pd.DataFrame | None = None) -> list[EditProposal]:
    """
    Constructs the master list of `EditProposal` rules for the generation engine.
    
    It extracts the functional hotspots from the design context. For each hotspot, 
    it calculates the evolutionary family frequencies and the specific target-host 
    frequencies. It merges these to generate a strict, biologically verified list 
    of `allowed_aas` for that position, returning the finalized dataclass objects.
    
    Example:
        build_edit_proposals_from_context(context, strict_df)
        -> [EditProposal(position=42, allowed_aas=['A', 'S'], ...), ...]
    """
    # Extracts the seed sequence, target host, position features (i.e. seed_frequency, family entropy etc), and structured windows from the context.json of Stage 07.
    seed_sequence = str(context["selected_seed"]["seed_sequence"])                                                          # Extracts the parental amino acid sequence from the JSON context.
    target_host = str(context["target_host"])                                                                               # Extracts the name of the bacteria the project aims to infect.
    position_features = list(context["editable_region"].get("position_features", []))                                       # Extracts the pre-calculated importance scores for allowable positions. 
    structured_windows = list(context["editable_region"].get("structured_windows", []))                                     # Extracts grouped position windows to maintain 3D structural logic.
    
    # Mapping residue positions to their parent structural window as defined in the context.json, with overlapping positions are mapped to the strongest parent window.
    # For example: position_to_region = {391: "window_391_414", ..., 405: "window_391_414", ..., 415: "window_404_427"}
    position_to_region: dict[int, str] = {}                                                                                 # Initializes a lookup mapping individual residue indices to their parent structural window.
    for window in structured_windows:                                                                                       # Iterates through the predefined structural groupings.
        region_name = str(window.get("name", f"window_{window.get('window_start', 'na')}_{window.get('window_end', 'na')}")) # Constructs a distinct label for the region, falling back to indices if unnamed.
        for pos in window.get("positions", []):                                                                             # Iterates through every specific integer position contained in this window.
            pos = int(pos)                                                                                                  # Casts the position string/number to an integer for strict mapping.
            if pos not in position_to_region:                                                                               # Prevents overwriting if a position was somehow claimed by a stronger earlier window.
                position_to_region[pos] = region_name                                                                       # Assigns the region label to the position index.

    # Recovers the evolutionary family sequences from the context.json and converts them into a list.
    family_rows = pd.DataFrame(context["family_context"].get("family_rows", []))                                            # Reconstructs the dataframe of evolutionarily related sequences from the context.
    family_sequences = family_rows.get("aa_sequence", pd.Series(dtype=str)).astype(str).tolist() if not family_rows.empty else [] # Extracts the column of sequence strings into a flat Python list.
    # Creates a simple (1-based numbered) list of all the integer positions allowed to be edited.
    positions_1based = [int(row["position"]) for row in position_features]                                                  # Creates a clean integer list of all targetable biological positions.
    # Scans the Strict Dataframe for viruses infecting the specific target host, and returns based on them a dictionary of probabilities of which amino acids are used in each of the allowed to be mutated positions (e.g., {42: {'S': 0.8, 'T': 0.2}}).
    target_preferences = infer_target_residue_preferences(strict_df, target_host, positions_1based) if strict_df is not None else {int(pos): {} for pos in positions_1based} # Fetches host-specific biases if database provided.
    
    # Initializes the final output list, looping over every allowed position.
    proposals: list[EditProposal] = []                                                                                      # Initializes the master list to hold the generated EditProposal objects.
    for row in position_features:                                                                                           # Loops over every allowable position evaluated by the prior context script.
        # Extracts the position integer, translates it to 0-based.
        pos = int(row["position"])                                                                                          # Extracts the 1-based biological coordinate.
        idx = pos - 1                                                                                                       # Translates the coordinate to a 0-based Python string index.
        # What is the wild-type amino acid currently at this position of the seed sequence.
        seed_aa = seed_sequence[idx]                                                                                        # Identifies what the wild-type amino acid currently is at this spot.

        # Looks at every single evolutionary family sequence and grabs the amino acid at that exact position. 
        # If it's a valid character (not a gap or junk), it adds it to the residues list.
        residues = []                                                                                                       # Initializes a temporary list to track amino acids at this position across the family.
        for seq in family_sequences:                                                                                        # Iterates through every evolutionarily related sequence.
            if idx < len(seq):                                                                                              # Confirms the current sequence isn't too short to possess this index.
                aa = seq[idx]                                                                                               # Extracts the specific character.
                if aa in VALID_AA:                                                                                          # Confirms the character isn't a placeholder, gap, or invalid symbol.
                    residues.append(aa)                                                                                     # Keeps the valid amino acid for frequency counting.
        # Calculates the frequency of each amino acid in the family, and converts them to a dictionary of probabilities for each. 
        family_counts = pd.Series(residues).value_counts(normalize=True) if residues else pd.Series(dtype=float)            # Calculates the fractional distribution of amino acids based on evolution.
        family_preference = {str(aa): float(freq) for aa, freq in family_counts.items()}                                    # Casts the pandas distribution into a native Python dictionary.

        # Creates an "allowed" mathematical set (no duplicates) by taking the top 4 most common amino acids 
        # from the evolutionary family and for the target host, and adds them to it.
        allowed = set()                                                                                                     # Starts an empty set to ensure uniqueness of allowed replacement candidates.
        allowed.update(list(family_preference.keys())[:4])                                                                  # Adds up to the top 4 most common amino acids found in the general protein family.
        allowed.update(list(target_preferences.get(pos, {}).keys())[:4])                                                    # Adds up to the top 4 most common amino acids found specifically infecting the target host.
        # Removes the amino acid of the seed sequence from the allowed set, and performs a final safety sweep on the set.
        allowed.discard(seed_aa)                                                                                            # Removes the original amino acid since mutating to the same thing does nothing.
        allowed = {aa for aa in allowed if aa in VALID_AA}                                                                  # Performs a final safety purge to ensure no garbage characters sneaked in.

        # If there are no plausible mutations remaining, abandon this position entirely and skip to the next loop iteration.
        if not allowed:                                                                                                     # Checks if no viable mutations could be identified.
            continue                                                                                                        # Skips generating a proposal for this rigid position.
        
        # If valid mutations exist, it packages all this math into the EditProposal object.
        proposals.append(                                                                                                   # Begins instantiating and saving the proposal object.
            EditProposal(                                                                                                   # Calls the dataclass constructor.
                position=pos,                                                                                               # Injects the biological index.
                seed_aa=seed_aa,                                                                                            # Injects the starting amino acid state.
                allowed_aas=sorted(allowed),                                                                                # Alphabetizes and injects the curated list of possible mutations.
                target_preference=target_preferences.get(pos, {}),                                                          # Attaches the host-specific probability weights.
                family_preference=family_preference,                                                                        # Attaches the general evolutionary probability weights.
                functional_weight=float(row.get("functional_weight", 0.0)),                                                 # Attaches the pre-calculated importance of editing this spot.
                conservation_penalty=float(row.get("seed_freq", 0.0)),                                                      # Attaches the evolutionary resistance to changing this spot.
                region_name=position_to_region.get(pos, "ungrouped"),                                                       # Attaches the structural tag, defaulting if unmapped.
            )                                                                                                               # Closes the constructor.
        )                                                                                                                   # Closes the append call.
    return proposals                                                                                                        # Returns the complete list of position editing rules.



def choose_editable_positions(proposals: list[EditProposal], max_positions: int, seed: int) -> list[int]:
    """
    Selects a highly optimized subset of editable positions to define the mutation space.
    
    To prevent geometric collapse, it groups the proposals by their structural regions. 
    It greedily picks the most biologically important position from each region to ensure 
    mutation diversity across the scaffold, then fills the remaining budget with the 
    highest priority global positions.
    
    Example:
        choose_editable_positions(all_proposals, max_positions=12, seed=42)
        -> [10, 15, 42, 105...] (A curated list of 12 highly impactful positions)
    """
    if len(proposals) <= max_positions:                                                                                     # Checks if the total available proposals is under the hardware generation limit.
        return [int(item.position) for item in proposals]                                                                   # Just returns all available positions if limits aren't breached.

    # Initialize the random number generator (for reproducibility)
    rng = random.Random(seed)                                                                                               # Instantiates an isolated random generator scoped only to this execution using the provided seed.
    # Create an empty dictionary to group proposals by region
    grouped: dict[str, list[EditProposal]] = {}                                                                             # Prepares a dictionary to cluster the proposals based on their 3D structural regions.
    # Iterates through all the proposals and groups them into buckets (lists) based on their region_name ("structured_window" name for mutations in the context.json, and "ungrouped" for the others).
    # For example: grouped = {'window_391_414': [P1, P2, P3, P4], '"window_323_352"': [P5, P6], ...}
    for item in proposals:                                                                                                  # Iterates through all valid edit rules.
        grouped.setdefault(item.region_name, []).append(item)                                                               # Places each rule into the list belonging to its specific region.
    # Sort the proposals within each region from best to worst primarily by functional_weight (host-binding importance), secondarily by conservation_penalty, and finally by position.
    for items in grouped.values():                                                                                          # Iterates over each region's specific list of proposals.
        items.sort(key=lambda x: (x.functional_weight, -x.conservation_penalty, x.position), reverse=True)                  # Sorts the region's proposals by highest biological importance and lowest risk.

    chosen: list[int] = []                                                                                                  # Initializes the final tracking list for selected integers.

    # Sorts the buckets based on which one contains the absolute highest-scoring proposal. 
    # Then, it takes exactly the #1 best proposal from each bucket in their sorted order.
    for region_name in sorted(grouped, key=lambda name: max(item.functional_weight for item in grouped[name]), reverse=True): # Loops through regions, prioritizing regions that contain highly important spots.
        chosen.append(int(grouped[region_name][0].position))                                                                # Immediately grabs the absolute best position from the current region.
        if len(chosen) >= max_positions:                                                                                    # Verifies if the quota was hit just by picking the top spots.
            return sorted(chosen[:max_positions])                                                                           # Truncates and sorts if quota is filled.

    # If one proposal has been selected from every region but not yet filled the mutation budget (max_positions), 
    # It gathers all the leftover proposals that weren't picked in the regional draft (remaining), sorts this entire leftover pile strictly by functional_weight.
    remaining = [item for item in proposals if int(item.position) not in chosen]
    remaining.sort(key=lambda x: (x.functional_weight, -x.conservation_penalty, x.position), reverse=True)
    
    slots_to_fill = max_positions - len(chosen)
    
    if slots_to_fill > 0:
        # Take half of the remaining slots greedily (best functional_weight)...
        greedy_slots = slots_to_fill // 2
        chosen.extend(int(item.position) for item in remaining[:greedy_slots])
        
        # ...and randomly sample the rest from a slightly wider pool (e.g., the next 20 best options)
        random_slots = slots_to_fill - greedy_slots
        pool_to_shuffle = remaining[greedy_slots :]
        rng.shuffle(pool_to_shuffle)
        
        chosen.extend(int(item.position) for item in pool_to_shuffle[:random_slots])

    return sorted(chosen[:max_positions])                                                                                   # Cleans up the array into numerical order and ensures the length cap is respected.


def substitution_priority(item: EditProposal, aa: str) -> float:
    """
    Computes a mathematical priority score for substituting a specific amino acid.
    
    It blends three factors using weighted addition: the frequency of the amino 
    acid in the target host bacteria (50%), the frequency in the broad evolutionary 
    family (30%), and the general biological importance of mutating that position (20%).
    
    Example:
        substitution_priority(proposal, "A")
        -> 0.82 (High priority: this amino acid is highly prevalent in the target host)
    """
    target = float(item.target_preference.get(aa, 0.0))                                                                     # Extracts the historical frequency of this specific amino acid in the target bacteria.
    family = float(item.family_preference.get(aa, 0.0))                                                                     # Extracts the historical frequency of this specific amino acid in general evolution.
    return float(0.50 * target + 0.30 * family + 0.20 * item.functional_weight)                                             # Blends the metrics together using a tuned weighting formula to output a final priority score.


# ----------------------------- Predictor and surrogate helpers ----------------------------- #


def _patch_logistic_regression_compat(model) -> object:
    """
    Restores missing metadata to scikit-learn models saved in older versions.
    
    When a serialized `.joblib` model is loaded into a newer scikit-learn environment, 
    internal attributes (like `n_features_in_`) might be missing. This function detects 
    those missing attributes and recalculates/injects them to prevent inference crashes.
    
    Example:
        model = _patch_logistic_regression_compat(loaded_model)
        -> (Returns the model safely patched for modern inference)
    """
    if not isinstance(model, LogisticRegression):                                                                           # Verifies the loaded object is actually a LogisticRegression model.
        return model                                                                                                        # Skips patching if it's a different type of object.
    if not hasattr(model, "multi_class"):                                                                                   # Checks if the multi_class attribute was lost during a scikit-learn version upgrade.
        model.multi_class = "auto"                                                                                          # Reinstates the default behavior manually.
    if not hasattr(model, "n_features_in_") and hasattr(model, "coef_"):                                                    # Checks if the feature count attribute was lost but coefficients are present.
        model.n_features_in_ = int(model.coef_.shape[1])                                                                    # Derives the missing feature count directly from the shape of the coefficient matrix.
    return model                                                                                                            # Returns the safely patched model ready for inference.



def load_target_predictor(model_path: str | Path, label_classes_path: str | Path | None) -> tuple[object, list[str]]:
    """
    Deserializes the Target Host Predictor model and its class mapping.
    
    It loads the `.joblib` model from disk, passes it through the version patcher, 
    and loads the corresponding JSON file containing the output class labels (e.g., 
    which column index maps to "Acinetobacter").
    
    Example:
        model, classes = load_target_predictor("model.joblib", "labels.json")
        -> (LogisticRegressionObject, ["Acinetobacter", "Klebsiella", ...])
    """
    # Deserialize the binary model file and immediately passes it through the version patcher.
    model = _patch_logistic_regression_compat(joblib.load(model_path))                                                      
    # If an explicit class order file was not provided, tries to rip the class strings directly out of the model's internal metadata.
    if label_classes_path is None:                                                                                          
        classes = [str(x) for x in getattr(model, "classes_", [])]                                                          
    # If an explicit class order .json file was provided, attempts to load it, and converts the labels to strings.
    else:                                                                                                                   
        with open(label_classes_path, "r", encoding="utf-8") as handle:                                                     # Opens the external class mapping JSON securely.
            classes = [str(x) for x in json.load(handle)]                                                                   # Loads and converts the labels to strings.
    if not classes:                                                                                                         # Validates that the class list was successfully generated via either method.
        raise ValueError("Could not recover label classes for the target predictor.")                                       # Crashes the script safely since probability routing is impossible without labels.
    return model, classes                                                                                                   # Returns the ready model and its corresponding output map.



def predict_target_probability(model, label_classes: list[str], target_host: str, embeddings: np.ndarray) -> np.ndarray:
    """
    Calculates the probability that a protein sequence targets a specific host bacteria.
    
    It accepts a batch of heavy ESM-2 embeddings, passes them through the loaded 
    Logistic Regression model, and slices out only the probability column that corresponds 
    to the designated `target_host`.
    
    Example:
        predict_target_probability(model, classes, "Acinetobacter", embeddings_array)
        -> np.array([0.98, 0.45, 0.88...])
    """
    if target_host not in label_classes:                                                                                    # Safeguards against asking the model to predict a bacteria it wasn't trained on.
        raise ValueError(f"Target host '{target_host}' not found in predictor label order: {label_classes}")                # Throws a descriptive error showing the user what the model actually knows.
    embeddings = np.asarray(embeddings, dtype=np.float32)                                                                   # Casts the incoming feature matrix to a memory-efficient 32-bit numpy array.
    # Verifies the dimensionality of the incoming data matches the model's architecture.
    if hasattr(model, "n_features_in_") and embeddings.shape[1] != int(model.n_features_in_):                               
        raise ValueError(f"Predictor expects {int(model.n_features_in_)} embedding features...")                            # Halts execution to prevent garbage vector inference.
    probs = model.predict_proba(embeddings)                                                                                 # Pushes the embeddings through the model to get a matrix of probabilities for all known hosts.
    target_idx = label_classes.index(target_host)                                                                           # Finds the exact column index corresponding to the bacteria we actually care about.
    return np.asarray(probs[:, target_idx], dtype=np.float32)                                                               # Slices out only the relevant column, returning a 1D array of success probabilities.



def maybe_load_surrogate(path: str | Path | None) -> dict | None:
    """
    Safely attempts to load the Structural Surrogate model bundle from disk.
    
    If the path is None (e.g., the user opted to run without a trained surrogate), 
    it returns None safely. Otherwise, it verifies the file exists and deserializes 
    the Random Forest models.
    
    Example:
        bundle = maybe_load_surrogate("surrogate_model.joblib")
        -> {"pass_model": RandomForestClassifier(...), "plddt_model": ...}
    """
    if path is None:                                                                                                        # Checks if the user opted to bypass surrogate loading.
        return None                                                                                                         # Exits the function cleanly without attempting to load.
    path = Path(path)                                                                                                       # Standardizes the string path to a pathlib object.
    if not path.exists():                                                                                                   # Verifies the file actually exists on the filesystem.
        raise FileNotFoundError(f"Missing surrogate model bundle: {path}")                                                  # Throws a hard error because a specified file is completely missing.
    return joblib.load(path)                                                                                                # Deserializes the full binary dictionary of surrogate machine learning models.



def surrogate_structural_risk(bundle: dict | None, feature_frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Estimates the physical folding risk of a sequence without running expensive ESMFold.
    
    If a trained Random Forest bundle is provided, it extracts the relevant features 
    from the dataframe and predicts the structural stability, mean pLDDT, and RMSD. 
    If no bundle is provided, it falls back to a mathematical heuristic equation 
    that penalizes sequence pathologies (like low complexity or extreme drift).
    
    Example:
        risk, plddt, rmsd = surrogate_structural_risk(bundle, df)
        -> ([0.15], [82.5], [1.8]) (Indicating low risk, good folding confidence, low drift)
    """
    df = feature_frame.copy()                                                                                               # Creates a protective copy of the dataframe to prevent modifying user variables.
    if bundle is None:                                                                                                      # Detects if the pipeline is running in heuristic mode without trained structural proxies.
        risk = (                                                                                                            # Begins calculating a mathematical penalty based strictly on biological rules of thumb.
            0.40 * (1.0 - np.clip(df["seed_identity"].to_numpy(dtype=float), 0.0, 1.0))                                     # Heavily penalizes candidates that have drifted too far from the parent sequence identity.
            + 0.15 * np.clip(df["outside_editable_fraction"].to_numpy(dtype=float), 0.0, 1.0)                               # Adds penalty if mutations occurred outside the allowed structural bounds.
            + 0.15 * np.clip(df["low_complexity_fraction"].to_numpy(dtype=float), 0.0, 1.0)                                 # Adds penalty for repetitive, unstructured garbage regions.
            + 0.10 * np.clip((df["longest_homopolymer_run"].to_numpy(dtype=float) - 4.0) / 4.0, 0.0, 1.0)                   # Adds penalty for long consecutive strings of the exact same amino acid.
            + 0.20 * np.clip(np.abs(df["mutation_count"].to_numpy(dtype=float) - 6.0) / 6.0, 0.0, 1.0)                      # Adds penalty if the number of edits strays too far from the optimal heuristic target of 6.
        )                                                                                                                   # Closes the heuristic risk equation.
        pred_plddt = 90.0 - 35.0 * risk                                                                                     # Estimates AI structure confidence by subtracting a multiplier of risk from an ideal 90 score.
        pred_rmsd = 1.0 + 6.0 * risk                                                                                        # Estimates structural drift by adding a multiplier of risk to an ideal 1.0 baseline.
        return np.asarray(risk, dtype=np.float32), np.asarray(pred_plddt, dtype=np.float32), np.asarray(pred_rmsd, dtype=np.float32) # Casts the math outputs to strict arrays and returns.

    feature_cols = list(bundle["feature_columns"])                                                                          # Extracts the exact column order the trained surrogate models expect.
    X = df.reindex(columns=feature_cols, fill_value=0.0).to_numpy(dtype=np.float32)                                         # Restructures the dataframe to match training, zero-filling missing data, and casts to matrix.

    pass_model = bundle.get("pass_model")                                                                                   # Retrieves the binary classifier predicting if the structure will break entirely.
    plddt_model = bundle.get("plddt_model")                                                                                 # Retrieves the regressor predicting the exact ESMFold confidence score.
    rmsd_model = bundle.get("rmsd_model")                                                                                   # Retrieves the regressor predicting the exact angstrom drift from the parent.

    pass_prob = pass_model.predict_proba(X)[:, 1] if pass_model is not None else np.full(len(df), 0.5, dtype=np.float32)     # Executes probability prediction for structure stability, defaulting to a coin flip if missing.
    pred_plddt = plddt_model.predict(X) if plddt_model is not None else np.full(len(df), 65.0, dtype=np.float32)            # Executes pLDDT inference, defaulting to a mediocre score if missing.
    pred_rmsd = rmsd_model.predict(X) if rmsd_model is not None else np.full(len(df), 3.5, dtype=np.float32)                # Executes RMSD inference, defaulting to a highly suspicious drift if missing.
    risk = 1.0 - np.clip(pass_prob, 0.0, 1.0)                                                                               # Inverts the stability probability to calculate the final structural risk factor.
    return np.asarray(risk, dtype=np.float32), np.asarray(pred_plddt, dtype=np.float32), np.asarray(pred_rmsd, dtype=np.float32) # Packs the three surrogate outputs into numpy arrays and returns.


# ----------------------------- Diversity and selection helpers ----------------------------- #


def greedy_diverse_pick(embeddings: np.ndarray, scores: np.ndarray, top_k: int, penalty_weight: float = 0.25) -> tuple[list[int], np.ndarray]:
    """
    Selects a highly diverse top-k candidate panel to prevent Mode Collapse.
    
    It iteratively selects the highest-scoring sequence, adds it to the "chosen" list, 
    and then mathematically penalizes all remaining candidates based on their cosine 
    similarity (closeness in embedding space) to the already chosen sequences. This 
    forces the algorithm to pick structurally distinct solutions rather than 50 clones 
    of the same local optimum.
    
    Example:
        chosen_indices, penalties = greedy_diverse_pick(embeddings, scores, top_k=5)
        -> ([42, 12, 88, 1, 105], [0.0, 0.12, 0.45, 0.15, 0.80])
    """
    if len(scores) == 0:                                                                                                    # Immediately checks for an empty pool of candidates.
        return [], np.asarray([], dtype=np.float32)                                                                         # Returns empty structures to prevent processing crashes downstream.
    order: list[int] = []                                                                                                   # Initializes the final ranked list of chosen candidate indices.
    penalties = np.zeros(len(scores), dtype=np.float32)                                                                     # Creates an empty array tracking how much mathematical penalty was applied to each item.
    normed = normalize_rows(np.asarray(embeddings, dtype=np.float32))                                                       # Safely normalizes the high-dimensional vectors to prepare for cosine similarity math.
    remaining = list(range(len(scores)))                                                                                    # Creates a pool of available indices pointing to all candidates.

    # Seed the diverse panel with the strongest-scoring candidate.
    first = int(np.argmax(scores))                                                                                          # Finds the index of the candidate with the highest raw performance score.
    order.append(first)                                                                                                     # Adds the best candidate to the final output list right away.
    remaining.remove(first)                                                                                                 # Removes the best candidate from the available pool so it can't be picked twice.

    # Add candidates greedily while penalizing those too close to the already selected ones.
    while remaining and len(order) < top_k:                                                                                 # Loops continuously until the target quota is hit or the candidate pool runs completely dry.
        best_idx = None                                                                                                     # Prepares to track the winner of the current algorithmic round.
        best_value = -math.inf                                                                                              # Prepares the benchmark to beat, set as low as possible.
        # Iterates through every candidate left in the pool.
        for idx in remaining:                                                                                               
            # Calculates the cosine similarity to all chosen items and finds the closest match (highest similarity).
            nearest = max(float(np.dot(normed[idx], normed[chosen])) for chosen in order)                                   
            # Subtracts a weighted fraction of that similarity from the candidate's base score.
            penalized = float(scores[idx] - penalty_weight * nearest)                                                       
            # Checks if this candidate's adjusted score is the best one seen so far this round.
            if penalized > best_value:                                                                                      
                # Updates the high score to beat and the tracking index to the current round leader.
                best_value = penalized                                                                                      # Updates the high score to beat.
                best_idx = idx                                                                                              # Updates the tracking index to the current round leader.
                penalties[idx] = nearest                                                                                    # Logs the exact similarity penalty that was applied for analytical transparency.
        # Accepts the round winner and removes it from the available pool.
        order.append(int(best_idx))                                                                                         # Formally accepts the round winner into the chosen list.
        remaining.remove(int(best_idx))                                                                                     # Removes the newly chosen candidate from the available pool.
    return order, penalties                                                                                                 # Returns the final diverse ordering array and the debug penalty metrics.