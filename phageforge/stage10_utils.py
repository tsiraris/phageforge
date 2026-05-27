"""Shared helpers for Stage 10 structure-conditioned redesign.

This module upgrades the sequence-first Stage 09 search into a true
backbone-conditioned redesign loop. The core idea is simple:
- keep a fixed seed scaffold,
- restrict edits to a compact local region,
- score each candidate against the fixed backbone with an inverse-folding model,
- and keep only candidates that remain both target-seeking and scaffold-compatible.
"""

from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from phageforge.stage07_utils import cosine_similarity, embed_sequences, normalize_rows, read_json as stage07_read_json, write_json as stage07_write_json
from phageforge.stage09_utils import (
    EditProposal,
    build_edit_proposals_from_context,
    choose_editable_positions,
    load_target_predictor,
    predict_target_probability,
)


@dataclass
class Stage10Candidate:
    """
    Compact container for one Stage 10 redesign candidate and its metadata.
    
    Holds the state of a single sequence as it moves through the inverse-folding Beam Search.
    
    Uses Python's `@dataclass` to automatically generate initialization and representation 
    methods, storing the sequence string, its parent's sequence, the specific mutations applied, 
    and the round in which it was generated.
    
    Example: cand = Stage10Candidate(candidate_sequence="MKA...", parent_sequence="MGA...", mutations=["2:G→K"], mutated_positions=[2], proposal_trace=["target_prior"], round_index=1)
    """

    # --- Defines the structural attributes of the candidate object ---
    candidate_sequence: str                                                                                                 # Stores the actual modified amino acid string of the candidate
    parent_sequence: str                                                                                                    # Stores the original amino acid string from which this candidate was derived
    mutations: list[str]                                                                                                    # Stores human-readable string representations of the edits (e.g., '14:A→S')
    mutated_positions: list[int]                                                                                            # Stores the exact 1-based integer indices where edits occurred
    proposal_trace: list[str]                                                                                               # Tracks the algorithmic origin or priority rules that suggested these edits
    round_index: int                                                                                                        # Tracks the specific iteration depth of the Beam Search that spawned this candidate


# ----------------------------- JSON helpers ----------------------------- #


def read_json(path: str | Path) -> dict:
    """
    Read a JSON file using the stable project helper.
    
    Safely loads a JSON file into a Python dictionary.
    
    Acts as a passthrough wrapper to `stage07_read_json`, maintaining consistent 
    file-handling behavior across the entire PhageForge pipeline.
    
    Example:
        data = read_json("results/stage10_context.json")
    """
    # --- Executes the actual file read operation ---
    return stage07_read_json(path)                                                                                          # Delegates the reading execution to the centralized Stage 07 JSON parser


def write_json(obj: dict, path: str | Path) -> None:
    """
    Write a JSON file using the stable project helper.
    
    Safely serializes a Python dictionary into a formatted JSON file.
    
    Acts as a passthrough wrapper to `stage07_write_json`, ensuring uniform 
    indentation and encoding standards.
    
    Example:
        write_json({"status": "ok"}, "results/output.json")
    """
    # --- Executes the actual file write operation ---
    stage07_write_json(obj, path)                                                                                           # Delegates the writing execution to the centralized Stage 07 JSON serializer


# ----------------------------- Path and PDB helpers ----------------------------- #


def resolve_seed_pdb(seed_pdb: str | Path | None, validation_dir: str | Path | None) -> Path:
    """
    Locates the exact physical 3D coordinate file (.pdb) that will act as the fixed scaffold.
    
    If a direct file path is provided, it validates it. If only a directory is provided, 
    it sequentially searches through a hardcoded list of known legacy Stage 08 filenames. 
    If those fail, it performs a recursive glob search for any file containing "seed" and "pdb".
    
    Example:
        path = resolve_seed_pdb(None, "results/stage08/")
        Returns: Path("results/stage08/selected_seed_structure.pdb")
    """
    # --- Checks for explicit user-provided seed path ---
    if seed_pdb is not None:                                                                                                # Evaluates if the user directly supplied a specific PDB file location
        path = Path(seed_pdb)                                                                                               # Converts the raw string input into a robust Path object
        if not path.exists():                                                                                               # Verifies that the specified file physically exists on the disk
            raise FileNotFoundError(f"Missing seed PDB: {path}")                                                            # Aborts execution to prevent downstream inverse-folding crashes
        return path                                                                                                         # Returns the verified explicit path immediately

    # --- Ensures at least one valid search parameter was provided ---
    if validation_dir is None:                                                                                              # Evaluates if the user failed to provide either a file or a fallback directory
        raise ValueError("Provide either --seed_pdb or --validation_dir so Stage 10 can anchor to a scaffold.")             # Aborts execution because Stage 10 cannot function without a 3D anchor

    # --- Constructs a prioritized list of expected legacy filenames ---
    validation_dir = Path(validation_dir)                                                                                   # Converts the fallback directory string into a robust Path object
    candidates = [                                                                                                          # Initializes a list of the most probable filenames from previous pipeline stages
        validation_dir / "selected_seed_structure.pdb",                                                                     # Appends the standard Stage 08 primary seed filename
        validation_dir / "seed_structure.pdb",                                                                              # Appends an alternative generic seed filename
        validation_dir / "seed.pdb",                                                                                        # Appends a highly simplified seed filename
        validation_dir / "structures" / "selected_seed_structure.pdb",                                                      # Appends the primary filename nested within a likely subdirectory
        validation_dir / "structures" / "seed_structure.pdb",                                                               # Appends the generic filename nested within a likely subdirectory
    ]                                                                                                                       # Closes the candidate path array
    
    # --- Checks the prioritized candidates sequentially ---
    for path in candidates:                                                                                                 # Iterates sequentially through the list of highly probable file locations
        if path.exists():                                                                                                   # Checks if the file actually exists at the current evaluated path
            return path                                                                                                     # Returns the first successfully located PDB path

    # --- Performs a desperate recursive search if standard locations fail ---
    matches = sorted(validation_dir.rglob("*seed*.pdb"))                                                                    # Executes a recursive wildcard search for any PDB file containing 'seed' and sorts the results
    if matches:                                                                                                             # Evaluates if the wildcard search yielded at least one valid file
        return matches[0]                                                                                                   # Returns the top alphabetically sorted match from the recursive search
    raise FileNotFoundError(f"Could not locate a seed structure PDB under {validation_dir}")                                # Aborts execution completely if no seed PDB could be found anywhere


def load_inverse_folding_structure(pdb_path: str | Path, chain_id: str = "A") -> tuple[np.ndarray, str]:
    """
    Load backbone coordinates and native sequence from a structure-conditioned IF backend.

    Parses the 3D atomic coordinates from a PDB file into a mathematical matrix required by ESM-IF1.
    
    It lazy-loads the `fair-esm` utilities, extracts the structure, and isolates the specific 
    3D coordinates and amino acid string for the requested protein chain.
    
    Example:
        coords, seq = load_inverse_folding_structure("seed.pdb", "A")
        Returns: (array([[[...]]]), "MGFYAG...")
    """
    # --- Attempts to safely import heavy FAIR-ESM dependencies ---
    try:                                                                                                                    # Initiates a try-block to gracefully catch missing specialized environment dependencies
        from esm.inverse_folding.util import extract_coords_from_structure, load_structure                                  # Imports the specific extraction modules directly from the FAIR-ESM library
    except Exception as exc:  # pragma: no cover - exercised only in target runtime                                         # Catches any import errors if the library is missing
        raise ImportError(                                                                                                  # Raises a clear, instructive error message for the user
            "Stage 10 requires the fair-esm inverse-folding utilities. Install `fair-esm` in SageMaker before running."     # Details exactly how to fix the missing dependency environment issue
        ) from exc                                                                                                          # Chains the original exception for deeper debugging tracebacks

    # --- Executes the atomic coordinate extraction ---
    # Parses the raw PDB file and isolates the specified biological chain
    structure = load_structure(str(pdb_path), chain_id)                                                                     
    # Strips out side-chains, returning only the backbone coordinates and the raw sequence
    coords, native_sequence = extract_coords_from_structure(structure)                                                      
    return np.asarray(coords, dtype=np.float32), str(native_sequence)                                                       # Casts the coordinates to a standard 32-bit float matrix and returns alongside the sequence


def load_inverse_folding_model(device: str = "cuda"):
    """
    Load the ESM-IF1 backbone-conditioned sequence model on the requested device.

    Instantiates the heavy neural network used to score sequences against 3D shapes.
    
    Lazy-loads PyTorch and ESM, downloads/loads the `esm_if1_gvp4` pretrained weights, 
    sets the model to evaluation mode, and pushes it to the GPU/CPU.
    
    Example:
        torch_lib, model, alphabet = load_inverse_folding_model("cuda")
    """
    # --- Attempts to safely import foundational ML libraries ---
    try:                                                                                                                    # Initiates a try-block to gracefully catch missing heavy ML dependencies
        import esm                                                                                                          # Imports the core FAIR-ESM neural network library
        import torch                                                                                                        # Imports the PyTorch deep learning framework
    except Exception as exc:  # pragma: no cover - exercised only in target runtime                                         # Catches any import errors if the environment is misconfigured
        raise ImportError("Stage 10 requires `fair-esm` and `torch` to load the inverse-folding model.") from exc           # Instructs the user that PyTorch and ESM are absolutely required

    # --- Instantiates and configures the Inverse-Folding Neural Network ---
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()                                                           # Downloads or loads the specific 142-million parameter Inverse Folding model and its token alphabet
    model = model.eval().to(device)                                                                                         # Freezes the model weights (eval mode) to prevent training updates and pushes it to the designated compute device
    return torch, model, alphabet                                                                                           # Returns the torch library reference, the loaded model, and the alphabet dictionary


def inverse_folding_log_likelihood(model, alphabet, coords: np.ndarray, sequence: str) -> float:
    """
    Score one sequence against a fixed scaffold using the inverse-folding likelihood.

    Mathematically evaluates how perfectly a 1D amino acid sequence fits inside a 3D shape.
    
    It uses the `score_sequence` utility from ESM-IF1, which feeds the 3D coordinates and 
    the text sequence into the Graph Neural Network. It returns the log-likelihood (a 
    negative float where closer to zero is better).
    
    Example:
        score = inverse_folding_log_likelihood(model, alpha, coords, "MGFY...")
        Returns: -1.25
    """
    # --- Attempts to safely import the specific scoring utility ---
    try:                                                                                                                    # Initiates a try-block for the specific sub-module import
        from esm.inverse_folding.util import score_sequence                                                                 # Imports the explicit scoring function from the FAIR-ESM utilities
    except Exception as exc:  # pragma: no cover - exercised only in target runtime                                         # Catches any import errors
        raise ImportError("Could not import `score_sequence` from fair-esm inverse folding utilities.") from exc            # Raises a descriptive error if the specific utility is missing

    # --- Executes the physical compatibility evaluation ---
    ll_full, _ = score_sequence(model, alphabet, coords, str(sequence))                                                     # Runs the neural network to calculate the full log-likelihood of the sequence fitting the coordinates
    return float(ll_full)                                                                                                   # Casts the resulting PyTorch tensor to a standard Python float and returns it


# ----------------------------- Sequence redesign helpers ----------------------------- #

def mutation_list(seed_sequence: str, candidate_sequence: str) -> list[str]:
    """
    Return 1-based human-readable mutation annotations such as 197:W→R.

    Identifies all the exact changes between an original seed and a new candidate sequence.
    
    It iterates over both strings simultaneously. Whenever the characters differ, it 
    formats a string indicating the position, the old amino acid, and the new one.
    
    Example:
        muts = mutation_list("MKA", "MGA")
        Returns: ["2:K→G"]
    """
    # --- Iterates through sequences to identify differences ---
    rows: list[str] = []                                                                                                    # Initializes an empty list to collect the formatted mutation strings
    for idx, (seed_aa, cand_aa) in enumerate(zip(seed_sequence, candidate_sequence), start=1):                              # Pairs up the characters from both sequences and loops through them, starting the index count at biological position 1
        if seed_aa != cand_aa:                                                                                              # Evaluates if the parental amino acid differs from the candidate's newly proposed amino acid
            rows.append(f"{idx}:{seed_aa}→{cand_aa}")                                                                       # Formats the discrepancy into standard biological notation and appends it to the list
    return rows                                                                                                             # Returns the complete list of identified mutations


def apply_mutation(sequence: str, position_1based: int, aa: str) -> str:
    """
    Return a new sequence with one 1-based substitution applied.

    Safely mutates a single character inside a string based on biological (1-based) indexing.
    
    Translates the 1-based index to a 0-based Python index, performs bounds checking, 
    converts the string to a mutable list, replaces the character, and rejoins it.
    
    Example:
        seq = apply_mutation("MKA", 2, "G")
        Returns: "MGA"
    """
    # --- Validates and executes a safe single-character replacement ---
    idx = int(position_1based) - 1                                                                                          # Translates the user-provided 1-based biological coordinate into a standard 0-based Python list index
    if idx < 0 or idx >= len(sequence):                                                                                     # Verifies that the calculated index strictly falls within the physical boundaries of the sequence length
        raise IndexError(f"Position {position_1based} is outside the sequence length {len(sequence)}")                      # Aborts execution with a clear error if the algorithm attempts to mutate non-existent residues
    if sequence[idx] == aa:                                                                                                 # Checks if the target position already contains the desired amino acid
        return sequence                                                                                                     # Skips the mutation process entirely and returns the unmodified string to save compute time
    chars = list(sequence)                                                                                                  # Casts the immutable sequence string into a mutable Python list of characters
    chars[idx] = aa                                                                                                         # Overwrites the target index with the newly specified amino acid character
    return "".join(chars)                                                                                                   # Recompiles the character list back into a continuous string and returns it


def choose_top_substitutions(proposal: EditProposal, top_k: int = 3) -> list[str]:
    """
    Return the strongest per-position substitution options from target and family priors.

    Ranks the allowed amino acids for a specific position based on their biological priority.
    
    Iterates over all allowed amino acids in the `EditProposal`. It calculates a blended 
    score heavily favoring target-host frequency (0.55), followed by family frequency (0.30), 
    and general functional weight (0.15). Returns the top-k highest scoring amino acids.
    
    Example:
        top_aas = choose_top_substitutions(proposal, top_k=2)
        Returns: ['S', 'T']
    """
    # --- Evaluates and ranks the provided biological substitution options ---
    aa_scores: list[tuple[str, float]] = []                                                                                 # Initializes an empty list to track tuples containing the amino acid and its computed priority score
    for aa in proposal.allowed_aas:                                                                                         # Iterates through every single legally permitted amino acid substitution defined in the proposal
        target = float(proposal.target_preference.get(aa, 0.0))                                                             # Extracts the specific probability that this amino acid successfully infects the target host
        family = float(proposal.family_preference.get(aa, 0.0))                                                             # Extracts the specific probability that this amino acid is naturally utilized by the evolutionary family
        score = 0.55 * target + 0.30 * family + 0.15 * float(proposal.functional_weight)                                    # Calculates a mathematically blended priority score weighted heavily toward target-host success
        aa_scores.append((aa, score))                                                                                       # Bundles the amino acid character and its calculated score into a tuple and stores it
    aa_scores.sort(key=lambda item: (item[1], item[0]), reverse=True)                                                       # Sorts the collected options primarily by descending score, and secondarily alphabetically for deterministic tie-breaking
    return [aa for aa, _ in aa_scores[:top_k]]                                                                              # Slices the list to keep only the requested top-k items and strips away the numerical scores, returning just the strings


def build_stage10_context(
    context: dict,
    seed_pdb_path: str | Path,
    strict_df: pd.DataFrame | None,
    max_edit_positions: int,
    soft_positions: int,
    min_mutations: int,
    max_mutations: int,
    seed: int,
) -> dict:
    """
    Construct the compact Stage 10 redesign context from Stage 07 artifacts.

    Builds the ultimate structural constraint blueprint. It shrinks the Stage 09 edit space 
    further to ensure inverse-folding stability.
    
    It re-calculates the EditProposals. It uses the geographic `choose_editable_positions` 
    to lock in the "hard" spots, assigns the next best to the "soft" spots, freezes everything 
    else, and packages this along with the strict physical PDB path into a master JSON dictionary.
    
    Example:
        stage10_ctx = build_stage10_context(ctx, "seed.pdb", df, 8, 4, 2, 6, 42)
        Returns: {"stage": "10a", "target_host": "Acinetobacter", ...}
    """
    # --- Rebuilds and rigidly filters the baseline substitution rulebook ---
    proposals = build_edit_proposals_from_context(context=context, strict_df=strict_df)                                     # Generates the foundational list of all biologically viable edit proposals using Stage 09 utilities
    if not proposals:                                                                                                       # Evaluates if the generator completely failed to find any mathematically viable mutations
        raise ValueError("No Stage 10 edit proposals could be built from the Stage 07 context.")                            # Aborts execution because the pipeline cannot redesign a completely locked protein

    hard_positions = choose_editable_positions(proposals=proposals, max_positions=max_edit_positions, seed=seed)            # Executes a geographically diverse selection algorithm to lock in the primary high-priority mutation sites
    sorted_proposals = sorted(                                                                                              # Begins sorting the remaining pool of proposals
        proposals,                                                                                                          # Passes the complete proposal list to the sorter
        key=lambda item: (item.functional_weight, -item.conservation_penalty, item.position),                               # Defines the sorting hierarchy: maximum functional reward, minimum evolutionary risk, then positional index
        reverse=True,                                                                                                       # Enforces descending order so the absolute best biological proposals sit at the top
    )                                                                                                                       # Closes the sorting execution block
    
    # --- Populates the flexible soft (secondary) editing buffer ---
    soft: list[int] = []                                                                                                    # Initializes an empty list to track the secondary backup mutation positions
    for item in sorted_proposals:                                                                                           # Iterates systematically down the perfectly sorted list of all proposals
        if int(item.position) in hard_positions:                                                                            # Checks if the current highest-scoring position was already drafted into the hard position list
            continue                                                                                                        # Skips the position to prevent duplicates across the edit tiers
        soft.append(int(item.position))                                                                                     # Assigns the next-best available position to the soft buffer pool
        if len(soft) >= soft_positions:                                                                                     # Checks if the soft buffer has reached its pre-defined capacity limit
            break                                                                                                           # Terminates the assignment loop immediately to preserve the strict mutation constraints

    # --- Defines the absolute frozen architectural boundaries ---
    seed_sequence = str(context["selected_seed"]["seed_sequence"])                                                          # Extracts the raw parental amino acid sequence string from the original blueprint
    total_editable_positions = sorted(set(hard_positions) | set(soft))                                                      # Mathematically unions the hard and soft sets to create the total authorized edit boundary
    frozen_positions = sorted(set(range(1, len(seed_sequence) + 1)) - set(total_editable_positions))                        # Subtracts the editable boundary from all possible positions, completely locking the rest of the protein chassis

    # --- Serializes the strict proposal rules for export ---
    proposal_rows: list[dict] = []                                                                                          # Initializes a payload list to hold the serialized dictionary versions of the EditProposals
    for item in proposals:                                                                                                  # Iterates through every single generated proposal object
        proposal_rows.append(                                                                                               # Begins appending a newly structured dictionary to the payload list
            {                                                                                                               # Opens the dictionary mapping for JSON serialization
                "position": int(item.position),                                                                             # Standardizes the biological sequence index as a pure integer
                "seed_aa": str(item.seed_aa),                                                                               # Standardizes the wild-type amino acid character as a string
                "allowed_aas": list(item.allowed_aas),                                                                      # Casts the unique set of mathematically permitted substitutions back into a standard JSON list
                "target_preference": dict(item.target_preference),                                                          # Ensures the target probability mapping is a native dictionary
                "family_preference": dict(item.family_preference),                                                          # Ensures the family probability mapping is a native dictionary
                "functional_weight": float(item.functional_weight),                                                         # Standardizes the calculated biological necessity score as a float
                "conservation_penalty": float(item.conservation_penalty),                                                   # Standardizes the calculated evolutionary risk score as a float
                "region_name": str(item.region_name),                                                                       # Standardizes the structural window grouping name as a string
                "edit_tier": "hard" if int(item.position) in hard_positions else ("soft" if int(item.position) in soft else "frozen"), # Dynamically evaluates and assigns the strict hierarchical constraint tier
            }                                                                                                               # Closes the dictionary mapping
        )                                                                                                                   # Closes the append execution

    # --- Compiles the ultimate overarching Stage 10 blueprint ---
    return {                                                                                                                # Begins returning the final master configuration dictionary
        "stage": "10a",                                                                                                     # Tags the artifact with its specific originating pipeline stage
        "target_host": str(context["target_host"]),                                                                         # Preserves the specific biological bacteria the project aims to infect
        "seed_pdb_path": str(seed_pdb_path),                                                                                # Embeds the absolutely critical absolute filepath to the physical 3D scaffold
        "selected_seed": dict(context["selected_seed"]),                                                                    # Carries forward the core identification metadata for the wild-type protein
        "family_context": dict(context.get("family_context", {})),                                                          # Carries forward the broad evolutionary context necessary for downstream metric tracking
        "editable_region": {                                                                                                # Opens a nested dictionary dedicated entirely to mutation restrictions
            "hard_positions": hard_positions,                                                                               # Saves the primary authorized mutation indices
            "soft_positions": soft,                                                                                         # Saves the secondary backup mutation indices
            "editable_positions": total_editable_positions,                                                                 # Saves the aggregate permitted boundary
            "frozen_positions": frozen_positions,                                                                           # Saves the massive array of structurally locked chassis indices
            "min_mutations": int(min_mutations),                                                                            # Enforces the absolute lower bound of sequence edits required per candidate
            "max_mutations": int(max_mutations),                                                                            # Enforces the absolute upper bound of sequence edits permitted per candidate
            "proposal_rows": proposal_rows,                                                                                 # Embeds the massive, strictly serialized substitution rulebook
        },                                                                                                                  # Closes the mutation restriction dictionary
        "source_stage07_context_json": str(context.get("upstream_artifacts", {}).get("phaseA_plan_json", "")),              # Provides an unbroken data provenance trail back to the initial experimental plan
    }                                                                                                                       # Closes the master configuration dictionary


# ----------------------------- Embedding-backbone helpers ----------------------------- #


def load_embedding_backend(model_name: str, device: str | None = None):
    """
    Load the ESM embedding backbone once so Stage 10 can reuse it across many rounds.

    Pre-loads the HuggingFace transformer model into GPU memory and keeps it there.
    
    Uses the `transformers` library to fetch the tokenizer and model. This avoids the 
    massive computational overhead of reloading the 650M parameter model from disk 
    on every single evaluation step of the Beam Search.
    
    Example:
        torch_lib, tok, mod, dev = load_embedding_backend("facebook/esm2_t33_650M_UR50D")
    """
    # --- Attempts to securely load the heavy Transformer infrastructure ---
    import torch                                                                                                            # Imports the underlying PyTorch tensor computation library
    from transformers import AutoModel, AutoTokenizer                                                                       # Imports the generalized HuggingFace loaders required for the ESM language model

    resolved_device = device or ("cuda" if torch.cuda.is_available() else "cpu")                                            # Dynamically detects GPU availability if no specific hardware was explicitly requested
    tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)                                              # Downloads or loads the specific text tokenizer, strictly preserving biological uppercase formatting
    model = AutoModel.from_pretrained(model_name).to(resolved_device).eval()                                                # Loads the massive transformer weights, maps them to the correct hardware, and freezes them for inference
    return torch, tokenizer, model, resolved_device                                                                         # Returns the live toolkit bundle back to the search engine for persistent reuse


def embed_sequences_with_backend(
    sequences: list[str],
    torch,
    tokenizer,
    model,
    device: str,
    batch_size: int = 4,
    max_length: int = 2048,
) -> np.ndarray:
    """
    Embed protein sequences with a preloaded ESM backend using masked mean pooling.

    Translates raw amino acid strings into high-dimensional mathematical vectors.
    
    It batches the sequences, tokenizes them, runs them through the Transformer `model`, 
    extracts the final hidden state, and calculates a mean-pooled vector (ignoring padding 
    tokens using the attention mask).
    
    Example:
        embs = embed_sequences_with_backend(["MKA", "MGA"], torch, tok, mod, "cuda")
        Returns: array([[ 0.01, -0.05, ... ], [ 0.02, -0.03, ... ]])
    """
    # --- Executes high-performance batched inference through the Transformer ---
    rows: list[np.ndarray] = []                                                                                             # Initializes an empty list to gather the resulting sequence vectors
    with torch.no_grad():                                                                                                   # Temporarily disables gradient tracking to drastically reduce memory consumption during pure inference
        for start in range(0, len(sequences), batch_size):                                                                  # Slices the massive candidate list into manageable hardware-appropriate chunks
            batch = sequences[start : start + batch_size]                                                                   # Extracts the specific subset of sequence strings for this processing loop
            toks = tokenizer(batch, padding=True, truncation=True, max_length=max_length, return_tensors="pt")              # Converts the biological text into standardized numerical tensors, padding short sequences with zeros
            toks = {key: value.to(device) for key, value in toks.items()}                                                   # Rapidly transfers all generated tensors from system RAM into the GPU VRAM
            hidden = model(**toks).last_hidden_state                                                                        # Pushes the tensors through the network, extracting the complex mathematical representations from the final layer
            mask = toks["attention_mask"].unsqueeze(-1)                                                                     # Isolates the attention mask and expands its dimensions to perfectly align with the hidden state structure
            pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1)                                              # Multiplies the hidden state by the mask to delete padding noise, then averages the remaining valid token vectors
            rows.append(pooled.detach().cpu().numpy())                                                                      # Severs the tensor from the compute graph, pulls it back to system RAM, converts it to a standard numpy array, and saves it
    return np.vstack(rows).astype(np.float32) if rows else np.zeros((0, 0), dtype=np.float32)                               # Vertically stacks all the processed batches into one massive matrix, returning a safe zero-array if the input was empty

# ----------------------------- Candidate evaluation helpers ----------------------------- #


def compute_family_centroid(
    family_sequences: list[str],
    embedding_model: str | None = None,
    batch_size: int = 4,
    torch=None,
    tokenizer=None,
    model=None,
    device: str | None = None,
) -> np.ndarray:
    """
    Embed family sequences and return their centroid under the predictor backbone.

    Calculates the mathematical "center of mass" for the entire evolutionary protein family.
    
    It embeds all provided natural family sequences into vectors, then averages them 
    along axis 0 to find the generalized evolutionary centroid. It uses the pre-loaded 
    backend if provided, or instantiates a fresh one if missing.
    
    Example:
        centroid = compute_family_centroid(["MKA", "MGA", "MTA"], tok=tok, mod=mod)
        Returns: array([ 0.015, -0.04, ... ])
    """
    # --- Computes the gravitational center of the evolutionary manifold ---
    if not family_sequences:                                                                                                # Evaluates if the user failed to provide any evolutionary cousins
        return np.zeros((1,), dtype=np.float32)                                                                             # Safely aborts by returning a dummy zero-vector to prevent downstream math crashes
    if tokenizer is not None and model is not None and torch is not None and device is not None:                            # Checks if the persistent high-performance backend tools were successfully injected
        embs = embed_sequences_with_backend(family_sequences, torch, tokenizer, model, device=device, batch_size=batch_size)# Executes embedding generation utilizing the fast pre-loaded environment
    else:                                                                                                                   # Engages fallback logic if the persistent tools were missing
        embs = embed_sequences(family_sequences, model_name=str(embedding_model), batch_size=batch_size)                    # Executes a slower, from-scratch embedding process using the legacy Stage 07 utility
    return np.asarray(embs.mean(axis=0), dtype=np.float32)                                                                  # Averages all the generated vectors mathematically to find the absolute center point, returning it as a strict 32-bit float array


def evaluate_candidate_table(
    sequences: list[str],
    target_host: str,
    predictor_model_path: str | Path | None,
    predictor_label_classes_path: str | Path | None,
    embedding_model: str,
    family_centroid: np.ndarray | None,
    coords: np.ndarray,
    if_model,
    if_alphabet,
    batch_size: int = 4,
    predictor=None,
    label_classes: list[str] | None = None,
    torch=None,
    tokenizer=None,
    model=None,
    device: str | None = None,
) -> pd.DataFrame:
    """
    Compute target and scaffold-conditioned scores for a batch of candidate sequences.

    The ultimate evaluation gauntlet. It grades every sequence on three critical metrics:
    Will it kill the target? Will it fold into the 3D scaffold? Is it evolutionarily realistic?
    
    It specific, the script embeds the sequences, runs the Target Predictor to get host probability, 
    runs the Inverse-Folding to get structural log-likelihood against the `coords`,
    calculates the Cosine Similarity against the `family_centroid`, and returns a DataFrame with all these scores.
    
    Example:
        df = evaluate_candidate_table(["MKA"], "Klebsiella", ..., coords, if_mod, if_alph)
        Returns: DataFrame with 'target_probability': 0.68, 'if1_log_likelihood': 0.12, etc.
    """
    # --- Orchestrates the massive multi-modal evaluation pass for all generated candidates ---
    if not sequences:                                                                                                       # First strictly verifies that there are actually candidate sequences to process
        return pd.DataFrame(columns=["candidate_sequence", "target_probability", "if1_log_likelihood", "family_cosine"])    # Gracefully returns an empty but correctly formatted dataframe structure to prevent downstream loop errors

    if predictor is None or label_classes is None:                                                                          # Evaluates if the predictive machine learning models were not persistently cached in memory
        predictor, label_classes = load_target_predictor(predictor_model_path, predictor_label_classes_path)                # Incurs the disk IO penalty to load the logistic regression models from scratch
    if tokenizer is not None and model is not None and torch is not None and device is not None:                            # Checks if the heavy persistent transformer backend is available
        embeddings = embed_sequences_with_backend(sequences, torch, tokenizer, model, device=device, batch_size=batch_size) # Rapidly generates high-dimensional vectors for all candidates using the persistent backend
    else:                                                                                                                   # Engages fallback logic for isolated execution environments
        embeddings = embed_sequences(sequences, model_name=embedding_model, batch_size=batch_size)                          # Lazily loads the transformer and generates embeddings using the older utility
    target_prob = predict_target_probability(predictor, label_classes, target_host, embeddings)                             # Executes the logistic regression inference to determine exactly how likely the proteins are to infect the designated bacteria

    if family_centroid is not None and family_centroid.size == embeddings.shape[1]:                                         # Safely verifies that the evolutionary center vector exists and mathematically aligns with the embedding dimensions
        family_cos = np.array([cosine_similarity(row, family_centroid) for row in embeddings], dtype=np.float32)            # Calculates the exact angular distance between each candidate and the evolutionary norm
    else:                                                                                                                   # Engages fallback logic if the family math is corrupt or missing
        family_cos = np.zeros(len(sequences), dtype=np.float32)                                                             # Injects a passive zero-array to ensure the pipeline doesn't crash during final scoring

    # --- Grade the physical 3D stability of every single sequence
    if_scores = np.array([inverse_folding_log_likelihood(if_model, if_alphabet, coords, seq) for seq in sequences], dtype=np.float32) 
    
    # --- Assembles the scores into a unified data structure ---
    frame = pd.DataFrame(                                                                                                   # Initiates the construction of the final evaluation table
        {                                                                                                                   # Opens the dictionary mapping columns to data arrays
            "candidate_sequence": sequences,                                                                                # Attaches the raw biological text string
            "target_probability": target_prob,                                                                              # Attaches the predictive infectivity scores
            "if1_log_likelihood": if_scores,                                                                                # Attaches the predictive 3D structural stability scores
            "family_cosine": family_cos,                                                                                    # Attaches the evolutionary realism scores
        }                                                                                                                   # Closes the dictionary mapping
    )                                                                                                                       # Closes the dataframe construction
    return frame                                                                                                            # Returns the fully populated intelligence table back to the search engine


def robust_minmax(values: Iterable[float]) -> np.ndarray:
    """
    Normalize a numeric vector to [0, 1] while guarding against degenerate ranges.

    Scales a list of numbers so the lowest becomes 0.0 and the highest becomes 1.0.
    
    Calculates standard (val - min) / (max - min). If all values are identical (max == min), 
    it avoids a divide-by-zero error by returning an array of 0.5.
    
    Example:
        normed = robust_minmax([10, 20, 30])
        Returns: array([0.0, 0.5, 1.0])
    """
    # --- Executes safe mathematical data scaling ---
    arr = np.asarray(list(values), dtype=np.float32)                                                                        # Casts the flexible input iterable into a strictly formatted, computationally efficient float array
    if arr.size == 0:                                                                                                       # Checks if the input array is completely devoid of data
        return arr                                                                                                          # Instantly returns the empty array to prevent downstream math errors
    lo = float(arr.min())                                                                                                   # Extracts the absolute lowest numerical value present in the entire dataset
    hi = float(arr.max())                                                                                                   # Extracts the absolute highest numerical value present in the entire dataset
    if hi <= lo:                                                                                                            # Evaluates the critical edge case where every single number in the dataset is identical
        return np.ones_like(arr) * 0.5                                                                                      # Circumvents an illegal divide-by-zero error by assigning a neutral median score to all items
    return (arr - lo) / (hi - lo)                                                                                           # Executes the standard feature scaling formula to force all data proportionately between zero and one


def composite_stage10_score(
    target_probability: np.ndarray,
    if1_log_likelihood: np.ndarray,
    family_cosine: np.ndarray,
    seed_identity: np.ndarray,
    mutation_count: np.ndarray,
    w_target: float = 0.30,                                                                                                 # Tunable weight; default reproduces the original hardcoded Stage 10 formula exactly
    w_if1: float = 0.45,                                                                                                    # Tunable weight; default reproduces the original hardcoded Stage 10 formula exactly
    w_family: float = 0.15,                                                                                                 # Tunable weight; default reproduces the original hardcoded Stage 10 formula exactly
    w_identity: float = 0.10,                                                                                               # Tunable weight; default reproduces the original hardcoded Stage 10 formula exactly
    w_mut_penalty: float = 0.10,                                                                                            # Tunable weight; default reproduces the original hardcoded Stage 10 formula exactly
) -> np.ndarray:
    """
    Combine Stage 10 signals into one ranking score.

    Calculates the ultimate fitness function determining which candidates survive the Beam Search.
    
    It normalizes all inputs to a 0.0-1.0 scale using `robust_minmax`. It then computes 
    a weighted sum where Inverse-Folding (0.45) is the dominant factor, followed by 
    Target Probability (0.30), penalizing excessive mutations.
    
    Example:
        score = composite_stage10_score(targets, if_scores, family, identity, muts)
        Returns: array([0.85, 0.42, ...])
    """
    # --- Normalizes disparate metrics into uniform comparable scales ---
    target_norm = robust_minmax(target_probability)                                                                         # Scales the biological infectivity predictions cleanly between zero and one
    if1_norm = robust_minmax(if1_log_likelihood)                                                                            # Scales the structural physics evaluations cleanly between zero and one
    family_norm = robust_minmax(family_cosine)                                                                              # Scales the evolutionary realism geometry cleanly between zero and one
    identity_norm = robust_minmax(seed_identity)                                                                            # Scales the wild-type preservation metrics cleanly between zero and one
    mut_penalty = robust_minmax(mutation_count)                                                                             # Scales the edit burden volume cleanly between zero and one
    
    # --- Computes the heavily-weighted structural survival calculus ---
    return w_target * target_norm + w_if1 * if1_norm + w_family * family_norm + w_identity * identity_norm - w_mut_penalty * mut_penalty   # Weighted blend (weights now tunable; defaults grant supreme authority to the 3D physics engine while punishing excessive deviation)


def greedy_diverse_subset(embeddings: np.ndarray, scores: np.ndarray, top_k: int, diversity_penalty_weight: float = 0.20) -> list[int]:
    """
    Choose a diversity-aware top-k subset from candidate embeddings.

    Selects the best candidates while preventing "Mode Collapse" 
    (for example the model submits 10 identical clones of the same protein).
    
    It iteratively selects the absolute highest-scoring candidate. For the remaining candidates, 
    it mathematically penalizes their score based on their cosine similarity to the already-chosen 
    sequences. This forces the selection panel to jump to different structural clusters.
    
    Example:
        chosen_indices = greedy_diverse_subset(embeddings, scores, top_k=3)
        Returns: [42, 12, 105]
    """
    # --- Validates constraints before initiating the complex spatial selection logic ---
    embeddings = np.asarray(embeddings, dtype=np.float32)                                                                   # Secures the high-dimensional spatial coordinates into a strict numeric matrix
    scores = np.asarray(scores, dtype=np.float32)                                                                           # Secures the composite fitness evaluations into a strict numeric array
    if len(scores) == 0:                                                                                                    # Verifies that there are actually candidates waiting to be evaluated
        return []                                                                                                           # Aborts the selection process harmlessly if the pool is empty
    if len(scores) <= top_k:                                                                                                # Evaluates if the available pool is already smaller than the requested output quota
        return list(np.argsort(-scores))                                                                                    # Bypasses the expensive spatial calculus entirely and simply returns a standard descending sort

    # --- Prepares the mathematical environment for the iterative distance penalties ---
    normed = normalize_rows(embeddings)                                                                                     # Converts raw vectors into unit vectors to drastically speed up the subsequent cosine similarity dot products
    remaining = list(range(len(scores)))                                                                                    # Initializes a tracking list containing the indices of all candidates still eligible for selection
    chosen = [int(np.argmax(scores))]                                                                                       # Bypasses the spatial check for the very first pick, directly granting victory to the absolute highest-scoring candidate
    remaining.remove(chosen[0])                                                                                             # Exiles the inaugural victor from the pool of eligible competitors

    # --- Executes the iterative, penalty-driven selection gauntlet ---
    while remaining and len(chosen) < top_k:                                                                                # Loops continuously until either the requested quota is filled or the entire candidate pool is exhausted
        best_idx = None                                                                                                     # Prepares a temporary variable to track the index of the winner for the current round
        best_value = -math.inf                                                                                              # Prepares an infinitely negative baseline score that the candidates must beat to take the lead
        for idx in remaining:                                                                                               # Iterates sequentially through every single candidate still waiting in the eligible pool
            nearest = max(float(np.dot(normed[idx], normed[picked])) for picked in chosen)                                  # Calculates the spatial distance to all previously selected victors and isolates the most similar one (the highest proximity threat)
            value = float(scores[idx]) - diversity_penalty_weight * nearest                                                                     # Radically slashes the candidate's fitness score based on how closely it mirrors the sequences already drafted to the panel
            if value > best_value:                                                                                          # Checks if this heavily penalized score is still somehow the highest number evaluated during this specific round
                best_value = value                                                                                          # Overwrites the benchmark with the new leading score
                best_idx = idx                                                                                              # Transfers the temporary crown to the current candidate's index
        chosen.append(int(best_idx))                                                                                        # Formally drafts the round's survivor into the elite final panel
        remaining.remove(best_idx)                                                                                          # Exiles the newly drafted candidate from the eligible pool to prevent duplicate selections
    return chosen                                                                                                           # Returns the fully assembled list of mathematically diverse indices


def seed_everything(seed: int) -> None:
    """
    Set deterministic random seeds for Python and NumPy for reproducible search order.

    Forces all random number generators to follow a predictable, repeatable pattern.
    
    Sets the manual seed for both the standard Python `random` module and the `numpy.random` 
    module.
    
    Example:
        seed_everything(42)
    """
    # --- Locks internal state engines to guarantee scientific reproducibility ---
    random.seed(seed)                                                                                                       # Hardcodes the entropy source for the standard Python standard library randomizer
    np.random.seed(seed)                                                                                                    # Hardcodes the entropy source for the complex NumPy mathematical array randomizer