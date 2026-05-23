"""Stage 11 utilities: universal, self-contained structure-conditioned RBP redesign.

Stage 11 differs from Stage 10 in exactly one architectural commitment: every
artifact required to redesign the seed (family centroid, target centroid,
edit space, hotspot priors) is recomputed inside the stage from the strict
RBP CSV + the cached ESM-2 embeddings + the trained host probe. No JSON
context, edit space, or surrogate model from Stages 06/07/09/10 is read.

The module is deliberately small. Heavy lifting (ESM-IF1 scoring, embedding
backbones, predictor loading, diversity reranking, composite scoring) is
reused from stage07_utils / stage09_utils / stage10_utils via explicit
re-exports so that the Stage 11 scripts can import everything from one
place. The new logic owned by this module covers:

1. Strict CSV loading and seed-row selection
2. ESMFold of a single sequence + Baseline Qualification gate
3. Family/target reference-set construction from cached embeddings
4. Per-position Shannon entropy + amino-acid preference computation
5. From-scratch EditProposal construction (no Stage 06 hotspots)
6. Hard/soft/frozen editable position selection with regional spreading
7. Run-metadata + run-name helpers (UTC ISO timestamps, GPU info, etc.)
"""

from __future__ import annotations

import json
import math
import os
import platform
import random
import socket
import subprocess
import sys
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

# --- Re-exports from earlier-stage utilities so 11a–11e can import from one place. -- #
from phageforge.stage07_utils import (                                                                                  # Re-export Stage 07 helpers Stage 11 scripts also need.
    AMINO_ACIDS,                                                                                                        # Canonical 20-letter alphabet used by the family/target preference math.
    cosine_similarity,                                                                                                  # Vector-vs-vector cosine, used to score family-cosine of candidates.
    embed_sequences,                                                                                                    # Fallback embedder when no persistent backend is loaded.
    normalize_rows,                                                                                                     # L2 row normalizer used by diversity reranking and centroid math.
    read_json as _read_json,                                                                                            # Stable JSON reader; aliased so the Stage 11 namespace keeps a clean name.
    write_json as _write_json,                                                                                          # Stable JSON writer; aliased so the Stage 11 namespace keeps a clean name.
    write_fasta,                                                                                                        # FASTA writer used by 11e for the handoff bundle.
    seed_everything,                                                                                                    # Master deterministic seeder for python/numpy/torch.
)
from phageforge.stage09_utils import (                                                                                  # Re-export the structural Stage 09 helpers Stage 11 reuses verbatim.
    EditProposal,                                                                                                       # Dataclass describing the allowed substitutions for one position.
    VALID_AA,                                                                                                           # Set of valid amino acids (fast membership testing).
    build_basic_sequence_features,                                                                                      # Cheap diagnostic features (used by 11e summaries).
    load_target_predictor,                                                                                              # Loader for the trained logistic regression host probe.
    parse_mutation_positions,                                                                                           # Robust integer extractor from mutation strings like "12:A→T".
    predict_target_probability,                                                                                         # Predictor inference, used by 11b's evaluation pass.
    sequence_identity,                                                                                                  # Positional identity used in composite scoring.
)
from phageforge.stage10_utils import (                                                                                  # Re-export Stage 10 core inverse-folding helpers used by 11b and 11c.
    Stage10Candidate,                                                                                                   # Candidate dataclass; Stage 11 reuses the same shape.
    apply_mutation,                                                                                                     # Single-position substitution helper.
    choose_top_substitutions,                                                                                           # Per-position substitution ranking.
    composite_stage10_score,                                                                                            # Locked composite scoring formula (target/IF1/family/identity/mut_count).
    compute_family_centroid,                                                                                            # ESM-2 family-centroid computation backed by the persistent backend.
    embed_sequences_with_backend,                                                                                       # Persistent-backend embedder used during search/prefilter.
    evaluate_candidate_table,                                                                                           # Multi-modal evaluation gauntlet (target prob + IF1 ll + family cos).
    greedy_diverse_subset,                                                                                              # ESM-2 embedding-space diversity reranking.
    inverse_folding_log_likelihood,                                                                                     # IF1 score for a single (coords, sequence) pair.
    load_embedding_backend,                                                                                             # Persistent ESM-2 backend loader (tokenizer + model).
    load_inverse_folding_model,                                                                                         # ESM-IF1 model loader.
    load_inverse_folding_structure,                                                                                     # PDB → (coords, native_seq) parser.
    mutation_list,                                                                                                      # Per-position mutation annotation (e.g. "197:W→R").
    robust_minmax,                                                                                                      # Safe min-max scaler used by composite scoring.
)

# Stage 11 reuses the Stage 10 candidate dataclass unchanged; expose it under the Stage 11 name
# so 11b's `from phageforge.stage11_utils import Stage11Candidate` resolves cleanly.
Stage11Candidate = Stage10Candidate                                                                                     # Explicit alias (identical shape; avoids an ImportError in 11b).


# Stage 11 exit codes — tightly scoped so wrapping shells can dispatch on outcome.
EXIT_OK = 0
EXIT_INPUT_ERROR = 1                                                                                                    # Bad CLI / missing files / malformed CSV.
EXIT_BASELINE_QUALIFICATION_FAILED = 2                                                                                  # Seed itself folds badly — pick a different seed and rerun.
EXIT_INFERENCE_ERROR = 3                                                                                                # ESMFold / IF1 / embedding OOM or other runtime crash.


# ----------------------------- JSON I/O passthroughs ----------------------------- #


def read_json(path: str | Path) -> dict:
    """Stable JSON reader (delegates to the canonical Stage 07 helper)."""
    return _read_json(path)                                                                                             # One single source of truth for read semantics across the project.


def write_json(obj: dict, path: str | Path) -> None:
    """Stable JSON writer (delegates to the canonical Stage 07 helper)."""
    _write_json(obj, path)                                                                                              # One single source of truth for write semantics across the project.


# ----------------------------- Time / naming helpers ----------------------------- #


def utc_timestamp() -> str:
    """Filesystem-safe UTC ISO-8601 timestamp like 20260520T143000Z."""
    return datetime.now(timezone.utc).strftime("%Y%m%d")                                                                # Date-only token so reruns on the same day reuse the same run folder (no timestamp churn).


def default_run_name(source_host: str, target_host: str, seed_protein_id: str, seed: int) -> str:
    """Run-name convention: <source>_to_<target>_<seed_id>_seed<N>_<UTC>.

    The protein_id is sanitized (dots → underscores) so the run name remains a
    safe directory segment on every filesystem we ship to.
    """
    sanitized = str(seed_protein_id).replace(".", "_").replace("/", "_").replace(" ", "_")                              # Keep the suffix filesystem-safe across Windows/Linux/macOS.
    return f"{source_host}_to_{target_host}_{sanitized}_seed{int(seed)}_{utc_timestamp()}"                              # Deterministic-up-to-UTC-second human-readable run name.


# ----------------------------- Strict CSV / seed selection ----------------------------- #


REQUIRED_STRICT_COLUMNS = ("virus_accession", "host_genus", "protein_id", "product", "aa_sequence")                     # Hard contract the strict CSV must satisfy.


def load_strict_dataset(strict_csv: str | Path) -> pd.DataFrame:
    """Load the strict RBP CSV and validate the required column contract.

    Returns a copy with `aa_sequence` and `host_genus` cast to clean strings, the
    1-based row-of-arrival preserved as `_row_id`, and any whitespace-only rows
    discarded.
    """
    path = Path(strict_csv)                                                                                             # Standardize into a robust pathlib object.
    if not path.exists():                                                                                               # Defensive existence check before pandas fails with a less-helpful error.
        raise FileNotFoundError(f"Missing strict CSV: {path}")                                                          # Surface the missing path directly so the operator sees it.
    df = pd.read_csv(path)                                                                                              # Eager-load the curated RBP bank into memory.
    missing = [col for col in REQUIRED_STRICT_COLUMNS if col not in df.columns]                                         # Identify any required columns that are absent.
    if missing:                                                                                                         # Refuse to proceed on a malformed CSV.
        raise ValueError(f"Strict CSV at {path} is missing required columns: {missing}")                                # Loud error so the operator can fix the input.
    df = df.copy()                                                                                                      # Avoid mutating the in-memory CSV object the caller may still hold.
    df["host_genus"] = df["host_genus"].astype(str).str.strip()                                                         # Normalize host genus to clean trimmed strings.
    df["aa_sequence"] = df["aa_sequence"].astype(str).str.strip()                                                       # Normalize sequence text to clean trimmed strings.
    df["protein_id"] = df["protein_id"].astype(str).str.strip()                                                         # Normalize protein_id to clean trimmed strings.
    df["_row_id"] = np.arange(len(df))                                                                                  # Stable integer row index for embedding lookups.
    return df                                                                                                           # Hand the validated frame back to the caller.


def select_seed_row(strict_df: pd.DataFrame, seed_protein_id: str, expected_source_host: str) -> pd.Series:
    """Return the strict-CSV row matching `seed_protein_id` and validate its host."""
    seed_protein_id = str(seed_protein_id).strip()                                                                      # Normalize the user's input to match CSV normalization.
    expected_source_host = str(expected_source_host).strip()                                                            # Same normalization for the source-host filter.
    hits = strict_df.loc[strict_df["protein_id"] == seed_protein_id]                                                    # Locate every row whose protein_id matches the seed.
    if hits.empty:                                                                                                      # Bail out clearly if the operator typo'd or chose a non-strict ID.
        raise ValueError(f"Seed protein_id '{seed_protein_id}' not found in strict CSV.")                               # Loud error before any compute is spent.
    if len(hits) > 1:                                                                                                   # Defensive: the strict CSV should not contain duplicates.
        raise ValueError(f"Seed protein_id '{seed_protein_id}' has multiple rows in strict CSV ({len(hits)}).")         # Surface the ambiguity rather than silently picking one.
    row = hits.iloc[0]                                                                                                  # The single canonical row for the seed.
    if str(row["host_genus"]) != expected_source_host:                                                                  # Enforce the operator's --source_host claim.
        raise ValueError(                                                                                               # Tell the operator exactly what went wrong.
            f"Seed protein_id '{seed_protein_id}' has host_genus='{row['host_genus']}', "
            f"but --source_host was '{expected_source_host}'. Re-check the strict CSV row."
        )
    seq = str(row["aa_sequence"]).strip().upper()                                                                       # Normalize the seed sequence to canonical uppercase.
    invalid = sorted(set(seq) - VALID_AA)                                                                               # Detect any non-canonical amino-acid characters.
    if invalid:                                                                                                         # Refuse to fold a sequence with stop codons or wildcards.
        raise ValueError(f"Seed sequence for '{seed_protein_id}' contains non-canonical residues: {invalid}")           # Surface every offending character.
    if "*" in seq:                                                                                                      # Defensive double-check that no in-frame stop codon snuck through.
        raise ValueError(f"Seed sequence for '{seed_protein_id}' contains an internal stop codon ('*').")               # Refuse to proceed on a corrupted sequence.
    if not (200 <= len(seq) <= 1500):                                                                                   # Soft length sanity check for RBPs (strict CSV ranges ~250–900 AA).
        raise ValueError(f"Seed sequence length {len(seq)} for '{seed_protein_id}' falls outside the 200..1500 range.") # Refuse pathologically short/long seeds.
    row = row.copy()                                                                                                    # Detach from the parent DataFrame so callers can edit safely.
    row["aa_sequence"] = seq                                                                                            # Persist the canonicalized sequence back into the row.
    return row                                                                                                          # Return the validated, normalized seed row.


# ----------------------------- ESMFold + Baseline Qualification gate ----------------------------- #


def _import_esmfold():
    """Import the ESMFold model lazily with a clear error if missing."""
    try:                                                                                                                # ESMFold pulls in ~8.4 GB of weights; defer the cost until we actually need it.
        import torch                                                                                                    # Torch is required for tensor operations.
        from transformers import AutoTokenizer, EsmForProteinFolding                                                    # ESMFold is exposed via the HuggingFace transformers package.
    except Exception as exc:                                                                                            # pragma: no cover - exercised only in the deploy environment.
        raise ImportError(                                                                                              # Surface a helpful install hint instead of a cryptic ImportError.
            "Stage 11 requires `torch` and `transformers` to load ESMFold. "
            "Install them before running 11a (`pip install torch transformers`)."
        ) from exc
    return torch, AutoTokenizer, EsmForProteinFolding                                                                   # Hand the imports back to the caller in a stable shape.


def _choose_device(torch_mod, requested: str) -> str:
    """Resolve the requested device, falling back to CPU if CUDA is unavailable."""
    requested = (requested or "auto").lower()                                                                           # Normalize the input string.
    if requested in {"auto", ""}:                                                                                       # Auto means "pick the best available device".
        return "cuda" if torch_mod.cuda.is_available() else "cpu"                                                       # Use CUDA if any GPU is visible, else CPU.
    if requested == "cuda" and not torch_mod.cuda.is_available():                                                       # The operator asked for CUDA but no GPU is present.
        print("[WARN] Requested cuda but no GPU is available; falling back to cpu.", flush=True)                        # Warn loudly so they don't blame the script for sluggish folding.
        return "cpu"                                                                                                    # Down-grade gracefully.
    return requested                                                                                                    # Trust the operator's choice when the device is available.


def _parse_pdb_plddt(pdb_text: str) -> tuple[list[float], int]:
    """Pull per-residue pLDDT (Cα B-factor) values out of an ESMFold PDB string."""
    per_res: list[float] = []                                                                                           # Accumulator for the per-residue confidence values.
    seen_residues: set[tuple[str, int]] = set()                                                                         # Dedup tracker (chain_id, residue_number) so we only count Cα once.
    for line in pdb_text.splitlines():                                                                                  # Walk the PDB text line by line.
        if not (line.startswith("ATOM") or line.startswith("HETATM")):                                                  # Only ATOM/HETATM records carry coordinates and B-factors.
            continue                                                                                                    # Skip headers, REMARK lines, etc.
        atom_name = line[12:16].strip()                                                                                 # The PDB column for atom name (e.g. "CA").
        if atom_name != "CA":                                                                                           # Only count the Cα atom per residue (one pLDDT per residue).
            continue                                                                                                    # Skip side-chain atoms.
        chain_id = line[21]                                                                                             # Chain identifier column.
        try:                                                                                                            # The residue sequence number may be malformed in rare cases.
            res_num = int(line[22:26].strip())                                                                          # 4-character residue number column.
            b_factor = float(line[60:66].strip())                                                                       # B-factor column where ESMFold writes pLDDT.
        except (ValueError, IndexError):                                                                                # Malformed line — skip rather than crash.
            continue                                                                                                    # Defensive parser.
        key = (chain_id, res_num)                                                                                       # Dedup key for the residue.
        if key in seen_residues:                                                                                        # We've already recorded this residue's Cα pLDDT.
            continue                                                                                                    # Skip the duplicate.
        seen_residues.add(key)                                                                                          # Mark this residue as seen.
        per_res.append(b_factor*100)                                                                                        # Collect the pLDDT value.
    return per_res, len(per_res)                                                                                        # Return the values and the residue count.


def esmfold_single_sequence(
    sequence: str,
    out_pdb_path: str | Path,
    device: str = "cuda",
    chunk_size: int = 128,
    num_recycles: int = 1,
) -> dict:
    """ESMFold a single sequence with facebook/esmfold_v1 and persist a PDB.

    Returns a dict with: mean_plddt (in 0-100 scale to match Stage 08 convention),
    per_residue_plddt list, n_residues, elapsed_seconds, device_used. pLDDT is
    parsed back from the PDB B-factor column rather than re-fetched from the
    PyTorch output object — this is the same convention Stage 08 uses, so the
    Baseline Qualification gate's threshold (default 70.0) is directly
    comparable to the Stage 08 oracle.
    """
    import time                                                                                                         # Defer the timing import to keep top-of-file imports lean.

    torch, AutoTokenizer, EsmForProteinFolding = _import_esmfold()                                                      # Lazy-load the heavy ESMFold dependencies.
    resolved_device = _choose_device(torch, device)                                                                     # Pick CUDA if available, otherwise CPU.

    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")                                                    # Pull the ESMFold tokenizer.
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)                         # Pull the ESMFold model with low-memory loading.
    model = model.eval().to(resolved_device)                                                                            # Move to the right device and lock into inference mode.
    try:                                                                                                                # The chunk/recycles knobs may not exist on every version.
        model.trunk.set_chunk_size(int(chunk_size))                                                                     # Smaller chunk_size = lower VRAM footprint at the cost of throughput.
    except Exception:                                                                                                   # Older transformers versions ignore the call gracefully.
        pass                                                                                                            # Continue without chunking if the API isn't present.

    inputs = tokenizer([str(sequence)], return_tensors="pt", add_special_tokens=False)["input_ids"].to(resolved_device) # Tokenize the seed (no special tokens for ESMFold).
    start_time = time.time()                                                                                            # Record wall-clock start for downstream reporting.
    try:                                                                                                                # Wrap the heavy inference so we can attach a hint if OOM strikes.
        with torch.no_grad():                                                                                           # Inference path — no gradients needed.
            output = model.infer_pdb(str(sequence))                                                                     # Run ESMFold; output is a PDB string.
    except RuntimeError as exc:                                                                                         # OOM on GPU surfaces here as a RuntimeError.
        if "out of memory" in str(exc).lower():                                                                         # Surface a targeted hint when this is the specific failure.
            raise RuntimeError(                                                                                         # Re-raise with the hint attached.
                f"ESMFold OOM while folding the seed. Reduce --esmfold_chunk_size "
                f"(current={chunk_size}) or use --esmfold_device cpu. Original: {exc}"
            ) from exc
        raise                                                                                                           # Re-raise unrelated errors unchanged.
    elapsed = float(time.time() - start_time)                                                                           # Record wall-clock end.

    out_pdb_path = Path(out_pdb_path)                                                                                   # Ensure the destination is a Path object.
    out_pdb_path.parent.mkdir(parents=True, exist_ok=True)                                                              # Create the destination directory tree if needed.
    out_pdb_path.write_text(output)                                                                                     # Persist the PDB text exactly as ESMFold emitted it.

    per_res, n_res = _parse_pdb_plddt(output)                                                                           # Read pLDDT back from the PDB so the value is comparable to Stage 08.
    if n_res == 0:                                                                                                      # Defensive: ESMFold should always emit at least one Cα.
        raise RuntimeError(f"ESMFold produced no Cα atoms in the PDB at {out_pdb_path}; cannot compute pLDDT.")         # Refuse to continue without confidence values.
    mean_plddt = float(np.mean(per_res))                                                                                # The headline metric used by the Baseline Qualification gate.

    plddt_arr = np.asarray(per_res, dtype=np.float32)                                                                   # Cast to numpy for vectorized stats.
    summary = {                                                                                                         # Compact stats consumed by the report layer.
        "min": float(plddt_arr.min()),
        "p10": float(np.percentile(plddt_arr, 10)),
        "median": float(np.percentile(plddt_arr, 50)),
        "max": float(plddt_arr.max()),
    }

    return {                                                                                                            # Bundle every output the gate / context need.
        "mean_plddt": mean_plddt,                                                                                       # Headline metric used by the gate.
        "per_residue_plddt": [float(v) for v in per_res],                                                               # Cast each value to a plain float for JSON serialization.
        "per_residue_plddt_summary": summary,                                                                           # Pre-computed summary for the report.
        "n_residues": int(n_res),                                                                                       # Total residues in the folded structure.
        "elapsed_seconds": elapsed,                                                                                     # Wall-clock cost (useful for SageMaker cost accounting).
        "device_used": resolved_device,                                                                                 # Record where the work happened (cpu vs cuda).
        "pdb_path": str(out_pdb_path.resolve()),                                                                        # Absolute path so downstream stages need no relative-path guessing.
    }


def baseline_qualification_gate(seed_metrics: dict, min_plddt: float) -> tuple[bool, str]:
    """The single most important guard in Stage 11.

    Returns (passed, reason). Pass means the wild-type seed folded above the
    pLDDT threshold and is structurally worth redesigning. Fail means the seed
    chassis is broken and Stage 11 must abort with EXIT_BASELINE_QUALIFICATION_FAILED
    so the operator picks a different `--seed_protein_id` rather than wasting
    compute on a corrupted scaffold (the entire reason Stage 10 failed).
    """
    actual = float(seed_metrics.get("mean_plddt", float("nan")))                                                        # Pull the headline confidence from the ESMFold output dict.
    threshold = float(min_plddt)                                                                                        # Cast the threshold to a clean float.
    if not math.isfinite(actual):                                                                                       # NaN/inf means folding silently produced no usable confidence.
        return False, f"Seed mean_plddt is not finite (got {actual!r}); refusing to redesign."                          # Defensive failure with a clear reason.
    if actual < threshold:                                                                                              # The core gate check.
        return False, (                                                                                                 # Compose a precise rejection message naming both numbers.
            f"Seed mean pLDDT={actual:.3f} is below the Baseline Qualification "
            f"threshold ({threshold:.3f}). Pick a different seed and rerun."
        )
    return True, f"Seed mean pLDDT={actual:.3f} passes threshold {threshold:.3f}."                                      # Acceptance message with both numbers.


# ----------------------------- Family / target context builders ----------------------------- #


def load_strict_embeddings(embeddings_pt: str | Path, index_csv: str | Path) -> tuple[np.ndarray, pd.DataFrame]:
    """Load cached ESM-2 strict embeddings + their protein_id index.

    The .pt file is a torch tensor of shape [N, D] written by scripts/02_embed_rbps.py.
    The index CSV maps row_id → protein_id. Returns (embeddings, index_df) where
    `embeddings` is a numpy float32 array of shape [N, D] and `index_df` has at
    least the columns: virus_accession, host_genus, protein_id, product, row_id.
    """
    try:                                                                                                                # Torch is required to deserialize the .pt file.
        import torch                                                                                                    # Local import keeps the heavy dependency lazy.
    except Exception as exc:                                                                                            # Pragma: this branch only fires in misconfigured environments.
        raise ImportError("Stage 11 requires `torch` to read the cached ESM-2 embeddings.") from exc                    # Surface a clear install hint.
    embeddings_pt = Path(embeddings_pt)                                                                                 # Standardize the path type.
    index_csv = Path(index_csv)                                                                                         # Standardize the path type.
    if not embeddings_pt.exists():                                                                                      # Refuse to proceed without the embeddings.
        raise FileNotFoundError(f"Missing strict embeddings tensor: {embeddings_pt}")                                   # Loud error.
    if not index_csv.exists():                                                                                          # Refuse to proceed without the index.
        raise FileNotFoundError(f"Missing strict embeddings index CSV: {index_csv}")                                    # Loud error.
    try:                                                                                                                # torch.load on the deploy box may or may not support weights_only.
        tensor = torch.load(embeddings_pt, map_location="cpu", weights_only=False)                                      # Force CPU load so we don't hold GPU memory we don't need.
    except TypeError:                                                                                                   # Older torch versions don't accept weights_only.
        tensor = torch.load(embeddings_pt, map_location="cpu")                                                          # Fall back to the legacy signature.
    embeddings = np.asarray(tensor.detach().cpu().numpy(), dtype=np.float32)                                            # Move to numpy float32 for downstream math.
    index_df = pd.read_csv(index_csv)                                                                                   # Read the row_id → metadata mapping.
    if len(embeddings) != len(index_df):                                                                                # Refuse mismatched lengths — that means the cache is stale.
        raise ValueError(                                                                                               # Surface the exact mismatch.
            f"Embeddings ({len(embeddings)} rows) and index ({len(index_df)} rows) "
            f"have different lengths; re-run scripts/02_embed_rbps.py."
        )
    return embeddings, index_df                                                                                         # Hand the aligned (matrix, dataframe) pair to the caller.


def embed_strict_set_fresh(
    strict_df: pd.DataFrame,
    embedding_model: str,
    batch_size: int = 4,
    torch=None,
    tokenizer=None,
    model=None,
    device: str | None = None,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Re-embed the entire strict CSV from scratch using the persistent ESM-2 backend.

    Exposed for users who want strict purity ("compute everything, inherit
    nothing"). Returns the same (embeddings, index_df) shape as
    `load_strict_embeddings`, so callers can swap one for the other without
    changing downstream code paths.
    """
    sequences = strict_df["aa_sequence"].astype(str).tolist()                                                           # Pull every strict-CSV sequence into a flat list.
    if tokenizer is not None and model is not None and torch is not None and device is not None:                        # Prefer the persistent backend if 11a already loaded it.
        embeddings = embed_sequences_with_backend(                                                                      # Use the high-performance path.
            sequences, torch, tokenizer, model, device=device, batch_size=batch_size,                                   # Pass the live components.
        )
    else:                                                                                                               # Fall back to the legacy from-scratch embedder.
        embeddings = embed_sequences(sequences, model_name=embedding_model, batch_size=batch_size)                      # Slower but self-contained.
    index_df = strict_df[["virus_accession", "host_genus", "protein_id", "product"]].copy()                             # Mirror the column shape of the cached index.
    index_df["row_id"] = np.arange(len(index_df))                                                                       # Stable row index for downstream lookups.
    return np.asarray(embeddings, dtype=np.float32), index_df                                                           # Return the freshly computed pair.


def _length_tolerance_mask(length_seed: int, lengths: np.ndarray, tolerance: float) -> np.ndarray:
    """Boolean mask: True for rows whose length differs from `length_seed` by ≤ `tolerance` × length_seed."""
    diffs = np.abs(lengths - length_seed) / max(length_seed, 1)                                                         # Relative-length deviation per row.
    return diffs <= float(tolerance)                                                                                    # Keep rows within tolerance.


def find_family_members(
    seed_embedding: np.ndarray,
    seed_protein_id: str,
    strict_embeddings: np.ndarray,
    strict_index: pd.DataFrame,
    strict_df: pd.DataFrame,
    top_n: int = 32,
    cosine_floor: float = 0.85,
    length_tolerance: float = 0.05,
) -> pd.DataFrame:
    """Top-N family members by ESM-2 cosine to the seed, gated by length similarity.

    Includes the seed itself as member rank 1 with cosine 1.0. Returns a DataFrame
    with columns: protein_id, host_genus, product, virus_accession, aa_sequence,
    cosine_to_seed, length_diff_frac, family_rank.
    """
    # --- Cosine to every strict row ---
    seed_norm = seed_embedding / max(float(np.linalg.norm(seed_embedding)), 1e-12)                                      # L2-normalize the seed vector for cosine math.
    matrix_norm = strict_embeddings / np.clip(np.linalg.norm(strict_embeddings, axis=1, keepdims=True), 1e-12, None)    # L2-normalize every row of the strict matrix.
    cosines = matrix_norm @ seed_norm                                                                                   # Dot product = cosine because both sides are unit-norm.

    # --- Pull sequence + length features by joining the index with the strict CSV on protein_id ---
    annotated = strict_index.merge(                                                                                     # Join the index dataframe with the strict CSV.
        strict_df[["protein_id", "aa_sequence"]],                                                                       # Bring in the sequence column.
        on="protein_id",                                                                                                # Match on protein_id.
        how="left",                                                                                                     # Keep every row in the index.
    ).copy()                                                                                                            # Detach for mutation safety.
    annotated["cosine_to_seed"] = cosines                                                                               # Attach the cosine column.
    annotated["aa_length"] = annotated["aa_sequence"].astype(str).str.len()                                             # Compute per-row sequence length.
    seed_length = int(annotated.loc[annotated["protein_id"] == seed_protein_id, "aa_length"].iloc[0])                   # Pull the seed's length out of the joined frame.
    annotated["length_diff_frac"] = (annotated["aa_length"] - seed_length).abs() / max(seed_length, 1)                  # Relative length deviation.

    # --- Apply the length filter and the cosine floor ---
    length_mask = annotated["length_diff_frac"].to_numpy() <= float(length_tolerance)                                   # True where length is within tolerance.
    cosine_mask = annotated["cosine_to_seed"].to_numpy() >= float(cosine_floor)                                         # True where cosine meets the floor.
    is_seed = annotated["protein_id"].to_numpy() == str(seed_protein_id)                                                # Always retain the seed itself regardless of filters.
    keep_mask = (length_mask & cosine_mask) | is_seed                                                                   # Survivors = (length OK AND cosine OK) OR the seed.
    pool = annotated.loc[keep_mask].copy()                                                                              # Build the surviving pool.

    # --- Rank survivors by cosine descending and pick top-N (seed pinned at rank 1) ---
    pool["is_seed"] = pool["protein_id"] == str(seed_protein_id)                                                        # Boolean marker so we can sort the seed to the top.
    pool = pool.sort_values(["is_seed", "cosine_to_seed"], ascending=[False, False]).reset_index(drop=True)             # Seed first, then highest cosine.
    pool = pool.head(int(top_n)).copy()                                                                                 # Truncate to the requested size.
    pool["family_rank"] = np.arange(1, len(pool) + 1)                                                                   # Re-index 1..N for human readability.
    return pool[["protein_id", "host_genus", "product", "virus_accession",                                              # Return a tidy, predictable column order.
                  "aa_sequence", "cosine_to_seed", "length_diff_frac", "family_rank"]]


def find_target_members(
    seed_embedding: np.ndarray,
    seed_protein_id: str,
    target_host: str,
    strict_embeddings: np.ndarray,
    strict_index: pd.DataFrame,
    strict_df: pd.DataFrame,
    top_m: int = 8,
    length_tolerance: float = 0.05,
) -> pd.DataFrame:
    """Top-M target-host RBPs by ESM-2 cosine to the seed, length-matched.

    Returns a DataFrame with the same columns as `find_family_members`. The seed
    is **not** included (we want the target manifold, not the source manifold).
    """
    # --- Same cosine math as the family helper ---
    seed_norm = seed_embedding / max(float(np.linalg.norm(seed_embedding)), 1e-12)                                      # L2-normalize the seed vector.
    matrix_norm = strict_embeddings / np.clip(np.linalg.norm(strict_embeddings, axis=1, keepdims=True), 1e-12, None)    # L2-normalize every row of the matrix.
    cosines = matrix_norm @ seed_norm                                                                                   # Cosine via unit-norm dot product.

    annotated = strict_index.merge(                                                                                     # Join with the strict CSV to pull sequences.
        strict_df[["protein_id", "aa_sequence"]],                                                                       # Only the columns we actually need.
        on="protein_id",                                                                                                # Match by protein_id.
        how="left",                                                                                                     # Keep all index rows.
    ).copy()                                                                                                            # Detach for mutation safety.
    annotated["cosine_to_seed"] = cosines                                                                               # Attach cosine values.
    annotated["aa_length"] = annotated["aa_sequence"].astype(str).str.len()                                             # Per-row length.
    seed_length = int(annotated.loc[annotated["protein_id"] == seed_protein_id, "aa_length"].iloc[0])                   # Seed length.
    annotated["length_diff_frac"] = (annotated["aa_length"] - seed_length).abs() / max(seed_length, 1)                  # Relative length deviation.

    # --- Restrict to target-host rows, length-matched ---
    target_mask = annotated["host_genus"].astype(str) == str(target_host)                                               # Keep only rows matching the target host.
    length_mask = annotated["length_diff_frac"].to_numpy() <= float(length_tolerance)                                   # Apply length tolerance.
    pool = annotated.loc[target_mask & length_mask].copy()                                                              # Build the surviving pool.

    if len(pool) < int(top_m):                                                                                          # Soft warning if length filtering shrank the pool below quota.
        print(                                                                                                          # Loud heads-up but not a fatal error.
            f"[WARN] Target-host '{target_host}' has only {len(pool)} length-matched rows "
            f"(requested top_m={int(top_m)}). Consider widening --length_tolerance.",
            flush=True,
        )
    pool = pool.sort_values("cosine_to_seed", ascending=False).reset_index(drop=True)                                   # Rank by cosine descending.
    pool = pool.head(int(top_m)).copy()                                                                                 # Truncate to top-M.
    pool["family_rank"] = np.arange(1, len(pool) + 1)                                                                   # Re-index for readability.
    return pool[["protein_id", "host_genus", "product", "virus_accession",                                              # Same column shape as the family helper.
                  "aa_sequence", "cosine_to_seed", "length_diff_frac", "family_rank"]]


def compute_centroid(embeddings: np.ndarray) -> np.ndarray:
    """L2-normalize rows then mean. Returns a 1-D float32 vector of shape (D,)."""
    if len(embeddings) == 0:                                                                                            # Defensive: empty input → zero-length vector.
        return np.zeros((0,), dtype=np.float32)                                                                         # Surface the empty case cleanly.
    normed = normalize_rows(np.asarray(embeddings, dtype=np.float32))                                                   # Unit-norm every row before averaging.
    return np.asarray(normed.mean(axis=0), dtype=np.float32)                                                            # Centroid in unit-vector space.


# ----------------------------- Per-position frequencies + entropy ----------------------------- #


def position_frequencies(sequences: Sequence[str]) -> list[dict[str, float]]:
    """Per-position AA frequencies across a list of equal-length sequences.

    Length is inferred from the first sequence; rows shorter than that length
    contribute their available residues only. Non-canonical characters are
    excluded from the count (the frequencies are computed over valid AAs only).
    """
    if not sequences:                                                                                                   # No data → empty list.
        return []                                                                                                       # Safe early return.
    target_length = len(sequences[0])                                                                                   # Take the seed length as the reference.
    freqs: list[dict[str, float]] = []                                                                                  # Result accumulator.
    for i in range(target_length):                                                                                      # Iterate over each column.
        counts: dict[str, int] = {}                                                                                     # Reset the counter for this column.
        total = 0                                                                                                       # Count valid observations at this column.
        for seq in sequences:                                                                                           # Walk every aligned sequence.
            if i < len(seq):                                                                                            # Skip rows that are too short here.
                aa = seq[i]                                                                                             # Pull the AA at column i.
                if aa in VALID_AA:                                                                                      # Ignore gaps, '*', X, etc.
                    counts[aa] = counts.get(aa, 0) + 1                                                                  # Increment the AA's count.
                    total += 1                                                                                          # Increment the valid-observation counter.
        if total == 0:                                                                                                  # No valid observations at this column.
            freqs.append({})                                                                                            # Empty dict → entropy-zero, no preferred AAs.
            continue                                                                                                    # Move on to the next column.
        freqs.append({aa: count / total for aa, count in counts.items()})                                               # Normalize to frequencies.
    return freqs                                                                                                        # One dict per column.


def normalized_shannon_entropy_per_position(freqs: list[dict[str, float]]) -> np.ndarray:
    """Per-position Shannon entropy normalized to [0, 1] (max = log2(20))."""
    if not freqs:                                                                                                       # Empty input → empty output.
        return np.zeros((0,), dtype=np.float32)                                                                         # Safe early return.
    max_entropy = math.log2(len(VALID_AA))                                                                              # Maximum possible entropy for 20 amino acids.
    entropy = np.zeros(len(freqs), dtype=np.float32)                                                                    # Result vector.
    for i, dist in enumerate(freqs):                                                                                    # Walk every column.
        if not dist:                                                                                                    # Empty distribution → entropy 0.
            continue                                                                                                    # Already zero, leave it.
        probs = np.asarray(list(dist.values()), dtype=np.float32)                                                       # Probability vector.
        probs = probs / max(probs.sum(), 1e-12)                                                                         # Defensive renormalization in case of tiny drift.
        h = float(-(probs * np.log2(probs + 1e-12)).sum())                                                              # Standard Shannon entropy in bits.
        entropy[i] = min(h / max_entropy, 1.0)                                                                          # Normalize to [0, 1], clamp top.
    return entropy                                                                                                      # One float per column.


def top_k_preferences(
    freqs: list[dict[str, float]],
    k: int = 4,
    exclude_aa: list[str] | None = None,
) -> list[list[str]]:
    """Top-K AAs per position by frequency, with optional per-position exclusions.

    If `exclude_aa` is provided it must align 1-to-1 with `freqs`; that AA is
    removed from the candidate list at the corresponding position (typically
    used to exclude the seed AA itself so we don't propose a no-op mutation).
    """
    result: list[list[str]] = []                                                                                        # Result accumulator.
    for i, dist in enumerate(freqs):                                                                                    # Walk every column.
        if not dist:                                                                                                    # Empty distribution → no preferences.
            result.append([])                                                                                           # Empty list at this position.
            continue                                                                                                    # Move on.
        excluded = (exclude_aa[i] if exclude_aa and i < len(exclude_aa) else None)                                      # Per-position exclusion (typically the seed AA).
        ranked = sorted(dist.items(), key=lambda item: (-item[1], item[0]))                                             # Sort by frequency desc, AA asc (deterministic tie-break).
        filtered = [aa for aa, _ in ranked if aa != excluded]                                                           # Drop the excluded AA.
        result.append(filtered[: int(k)])                                                                               # Keep top-K.
    return result                                                                                                       # One list per position.


# ----------------------------- Edit-space construction ----------------------------- #


def build_stage11_edit_proposals(
    seed_sequence: str,
    family_aligned_seqs: list[str],
    target_aligned_seqs: list[str],
    entropy_floor: float = 0.20,
    family_top_k: int = 4,
    target_top_k: int = 4,
    max_allowed_aas_per_pos: int = 6,
    region_block: int = 50,
) -> list[EditProposal]:
    """Build the from-scratch EditProposal list. No Stage 06/07 context required.

    Algorithm (mirroring §10 of the design doc):

    1. Compute per-position family + target AA frequencies.
    2. Compute normalized Shannon entropy from family frequencies.
    3. For every position with entropy >= entropy_floor:
         - allowed_aas = union(top-K family AAs, top-K target AAs) minus seed_aa,
           capped at max_allowed_aas_per_pos.
         - functional_weight = entropy[i] * (1 - family_freqs[i].get(seed_aa, 0)).
         - conservation_penalty = 1 - entropy[i].
         - region_name = stage11_window_<i // region_block>.
    """
    if not family_aligned_seqs:                                                                                         # Refuse to build an edit space with no family signal.
        raise ValueError("Cannot build Stage 11 edit space: family_aligned_seqs is empty.")                             # Loud error.
    L = len(seed_sequence)                                                                                              # Convenience: seed length.

    family_aligned = [seq[:L] if len(seq) >= L else seq.ljust(L, "X") for seq in family_aligned_seqs]                   # Truncate or right-pad to seed length so columns align.
    target_aligned = [seq[:L] if len(seq) >= L else seq.ljust(L, "X") for seq in target_aligned_seqs]                   # Same handling for the target pool.

    family_freqs = position_frequencies(family_aligned)                                                                 # Per-position family AA frequencies.
    target_freqs = position_frequencies(target_aligned)                                                                 # Per-position target-host AA frequencies.
    entropy = normalized_shannon_entropy_per_position(family_freqs)                                                     # Per-position entropy from family signal.

    seed_aas = list(seed_sequence)                                                                                      # Per-position exclusion list (the seed AA itself).
    family_top = top_k_preferences(family_freqs, k=int(family_top_k), exclude_aa=seed_aas)                              # Top-K family AAs per position, excluding seed AA.
    target_top = top_k_preferences(target_freqs, k=int(target_top_k), exclude_aa=seed_aas)                              # Top-K target AAs per position, excluding seed AA.

    proposals: list[EditProposal] = []                                                                                  # Result accumulator.
    for i in range(L):                                                                                                  # 0-based iteration; we'll convert to 1-based positions on the fly.
        if entropy[i] < float(entropy_floor):                                                                           # Too conserved to safely edit.
            continue                                                                                                    # Skip.
        seed_aa = seed_sequence[i]                                                                                      # The current wild-type residue.
        if seed_aa not in VALID_AA:                                                                                     # Defensive: skip non-canonical positions in the seed itself.
            continue                                                                                                    # Skip.
        allowed = []                                                                                                    # Build the ordered allowed-AA list.
        for aa in family_top[i] + target_top[i]:                                                                        # Walk family preferences then target preferences.
            if aa in VALID_AA and aa != seed_aa and aa not in allowed:                                                  # Skip seed AA, duplicates, and non-canonical letters.
                allowed.append(aa)                                                                                      # Append in priority order.
            if len(allowed) >= int(max_allowed_aas_per_pos):                                                            # Cap branching per position.
                break                                                                                                   # Stop adding AAs.
        if not allowed:                                                                                                 # No legal substitutions at this position.
            continue                                                                                                    # Skip rather than emit an empty proposal.
        family_pref = {aa: float(family_freqs[i].get(aa, 0.0)) for aa in family_top[i]}                                 # Family preferences dict for the proposal.
        target_pref = {aa: float(target_freqs[i].get(aa, 0.0)) if i < len(target_freqs) else 0.0                        # Target preferences dict for the proposal.
                        for aa in target_top[i]}
        functional_weight = float(entropy[i]) * (1.0 - float(family_freqs[i].get(seed_aa, 0.0)))                        # Variable position AND seed is an outlier here.
        conservation_penalty = 1.0 - float(entropy[i])                                                                  # Inverse of entropy.
        region_name = f"stage11_window_{i // int(region_block)}"                                                        # Coarse grouping for spreading.

        proposals.append(                                                                                               # Build and append the EditProposal.
            EditProposal(
                position=i + 1,                                                                                         # Convert to 1-based biology position.
                seed_aa=seed_aa,                                                                                        # Wild-type AA at this position.
                allowed_aas=allowed,                                                                                    # The capped list of legal substitutions.
                target_preference=target_pref,                                                                          # Target-host AA preferences for blending.
                family_preference=family_pref,                                                                          # Family AA preferences for blending.
                functional_weight=functional_weight,                                                                    # Composite priority signal.
                conservation_penalty=conservation_penalty,                                                              # Risk signal (high for conserved positions).
                region_name=region_name,                                                                                # Geographic bucket for spreading.
            )
        )
    return proposals                                                                                                    # The full set of legal edit proposals.


def stage11_choose_editable_positions(
    proposals: list[EditProposal],
    max_hard: int,
    max_soft: int,
    seed: int,
) -> tuple[list[int], list[int], list[int]]:
    """Select (hard, soft, frozen_remaining) editable positions with regional spreading.

    Hard positions are the highest-functional_weight proposal in each region, up
    to `max_hard` positions. Soft positions are the next-best by functional_weight
    overall, skipping anything already chosen as hard, up to `max_soft`. The
    `frozen_remaining` list is purely informational at this layer — the caller
    is responsible for unioning into seed length and writing the final frozen list.
    """
    if not proposals:                                                                                                   # Defensive: refuse empty input.
        return [], [], []                                                                                               # Empty triple.
    rng = random.Random(int(seed))                                                                                      # Reproducible deterministic RNG (only used for tie-breaking).
    by_region: dict[str, list[EditProposal]] = {}                                                                       # Bucket the proposals by region.
    for item in proposals:                                                                                              # Walk every proposal once.
        by_region.setdefault(item.region_name, []).append(item)                                                         # Append to its region's bucket.
    for items in by_region.values():                                                                                    # Sort each bucket by descending priority.
        items.sort(key=lambda x: (x.functional_weight, -x.conservation_penalty, x.position), reverse=True)              # Highest functional_weight, lowest conservation penalty, then position.

    hard: list[int] = []                                                                                                # Selected hard positions.
    region_order = sorted(by_region, key=lambda name: max(it.functional_weight for it in by_region[name]), reverse=True)# Sort regions by the strength of their best member.
    for region_name in region_order:                                                                                    # Walk regions in priority order.
        bucket = by_region[region_name]                                                                                 # The proposals in this region.
        for item in bucket:                                                                                             # Iterate inside the region.
            pos = int(item.position)                                                                                    # Position as int.
            if pos in hard:                                                                                             # Skip duplicates if any.
                continue                                                                                                # Already chosen.
            hard.append(pos)                                                                                            # Lock in the region's best position.
            break                                                                                                       # One position per region in this pass.
        if len(hard) >= int(max_hard):                                                                                  # Cap reached.
            break                                                                                                       # Stop after the cap is full.

    if len(hard) < int(max_hard):                                                                                       # If we still have room, fill greedily by global priority.
        leftover = [it for it in proposals if int(it.position) not in hard]                                             # Pool of unchosen proposals.
        leftover.sort(key=lambda x: (x.functional_weight, -x.conservation_penalty, x.position), reverse=True)           # Global priority order.
        for it in leftover:                                                                                             # Walk top-priority leftovers.
            hard.append(int(it.position))                                                                               # Promote to hard.
            if len(hard) >= int(max_hard):                                                                              # Cap reached.
                break                                                                                                   # Done.

    hard_set = set(hard)                                                                                                # Fast lookup for soft selection.
    soft_pool = [it for it in proposals if int(it.position) not in hard_set]                                            # Remaining proposals.
    soft_pool.sort(key=lambda x: (x.functional_weight, -x.conservation_penalty, x.position), reverse=True)              # Sort by priority.
    soft: list[int] = [int(it.position) for it in soft_pool[: int(max_soft)]]                                           # Top-N go to soft.

    hard = sorted(set(hard))                                                                                            # Final hard list — sorted, deduped.
    soft = sorted(set(soft))                                                                                            # Final soft list — sorted, deduped.
    return hard, soft, []                                                                                               # Frozen is computed by the caller against seed length.


def serialize_edit_proposals(
    proposals: list[EditProposal],
    hard_positions: list[int],
    soft_positions: list[int],
) -> list[dict]:
    """Convert EditProposal objects into JSON-safe dicts with an `edit_tier` tag."""
    hard_set = set(int(p) for p in hard_positions)                                                                      # Fast lookup.
    soft_set = set(int(p) for p in soft_positions)                                                                      # Fast lookup.
    rows: list[dict] = []                                                                                               # Result accumulator.
    for item in proposals:                                                                                              # Walk every proposal.
        pos = int(item.position)                                                                                        # 1-based position as int.
        if pos in hard_set:                                                                                             # Determine the tier.
            tier = "hard"                                                                                                # Primary edit position.
        elif pos in soft_set:                                                                                           # Soft buffer.
            tier = "soft"                                                                                                # Secondary edit position.
        else:                                                                                                           # Discarded.
            tier = "frozen"                                                                                              # Excluded from editing.
        rows.append(                                                                                                    # Serialize to a JSON-safe dict.
            {
                "position": pos,                                                                                        # 1-based position.
                "seed_aa": str(item.seed_aa),                                                                           # Wild-type AA.
                "allowed_aas": list(item.allowed_aas),                                                                  # Legal substitutions.
                "target_preference": {str(k): float(v) for k, v in item.target_preference.items()},                     # Target preferences dict.
                "family_preference": {str(k): float(v) for k, v in item.family_preference.items()},                     # Family preferences dict.
                "functional_weight": float(item.functional_weight),                                                     # Composite priority signal.
                "conservation_penalty": float(item.conservation_penalty),                                               # Risk signal.
                "region_name": str(item.region_name),                                                                   # Region bucket name.
                "edit_tier": tier,                                                                                      # hard / soft / frozen.
            }
        )
    return rows                                                                                                         # Done.


# ----------------------------- Run metadata writer ----------------------------- #


def write_run_metadata(
    out_path: str | Path,
    cli_args: dict,
    extras: dict | None = None,
) -> None:
    """Capture every CLI arg, Python info, GPU info, package versions, UTC timestamp.

    Best-effort: any subcomponent (git, pynvml, torch) is wrapped in a try/except
    so the metadata write succeeds even on minimal environments.
    """
    info: dict = {                                                                                                      # Build the metadata dict incrementally.
        "stage": "11a",                                                                                                 # Tag this artifact.
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),               # ISO-8601 UTC stamp.
        "cli_args": {str(k): _jsonable(v) for k, v in cli_args.items()},                                                # Best-effort JSON-safe rendering of CLI args.
        "python": {                                                                                                     # Interpreter info.
            "version": sys.version,                                                                                     # Full Python version string.
            "executable": sys.executable,                                                                               # Absolute path to the interpreter.
            "platform": platform.platform(),                                                                            # Human-readable platform string.
        },
        "host": {                                                                                                       # Host machine info.
            "hostname": socket.gethostname(),                                                                           # Hostname (useful for SageMaker job IDs).
        },
    }
    try:                                                                                                                # Pull GPU info if available.
        import torch                                                                                                    # Local import — torch is optional for metadata writing.
        gpu = {                                                                                                         # GPU sub-dict.
            "cuda_available": bool(torch.cuda.is_available()),                                                          # Boolean.
            "torch_version": str(torch.__version__),                                                                    # Torch version.
        }
        if torch.cuda.is_available():                                                                                   # Add GPU name if a GPU is present.
            gpu["gpu_name"] = str(torch.cuda.get_device_name(0))                                                        # First device name.
            gpu["gpu_count"] = int(torch.cuda.device_count())                                                           # Total visible devices.
        info["gpu"] = gpu                                                                                                # Attach to the metadata dict.
    except Exception:                                                                                                   # Torch isn't installed — skip silently.
        info["gpu"] = {"cuda_available": False, "note": "torch not importable"}                                          # Note the missing torch.

    info["packages"] = _collect_package_versions(["torch", "transformers", "fair-esm",                                  # Capture pinned versions for reproducibility.
                                                   "scikit-learn", "biopython", "numpy", "pandas"])

    info["git_commit"] = _try_git_commit()                                                                              # Best-effort git rev.

    if extras:                                                                                                          # Merge any caller-provided extras.
        info.update(extras)                                                                                              # Last-write wins so callers can override.

    write_json(info, out_path)                                                                                          # Persist as JSON.


def _jsonable(value):
    """Best-effort coercion of arbitrary values into JSON-safe primitives."""
    if isinstance(value, (str, int, float, bool)) or value is None:                                                     # Primitives pass through.
        return value                                                                                                    # No change.
    if isinstance(value, Path):                                                                                         # Paths to strings.
        return str(value)                                                                                               # Native string representation.
    if isinstance(value, (list, tuple)):                                                                                # Lists/tuples recurse.
        return [_jsonable(v) for v in value]                                                                            # Recurse.
    if isinstance(value, dict):                                                                                         # Dicts recurse.
        return {str(k): _jsonable(v) for k, v in value.items()}                                                         # Recurse.
    return str(value)                                                                                                   # Fallback: repr-ish string.


def _collect_package_versions(names: list[str]) -> dict[str, str]:
    """Best-effort `package: version` mapping. Missing packages map to 'not_installed'."""
    try:                                                                                                                # importlib.metadata is the modern entry point.
        from importlib.metadata import version, PackageNotFoundError                                                    # Python 3.8+ standard library.
    except Exception:                                                                                                   # Pragma: extremely old Pythons only.
        return {name: "unknown" for name in names}                                                                      # Bail gracefully.
    out: dict[str, str] = {}                                                                                            # Result accumulator.
    for name in names:                                                                                                  # Walk every requested package.
        try:                                                                                                            # Each package may or may not be installed.
            out[name] = str(version(name))                                                                              # Successfully resolved.
        except PackageNotFoundError:                                                                                    # Package isn't installed.
            out[name] = "not_installed"                                                                                  # Sentinel for the report.
        except Exception:                                                                                               # Catch-all for distro weirdness.
            out[name] = "unknown"                                                                                        # Mark unknown rather than crash.
    return out                                                                                                          # Done.


def _try_git_commit() -> str:
    """Return the current short git commit or 'unknown' if git is unavailable."""
    try:                                                                                                                # Best-effort only — no project should require git to be installed.
        out = subprocess.run(                                                                                            # Capture git rev-parse output.
            ["git", "rev-parse", "--short", "HEAD"],                                                                    # Standard invocation.
            stdout=subprocess.PIPE,                                                                                     # Capture stdout.
            stderr=subprocess.DEVNULL,                                                                                  # Swallow stderr noise.
            timeout=2,                                                                                                   # Don't hang if git is broken.
        )
        return out.stdout.decode().strip() or "unknown"                                                                 # Clean up and return.
    except Exception:                                                                                                   # Pragma: only fires when git is missing or fails.
        return "unknown"                                                                                                # Sentinel for downstream tooling.


# ----------------------------- Final context assembly ----------------------------- #


def build_stage11_context(
    *,
    seed_row: pd.Series,
    target_host: str,
    source_host: str,
    seed_pdb_path: str | Path,
    baseline_qualification: dict,
    family_df: pd.DataFrame,
    target_df: pd.DataFrame,
    family_centroid: np.ndarray,
    target_centroid: np.ndarray,
    proposals: list[EditProposal],
    hard_positions: list[int],
    soft_positions: list[int],
    min_mutations: int,
    max_mutations: int,
    cli_args: dict,
    embedding_model: str,
    predictor_model_path: str | Path,
    predictor_label_classes_path: str | Path,
    strict_csv_path: str | Path,
    strict_embeddings_path: str | Path | None,
) -> dict:
    """Assemble the complete `stage11_context.json` dictionary.

    The outer shape mirrors `stage10_context.json` so 11d can reuse the
    unmodified 08a validator. The Stage 11 additions are:

    - `baseline_qualification`: gate metrics + threshold + pass flag.
    - `target_context`: top-M target-host RBPs + their centroid.
    - `stage11_provenance`: CLI args, model paths, timestamp, git commit.

    The `editable_region.hotspot_positions` alias is included so 08a's legacy
    field read does not break.
    """
    seed_sequence = str(seed_row["aa_sequence"])                                                                        # Wild-type seed sequence (already validated).
    L = len(seed_sequence)                                                                                              # Seed length for frozen-position math.

    editable_positions = sorted(set(int(p) for p in hard_positions) | set(int(p) for p in soft_positions))              # Union of hard + soft, sorted.
    frozen_positions = sorted(set(range(1, L + 1)) - set(editable_positions))                                           # Everything else is frozen.

    proposal_rows = serialize_edit_proposals(proposals, hard_positions, soft_positions)                                 # JSON-safe edit proposals.

    family_member_ids = family_df["protein_id"].astype(str).tolist() if len(family_df) else []                          # IDs for family table.
    family_member_cosines = [float(x) for x in family_df["cosine_to_seed"].tolist()] if len(family_df) else []          # Cosines for family table.

    family_rows_payload = []                                                                                            # Embed the family rows so downstream stages don't need to re-merge.
    for _, row in family_df.iterrows():                                                                                 # Walk every family member.
        family_rows_payload.append(                                                                                     # Append a JSON-safe dict per row.
            {
                "protein_id": str(row["protein_id"]),                                                                   # Member ID.
                "host_genus": str(row["host_genus"]),                                                                   # Member host (mostly source-host for the family).
                "product": str(row["product"]),                                                                         # Annotation string.
                "virus_accession": str(row["virus_accession"]),                                                         # Phage accession.
                "aa_sequence": str(row["aa_sequence"]),                                                                 # Full sequence (needed for entropy + edit-space rebuilds).
                "cosine_to_seed": float(row["cosine_to_seed"]),                                                         # Distance metric.
                "length_diff_frac": float(row["length_diff_frac"]),                                                     # Length filter outcome.
                "family_rank": int(row["family_rank"]),                                                                 # 1-based rank in the family.
            }
        )

    target_member_ids = target_df["protein_id"].astype(str).tolist() if len(target_df) else []                          # IDs for target table.
    target_member_cosines = [float(x) for x in target_df["cosine_to_seed"].tolist()] if len(target_df) else []          # Cosines for target table.
    target_rows_payload = []                                                                                            # Embed the target rows likewise.
    for _, row in target_df.iterrows():                                                                                 # Walk every target member.
        target_rows_payload.append(                                                                                     # Append a JSON-safe dict per row.
            {
                "protein_id": str(row["protein_id"]),                                                                   # Member ID.
                "host_genus": str(row["host_genus"]),                                                                   # Member host (== target_host).
                "product": str(row["product"]),                                                                         # Annotation.
                "virus_accession": str(row["virus_accession"]),                                                         # Phage accession.
                "aa_sequence": str(row["aa_sequence"]),                                                                 # Full sequence.
                "cosine_to_seed": float(row["cosine_to_seed"]),                                                         # Distance metric.
                "length_diff_frac": float(row["length_diff_frac"]),                                                     # Length filter outcome.
                "family_rank": int(row["family_rank"]),                                                                 # 1-based rank within the target pool.
            }
        )

    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")                         # UTC ISO-8601 stamp.

    context = {                                                                                                         # The master context dict.
        "stage": "11a",                                                                                                 # Origin stage tag.
        "target_host": str(target_host),                                                                                # The host we are flipping toward.
        "source_host": str(source_host),                                                                                # The seed's native host.
        "seed_pdb_path": str(Path(seed_pdb_path).resolve()),                                                            # Absolute path to the freshly folded seed PDB.

        "selected_seed": {                                                                                              # Seed identification block.
            "seed_protein_id": str(seed_row["protein_id"]),                                                             # Strict-CSV protein_id.
            "seed_identifier_hint": str(seed_row["protein_id"]),                                                        # Hint duplicates the ID for legacy readers.
            "virus_accession": str(seed_row["virus_accession"]),                                                        # Phage accession.
            "source_host": str(source_host),                                                                            # Native host.
            "product": str(seed_row["product"]),                                                                        # Annotation.
            "seed_sequence": seed_sequence,                                                                             # Full sequence — required by 08a.
            "sequence_length": int(L),                                                                                  # Convenience for downstream readers.
            "seed_source_kind": "strict_csv_direct",                                                                    # Provenance tag.
            "seed_source_desc": f"{Path(strict_csv_path).name} row protein_id={seed_row['protein_id']}",                # Human-readable provenance string.
        },

        "baseline_qualification": dict(baseline_qualification),                                                         # Gate metrics + pass flag.

        "family_context": {                                                                                             # Family signal block.
            "family_member_count": int(len(family_df)),                                                                 # Pool size.
            "family_member_ids": family_member_ids,                                                                     # IDs only (compact).
            "family_member_cosines": family_member_cosines,                                                             # Cosines only (compact).
            "family_centroid": [float(x) for x in family_centroid.tolist()],                                            # Full centroid vector.
            "family_rows": family_rows_payload,                                                                         # Full rows for downstream stages.
        },

        "target_context": {                                                                                             # Target-host signal block (new in Stage 11).
            "target_host": str(target_host),                                                                            # Target host name.
            "target_member_count": int(len(target_df)),                                                                 # Pool size.
            "target_member_ids": target_member_ids,                                                                     # IDs only.
            "target_member_cosines": target_member_cosines,                                                             # Cosines only.
            "target_centroid": [float(x) for x in target_centroid.tolist()],                                            # Centroid vector.
            "target_rows": target_rows_payload,                                                                         # Full rows.
        },

        "editable_region": {                                                                                            # Edit-space block (08a-compatible).
            "hard_positions": [int(p) for p in hard_positions],                                                         # Primary edit sites.
            "soft_positions": [int(p) for p in soft_positions],                                                         # Secondary edit sites.
            "editable_positions": editable_positions,                                                                   # Union (= 08a-compatible field).
            "hotspot_positions": editable_positions,                                                                    # Legacy alias the 08a validator reads (line 601).
            "frozen_positions": frozen_positions,                                                                       # Everything else (informational).
            "min_mutations": int(min_mutations),                                                                        # Lower mutation budget.
            "max_mutations": int(max_mutations),                                                                        # Upper mutation budget.
            "proposal_rows": proposal_rows,                                                                             # The full from-scratch edit proposals.
        },

        "stage11_provenance": {                                                                                         # Provenance block.
            "strict_csv_path": str(Path(strict_csv_path).resolve()),                                                    # Where the dataset came from.
            "strict_embeddings_path": (str(Path(strict_embeddings_path).resolve())                                      # Where the cached embeddings came from (or None).
                                       if strict_embeddings_path else None),
            "embedding_model": str(embedding_model),                                                                    # The HuggingFace ESM-2 id.
            "predictor_model_path": str(Path(predictor_model_path).resolve()),                                          # Where the host probe came from.
            "predictor_label_classes_path": str(Path(predictor_label_classes_path).resolve()),                          # Where the label mapping came from.
            "cli_args": {str(k): _jsonable(v) for k, v in cli_args.items()},                                            # Every CLI knob the operator passed.
            "git_commit": _try_git_commit(),                                                                            # Best-effort repo state.
            "timestamp_utc": timestamp,                                                                                 # When the context was built.
        },
    }
    return context                                                                                                      # Done.
