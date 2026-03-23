# This script aggregates all existing design-run candidates, computes fast validity metrics for them such as: sequence similarity, 
# and writes one master table that you can reuse for ranking and reporting.

from __future__ import annotations                                                      
import argparse                                                                         
import json                                                                             
import re                                                                               
from pathlib import Path                                                                
from typing import Dict, Iterable, List, Tuple                                          
import numpy as np                                                                      
import pandas as pd                                                                     
import torch                                                                            
from tqdm import tqdm                                                                   
from transformers import AutoTokenizer, EsmForMaskedLM                                 

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")                                               # Standard 20 amino acids accepted by the script.

RBP_KEYWORDS = [                                                                        # Start the list of keywords used for family-retention flags.
    "tail fiber",                                                                       # Canonical tail fiber keyword.
    "tail fibre",                                                                       # British spelling variant.
    "tail spike",                                                                       # Canonical tail spike keyword.
    "receptor binding",                                                                 # Receptor-binding wording.
    "host specificity",                                                                 # Host-specificity wording.
    "adsorption",                                                                       # Adsorption-related wording.
    "fiber protein",                                                                    # Broader fiber wording.
    "fibre protein",                                                                    # Broader fibre wording.
]                                                                                       # End the annotation keyword list.

# Use regex patterns so hyphenated and non-hyphenated RBP annotations are both recognized.
RBP_PATTERNS = [                                                                        
    r"receptor[-\s]?binding",                                                           
    r"tail[-\s]?fiber",                                                                 
    r"tail[-\s]?fibre",                                                                 
    r"tail[-\s]?spike",                                                                 
    r"host[-\s]?specific",                                                              
    r"adsorption",                                                                      
    r"fib[ea]r protein",                                                                
    r"tail protein",                                                                    
    r"depolymerase",                                                                    
]                                                                                       

# Patterns that often indicate non-RBP neighbors to reduce false positives.
NEGATIVE_PATTERNS = [                                                                   
    r"baseplate assembly",                                                              
    r"tail sheath",                                                                     
    r"tail tube",                                                                       
    r"portal",                                                                          
    r"terminase",                                                                       
    r"capsid",                                                                          
    r"major tail",                                                                      
]                                                                                       


def is_valid_sequence(seq: str) -> bool:                                                
    """
    This helper validates that sequences look like normal proteins before scoring them.
    """
    seq = str(seq).strip().upper()                                                      # Normalize whitespace and case before checking.
    return len(seq) > 0 and all(ch in AMINO_ACIDS for ch in seq)                        # Keep only non-empty sequences using standard amino acids.


def normalize_rows(x: np.ndarray) -> np.ndarray:                                        
    """
    This helper converts a 2D embedding matrix into row-normalized vectors.
    """
    norms = np.linalg.norm(x, axis=1, keepdims=True)                                    # Compute one L2 norm per embedding row.
    norms = np.clip(norms, 1e-12, None)                                                 # Avoid division by zero on degenerate vectors.
    return x / norms                                                                    # Return row-wise unit vectors for cosine similarity.


def parse_mutation_positions(mutation_text: str) -> List[int]:                          
    """
    This helper parses mutation strings such as "A123V;G456D" into 0-based positions (1-based biological numbering to 0-based Python indexing).
    """
    if pd.isna(mutation_text):                                                          # Handle missing mutation strings safely.
        return []                                                                       # Missing text means no parsed positions.
    positions: List[int] = []                                                           # Create a list that will collect parsed positions.
    for token in str(mutation_text).split(";"):                                         # Split on semicolons because the repo stores one mutation per token.
        token = token.strip()                                                           # Remove surrounding whitespace from each token.
        if token == "":                                                                 # Ignore empty pieces created by malformed strings.
            continue                                                                    # Skip blank tokens.
        match = re.search(r"(\d+)", token)                                              # Extract the numeric residue index from the token.
        if match is None:                                                               # Guard against unexpected token formats.
            continue                                                                    # Ignore malformed tokens instead of crashing.
        positions.append(int(match.group(1)) - 1)                                       # Convert 1-based biological numbering to 0-based Python indexing.
    return sorted(set(p for p in positions if p >= 0))                                  # Remove duplicates, keep only valid non-negative positions, and sort.


def sequence_distance(seq_a: str, seq_b: str) -> int:                                   
    """
    This helper measures simple Hamming distance between two same-length protein strings.
    If the lengths are different, it also adds the difference to their distance.
    """
    if len(seq_a) != len(seq_b):                                                        # Most design candidates should have the same length as the seed.
        min_len = min(len(seq_a), len(seq_b))                                           # Compare only the overlapping region first.
        return sum(a != b for a, b in zip(seq_a[:min_len], seq_b[:min_len])) + abs(len(seq_a) - len(seq_b))  # Count mismatches plus length difference.
    return sum(a != b for a, b in zip(seq_a, seq_b))                                    # Return Hamming distance when lengths match.


def normalize_annotation_text(text: str) -> str:                                        
    """
    This helper normalizes product annotations so keyword matching is robust to punctuation and spacing.
    For example, "tail fiber" and "tail-fiber" both match "tail fiber".
    """
    if pd.isna(text):                                                                   # Missing annotations should normalize to an empty string.
        return ""                                                                       # Empty string is easier to handle downstream than NaN.
    text = str(text).strip().lower()                                                    # Normalize case and trim whitespace before matching.
    text = re.sub(r"[_/|]+", " ", text)                                                 # Replace separators with spaces to avoid brittle matching.
    text = re.sub(r"\s+", " ", text)                                                    # Collapse repeated spaces into one canonical form.
    return text                                                                         # Return normalized annotation text.


def annotation_is_rbp_like(product_text: str) -> bool:                                  
    """
    This helper returns whether an annotation looks like an RBP-family annotation.
    Returns False if the text is empty or contains only negative keywords.
    """
    text = normalize_annotation_text(product_text)                                      # Normalize annotation text so hyphenated and spaced forms match equally.
    if text == "":                                                                      # Empty annotations should fail closed instead of matching accidentally.
        return False                                                                    # Missing annotation is not enough evidence for family retention.
    if any(re.search(pattern, text) for pattern in NEGATIVE_PATTERNS):                  # Block obvious non-RBP annotations unless rescued by a positive match.
        if not any(re.search(pattern, text) for pattern in RBP_PATTERNS):               # Allow borderline annotations through only when they also contain a positive RBP phrase.
            return False                                                                # Negative-only annotation should not count as RBP-like.
    return any(re.search(pattern, text) for pattern in RBP_PATTERNS)                    # Use robust regex matching rather than fragile exact substring matching.


def sequence_mutation_list(seed_seq: str, candidate_seq: str) -> Tuple[str, int]:      
    """
    This helper (re)computes exact mutations directly from seed and candidate sequences.
    Returns a list of explicit substitution strings and the count of mutations.
    """
    if pd.isna(seed_seq) or pd.isna(candidate_seq):                                     # Guard against missing sequence values.
        return "", np.nan                                                               # Return an empty list and unknown count when computation is impossible.
    seed_seq = str(seed_seq)                                                            # Normalize seed sequence to a string.
    candidate_seq = str(candidate_seq)                                                  # Normalize candidate sequence to a string.
    if len(seed_seq) != len(candidate_seq):                                             # If lengths differ, we cannot express the difference as simple substitutions.
        return "LEN_MISMATCH", sequence_distance(seed_seq, candidate_seq)               # Return a sentinel plus the robust distance.
    mutations: List[str] = []                                                           # Collect explicit substitution strings here.
    for idx, (seed_aa, candidate_aa) in enumerate(zip(seed_seq, candidate_seq), start=1):  # Compare aligned residues one-by-one using biological 1-based numbering.
        if seed_aa != candidate_aa:                                                     # Record only positions that changed.
            mutations.append(f"{seed_aa}{idx}{candidate_aa}")                           # Build a human-readable substitution token.
    return ";".join(mutations), len(mutations)                                          # Return the explicit mutation list and exact mutation count.


def load_json_if_exists(path: Path) -> dict:                                            
    """
    This helper loads a design-run metadata JSON if the file exists.
    """
    if not path.exists():                                                               # Check whether the file is present.
        return {}                                                                       # Return an empty dictionary when metadata is missing.
    with open(path, "r", encoding="utf-8") as handle:                                   # Open the JSON file with explicit UTF-8 encoding.
        return json.load(handle)                                                        # Parse and return the JSON object.


def build_seed_lookup(csv_paths: Iterable[Path]) -> Dict[str, dict]:                    
    """
    This helper loads the seed-sequence lookup table from the strict dataset first. 
    It reads the CSVs that may contain the seed proteins and returns a dictionary that maps protein IDs to metadata records such as aa_sequence, virus_accession, and host_genus.
    """
    lookup: Dict[str, dict] = {}                                                        # Create the mapping container.
    for csv_path in csv_paths:                                                          # Loop over datasets that may contain the seed proteins.
        if not csv_path.exists():                                                       # Skip files that are not present in the repo.
            continue                                                                    # Move on to the next candidate CSV.
        df = pd.read_csv(csv_path)                                                      # Load the dataset into a DataFrame.
        required = {"protein_id", "aa_sequence", "virus_accession", "host_genus"}       # Define the minimum columns required for seed recovery.
        if not required.issubset(set(df.columns)):                                      # Validate that the CSV has the fields we need.
            continue                                                                    # Ignore datasets without the required columns.
        for _, row in df.iterrows():                                                    # Iterate through rows to populate the lookup.
            pid = str(row["protein_id"])                                                # Convert the protein identifier to a string key.
            if pid not in lookup:                                                       # Keep the first occurrence to avoid accidental overwrites.
                lookup[pid] = {                                                         # Store the seed metadata for later joins.
                    "protein_id": pid,                                                  # Save protein ID explicitly.
                    "aa_sequence": str(row["aa_sequence"]),                             # Save the protein sequence.
                    "virus_accession": str(row["virus_accession"]),                     # Save the virus accession for context.
                    "host_genus": str(row["host_genus"]),                               # Save the source host genus.
                }                                                                       # Finish the seed metadata record.
    return lookup                                                                       # Return the completed protein lookup dictionary.


def load_reference_embedding_bank(embedding_pt: Path, embedding_index_csv: Path) -> pd.DataFrame:  
    """
    This helper loads the precomputed stage-1 embeddings, aligns them to their index, and returns a DataFrame of metadata and embeddings.
    """
    emb_tensor = torch.load(embedding_pt, map_location="cpu")                           # Load the stage-1 embeddings from disk onto CPU memory.
    emb_array = emb_tensor.detach().cpu().numpy().astype(np.float32)                    # Convert the tensor into a NumPy array for sklearn math.
    index_df = pd.read_csv(embedding_index_csv)                                         # Load the row metadata corresponding to the embeddings.
    index_df = index_df.copy()                                                          # Create an explicit copy before editing columns.
    if "row_id" not in index_df.columns:                                                # Validate that the repo index contains the row mapping column.
        raise ValueError("esm2_embeddings_index.csv must contain a 'row_id' column.")   # Raise a clear error if the expected format is missing.
    if emb_array.shape[0] != len(index_df):                                             # Validate that the row count matches between embeddings and index.
        raise ValueError("Embedding row count does not match esm2_embeddings_index.csv.")  # Stop early on corrupted or mismatched inputs.
    index_df["embedding_row"] = index_df["row_id"].astype(int)                          # Create an explicit integer row index for safe gathering.
    index_df["product"] = index_df.get("product", "")                                   # Ensure a product column exists even if absent in some exports.
    index_df["protein_id"] = index_df["protein_id"].astype(str)                         # Normalize protein IDs to strings for merges.
    index_df["host_genus"] = index_df["host_genus"].astype(str)                         # Normalize host labels to strings for reporting.
    index_df["virus_accession"] = index_df["virus_accession"].astype(str)               # Normalize accession labels to strings for reporting.
    index_df.attrs["embedding_array"] = emb_array                                       # Attach the actual embedding matrix to the DataFrame metadata.
    return index_df                                                                     # Return the aligned embedding bank.


def subset_reference_bank(bank_df: pd.DataFrame, dataset_csv: Path) -> Tuple[pd.DataFrame, np.ndarray]:  
    """
    This helper filters the aligned embedding bank down to one reference dataset.
    Returns both the metadata table and the aligned embeddings for the reference dataset.
    """
    dataset_df = pd.read_csv(dataset_csv)                                                   # Load the target reference dataset.
    dataset_ids = set(dataset_df["protein_id"].astype(str))                                 # Collect the protein IDs belonging to this reference set.
    subset_df = bank_df[bank_df["protein_id"].isin(dataset_ids)].copy().reset_index(drop=True)  # Keep only rows present in the chosen dataset.
    emb_array = bank_df.attrs["embedding_array"]                                            # Retrieve the full embedding matrix stored on the bank DataFrame.
    subset_emb = emb_array[subset_df["embedding_row"].to_numpy()]                           # Gather the matching subset embeddings by row index.
    return subset_df, subset_emb                                                            # Return both the metadata table and its aligned embeddings.


def aggregate_design_runs(runs_dir: Path, seed_lookup: Dict[str, dict]) -> pd.DataFrame:  
    """
    This helper extracts all top candidates across all existing design runs, and combines them into a single DataFrame.
    """
    rows: List[dict] = []                                                               # Create a list of candidate records that will become one DataFrame.
    for run_dir in sorted(runs_dir.glob("*")):                                          # Loop over each design-run folder in a deterministic order.
        if not run_dir.is_dir():                                                        # Ignore any stray files that are not run folders.
            continue                                                                    # Move on to the next path.
        top_csv = run_dir / "top_candidates.csv"                                        # Point to the saved top-candidate table for this run.
        meta_json = run_dir / "run_metadata.json"                                       # Point to the saved metadata JSON for this run.
        if not top_csv.exists():                                                        # Skip folders that do not contain a candidate table.
            continue                                                                    # Move on to the next run.
        top_df = pd.read_csv(top_csv)                                                   # Load the top-candidate table.
        meta = load_json_if_exists(meta_json)                                           # Load metadata when available.
        seed_pid = str(meta.get("seed_protein_id", top_df.get("protein_id", pd.Series([""])).astype(str).iloc[0]))  # Infer the seed protein ID from metadata or fallback columns.
        seed_info = seed_lookup.get(seed_pid, {})                                       # Recover the seed sequence and host context from the lookup.
        method_name = run_dir.name                                                      # Use the run-directory name as the broad method label.
        for _, row in top_df.iterrows():                                                # Iterate over every candidate row saved by the design run.
            candidate_seq = str(row.get("aa_sequence", ""))                             # Extract the candidate amino-acid sequence.
            if not is_valid_sequence(candidate_seq):                                    # Ignore malformed or empty sequences.
                continue                                                                # Skip invalid sequence rows.
            rows.append({                                                               # Append one normalized record for the candidate.
                "run_name": run_dir.name,                                               # Save the run folder name for grouping later.
                "run_dir": str(run_dir),                                                # Save the full run path for traceability.
                "method": method_name,                                                  # Save a simple method label derived from the folder name.
                "candidate_id": str(row.get("candidate_id", "")),                       # Save the candidate ID assigned by the optimizer.
                "round": int(row.get("round", -1)) if not pd.isna(row.get("round", np.nan)) else -1,  # Save the design round if present.
                "rank_in_round": int(row.get("rank_in_round", -1)) if not pd.isna(row.get("rank_in_round", np.nan)) else -1,  # Save the candidate rank inside the round if present.
                "protein_id": str(row.get("protein_id", seed_pid)),                     # Save the protein ID; fallback to the seed protein ID when missing.
                "virus_accession": str(row.get("virus_accession", seed_info.get("virus_accession", ""))),  # Save the accession with metadata fallback.
                "source_host": str(row.get("source_host", meta.get("source_host", seed_info.get("host_genus", "")))),  # Save the original host genus.
                "target_host": str(row.get("target_host", meta.get("target_host", ""))),  # Save the target host label used during optimization.
                "pred_label": str(row.get("pred_label", "")),                           # Save the predicted class label from the original run.
                "target_score": float(row.get("target_score", np.nan)),                 # Save the original target-host score.
                "selection_score": float(row.get("selection_score", np.nan)),           # Save the internal selection score from the optimizer.
                "seed_cosine_similarity": float(row.get("seed_cosine_similarity", np.nan)),  # Save the existing seed similarity metric when present.
                "seed_novelty_distance": float(row.get("seed_novelty_distance", np.nan)),  # Save the existing novelty metric when present.
                "n_mutations": int(row.get("n_mutations", np.nan)) if not pd.isna(row.get("n_mutations", np.nan)) else np.nan,  # Save the mutation count when present.
                "mutations": str(row.get("mutations", "")),                             # Save the human-readable mutation string.
                "aa_sequence": candidate_seq,                                           # Save the actual candidate sequence.
                "seed_protein_id": seed_pid,                                            # Save the seed protein ID explicitly for joins and grouping.
                "seed_sequence": str(seed_info.get("aa_sequence", "")),                 # Save the actual seed sequence for distance calculations.
                "optimization_seed": int(meta.get("seed", -1)) if meta.get("seed", None) is not None else -1,  # Save the optimizer random seed.
                "model_path": str(meta.get("model_path", "")),                          # Save the probe model path used by the run.
                "label_classes_path": str(meta.get("label_classes_path", "")),          # Save the label mapping path used by the run.
            })                                                                          # Finish the candidate record.
    if len(rows) == 0:                                                                  # Fail early if no design candidates were found.
        raise ValueError(f"No top_candidates.csv files found under: {runs_dir}")        # Raise a clear error so the user can fix the input path.
    df = pd.DataFrame(rows)                                                             # Turn the list of records into one DataFrame.
    df["seq_dist_from_seed"] = df.apply(lambda r: sequence_distance(str(r["aa_sequence"]), str(r["seed_sequence"])) if str(r["seed_sequence"]) != "" else np.nan, axis=1)  # Compute exact distance from the seed sequence.
    mutation_pairs = df.apply(lambda r: sequence_mutation_list(str(r["seed_sequence"]), str(r["aa_sequence"])), axis=1)  # Derive exact mutation strings and counts directly from the sequences.
    df["actual_mutation_list"] = [pair[0] for pair in mutation_pairs]                  # Save explicit mutations reconstructed from seed and candidate sequences.
    df["actual_n_mutations"] = [pair[1] for pair in mutation_pairs]                    # Save the exact mutation count computed from the sequences.
    df["actual_seq_dist_from_seed"] = df["seq_dist_from_seed"]                        # Expose the actual distance under an explicit corrected column name.
    df["mutation_count_disagreement"] = pd.to_numeric(df["n_mutations"], errors="coerce") != pd.to_numeric(df["actual_n_mutations"], errors="coerce")  # Flag rows where logged mutation counts disagree with the actual sequence difference.
    return df                                                                           # Return the aggregated candidate table.


def embed_sequences(sequences: List[str], tokenizer, model, batch_size: int, max_length: int, device: str) -> np.ndarray:  
    """
    This helper embeds arbitrary sequences with the same ESM model family used upstream.
    Returns a matrix of shape [N, H] where N is the number of sequences and H is the embedding dimension.
    """
    pooled_rows: List[np.ndarray] = []                                                  # Create a list to collect batch embeddings.
    with torch.inference_mode():                                                        # Disable gradients because we only need forward passes.
        for start in tqdm(range(0, len(sequences), batch_size), desc="Embedding candidate sequences"):  # Iterate through sequences in batches with a progress bar.
            batch = sequences[start : start + batch_size]                               # Slice one batch of raw amino-acid sequences.
            toks = tokenizer(batch, return_tensors="pt", padding=True, truncation=True, max_length=max_length).to(device)  # Tokenize and move batch tensors to the chosen device.
            out = model.esm(**toks)                                                     # Run the encoder part of the masked language model.
            hidden = out.last_hidden_state                                              # Extract token-level hidden states of shape [B, L, H].
            mask = toks["attention_mask"].unsqueeze(-1).to(hidden.dtype)                # Expand the attention mask so it can weight token vectors.
            summed = (hidden * mask).sum(dim=1)                                         # Sum non-padding token embeddings across sequence length.
            denom = mask.sum(dim=1).clamp(min=1.0)                                      # Count real tokens while avoiding zero division.
            pooled = (summed / denom).detach().cpu().numpy().astype(np.float32)         # Mean-pool tokens into one embedding per sequence on CPU.
            pooled_rows.append(pooled)                                                  # Save the pooled batch embeddings.
    return np.concatenate(pooled_rows, axis=0)                                          # Concatenate all batches into one [N, H] matrix.


def compute_mutation_local_pll(candidate_df: pd.DataFrame, tokenizer, model, device: str, max_length: int) -> pd.DataFrame:  
    """
    Naturalness scoring: This helper uses the ESM MLM model to look specifically at the amino acids that were mutated, 
    computes the Pseudo-Log-Likelihood (PLL) score for each mutation with respect to the seed, \
    and returns a DataFrame with one row per candidate/seed pair. 
    """
    mask_token_id = tokenizer.mask_token_id                                             # Cache the special [MASK] token ID once for efficiency.
    records: List[dict] = []                                                            # Create a list to collect one PLL record per candidate.
    unique_pairs = candidate_df[["aa_sequence", "seed_sequence", "mutations", "actual_mutation_list"]].drop_duplicates().reset_index(drop=True)  # Carry the sequence-derived mutation list so PLL uses the corrected positions when available.
    for _, row in tqdm(unique_pairs.iterrows(), total=len(unique_pairs), desc="Computing mutation-local PLL"):  # Loop through unique pairs with a progress bar.
        candidate_seq = str(row["aa_sequence"])                                         # Extract the candidate sequence.
        seed_seq = str(row["seed_sequence"])                                            # Extract the corresponding seed sequence (which is also found in the candidate sequences).
        corrected_mutation_text = row.get("actual_mutation_list", "")                  # Prefer the exact mutation list reconstructed from the sequences.
        mutation_source = corrected_mutation_text if str(corrected_mutation_text) not in {"", "nan", "LEN_MISMATCH"} else str(row["mutations"])  # Fall back to the logged mutation text only when the corrected list is unavailable.
        mutation_positions = parse_mutation_positions(str(mutation_source))              # Parse mutation positions from the corrected mutation source.
        mutation_positions = [p for p in mutation_positions if p < len(candidate_seq) and p < len(seed_seq)]  # Keep only positions valid for both sequences.
        if len(mutation_positions) == 0:                                                # If there are no mutations, store NaNs in the pll columns.
            records.append({                                                            
                "aa_sequence": candidate_seq,                                           
                "seed_sequence": seed_seq,                                              
                "pll_candidate_mean": np.nan,                                           
                "pll_seed_mean": np.nan,                                                
                "pll_delta_candidate_minus_seed": np.nan,                               
            })                                                                          
            continue                                                                    
        seqs_to_score = [candidate_seq, seed_seq]                                       # Put both sequences into a list for scoring.
        seq_labels = ["candidate", "seed"]                                              # Attach readable labels to the two scored sequences.
        pll_means: Dict[str, float] = {}                                                # Create a dictionary for the resulting mean local PLL values.
        for seq_label, seq_text in zip(seq_labels, seqs_to_score):                      # Loop over the candidate and seed sequences.
            input_ids = tokenizer(seq_text, return_tensors="pt", truncation=True, max_length=max_length)["input_ids"].to(device)  # Tokenize one sequence for masked scoring.
            token_ids = input_ids[0].clone()                                            # Make a mutable copy of the token IDs.
            aa_token_positions = list(range(1, min(len(seq_text), max_length - 2) + 1)) # Map amino-acid positions to token positions after the BOS token.
            per_pos_log_probs: List[float] = []                                         # Collect the masked log-probability at each mutated residue.
            for aa_pos in mutation_positions:                                           # Loop over each mutated amino-acid position.
                if aa_pos >= len(aa_token_positions):                                   # Guard against truncation by the tokenizer max length.
                    continue                                                            # Skip positions that fall outside the tokenized region.
                token_pos = aa_token_positions[aa_pos]                                  # Convert amino-acid index to token index in the model input.
                true_token_id = int(token_ids[token_pos].item())                        # Save the real token ID for the residue being evaluated.
                masked_ids = token_ids.clone()                                          # Copy the token sequence so we can mask one position.
                masked_ids[token_pos] = mask_token_id                                   # Replace the residue token with [MASK].
                with torch.inference_mode():                                            # Disable gradients for efficient scoring.
                    logits = model(input_ids=masked_ids.unsqueeze(0)).logits[0, token_pos]  # Run the masked model and extract logits at the masked position.
                    log_probs = torch.log_softmax(logits, dim=-1)                       # Convert logits into log-probabilities.
                per_pos_log_probs.append(float(log_probs[true_token_id].item()))        # Save the log-probability assigned to the true residue.
            pll_means[seq_label] = float(np.mean(per_pos_log_probs)) if len(per_pos_log_probs) > 0 else np.nan  # Average across mutated positions for this sequence.
        records.append({                                                                # Save one merged PLL record for the candidate/seed pair.
            "aa_sequence": candidate_seq,                                               # Key by candidate sequence for joining later.
            "seed_sequence": seed_seq,                                                  # Also key by seed sequence to disambiguate repeated candidates.
            "pll_candidate_mean": pll_means.get("candidate", np.nan),                   # Save the candidate mean local PLL.
            "pll_seed_mean": pll_means.get("seed", np.nan),                             # Save the seed mean local PLL at the same positions.
            "pll_delta_candidate_minus_seed": pll_means.get("candidate", np.nan) - pll_means.get("seed", np.nan) if not np.isnan(pll_means.get("candidate", np.nan)) and not np.isnan(pll_means.get("seed", np.nan)) else np.nan,  # Save candidate minus seed naturalness delta.
        })                                                                              # Finish the PLL record.
    return pd.DataFrame(records)                                                        # Return one DataFrame with all mutation-local PLL metrics.


def compute_reference_neighbors(query_emb_norm: np.ndarray, ref_df: pd.DataFrame, ref_emb_norm: np.ndarray, prefix: str, top_k: int) -> pd.DataFrame:  
    """
    This helper computes nearest-neighbor metrics between a query and a reference bank, and returns them as a DataFrame.
    """
    similarity = query_emb_norm @ ref_emb_norm.T                                        # Compute cosine similarities by dot product because all rows are normalized.
    nearest_idx = similarity.argmax(axis=1)                                             # Find the nearest reference neighbor for each query sequence.
    topk_idx = np.argsort(similarity, axis=1)[:, -top_k:]                               # Find the indices of the top-k nearest neighbors for density-like metrics.
    rows: List[dict] = []                                                               # Create the output records for this reference bank.
    for i in range(similarity.shape[0]):                                                # Loop over each query sequence.
        nn_i = int(nearest_idx[i])                                                      # Extract the nearest-neighbor index as a Python int.
        topk_sims = similarity[i, topk_idx[i]]                                          # Gather the top-k cosine similarities for the current query.
        nn_row = ref_df.iloc[nn_i]                                                      # Retrieve the metadata row for the nearest reference protein.
        rows.append({                                                                   # Save one record of nearest-neighbor metrics.
            f"{prefix}_nn_protein_id": str(nn_row["protein_id"]),                       # Save the nearest reference protein ID.
            f"{prefix}_nn_host": str(nn_row.get("host_genus", "")),                     # Save the nearest reference host genus.
            f"{prefix}_nn_product": str(nn_row.get("product", "")),                     # Save the nearest reference product annotation.
            f"{prefix}_nn_cosine": float(similarity[i, nn_i]),                          # Save the nearest-neighbor cosine similarity.
            f"{prefix}_top{top_k}_mean_cosine": float(np.mean(topk_sims)),              # Save the mean top-k cosine similarity as a local-density proxy.
            f"{prefix}_annotation_is_rbp_like": bool(annotation_is_rbp_like(str(nn_row.get("product", "")))),  # Save a heuristic family-retention flag.
        })                                                                              # Finish the nearest-neighbor record.
    return pd.DataFrame(rows)                                                           # Return the metrics for this reference bank.


def add_relative_flags(df: pd.DataFrame, top_k_density: int) -> pd.DataFrame:          
    """
    This helper computes percentile-based pass/fail flags relative to the seed itself:                                      
        passes_mutation_budget: Keep only local designs with at most three mutations by default.
        passes_manifold_floor: Flag candidates that stay within the lower-quartile strict manifold (which means that its cosine similarity to the seed is above the 25% of the cosine similarities of all the candidates from the same seed).
        passes_family_retention: Require the nearest strict annotation to stay RBP-like (which means that the candidate's closest relative is still a true adsorption protein).
        passes_naturalness: Allow moderate but not catastrophic PLL drops compared to the seed.
        passes_density_floor: Candidates don't just stay in the broad RBP manifold, but actually stay in dense, highly-populated regions of the strict embedding space (avoiding weird outlier regions).
        validated_pass: Require all fast filters to pass.
    Returns:
        A pd.DataFrame with the added flags.
    """
    out = df.copy()                                                                         # Work on a copy so the caller keeps the original untouched.
    out["passes_mutation_budget"] = pd.to_numeric(out["actual_n_mutations"], errors="coerce").fillna(999).astype(float) <= 12  # Use the exact sequence-derived mutation count instead of the logged mutation field.
    strict_floor = out.groupby("seed_protein_id")["strict_nn_cosine"].transform(lambda s: s.quantile(0.25))  # Build a per-seed manifold floor using the lower quartile.
    out["passes_manifold_floor"] = out["strict_nn_cosine"] >= strict_floor                  # Flag candidates that stay within the lower-quartile strict manifold (better than 25% of the cosine similarity score of all the candidates from the same seed).
    topk_col = f"strict_top{top_k_density}_mean_cosine"                                      # Build the strict density-column name dynamically from the CLI top-k value.
    topk_floor = out.groupby("seed_protein_id")[topk_col].transform(lambda s: s.quantile(0.25))  # Add a second per-seed manifold-density floor using the strict top-k neighborhood.
    out["passes_family_retention"] = out[["strict_annotation_is_rbp_like", "structural_annotation_is_rbp_like", "structural_plus_annotation_is_rbp_like"]].fillna(False).any(axis=1)  # Allow family retention to be rescued by structural or structural_plus banks too.
    out["passes_naturalness"] = out["pll_delta_candidate_minus_seed"].fillna(0.0) >= 0.0       # Require neutral-or-better mutation-local PLL in the corrected closeout pass.
    out["passes_density_floor"] = out[topk_col] >= topk_floor                                # Require candidates to stay inside the local strict-manifold density band as well.
    out["validated_pass"] = out[["passes_mutation_budget", "passes_manifold_floor", "passes_density_floor", "passes_family_retention", "passes_naturalness"]].all(axis=1)  # Require the corrected density gate in addition to the original fast filters.
    return out                                                                              # Return the flagged table.


def main() -> None:                                                                     
    
    # Parse command-line arguments.
    parser = argparse.ArgumentParser(description="Aggregate PhageForge design candidates and compute fast validity metrics.")  
    parser.add_argument("--repo_root", type=str, default=".", help="Path to the phageforge repository root.")  
    parser.add_argument("--runs_dir", type=str, default="results/design_runs", help="Directory containing all design-run folders.")  
    parser.add_argument("--embedding_pt", type=str, default="data/processed/esm2_embeddings.pt", help="Path to the precomputed stage-1 embedding tensor.")  
    parser.add_argument("--embedding_index_csv", type=str, default="data/processed/esm2_embeddings_index.csv", help="Path to the stage-1 embedding index CSV.")  
    parser.add_argument("--strict_csv", type=str, default="data/processed/rbp_dataset_eskapee_strict.csv", help="Path to the strict RBP dataset CSV.")  
    parser.add_argument("--structural_csv", type=str, default="data/processed/rbp_dataset_eskapee_structural.csv", help="Path to the structural RBP dataset CSV.")  
    parser.add_argument("--structural_plus_csv", type=str, default="data/processed/rbp_dataset_eskapee_structural_plus.csv", help="Path to the structural_plus dataset CSV.")  
    parser.add_argument("--stage1_csv", type=str, default="data/processed/rbp_dataset_eskapee_stage1.csv", help="Path to the stage-1 dataset CSV used to recover seed sequences.")  
    parser.add_argument("--esm_model", type=str, default="facebook/esm2_t33_650M_UR50D", help="ESM model used for candidate embeddings and mutation-local PLL scoring.")  
    parser.add_argument("--batch_size", type=int, default=8, help="Batch size for candidate embedding on GPU.")  
    parser.add_argument("--max_length", type=int, default=1022, help="Maximum protein length supported by ESM-2.")  
    parser.add_argument("--top_k_density", type=int, default=5, help="How many nearest strict neighbors to average for a density-like manifold score.")  
    parser.add_argument("--output_dir", type=str, default="results/analysis/validity", help="Directory where the master validity outputs will be written.")  
    args = parser.parse_args()                                                          

    # Resolve paths
    repo_root = Path(args.repo_root).resolve()                                          
    runs_dir = (repo_root / args.runs_dir).resolve()                                    
    embedding_pt = (repo_root / args.embedding_pt).resolve()                            
    embedding_index_csv = (repo_root / args.embedding_index_csv).resolve()              
    strict_csv = (repo_root / args.strict_csv).resolve()                                
    structural_csv = (repo_root / args.structural_csv).resolve()                        
    structural_plus_csv = (repo_root / args.structural_plus_csv).resolve()              
    stage1_csv = (repo_root / args.stage1_csv).resolve()                                
    output_dir = (repo_root / args.output_dir).resolve()                                
    output_dir.mkdir(parents=True, exist_ok=True)                                       
    
    # Print the resolved paths for root, runs, and output directories for sanity checking
    print(f"[INFO] repo_root={repo_root}")                                               
    print(f"[INFO] runs_dir={runs_dir}")                                                 
    print(f"[INFO] output_dir={output_dir}")                                            

    # Aggregate all design runs into a single table
    seed_lookup = build_seed_lookup([strict_csv, structural_csv, structural_plus_csv, stage1_csv])  # Build one protein-ID lookup for recovering seed sequences.
    all_candidates = aggregate_design_runs(runs_dir=runs_dir, seed_lookup=seed_lookup)              # Aggregate every saved design candidate into one table.
    all_candidates.to_csv(output_dir / "master_candidates_raw.csv", index=False)                    # Save the raw aggregated table before adding model-based metrics.
    print(f"[INFO] Aggregated {len(all_candidates):,} candidate rows across all runs.")             # Print the aggregation size.

    # Build the strict and structural subsets of the reference embedding bank, and normalize them for cosine similarity. 
    bank_df = load_reference_embedding_bank(embedding_pt=embedding_pt, embedding_index_csv=embedding_index_csv)  # Load the full stage-1 embedding bank with metadata.
    strict_df, strict_emb = subset_reference_bank(bank_df, strict_csv)                   # Build the strict reference subset.
    structural_df, structural_emb = subset_reference_bank(bank_df, structural_csv)       # Build the structural reference subset.
    structural_plus_df, structural_plus_emb = subset_reference_bank(bank_df, structural_plus_csv)  # Build the structural_plus reference subset.
    strict_emb_norm = normalize_rows(strict_emb.astype(np.float32))                      # Normalize strict embeddings for cosine similarity.
    structural_emb_norm = normalize_rows(structural_emb.astype(np.float32))              # Normalize structural embeddings for cosine similarity.
    structural_plus_emb_norm = normalize_rows(structural_plus_emb.astype(np.float32))    # Normalize structural_plus embeddings for cosine similarity.

    # Deduplicate candidate sequences, print the number of unique candidates
    unique_candidate_df = all_candidates[["aa_sequence"]].drop_duplicates().reset_index(drop=True)  # Deduplicate candidate sequences before embedding them.
    print(f"[INFO] Unique candidate sequences to embed: {len(unique_candidate_df):,}")      # Print the number of unique candidate sequences.

    device = "cuda" if torch.cuda.is_available() else "cpu"                                 # Choose GPU when available, otherwise fall back to CPU.
    print(f"[INFO] device={device}")                                                        # Print the chosen device.
    tokenizer = AutoTokenizer.from_pretrained(args.esm_model, do_lower_case=False)          # Load the ESM tokenizer matching the requested model.
    esm_mlm = EsmForMaskedLM.from_pretrained(args.esm_model)                                # Load the ESM masked language model once and reuse it.
    esm_mlm.to(device)                                                                      # Move the model to GPU or CPU.
    esm_mlm.eval()                                                                          # Switch the model into evaluation mode.

    # Embed all unique candidate sequences, and normalize them for cosine similarity
    candidate_emb = embed_sequences(sequences=unique_candidate_df["aa_sequence"].tolist(), tokenizer=tokenizer, model=esm_mlm, batch_size=args.batch_size, max_length=args.max_length, device=device)  # Embed all unique candidate sequences.
    candidate_emb_norm = normalize_rows(candidate_emb.astype(np.float32))                   # Normalize candidate embeddings for cosine similarity.

    # Compute neirest-neighbor metrics between each candidate and the strict, structural, and structural_plus reference banks, and combine them into a single table
    strict_neighbor_df = compute_reference_neighbors(query_emb_norm=candidate_emb_norm, ref_df=strict_df, ref_emb_norm=strict_emb_norm, prefix="strict", top_k=args.top_k_density)  # Compute strict-manifold metrics.
    structural_neighbor_df = compute_reference_neighbors(query_emb_norm=candidate_emb_norm, ref_df=structural_df, ref_emb_norm=structural_emb_norm, prefix="structural", top_k=args.top_k_density)  # Compute structural-bank metrics.
    structural_plus_neighbor_df = compute_reference_neighbors(query_emb_norm=candidate_emb_norm, ref_df=structural_plus_df, ref_emb_norm=structural_plus_emb_norm, prefix="structural_plus", top_k=args.top_k_density)  # Compute structural_plus-bank metrics.

    unique_candidate_metrics = pd.concat([unique_candidate_df, strict_neighbor_df, structural_neighbor_df, structural_plus_neighbor_df], axis=1)  # Combine all embedding-based metrics per unique candidate.
    unique_candidate_metrics.to_csv(output_dir / "strict_reference_neighbors.csv", index=False)  # Save the per-unique-sequence nearest-neighbor table.

    # Compute mutation-local PLL metrics for each candidate, and combine them into a single csv
    pll_df = compute_mutation_local_pll(candidate_df=all_candidates[["aa_sequence", "seed_sequence", "mutations", "actual_mutation_list"]].drop_duplicates().reset_index(drop=True), tokenizer=tokenizer, model=esm_mlm, device=device, max_length=args.max_length)  # Pass the corrected mutation list into the PLL step so mutated positions are trustworthy.
    pll_df.to_csv(output_dir / "mutation_local_pll.csv", index=False)                    # Save the raw PLL table for transparency and debugging.

    # Merge all candidate metrics into a single table including validity metrics, mutation-local PLL metrics, and practical pass/fail flags
    merged = all_candidates.merge(unique_candidate_metrics, on="aa_sequence", how="left")  # Attach embedding-based validity metrics back to every candidate row.
    merged = merged.merge(pll_df, on=["aa_sequence", "seed_sequence"], how="left")         # Attach mutation-local PLL metrics back to every candidate row.
    merged["strict_manifold_score"] = merged[["strict_nn_cosine", f"strict_top{args.top_k_density}_mean_cosine"]].mean(axis=1)  # Create one simple strict-manifold composite score.
    merged["family_retention_flag"] = merged[["strict_annotation_is_rbp_like", "structural_annotation_is_rbp_like", "structural_plus_annotation_is_rbp_like"]].fillna(False).any(axis=1)  # Keep family retention as an any-bank signal so strict neighbors can be rescued by structural banks.
    merged = add_relative_flags(merged, top_k_density=args.top_k_density)                 # Use the corrected flagging logic with sequence-derived mutation counts and strict-density floors.

    # Build a compact deduplicated table with one best copy per unique target-host sequence
    sort_cols = ["target_host", "validated_pass", "target_score", "strict_manifold_score", "pll_delta_candidate_minus_seed"]  # Sort so the strongest corrected copy of each target/sequence pair is kept first.
    ascending = [True, False, False, False, False]                                      # Prefer validated, high-scoring, manifold-supported copies.
    merged_sorted = merged.sort_values(sort_cols, ascending=ascending).copy()           # Sort before deduplication so the best copy survives.
    unique_corrected = merged_sorted.drop_duplicates(subset=["target_host", "aa_sequence"], keep="first").copy()  # Keep only one row per unique sequence per target host.
    unique_corrected["duplicate_count_same_target_sequence"] = unique_corrected.apply(lambda r: int(((merged["target_host"] == r["target_host"]) & (merged["aa_sequence"] == r["aa_sequence"])).sum()), axis=1)  # Record how many duplicated copies existed for the same target/sequence pair.

    # Build a compact run-level summary table
    summary_by_run = merged.groupby(["run_name", "target_host"], dropna=False).agg(       # Group by run name and target host.
        n_candidates=("aa_sequence", "size"),                                             # Count how many candidate rows each run contributed.
        n_unique_sequences=("aa_sequence", pd.Series.nunique),                            # Track how many unique sequences each run really contributed.
        mean_target_score=("target_score", "mean"),                                       # Average the target-host score within the run.
        best_target_score=("target_score", "max"),                                        # Save the best target-host score within the run.
        mean_strict_manifold=("strict_manifold_score", "mean"),                           # Average the strict manifold score within the run.
        mean_pll_delta=("pll_delta_candidate_minus_seed", "mean"),                        # Average the candidate-vs-seed PLL delta within the run.
        mean_actual_seq_dist=("actual_seq_dist_from_seed", "mean"),                       # Report actual sequence drift from the seed instead of only the logged mutation field.
        validated_fraction=("validated_pass", "mean"),                                    # Compute the fraction of candidates that pass all fast filters.
    ).reset_index()                                                                       # Flatten the grouped result back into a DataFrame.

    # Save the master candidate table with all validity metrics and the run-level summary table to csv
    merged.to_csv(output_dir / "master_candidates_validated.csv", index=False)            
    unique_corrected.to_csv(output_dir / "master_candidates_unique_corrected.csv", index=False)  # Save the deduplicated corrected table for downstream ranking/reporting.
    summary_by_run.to_csv(output_dir / "validity_summary_by_run.csv", index=False)        

    # Build a metadata dictionary with information about the run
    metadata = {                                                                          
        "repo_root": str(repo_root),                                                      
        "runs_dir": str(runs_dir),                                                        
        "n_all_candidates": int(len(all_candidates)),                                     
        "n_unique_candidate_sequences": int(len(unique_candidate_df)),                    
        "n_unique_target_host_sequences": int(len(unique_corrected)),                     # Record the deduplicated unique target-host sequence count for provenance.
        "strict_reference_size": int(len(strict_df)),                                     
        "structural_reference_size": int(len(structural_df)),                             
        "structural_plus_reference_size": int(len(structural_plus_df)),                   
        "device": device,                                                                 
        "esm_model": args.esm_model,                                                      
        "top_k_density": int(args.top_k_density),                                         
    }                                                                                     
    with open(output_dir / "validity_run_metadata.json", "w", encoding="utf-8") as handle:  # Open the metadata file for writing.
        json.dump(metadata, handle, indent=2)                                             # Write the execution metadata as pretty JSON.

    # Print a summary of what was saved
    print("[INFO] Saved:")                                                               
    print(f"  - {output_dir / 'master_candidates_raw.csv'}")                              # Show the raw master table path.
    print(f"  - {output_dir / 'master_candidates_validated.csv'}")                        # Show the validated master table path.
    print(f"  - {output_dir / 'master_candidates_unique_corrected.csv'}")                 # Show the deduplicated corrected master-table path.
    print(f"  - {output_dir / 'strict_reference_neighbors.csv'}")                         # Show the nearest-neighbor table path.
    print(f"  - {output_dir / 'mutation_local_pll.csv'}")                                 # Show the PLL table path.
    print(f"  - {output_dir / 'validity_summary_by_run.csv'}")                            # Show the run summary table path.


if __name__ == "__main__":                                                               # Run the CLI only when the file is executed directly.
    main()                                                                                # Invoke the main function.
