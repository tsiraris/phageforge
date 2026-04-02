
"""
==============================================================================
Stage 07c: Score generated candidates with sequence + structure-aware proxies.
==============================================================================

This script adds lightweight structure-aware reranking features to Stage 07.
It is intentionally practical rather than claiming full protein folding:
- embed seed and generated candidates with ESM
- measure seed similarity
- measure family-centroid similarity when reference embeddings are available
- otherwise fall back to the family / target centroids already stored in the Stage 07 context
- combine these terms into one structure-aware proxy score
"""

from __future__ import annotations                                                    # Enable postponed annotation evaluation for cleaner typing.
import argparse                                                                       # Parse command-line arguments.
import json                                                                           # Load the Stage 07 context JSON.
from pathlib import Path                                                              # Work with filesystem paths robustly.
import numpy as np                                                                    # Compute cosine similarities and numeric score combinations.
import pandas as pd                                                                   # Read and write candidate tables.
import torch                                                                          # Run embedding inference and load cached tensors.
from transformers import AutoTokenizer, EsmModel                                      # Use ESM embeddings as a practical structure-aware sequence proxy.


def parse_args() -> argparse.Namespace:
    """Return command-line arguments for Stage 07 structure-aware candidate scoring."""
    ap = argparse.ArgumentParser(description="Score Stage 07 generated candidates with structure-aware proxies.")  # Create the CLI parser.
    ap.add_argument("--context_json", type=str, required=True, help="Context JSON from 07a_prepare_stage07_design_context.py.")  # Point to the Stage 07 context JSON.
    ap.add_argument("--generated_csv", type=str, required=True, help="Generated candidate CSV from 07b.")  # Point to the generated-candidate table.
    ap.add_argument("--reference_embeddings_pt", type=str, default="", help="Optional broad embedding tensor from 02_embed_rbps.py.")  # Provide optional cached reference embeddings.
    ap.add_argument("--reference_index_csv", type=str, default="", help="Optional index CSV matching the reference embedding tensor.")  # Provide optional metadata aligned to reference embeddings.
    ap.add_argument("--scored_csv", type=str, required=True, help="Where to write the structure-scored candidate table.")  # Point to the final scored output CSV.
    ap.add_argument("--esm_model", type=str, default="auto", help="ESM checkpoint used for practical sequence embeddings. 'auto' chooses a dimension-compatible model.")  # Default to automatic model selection so context centroids can be used cleanly.
    ap.add_argument("--batch_size", type=int, default=8, help="Embedding batch size.")  # Control inference batch size.
    ap.add_argument("--max_aa", type=int, default=1022, help="Maximum amino-acid length for embedding.")  # Control the ESM truncation limit.
    return ap.parse_args()                                                            # Parse the CLI and return the namespace.


# Embedding helpers to normalize embeddings and run ESM inference.
def l2_normalize(x: np.ndarray) -> np.ndarray:
    """Return row-wise L2-normalized vectors suitable for cosine similarity."""
    norms = np.linalg.norm(x, axis=1, keepdims=True)                                  # Compute one vector norm per embedding row.
    norms = np.clip(norms, 1e-12, None)                                               # Prevent divide-by-zero for degenerate vectors.
    return x / norms                                                                  # Divide every row by its norm to obtain unit vectors.


def mean_pool(hidden: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    """Average token embeddings over valid positions only."""
    mask = mask.unsqueeze(-1).to(hidden.dtype)                                        # Expand the attention mask to match hidden-state shape.
    summed = (hidden * mask).sum(dim=1)                                               # Sum token embeddings only over valid positions.
    denom = mask.sum(dim=1).clamp(min=1.0)                                            # Count how many valid tokens contributed to each sum.
    return summed / denom                                                             # Return the mean pooled embedding per sequence.


def load_context(path: Path) -> dict:
    """Load the Stage 07 context JSON from disk."""
    with open(path, "r", encoding="utf-8") as handle:                                 # Open the JSON file in UTF-8 mode.
        return json.load(handle)                                                      # Parse and return the context dictionary.


def infer_target_centroid(context: dict) -> list[float]:
    """Return the best available target-anchor centroid from the Stage 07 context."""
    target_block = context.get("target_anchor_context", {})                           # Read the target-anchor block from the Stage 07 context.
    centroid = target_block.get("target_reference_centroid", [])                      # Prefer the already-selected target-host centroid when present.
    if centroid:                                                                      # Return the direct centroid when it exists.
        return centroid
    refs = target_block.get("target_anchor_references", [])                           # Fall back to raw reference rows when only those are present.
    if refs and isinstance(refs[0], dict):                                            # Guard against non-dictionary reference rows.
        for key in ["embedding", "embedding_vector", "centroid"]:                    # Check a few plausible embedding key names.
            vectors = [x.get(key) for x in refs if isinstance(x.get(key), list)]
            if vectors:
                return np.asarray(vectors, dtype=np.float32).mean(axis=0).tolist()    # Average reference vectors into one centroid if needed.
    return []                                                                         # Return an empty list when no target centroid is available.


def choose_esm_model(args_model: str, context: dict) -> str:
    """Choose an ESM model whose hidden size matches the stored Stage 07 centroids when possible."""
    if args_model != "auto":                                                          # Respect an explicit user override.
        return args_model
    family_centroid = context.get("family_context", {}).get("family_centroid", [])   # Read the stored family centroid from the context.
    target_centroid = infer_target_centroid(context)                                  # Read the stored target-anchor centroid from the context.
    candidate_len = len(family_centroid) or len(target_centroid)                      # Infer the intended embedding dimensionality from the context.
    if candidate_len == 1280:                                                         # Stage 06 centroids in this repo were built with ESM2 650M embeddings.
        return "facebook/esm2_t33_650M_UR50D"
    if candidate_len == 480:                                                          # Support smaller centroids when they exist.
        return "facebook/esm2_t12_35M_UR50D"
    return "facebook/esm2_t33_650M_UR50D"                                             # Prefer the richer default when the context does not specify a dimension.


def embed_sequences(sequences: list[str], model_name: str, batch_size: int, max_aa: int) -> np.ndarray:
    """
    Embed a list of amino-acid sequences with an ESM encoder and mean pooling.
    Returns a matrix of shape [N, H] where N is the number of sequences and H is the embedding dimension.
    """
    device = "cuda" if torch.cuda.is_available() else "cpu"                           # Prefer GPU for ESM inference when available.
    tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)        # Load the raw-sequence tokenizer for the chosen ESM checkpoint.
    model = EsmModel.from_pretrained(model_name).to(device).eval()                    # Load the ESM encoder and switch it to eval mode.
    rows = []                                                                         # Collect pooled embeddings batch by batch.

    with torch.no_grad():                                                             # Disable gradient tracking because this is inference only.
        for start in range(0, len(sequences), batch_size):                            # Iterate through the sequence list in mini-batches.
            batch = sequences[start:start + batch_size]                               # Slice one batch of raw sequences.
            toks = tokenizer(                                                         # Tokenize the batch into input IDs and attention masks.
                batch,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=max_aa,
            )
            toks = {k: v.to(device) for k, v in toks.items()}                         # Move the token tensors to the selected device.
            out = model(**toks)                                                       # Run the ESM encoder forward on the tokenized batch.
            pooled = mean_pool(out.last_hidden_state, toks["attention_mask"])         # Collapse token embeddings into one vector per sequence.
            rows.append(pooled.cpu())                                                 # Move pooled embeddings back to CPU and store them.

    return torch.cat(rows, dim=0).numpy()                                             # Concatenate all batches and return a NumPy array.


def add_reference_scores(gen: pd.DataFrame, context: dict, ref_emb: np.ndarray, ref_idx: pd.DataFrame) -> pd.DataFrame:
    """Attach family and target-anchor centroid similarities using cached reference embeddings."""
    out = gen.copy()                                                                  # Work on a copy so the caller keeps the input DataFrame unchanged.
    row_col = "row_id" if "row_id" in ref_idx.columns else ("idx" if "idx" in ref_idx.columns else None)  # Detect the embedding-row index column.
    if row_col is None:                                                               # Ensure the index table can point into the embedding tensor.
        raise ValueError("reference_index_csv must contain row_id or idx.")           # Fail with a clear error if no row-alignment column exists.
    ref_idx = ref_idx.copy()                                                          # Copy the index table before adding new columns.
    ref_idx["embedding_row"] = ref_idx[row_col].astype(int)                          # Normalize the reference embedding row selector to integers.

    family_member_ids = set(str(x) for x in context.get("family_context", {}).get("family_member_ids", []))  # Read the Stage 06 family member IDs.
    family_rows = ref_idx[ref_idx["protein_id"].astype(str).isin(family_member_ids)] if "protein_id" in ref_idx.columns else pd.DataFrame()  # Keep only family members when possible.
    if not family_rows.empty:                                                         # Compute a family centroid only if family members are available.
        family_matrix = ref_emb[family_rows["embedding_row"].to_numpy()]             # Gather family-member embeddings from the cached tensor.
        family_centroid = l2_normalize(family_matrix.mean(axis=0, keepdims=True)).reshape(-1)  # Compute one normalized family centroid.
        out["family_cosine"] = out["candidate_emb"].apply(lambda x: float(np.dot(x, family_centroid)))  # Score every candidate against the family centroid.
    else:
        out["family_cosine"] = np.nan                                                # Leave as NaN so a later fallback branch can repair it cleanly.

    anchor_refs = context.get("target_anchor_context", {}).get("target_anchor_references", [])  # Read the target-anchor records from the Stage 07 context.
    anchor_ids = {str(x.get("protein_id")) for x in anchor_refs if isinstance(x, dict) and x.get("protein_id") is not None}  # Extract the anchor protein IDs.
    anchor_rows = ref_idx[ref_idx["protein_id"].astype(str).isin(anchor_ids)] if anchor_ids and "protein_id" in ref_idx.columns else pd.DataFrame()  # Keep only matched anchor rows.
    if not anchor_rows.empty:                                                         # Compute the anchor centroid only if anchor proteins are available.
        anchor_matrix = ref_emb[anchor_rows["embedding_row"].to_numpy()]             # Gather anchor embeddings from the cached tensor.
        anchor_centroid = l2_normalize(anchor_matrix.mean(axis=0, keepdims=True)).reshape(-1)  # Compute one normalized anchor centroid.
        out["target_anchor_cosine"] = out["candidate_emb"].apply(lambda x: float(np.dot(x, anchor_centroid)))  # Score every candidate against the anchor centroid.
    else:
        out["target_anchor_cosine"] = np.nan                                         # Leave as NaN so a later fallback branch can repair it cleanly.

    return out                                                                        # Return the DataFrame with new family and anchor similarity columns.


def add_context_centroid_scores(gen: pd.DataFrame, context: dict) -> pd.DataFrame:
    """Attach family and target-anchor cosine scores directly from centroids stored in the Stage 07 context."""
    out = gen.copy()                                                                  # Work on a copy so the caller keeps the input table unchanged.
    emb_dim = len(out.iloc[0]["candidate_emb"]) if len(out) else 0                    # Infer candidate embedding dimensionality from the first embedding row.

    family_centroid = np.asarray(context.get("family_context", {}).get("family_centroid", []), dtype=np.float32)  # Read the stored family centroid.
    if family_centroid.size == emb_dim and emb_dim > 0:                               # Only use the centroid when dimensions match the candidate embeddings.
        family_centroid = l2_normalize(family_centroid.reshape(1, -1)).reshape(-1)
        out["family_cosine"] = out["candidate_emb"].apply(lambda x: float(np.dot(x, family_centroid)))
    elif "family_cosine" not in out.columns:
        out["family_cosine"] = 0.0                                                   # Fill with zeros only when no other family score exists.

    target_centroid = np.asarray(infer_target_centroid(context), dtype=np.float32)    # Read the stored target-anchor centroid.
    if target_centroid.size == emb_dim and emb_dim > 0:                               # Only use the centroid when dimensions match the candidate embeddings.
        target_centroid = l2_normalize(target_centroid.reshape(1, -1)).reshape(-1)
        out["target_anchor_cosine"] = out["candidate_emb"].apply(lambda x: float(np.dot(x, target_centroid)))
    elif "target_anchor_cosine" not in out.columns:
        out["target_anchor_cosine"] = 0.0                                            # Fill with zeros only when no other anchor score exists.

    return out                                                                        # Return the candidate table with repaired family / anchor scores.


def finalize_structure_scores(gen: pd.DataFrame) -> pd.DataFrame:
    """Fill missing score columns and compute one compact structure-aware proxy score."""
    out = gen.copy()                                                                  # Work on a copy so the caller keeps the original table unchanged.
    for col in ["family_cosine", "target_anchor_cosine"]:                            # Ensure both similarity columns exist before the final score is built.
        if col not in out.columns:
            out[col] = 0.0
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0.0)               # Replace NaNs from optional branches with neutral zeros.
    out["seed_cosine"] = pd.to_numeric(out["seed_cosine"], errors="coerce").fillna(0.0)  # Normalize seed similarity before scoring.
    out["structure_score"] = (                                                       # Combine the proxy terms into one pragmatic structure-aware score.
        0.45 * out["seed_cosine"].astype(float) +
        0.35 * out["family_cosine"].astype(float) +
        0.20 * out["target_anchor_cosine"].astype(float)
    )
    return out                                                                        # Return the candidate table with a final structure score.


# Main script: Reads Stage 07b generated candidates, embeds candidates and writes the scored CSV.
def main() -> None:
    args = parse_args()                                                               # Parse command-line arguments.
    context = load_context(Path(args.context_json))                                   # Load the Stage 07 context JSON from disk.
    gen = pd.read_csv(args.generated_csv)                                             # Read the generated candidate table from Stage 07b.
    if len(gen) == 0:                                                                 # Ensure the input table is non-empty before embedding.
        raise ValueError("generated_csv is empty.")

    esm_model = choose_esm_model(args.esm_model, context)                             # Choose a dimension-compatible embedding model automatically when requested.
    seed_seq = context["selected_seed"]["seed_sequence"]                              # Recover the seed sequence from the Stage 07 context.
    all_sequences = [seed_seq] + gen["candidate_sequence"].astype(str).tolist()       # Embed the seed first, followed by all generated candidates.
    emb = l2_normalize(embed_sequences(all_sequences, esm_model, args.batch_size, args.max_aa))  # Embed and normalize all sequences.
    seed_emb = emb[0]                                                                 # Keep the first embedding as the seed embedding.
    cand_emb = emb[1:]                                                                # Keep the remaining embeddings as candidate embeddings.

    gen = gen.copy()                                                                  # Work on a copy of the generated candidate table.
    gen["candidate_emb"] = list(cand_emb)                                             # Store candidate embeddings temporarily for later similarity calculations.
    gen["seed_cosine"] = [float(np.dot(x, seed_emb)) for x in cand_emb]               # Score each candidate by cosine similarity to the seed sequence.
    gen["esm_embedding_model"] = esm_model                                             # Record which embedding model produced the structure-aware proxies.

    if args.reference_embeddings_pt and args.reference_index_csv:                     # Prefer direct cached reference embeddings when both files were provided.
        ref_emb = torch.load(args.reference_embeddings_pt, map_location="cpu")        # Load the cached reference embedding tensor from disk.
        if isinstance(ref_emb, torch.Tensor):
            ref_emb = ref_emb.numpy()
        ref_emb = l2_normalize(np.asarray(ref_emb, dtype=np.float32))                 # Normalize cached reference embeddings for cosine similarity.
        ref_idx = pd.read_csv(args.reference_index_csv)                               # Read the metadata aligned to the cached reference tensor.
        gen = add_reference_scores(gen, context, ref_emb, ref_idx)                    # Attach family and target-anchor scores from the reference bank.

    gen = add_context_centroid_scores(gen, context)                                   # Repair missing family / anchor scores directly from the Stage 07 context centroids.
    gen = finalize_structure_scores(gen)                                              # Fill missing values and compute the final structure-aware proxy score.
    gen = gen.drop(columns=["candidate_emb"])                                         # Remove the temporary high-dimensional embedding column before writing CSV.

    scored_path = Path(args.scored_csv)                                               # Convert the scored output path into a Path object.
    scored_path.parent.mkdir(parents=True, exist_ok=True)                             # Create the output directory if needed.
    gen.to_csv(scored_path, index=False)                                              # Save the scored candidate table to disk.
    print(f"Wrote: {scored_path}")                                                    # Print the output path for quick confirmation.


if __name__ == "__main__":                                                            # Standard Python entrypoint guard.
    main()                                                                            # Execute the structure-aware scoring CLI.
