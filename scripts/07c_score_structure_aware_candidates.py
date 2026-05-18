"""Stage 07c: Score generated candidates against the seed, family manifold, and target centroid."""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from phageforge.stage07_utils import cosine_similarity, embed_sequences, read_json


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for structure-aware candidate scoring."""
    ap = argparse.ArgumentParser(description="Score Stage 07 candidates with ESM embeddings.")                                                  # Initializes the argument parser object with a description
    ap.add_argument("--context_json", type=str, required=True, help="Stage 07 context JSON produced by 07a_prepare_stage07_design_context.py.") # Adds required argument for the input context JSON path
    ap.add_argument("--generated_csv", type=str, required=True, help="Generated candidate CSV produced by 07b_generate_rbps_with_esm3.py.")     # Adds required argument for the input generated candidates CSV
    ap.add_argument("--scored_csv", type=str, required=True, help="Where to write the structure-aware scored CSV.")                             # Adds required argument for the final output CSV destination
    ap.add_argument("--esm_model", type=str, default="facebook/esm2_t33_650M_UR50D", help="ESM embedding model used for manifold scoring.")     # Adds optional argument to specify which HuggingFace model to use
    ap.add_argument("--batch_size", type=int, default=2, help="Batch size for embedding candidate sequences.")                                  # Adds optional argument to control GPU batching for memory limits
    ap.add_argument("--max_aa", type=int, default=2048, help="Maximum sequence length passed into the embedding model.")                        # Adds optional argument to truncate overly long protein sequences
    return ap.parse_args()                                                                                                                      # Parses the provided command-line inputs and returns the namespace


def main() -> None:
    # Read the Stage 07 context and the generated candidates to prepare embedding-based scoring.
    args = parse_args()                                                                                                                         # Extracts the command-line arguments into a usable object
    context = read_json(args.context_json)                                                                                                      # Loads the configuration rules and baseline vectors from the JSON file
    generated_df = pd.read_csv(args.generated_csv)                                                                                              # Loads the raw candidate strings produced by the generative model
    if generated_df.empty:                                                                                                                      # Checks if the loaded table has absolutely no rows
        raise ValueError(f"Generated CSV is empty: {args.generated_csv}")                                                                       # Halts the pipeline with an error if no inputs exist to score

    valid_df = generated_df.loc[generated_df["generation_status"].astype(str) == "ok"].copy().reset_index(drop=True)                            # Filters to keep only successfully generated proteins, resetting row indices
    if valid_df.empty:                                                                                                                          # Checks if the filtered dataframe of successful candidates is empty
        raise ValueError("No successful ESM3 candidates were found to score.")                                                                  # Halts execution because all previous generation attempts failed

    # Embed the selected seed and the candidate sequences using the requested ESM representation model.
    seed_sequence = str(context["selected_seed"]["seed_sequence"])                                                                              # Extracts the unmodified parent sequence from the configuration dictionary
    sequences = [seed_sequence] + valid_df["candidate_sequence"].astype(str).tolist()                                                           # Combines the parent and new variant strings into a single ordered list
    embeddings = embed_sequences(sequences, model_name=args.esm_model, batch_size=args.batch_size, max_length=args.max_aa)                      # Converts the text sequences into high-dimensional numerical vectors
    seed_embedding = embeddings[0]                                                                                                              # Isolates the parent vector using its known position at index 0
    candidate_embeddings = embeddings[1:]                                                                                                       # Slices the list to separate out the vectors of the new variants

    # Compare each candidate embedding to the seed, the family centroid, and the target-host centroid.
    family_centroid = np.asarray(context["family_context"].get("family_centroid", []), dtype=np.float32)                                        # Safely loads the average vector representing the broad protein family
    target_centroid = np.asarray(context["target_context"].get("target_centroid", []), dtype=np.float32)                                        # Safely loads the vector representing the new host bacterial target
    valid_df["seed_cosine"] = [cosine_similarity(emb, seed_embedding) for emb in candidate_embeddings]                                          # Computes structural similarity directly between variant and parent
    valid_df["esm_embedding_model"] = args.esm_model                                                                                            # Records the specific embedding model used for data provenance
    valid_df["family_cosine"] = [cosine_similarity(emb, family_centroid) if family_centroid.size else 0.0 for emb in candidate_embeddings]      # Computes structural drift from the family average, defaulting to 0 if missing
    valid_df["target_anchor_cosine"] = [cosine_similarity(emb, target_centroid) if target_centroid.size else 0.0 for emb in candidate_embeddings] # Computes closeness to the target host requirement, defaulting to 0 if missing
    valid_df["structure_score"] = (                                                                                                             # Begins computing the final weighted proxy score for plausibility
        0.25 * valid_df["seed_cosine"]                                                                                                          # Applies a 25% importance weight to maintaining the exact parent shape
        + 0.30 * valid_df["family_cosine"]                                                                                                      # Applies a 30% importance weight to staying inside the natural family manifold
        + 0.45 * valid_df["target_anchor_cosine"]                                                                                               # Applies a 45% importance weight to shifting binding toward the target host
    )                                                                                                                                           # Closes the arithmetic formulation for the overall score

    # Merge the scores back onto the original generated table so provenance rows are preserved as well.
    scored_df = generated_df.merge(                                                                                                             # Initiates a join to attach the new calculations back to the original dataset
        valid_df[["sample_id", "seed_cosine", "esm_embedding_model", "family_cosine", "target_anchor_cosine", "structure_score"]],              # Selects only the newly computed columns to avoid duplication
        on="sample_id",                                                                                                                         # Joins based on the unique identifier assigned in the previous stage
        how="left",                                                                                                                             # Uses a left join so that failed generation rows are still preserved in the final output
    )                                                                                                                                           # Closes the pandas merge function
    out_path = Path(args.scored_csv)                                                                                                            # Converts the string destination into a system path object
    out_path.parent.mkdir(parents=True, exist_ok=True)                                                                                          # Ensures the destination folder structure exists, ignoring if already there
    scored_df.to_csv(out_path, index=False)                                                                                                     # Saves the combined, fully scored table to disk without row numbers
    print(f"Wrote: {out_path}")                                                                                                                 # Prints a confirmation that the file was successfully saved


if __name__ == "__main__":                                                                                                                      # Checks if the file is being executed directly rather than imported
    main()                                                                                                                                      # Runs the main orchestration function