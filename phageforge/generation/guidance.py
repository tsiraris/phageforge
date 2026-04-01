"""Helper functions for Stage 07 seed encoding, hotspot masking, and edit sampling."""

from __future__ import annotations                                                   # Delay annotation evaluation for cleaner typing.
from typing import Iterable, List                                                    # Use explicit typing for token lists and iterables.
import torch                                                                         # Sample amino-acid replacements from model logits.
from .dataset import AMINO_ACIDS, BOS_ID, EOS_ID, ID_TO_TOKEN, MASK_ID, TOKEN_TO_ID  # Reuse the shared token vocabulary constants.


# Sequence conversion helpers to convert between strings and token IDs.
def encode_seed(seq: str) -> List[int]:
    """Convert a raw amino-acid seed sequence into token IDs with BOS/EOS markers."""
    seq = "".join([aa for aa in str(seq).strip().upper() if aa in AMINO_ACIDS])       # Normalize to uppercase and keep only canonical amino acids.
    return [BOS_ID] + [TOKEN_TO_ID[aa] for aa in seq] + [EOS_ID]                      # Wrap the sequence with BOS/EOS and map residues to IDs.


def decode_ids(ids: Iterable[int]) -> str:
    """Convert token IDs back into a plain amino-acid sequence string."""
    toks = [ID_TO_TOKEN[int(i)] for i in ids if int(i) in ID_TO_TOKEN]                # Convert every valid ID back into its token string.
    return "".join([tok for tok in toks if tok in AMINO_ACIDS])                        # Keep only amino-acid tokens and drop special symbols.


def parse_hotspots(text: str) -> List[int]:
    """Parse a comma-separated hotspot string into sorted unique 0-based positions."""
    if not text:                                                                      # Handle the empty-string case used by some CLI calls.
        return []                                                                     # Return an empty hotspot list when no string was provided.
    values = sorted({int(x.strip()) for x in str(text).split(",") if x.strip()})      # Parse integers, deduplicate them, and sort the result.
    return [v for v in values if v >= 0]                                              # Keep only non-negative positions.


# Masking and sampling helpers to define where and how edits are sampled.
def hotspot_masked_ids(seed_ids: List[int], hotspot_positions: List[int]) -> List[int]:
    """Replace hotspot token positions with mask tokens while leaving BOS/EOS intact."""
    ids = list(seed_ids)                                                              # Copy the seed token list so the original is not modified in place.
    for pos in hotspot_positions:                                                     # Iterate over requested 0-based residue positions.
        token_pos = pos + 1                                                           # Shift by one because BOS occupies token position 0.
        if 0 < token_pos < len(ids) - 1:                                              # Only mask true internal sequence positions.
            ids[token_pos] = MASK_ID                                                  # Replace the residue token with the mask token ID.
    return ids                                                                        # Return the masked token sequence.


def sample_hotspot_edits(
    logits: torch.Tensor,
    base_ids: List[int],
    hotspot_positions: List[int],
    temperature: float = 1.0,
    top_k: int = 5,
) -> List[int]:
    """
    Sample amino-acid replacements at hotspot positions from token logits.

    The function supports temperature scaling and optional top-k restriction so
    generation remains controllable and stays close to the learned sequence space.
    """
    sampled = list(base_ids)                                                          # Start from the base token sequence and edit selected positions only.
    for pos in hotspot_positions:                                                     # Visit each editable hotspot position.
        token_pos = pos + 1                                                           # Shift by one because BOS occupies the first token slot.
        if token_pos <= 0 or token_pos >= logits.shape[0] - 1:                        # Skip invalid or boundary positions.
            continue                                                                  # Continue to the next hotspot without sampling.
        row = logits[token_pos] / max(float(temperature), 1e-6)                       # Apply temperature scaling to the token logits.
        probs = torch.softmax(row, dim=-1)                                            # Convert the scaled logits into a probability distribution.
        if 0 < top_k < probs.numel():                                                 # Restrict sampling to the top-k tokens when requested.
            top_vals, top_idx = torch.topk(probs, k=top_k)                            # Keep only the top-k token probabilities and IDs.
            top_vals = top_vals / top_vals.sum()                                      # Renormalize the top-k probabilities so they sum to one.
            chosen = top_idx[torch.multinomial(top_vals, num_samples=1)[0]].item()    # Sample one token from the restricted top-k distribution.
        else:                                                                         # Otherwise sample from the full vocabulary distribution.
            chosen = torch.multinomial(probs, num_samples=1)[0].item()                # Draw one token ID from the full probability vector.
        sampled[token_pos] = int(chosen)                                              # Insert the sampled token into the editable position.
    return sampled                                                                    # Return the edited token sequence.


def mutation_burden(seed_seq: str, candidate_seq: str) -> int:
    """Count how many residue changes separate the candidate sequence from the seed."""
    n = min(len(seed_seq), len(candidate_seq))                                        # Compare only over the shared prefix length first.
    burden = sum(1 for i in range(n) if seed_seq[i] != candidate_seq[i])              # Count mismatched residues over the shared prefix.
    burden += abs(len(seed_seq) - len(candidate_seq))                                 # Add any extra insertions or deletions as additional burden.
    return burden                                                                     # Return the total mutation burden.


def approximate_edit_span(hotspots: List[int], flank: int, seq_len: int) -> tuple[int, int]:
    """Convert hotspot positions into one approximate editable window span."""
    if len(hotspots) == 0:                                                            # If no hotspots exist, expose the full sequence span.
        return 0, seq_len                                                             # Return a whole-sequence span.
    lo = max(0, min(hotspots) - flank)                                                # Extend slightly before the first hotspot while staying in bounds.
    hi = min(seq_len, max(hotspots) + flank + 1)                                      # Extend slightly after the last hotspot while staying in bounds.
    return lo, hi                                                                     # Return the inferred 0-based half-open span.
