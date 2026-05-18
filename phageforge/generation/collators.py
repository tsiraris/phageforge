"""
Batch collation helpers for Stage 07 masked-denoising training.
Builds masked denoising batches for the Stage 07 local conditional generator, hiding 15% of the each sequence.
"""

from __future__ import annotations                                                   # Delay type annotation evaluation for cleaner typing.
from typing import List                                                              # Express that the collator receives a list of dataset examples.
import torch                                                                         # Build padded batch tensors and random masks.
from .dataset import BOS_ID, EOS_ID, MASK_ID, PAD_ID                                 # Reuse shared special-token IDs from the dataset module.


# A class that pads sequences and creates masked labels.
class MaskedConditioningCollator:
    """
    Randomly mask internal amino-acid tokens for conditional denoising training.
    Return the full training batch expected by the generator.
    """

    def __init__(self, mask_prob: float = 0.15) -> None:
        """Store the probability of masking each non-special amino-acid position."""
        self.mask_prob = float(mask_prob)                                             # Convert the user-supplied masking probability to a float.

    def __call__(self, batch: List[dict]) -> dict:
        """Pad a batch of variable-length sequences and create masked-token labels."""
        max_len = max(item["input_ids"].shape[0] for item in batch)                   # Find the longest sequence length in the incoming batch.
        input_ids = torch.full((len(batch), max_len), PAD_ID, dtype=torch.long)       # Preallocate a padded token matrix filled with PAD tokens.
        labels = torch.full((len(batch), max_len), -100, dtype=torch.long)            # Preallocate label targets with ignore-index values for CE loss.
        attention_mask = torch.zeros((len(batch), max_len), dtype=torch.long)         # Preallocate the attention mask marking real versus padded tokens.
        host_ids = torch.stack([item["host_id"] for item in batch])                   # Stack per-example host condition IDs into one tensor.
        family_ids = torch.stack([item["family_id"] for item in batch])               # Stack per-example family condition IDs into one tensor.
        row_ids = torch.stack([item["row_id"] for item in batch])                     # Stack per-example row IDs for tracing back predictions.

        # Fill the padded tensors one example at a time and decide which positions to mask.
        for i, item in enumerate(batch):                                              # Iterate over examples in the incoming mini-batch.
            ids = item["input_ids"]                                                   # Read the raw token IDs for this example.
            n = ids.shape[0]                                                          # Measure the true sequence length including BOS/EOS.
            input_ids[i, :n] = ids                                                    # Copy the token IDs into the padded batch matrix.
            attention_mask[i, :n] = 1                                                 # Mark the copied positions as valid tokens.

            rand = torch.rand(n)                                                      # Sample one random value per token position.
            maskable = (ids != PAD_ID) & (ids != BOS_ID) & (ids != EOS_ID)            # Forbid masking padding and boundary tokens.
            chosen = (rand < self.mask_prob) & maskable                               # Choose internal amino-acid positions to mask.
            labels[i, chosen] = ids[chosen]                                           # Store the original tokens as prediction targets.
            input_ids[i, chosen] = MASK_ID                                            # Replace the chosen tokens with the mask token.

        return {                                                                      # Return the full training batch expected by the generator.
            "input_ids": input_ids,                                                  # Provide masked token IDs to the model.
            "labels": labels,                                                        # Provide original token labels for loss computation.
            "attention_mask": attention_mask,                                        # Provide the mask separating real tokens from padding.
            "host_ids": host_ids,                                                    # Provide host conditioning labels.
            "family_ids": family_ids,                                                # Provide family conditioning labels.
            "row_ids": row_ids,                                                      # Provide row IDs for debugging or future analysis.
        }
