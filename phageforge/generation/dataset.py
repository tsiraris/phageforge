"""Stage 07 dataset utilities for conditional RBP generation.
This 
"""

from __future__ import annotations                                                   # Delay type annotation evaluation for cleaner forward references.
from dataclasses import dataclass                                                    # Define a compact typed container for one training example.
from typing import Dict, List                                                        # Use explicit type hints for mappings and token lists.
import pandas as pd                                                                  # Read and normalize generator training tables.
import torch                                                                         # Build PyTorch tensors returned by the dataset.
from torch.utils.data import Dataset                                                 # Implement a standard PyTorch dataset class.


# Define the amino-acid vocabulary
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")                                           # Store the 20 canonical amino acids accepted by the tokenizer.
SPECIAL = ["<pad>", "<mask>", "<bos>", "<eos>"]                                      # Define padding, masking, and boundary tokens used by the model.
VOCAB = SPECIAL + AMINO_ACIDS                                                        # Concatenate special tokens and amino acids into one vocabulary.
TOKEN_TO_ID = {tok: i for i, tok in enumerate(VOCAB)}                                # Create a fast token-to-index lookup table.
ID_TO_TOKEN = {i: tok for tok, i in TOKEN_TO_ID.items()}                             # Create the inverse index-to-token lookup table.
PAD_ID = TOKEN_TO_ID["<pad>"]                                                        # Cache the padding token ID for collators and models.
MASK_ID = TOKEN_TO_ID["<mask>"]                                                      # Cache the mask token ID used during denoising training.
BOS_ID = TOKEN_TO_ID["<bos>"]                                                        # Cache the beginning-of-sequence token ID.
EOS_ID = TOKEN_TO_ID["<eos>"]                                                        # Cache the end-of-sequence token ID.


# Define the Stage 07 example type
@dataclass
class Stage07Example:
    """Store one Stage 07 training example in a typed, readable form."""

    sequence: str                                                                     # Keep the raw amino-acid sequence.
    host_genus: str                                                                   # Keep the conditioning host genus label.
    family_id: str                                                                    # Keep the conditioning family label.


# Convert Dataset rows into token IDs and condition labels, and return as a dictionary of PyTorch tensors consumed by the collator.
class Stage07SequenceDataset(Dataset):
    """PyTorch dataset that tokenizes Stage 07 sequences and conditioning labels."""

    def __init__(
        self,
        df: pd.DataFrame,
        host_to_id: Dict[str, int],
        family_to_id: Dict[str, int],
        max_len: int = 1022,
    ) -> None:
        """Store the Stage 07 table and the label vocabularies used for conditioning."""
        self.df = df.reset_index(drop=True)                                           # Reset the row index so dataset row IDs remain contiguous.
        self.host_to_id = host_to_id                                                  # Store the host-genus mapping used to build condition IDs.
        self.family_to_id = family_to_id                                              # Store the family mapping used to build condition IDs.
        self.max_len = int(max_len)                                                   # Store the maximum tokenized sequence length allowed.

    def __len__(self) -> int:
        """Return the number of rows available for training or validation."""
        return len(self.df)                                                           # Report dataset length directly from the stored DataFrame.

    def encode_seq(self, seq: str) -> List[int]:
        """Convert a raw amino-acid sequence into model token IDs.

        The function removes unsupported characters, truncates to the allowed length,
        and wraps the sequence with BOS and EOS boundary tokens.
        """
        seq = str(seq).strip().upper()                                                # Normalize the sequence to an uppercase stripped string.
        seq = "".join([aa for aa in seq if aa in AMINO_ACIDS])                        # Keep only canonical amino acids accepted by the vocabulary.
        seq = seq[: self.max_len - 2]                                                 # Truncate so BOS and EOS still fit inside the length limit.
        ids = [BOS_ID] + [TOKEN_TO_ID[a] for a in seq] + [EOS_ID]                     # Encode the sequence and add the boundary tokens.
        return ids                                                                    # Return the token ID list ready for batching.

    def __getitem__(self, idx: int) -> dict:
        """Return one dataset item containing token IDs and condition labels."""
        row = self.df.iloc[idx]                                                       # Select the requested DataFrame row by integer index.
        ids = self.encode_seq(row["sequence"])                                        # Tokenize the sequence column into model IDs.
        return {                                                                      # Return a dictionary consumed by the collator.
            "input_ids": torch.tensor(ids, dtype=torch.long),                         # Store sequence token IDs as a long tensor.
            "host_id": torch.tensor(self.host_to_id[str(row["host_genus"])], dtype=torch.long),     # Map the host label to an integer condition ID.
            "family_id": torch.tensor(self.family_to_id[str(row["family_id"])], dtype=torch.long),  # Map the family label to an integer condition ID.
            "row_id": torch.tensor(idx, dtype=torch.long),                            # Preserve the original row index for debugging and traceability.
        }
