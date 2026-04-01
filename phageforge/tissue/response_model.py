"""Optional neural heads for Stage 07 tissue-context encoding and compatibility scoring."""

from __future__ import annotations                                                   # Delay annotation evaluation for cleaner modern type hints.

import torch                                                                         # Type tensors passed through the optional tissue branch.
import torch.nn as nn                                                                # Build the MLP blocks used for tissue adaptation and scoring.


# ------------------------------ tissue adapter ------------------------------ # Explain that the next class projects raw tissue vectors into a learned latent space.
class TissueContextAdapter(nn.Module):
    """Small adapter that projects tissue/context vectors into a learned hidden space."""

    def __init__(self, tissue_dim: int, hidden_dim: int = 128) -> None:
        """Create a compact MLP used to transform precomputed tissue embeddings."""
        super().__init__()                                                            # Initialize the parent PyTorch module.
        self.net = nn.Sequential(                                                     # Stack simple linear, nonlinear, and dropout layers.
            nn.Linear(tissue_dim, hidden_dim),                                        # Map the raw tissue vector into the hidden space.
            nn.GELU(),                                                                # Apply a smooth nonlinear activation.
            nn.Dropout(0.1),                                                          # Add mild regularization.
            nn.Linear(hidden_dim, hidden_dim),                                        # Refine the tissue representation in the same hidden space.
            nn.GELU(),                                                                # Apply the second nonlinear activation.
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Return the adapted tissue embedding produced by the internal MLP."""
        return self.net(x)                                                            # Pass the raw tissue tensor through the adapter network.


# ------------------------------ compatibility head ------------------------------ # Explain that the next class combines protein and tissue signals into one scalar score.
class TissueCompatibilityHead(nn.Module):
    """Predict a scalar compatibility score from protein and optional tissue embeddings."""

    def __init__(self, protein_dim: int, tissue_dim: int, hidden_dim: int = 256) -> None:
        """Create the small MLP used to fuse protein and tissue representations."""
        super().__init__()                                                            # Initialize the parent PyTorch module.
        self.net = nn.Sequential(                                                     # Build the final compatibility-scoring network.
            nn.Linear(protein_dim + tissue_dim, hidden_dim),                          # Fuse the protein and tissue vectors into one hidden representation.
            nn.GELU(),                                                                # Apply a nonlinear activation after fusion.
            nn.Dropout(0.1),                                                          # Add light regularization.
            nn.Linear(hidden_dim, 1),                                                 # Predict one scalar compatibility score.
        )

    def forward(self, protein_emb: torch.Tensor, tissue_emb: torch.Tensor | None = None) -> torch.Tensor:
        """Return one compatibility score per protein, optionally conditioned on tissue."""
        if tissue_emb is None:                                                        # Handle the fallback mode where no tissue context is available.
            tissue_emb = torch.zeros(                                                 # Create an empty-width tensor that preserves the batch dimension.
                protein_emb.shape[0],                                                 # Match the number of protein examples in the batch.
                0,                                                                    # Use zero extra features when no tissue signal is present.
                device=protein_emb.device,                                            # Place the fallback tensor on the same device as the protein tensor.
                dtype=protein_emb.dtype,                                              # Use the same numeric type as the protein tensor.
            )
        x = torch.cat([protein_emb, tissue_emb], dim=-1)                              # Concatenate the protein and tissue embeddings along the feature axis.
        return self.net(x).squeeze(-1)                                                # Predict one scalar score per row and remove the singleton dimension.
