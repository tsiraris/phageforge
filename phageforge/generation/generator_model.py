"""Conditional masked generator (Transformer Encoder) used for Stage 07 hotspot-limited sequence editing."""

from __future__ import annotations                                                   # Delay type annotation evaluation for cleaner forward references.
import torch                                                                         # Build and type the input/output tensors passed through the network.
import torch.nn as nn                                                                # Assemble the transformer layers and embedding blocks.


# The main conditional sequence generator.
class ConditionalMaskedRBPGenerator(nn.Module):
    """Lightweight conditional transformer for masked amino-acid reconstruction."""

    def __init__(
        self,
        vocab_size: int,
        n_hosts: int,
        n_families: int,
        d_model: int = 512,
        n_heads: int = 8,
        n_layers: int = 6,
        dropout: float = 0.1,
        max_len: int = 1024,
    ) -> None:
        """Create token, position, host, and family embeddings plus a transformer encoder."""
        super().__init__()                                                            # Initialize the base PyTorch module.
        self.token_emb = nn.Embedding(vocab_size, d_model)                            # Learn an embedding vector for every vocabulary token.
        self.pos_emb = nn.Embedding(max_len, d_model)                                 # Learn a positional embedding for each sequence index.
        self.host_emb = nn.Embedding(n_hosts, d_model)                                # Learn one conditioning vector for each host genus.
        self.family_emb = nn.Embedding(n_families, d_model)                           # Learn one conditioning vector for each scaffold family.
        self.dropout = nn.Dropout(dropout)                                            # Regularize the model with dropout before the encoder.
        layer = nn.TransformerEncoderLayer(                                           # Build one reusable transformer encoder block.
            d_model=d_model,                                                          # Set the hidden dimensionality of the model.
            nhead=n_heads,                                                            # Set the number of self-attention heads.
            dim_feedforward=4 * d_model,                                              # Set the feed-forward width inside the transformer block.
            dropout=dropout,                                                          # Reuse the same dropout value inside the transformer block.
            activation="gelu",                                                        # Use GELU as the nonlinear activation.
            batch_first=True,                                                         # Keep tensors in [batch, length, hidden] format.
            norm_first=True,                                                          # Apply layer normalization before the sublayers for stability.
        )
        self.encoder = nn.TransformerEncoder(layer, num_layers=n_layers)              # Stack several encoder layers into the final sequence model.
        self.norm = nn.LayerNorm(d_model)                                             # Normalize hidden states before token prediction.
        self.head = nn.Linear(d_model, vocab_size)                                    # Predict a vocabulary distribution at every sequence position.

    def forward(
        self,
        input_ids: torch.Tensor,
        attention_mask: torch.Tensor,
        host_ids: torch.Tensor,
        family_ids: torch.Tensor,
    ) -> torch.Tensor:
        """Return token logits for a batch of masked sequences under host/family conditioning."""
        batch_size, seq_len = input_ids.shape                                         # Read the batch size and sequence length from the token matrix.
        positions = torch.arange(seq_len, device=input_ids.device).unsqueeze(0).expand(batch_size, seq_len)  # Build position indices for every example.
        x = self.token_emb(input_ids) + self.pos_emb(positions)                       # Add token and positional embeddings together.
        condition = self.host_emb(host_ids).unsqueeze(1) + self.family_emb(family_ids).unsqueeze(1)  # Create one broadcastable conditioning vector per example.
        x = self.dropout(x + condition)                                               # Inject the conditioning signal and regularize with dropout.
        key_padding_mask = attention_mask == 0                                        # Mark padded positions so attention ignores them.
        x = self.encoder(x, src_key_padding_mask=key_padding_mask)                    # Run the conditioned sequence through the transformer stack.
        x = self.norm(x)                                                              # Normalize the hidden states before projection.
        return self.head(x)                                                           # Return token logits used for masked-token prediction.
