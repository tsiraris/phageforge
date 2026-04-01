"""
===========================================================================
Stage 07b: Train the Stage 07 local conditional masked denoising generator.
===========================================================================

This script trains the (fallback) practical local generator model. 
Trains a model that given a masked protein sequence, plus the target host and family, learns to predict what amino acid should go in the masked positions.
The model learns to reconstruct masked amino acids while being conditioned on:
- host genus
- scaffold family
"""

from __future__ import annotations                                                    # Enable postponed annotation evaluation for cleaner typing.
import argparse                                                                       # Parse command-line arguments.
import json                                                                           # Save training history as JSON.
from pathlib import Path                                                              # Build output paths robustly.
import pandas as pd                                                                   # Load the generator dataset CSV.
import torch                                                                         # Train the PyTorch generator model.
import torch.nn.functional as F                                                      # Use cross-entropy loss for masked-token prediction.
from torch.utils.data import DataLoader                                              # Batch the dataset during training.
from phageforge.generation.collators import MaskedConditioningCollator               # Build masked denoising batches.
from phageforge.generation.dataset import Stage07SequenceDataset, VOCAB              # Tokenize sequences and expose the vocabulary size.
from phageforge.generation.generator_model import ConditionalMaskedRBPGenerator      # Instantiate the conditional transformer generator.


def parse_args() -> argparse.Namespace:
    """Return command-line arguments for Stage 07 generator training."""
    ap = argparse.ArgumentParser(description="Train the Stage 07 conditional generator.")  # Create the parser shown in `--help` output.
    ap.add_argument("--dataset_csv", type=str, required=True, help="CSV produced by Stage 07 dataset preparation.")  # Point to the prepared generator dataset.
    ap.add_argument("--out_dir", type=str, required=True, help="Directory where checkpoints and training history will be written.")  # Point to the output directory.
    ap.add_argument("--epochs", type=int, default=20, help="Number of full training epochs to run.")  # Control training duration.
    ap.add_argument("--batch_size", type=int, default=16, help="Mini-batch size used by the DataLoader.")  # Control batch size.
    ap.add_argument("--lr", type=float, default=3e-4, help="Learning rate for AdamW.")  # Control the optimizer learning rate.
    ap.add_argument("--mask_prob", type=float, default=0.15, help="Probability of masking each internal amino-acid token.")  # Control denoising difficulty.
    ap.add_argument("--max_len", type=int, default=1022, help="Maximum amino-acid length including BOS/EOS room.")  # Control the input truncation limit.
    return ap.parse_args()                                                            # Parse the CLI and return the arguments namespace.


def build_maps(df: pd.DataFrame) -> tuple[dict[str, int], dict[str, int]]:
    """
    Build integer label vocabularies for hosts and family IDs.
    Example:
        Acinetobacter -> 0, Klebsiella -> 1
        and 
        tail_fiber_family_A -> 0, tail_spike_family_B -> 1
    """
    hosts = sorted(df["host_genus"].astype(str).unique().tolist())                    # Collect and sort all unique host genera.
    families = sorted(df["family_id"].astype(str).unique().tolist())                  # Collect and sort all unique family IDs.
    return (                                                                          # Return two dictionaries used by the dataset and checkpoint.
        {h: i for i, h in enumerate(hosts)},                                          # Map each host genus string to a numeric ID.
        {f: i for i, f in enumerate(families)},                                       # Map each family ID string to a numeric ID.
    )


def run_epoch(model, loader, device, optimizer=None) -> float:
    """Run one full train or validation epoch and return mean masked-token loss."""
    # Determine whether this is a training or validation epoch by whether an optimizer was supplied.
    train = optimizer is not None                                                     
    model.train(train)                                                                # Switch the model to train or eval mode accordingly.
    total_loss = 0.0                                                                  # Accumulate loss over all masked tokens.
    total_tokens = 0                                                                  # Count how many masked tokens contributed to the loss.

    # Iterate over mini-batches, move tensors to the chosen device, and compute loss.
    for batch in loader:                                                              # Visit each mini-batch from the DataLoader.
        batch = {k: v.to(device) for k, v in batch.items()}                           # Move every tensor in the batch dictionary to the selected device.
        logits = model(                                                               # Run the conditioned generator forward on the masked sequences.
            input_ids=batch["input_ids"],                                             # Provide masked token IDs.
            attention_mask=batch["attention_mask"],                                   # Provide the real-token versus padding mask.
            host_ids=batch["host_ids"],                                               # Provide host conditioning IDs.
            family_ids=batch["family_ids"],                                           # Provide family conditioning IDs.
        )
        loss = F.cross_entropy(                                                       # Compute masked-token cross-entropy over vocabulary logits.
            logits.view(-1, logits.shape[-1]),                                        # Flatten batch and sequence positions into one axis.
            batch["labels"].view(-1),                                                 # Flatten target labels to match the logits.
            ignore_index=-100,                                                        # Ignore unmasked positions that should not contribute to the loss.
        )

        # Only backpropagate and update weights during training epochs.
        if train:                                                                     # Check whether the current epoch is a training epoch.
            optimizer.zero_grad()                                                     # Clear old gradients from the previous optimization step.
            loss.backward()                                                           # Backpropagate the masked-token loss through the model.
            optimizer.step()                                                          # Update the model weights using AdamW.

        valid = int((batch["labels"] != -100).sum().item())                           # Count how many masked positions were supervised in this batch.
        total_loss += float(loss.item()) * max(valid, 1)                              # Weight the batch loss by the number of supervised tokens.
        total_tokens += max(valid, 1)                                                 # Update the total supervised-token count.

    return total_loss / max(total_tokens, 1)                                          # Return average loss per supervised masked token.


def main() -> None:
    # Parse command-line arguments and load the dataset.
    args = parse_args()                                                               # Parse command-line arguments.
    df = pd.read_csv(args.dataset_csv)                                                # Read the prepared Stage 07 generator dataset from disk.
    if "split" not in df.columns:                                                     # Validate that the dataset already contains train/val/test splits.
        raise ValueError("dataset_csv must contain a split column from the Stage 07 dataset builder.")

    # Build the training and validation datasets and create data loaders.
    host_to_id, family_to_id = build_maps(df)                                         # Build numeric conditioning vocabularies from the dataset table.
    train_df = df[df["split"] == "train"].reset_index(drop=True)                      # Keep the training subset only.
    val_df = df[df["split"] == "val"].reset_index(drop=True)                          # Keep the validation subset only.
    if len(val_df) == 0:                                                              # Ensure early-stopping selection has a validation set.
        raise ValueError("No validation rows found. Please regenerate the dataset with a non-empty val split.")

    train_ds = Stage07SequenceDataset(train_df, host_to_id, family_to_id, max_len=args.max_len)          # Build the training dataset object.
    val_ds = Stage07SequenceDataset(val_df, host_to_id, family_to_id, max_len=args.max_len)              # Build the validation dataset object.
    collator = MaskedConditioningCollator(mask_prob=args.mask_prob)                                      # Create the masking collator used by both loaders.
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, collate_fn=collator)   # Create the shuffled training loader.
    val_loader = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False, collate_fn=collator)      # Create the deterministic validation loader.

    # Build the model, and the optimizer.
    device = "cuda" if torch.cuda.is_available() else "cpu"                           # Prefer GPU when available, otherwise fall back to CPU.
    model = ConditionalMaskedRBPGenerator(                                            # Build the conditional transformer generator.
        vocab_size=len(VOCAB),                                                        # Set the vocabulary size from the shared token bank.
        n_hosts=len(host_to_id),                                                      # Set the number of host condition labels.
        n_families=len(family_to_id),                                                 # Set the number of family condition labels.
    ).to(device)                                                                      # Move the generator to the selected device.
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr)                     # Create the AdamW optimizer over all model parameters.

    out_dir = Path(args.out_dir)                                                      # Convert the output directory string into a Path object.
    out_dir.mkdir(parents=True, exist_ok=True)                                        # Create the output directory if it does not exist.

    best_val = float("inf")                                                           # Track the best validation loss seen so far.
    history = []                                                                      # Store per-epoch metrics for plotting or later inspection.

    # Run the train/validation loop for several epochs, and keep only the best checkpoint.
    for epoch in range(1, args.epochs + 1):                                           # Iterate over the requested number of epochs.
        tr = run_epoch(model, train_loader, device, optimizer=optimizer)              # Run one full training epoch.
        va = run_epoch(model, val_loader, device, optimizer=None)                     # Run one full validation epoch without gradient updates.
        history.append({"epoch": epoch, "train_mask_loss": tr, "val_mask_loss": va})  # Record the epoch metrics in memory.
        print(f"epoch={epoch:02d} train_mask_loss={tr:.5f} val_mask_loss={va:.5f}")   # Print a concise progress line.

        # Save a checkpoint only when validation improves.
        if va < best_val:                                                             
            best_val = va                                                             # Update the best validation loss.
            torch.save(                                                               # Serialize the best checkpoint to disk.
                {
                    "model": model.state_dict(),                                      # Store the learned model parameters.
                    "host_to_id": host_to_id,                                         # Store the host vocabulary used during training.
                    "family_to_id": family_to_id,                                     # Store the family vocabulary used during training.
                    "history": history,                                               # Store the training history for provenance.
                },
                out_dir / "best_generator.pt",                                        # Write the checkpoint under a fixed filename.
            )
    # Save the training history to disk, and print the final checkpoint path.
    (out_dir / "training_history.json").write_text(json.dumps(history, indent=2), encoding="utf-8")  # Save the epoch history as JSON.
    print(f"Wrote: {out_dir / 'best_generator.pt'}")                                  # Print the final checkpoint path for quick confirmation.


if __name__ == "__main__":                                                            # Standard Python entrypoint guard.
    main()                                                                            # Execute the training CLI.
