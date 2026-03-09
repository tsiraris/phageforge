from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import classification_report, f1_score, accuracy_score
from sklearn.preprocessing import LabelEncoder
from tqdm import tqdm


def set_seed(seed: int):
    """ Set Python / NumPy / PyTorch random seeds for reproducibility. """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def split_by_virus(df: pd.DataFrame, train_frac=0.8, val_frac=0.1, seed=42):
    """
    Group split by virus_accession per class to reduce leakage:
    - for each host_genus: split its unique viruses into train/val/test
    - assign all proteins of a virus to the same split
    Returns a split array where each element is either "train", "val" or "test".
    """
    rng = np.random.default_rng(seed)
    split = np.full(len(df), "train", dtype=object)     # Create an array of length equal to the number of rows in df, filled with the string "train".

    for genus, g in df.groupby("host_genus"):               # For each unique value of host_genus in df, create a new group g
        viruses = g["virus_accession"].unique().tolist()    # Get a list of unique values of virus_accession in g
        rng.shuffle(viruses)                                # Shuffle the list

        n = len(viruses)                                    # Get the length of the list
        n_train = int(round(train_frac * n))                # Round the train_frac to the nearest integer
        n_val = int(round(val_frac * n))                    # Round the val_frac to the nearest integer
        train_set = set(viruses[:n_train])                  # Create a set of the first n_train elements of viruses
        val_set = set(viruses[n_train : n_train + n_val])   # Create a set of the elements from n_train to n_train + n_val
        test_set = set(viruses[n_train + n_val :])          # Create a set of the elements from n_train + n_val to the end

        idx = g.index                                       # Returns the index labels that the group g contains in the original df
        v = df.loc[idx, "virus_accession"]                  # Retrieves only the "virus_accession" column values for the rows of this current group g
        split[idx[v.isin(train_set)]] = "train"             # Find which of those viruses are assigned to train_set and updates the split array to "train" for the corresponding indices
        split[idx[v.isin(val_set)]] = "val"                 # Assign in the split array accordingly the values of v that are in val_set to "val"
        split[idx[v.isin(test_set)]] = "test"               # Assign the values of v that are in test_set to "test"

    return split            # Return a split array where each element is either "train", "val" or "test"


class MLP(nn.Module):
    """ MLP classifier: Takes as inputs the embeddings, comprises 3 hidden layers with dropout and returns the logits. """
    def __init__(self, d_in: int, n_classes: int):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(d_in, 512),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, n_classes),
        )

    def forward(self, x):
        return self.net(x)


def save_split_tables(df: pd.DataFrame, out_dir: Path):
    """ Save virus-level and protein-level split summaries as CSV files. """
    virus_counts = (
        df.groupby(["host_genus", "split"])["virus_accession"]
        .nunique()
        .unstack(fill_value=0)
        .reset_index()
    )
    virus_counts.to_csv(out_dir / "split_counts_viruses.csv", index=False)

    protein_counts = (
        df.groupby(["host_genus", "split"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    protein_counts.to_csv(out_dir / "split_counts_proteins.csv", index=False)


def main():
    # Create an ArgumentParser object and parse the command line arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb", type=str, default="data/processed/esm2_embeddings.pt")
    ap.add_argument("--idx", type=str, default="data/processed/esm2_embeddings_index.csv")
    ap.add_argument("--out_dir", type=str, default="data/processed/phi_mlp")
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--batch_size", type=int, default=64)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--dataset_name", type=str, default="broad")
    ap.add_argument("--use_class_weights", action="store_true")
    args = ap.parse_args()

    set_seed(args.seed)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Load the embeddings and initial dataframe
    emb = torch.load(args.emb, map_location="cpu")
    df = pd.read_csv(args.idx)

    # Labels
    le = LabelEncoder()                                 # Create a LabelEncoder object to convert string labels into integers
    y = le.fit_transform(df["host_genus"].astype(str))  # Scans host_genus column, finds all unique strings, sorts them alphabetically, assigns an integer to each and replaces every string in the dataframe with its corresponding integer.
    n_classes = len(le.classes_)                        # Get the number of classes; le.classes_ is an array of the unique strings the encoder found.

    # Split (by virus)
    df["split"] = split_by_virus(df, seed=args.seed)    # Add a new column "split" to the dataframe 
    
    # Print Viruses per genus per split 
    print("\nViruses per host genus per split:")
    print(
        df.groupby(["host_genus", "split"])["virus_accession"]
        .nunique()
        .unstack(fill_value=0)
    )

    # Print Proteins per host genus per split
    print("\nProteins per host genus per split:")
    print(
        df.groupby(["host_genus", "split"])
        .size()
        .unstack(fill_value=0)
    )

    # Print total proteins per host genus
    print("\nTotal proteins per host genus:")
    print(df["host_genus"].value_counts())

    save_split_tables(df, out_dir)
    
    split = df["split"].values                          # Get the values of the "split" column as a numpy array

    X = emb.float()                                     # Convert the embeddings to float
    y_t = torch.tensor(y, dtype=torch.long)             # Convert the labels to torch tensor

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[INFO] device={device} | classes={list(le.classes_)}")

    # Class target weights for the weighted loss function (to handle class imbalance)
    counts = np.bincount(y)                                                 # Get the counts of each class (how many times each integer appears in y); Shape: [1,n_classes]
    weights = (counts.sum() / np.maximum(counts, 1)).astype(np.float32)     # Compute the raw class weights for each class: Total/ClassCount
    weights = weights / weights.mean()                                      # Normalize the class weights to avoid gradient explosion
    class_w = torch.tensor(weights, dtype=torch.float32).to(device)         # Convert the weights into a 32-bit float PyTorch tensor and move it to the GPU

    model = MLP(d_in=X.shape[1], n_classes=n_classes).to(device)            
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr)                 # AdamW optimizer

    if args.use_class_weights:
        loss_fn = nn.CrossEntropyLoss(weight=class_w)                       # Weighted cross entropy loss
    else:
        loss_fn = nn.CrossEntropyLoss()                                     # Standard cross entropy loss

    history_rows = []

    def run_epoch(split_name: str, train: bool):
        """ Run a single epoch of training or evaluation. """
        idxs = np.where(split == split_name)[0]     # Get the indices of the "split" column that are equal to split_name (train, val, test)
        if len(idxs) == 0:
            return None

        if train:           # If train is True
            model.train()   # Set the model to training mode
        else:               # If train is False
            model.eval()    # Set the model to evaluation mode

        total_loss = 0.0
        preds_all, y_all = [], [] 

        for i in range(0, len(idxs), args.batch_size):  # Loop over the indices
            b = idxs[i : i + args.batch_size]           # Get a batch of indices
            xb = X[b].to(device)                        # Move the batch to the device
            yb = y_t[b].to(device)                      # Move the labels to the device

            with torch.set_grad_enabled(train):         # Enable/Disable gradient computation   
                logits = model(xb)                      # Get the logits
                loss = loss_fn(logits, yb)              # Compute the loss

                if train:                               # If train is True
                    opt.zero_grad()                     # Set previous batch gradients to zero
                    loss.backward()                     # Compute loss gradients with respect to model parameters
                    opt.step()                          # Update the parameters

            total_loss += loss.item() * len(b)          # Total sum loss for the batch
            preds_all.append(logits.argmax(dim=1).detach().cpu().numpy())   # Return the index of the highest logit (predicted class) for each sample, convert it to a numpy array and add it to the list
            y_all.append(yb.detach().cpu().numpy())                         # Add the ground truth labels as a numpy array to the list

        preds_all = np.concatenate(preds_all)           # Concatenate all the predictions
        y_all = np.concatenate(y_all)                   # Concatenate all the labels
        avg_loss = total_loss / len(idxs)               # Average loss for the epoch
        macro_f1 = f1_score(y_all, preds_all, average="macro")  # Macro F1 score for the epoch
        weighted_f1 = f1_score(y_all, preds_all, average="weighted")  # Weighted F1 score for the epoch
        acc = accuracy_score(y_all, preds_all)                  # Accuracy for the epoch
        return avg_loss, macro_f1, weighted_f1, acc

    best_val = -1.0
    best_path = out_dir / "best_model.pt"

    # Train and eval for args.epochs
    for epoch in range(1, args.epochs + 1):
        tr = run_epoch("train", train=True) # Run train epoch: tr = [loss, f1, acc]
        va = run_epoch("val", train=False)  # Run val epoch: va = [loss, f1, acc]   

        if va is None:
            print("[WARN] No val split created; check dataset.")
            break

        print(f"Epoch {epoch:02d} | train loss={tr[0]:.4f} f1={tr[1]:.3f} acc={tr[3]:.3f} "
              f"| val loss={va[0]:.4f} f1={va[1]:.3f} acc={va[3]:.3f}")

        history_rows.append(
            {
                "epoch": epoch,
                "train_loss": tr[0],
                "train_macro_f1": tr[1],
                "train_weighted_f1": tr[2],
                "train_acc": tr[3],
                "val_loss": va[0],
                "val_macro_f1": va[1],
                "val_weighted_f1": va[2],
                "val_acc": va[3],
            }
        )

        # Save the best model as a state dictionary and the label encoder as a list
        if va[1] > best_val:
            best_val = va[1]
            torch.save(
                {"model": model.state_dict(), "label_encoder": le.classes_.tolist()},
                best_path           # best_path = out_dir / "best_model.pt"
            )

    pd.DataFrame(history_rows).to_csv(out_dir / "history.csv", index=False)

    # Testing stage using the best model
    ckpt = torch.load(best_path, map_location=device)   # Load the best model
    model.load_state_dict(ckpt["model"])                # Load the model state dictionary, i.e. best model parameters
    model.eval()                                        # Set the model to evaluation mode

    test_idxs = np.where(split == "test")[0]            # Get the indices of the split array that are equal to "test"
    xb = X[test_idxs].to(device)                        # Move the test data to the device
    yb = y_t[test_idxs].cpu().numpy()                   # Move the test labels to the cpu and convert them to a numpy array
    with torch.no_grad():
        preds = model(xb).argmax(dim=1).cpu().numpy()   # Get the predictions

    # Save the test report
    report = classification_report(yb, preds, target_names=le.classes_, digits=3, zero_division=0)   # Get the sklearn classification report 
    (out_dir / "test_report.txt").write_text(report)                                # Write the report to a text file

    # Calculate test metrics: macro F1, weighted F1, and accuracy
    test_macro_f1 = f1_score(yb, preds, average="macro", zero_division=0)        # Computes the metric independently for each class, then average equally (all classes have the same weight/importance). 
    test_weighted_f1 = f1_score(yb, preds, average="weighted", zero_division=0)  # Averages per-class scores, but weights each class by its support (number of protein sequences per class in the test set)
    test_acc = accuracy_score(yb, preds)

    metrics = {
        "dataset_name": args.dataset_name,
        "model_name": "mlp",
        "seed": args.seed,
        "n_total": int(len(df)),
        "n_train": int((df["split"] == "train").sum()),
        "n_val": int((df["split"] == "val").sum()),
        "n_test": int((df["split"] == "test").sum()),
        "accuracy": float(test_acc),
        "macro_f1": float(test_macro_f1),
        "weighted_f1": float(test_weighted_f1),
        "classes": le.classes_.tolist(),
        "use_class_weights": bool(args.use_class_weights),
        "emb_path": args.emb,
        "idx_path": args.idx,
    }
    (out_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))

    predictions_df = pd.DataFrame(
        {
            "true_label": [le.classes_[i] for i in yb],
            "pred_label": [le.classes_[i] for i in preds],
        }
    )
    predictions_df.to_csv(out_dir / "test_predictions.csv", index=False)

    print("\n=== TEST REPORT ===\n")
    print(report)
    print(f"✅ Saved: {best_path}")
    print(f"✅ Saved: {out_dir / 'test_report.txt'}")
    print(f"✅ Saved: {out_dir / 'metrics.json'}")
    print(f"✅ Saved: {out_dir / 'history.csv'}")
    print(f"✅ Saved: {out_dir / 'split_counts_viruses.csv'}")
    print(f"✅ Saved: {out_dir / 'split_counts_proteins.csv'}")
    print(f"✅ Saved: {out_dir / 'test_predictions.csv'}")


if __name__ == "__main__":
    main()