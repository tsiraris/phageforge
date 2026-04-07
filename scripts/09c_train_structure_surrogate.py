#!/usr/bin/env python
"""Stage 09c: Train or configure a Stage 09 structural surrogate.

The surrogate is designed to be practical rather than fragile:
- if enough structurally labeled candidates exist, it fits lightweight sklearn models;
- if supervision is too sparse, it writes a rule-based bundle that Stage 09 can still use.

The output is always a single joblib bundle that later scripts can load uniformly.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score

from phageforge.stage09_utils import read_json



def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the Stage 09 structural-surrogate trainer."""
    ap = argparse.ArgumentParser(description="Train or configure the Stage 09 structural surrogate.")                                                                 # Initialize parser with program description
    ap.add_argument("--dataset_csv", type=str, required=True, help="Stage 09 structural-surrogate dataset CSV produced by 09b_build_structure_surrogate_dataset.py.") # Add argument for input dataset CSV path
    ap.add_argument("--summary_json", type=str, required=True, help="JSON summary written by 09b_build_structure_surrogate_dataset.py.")                              # Add argument for input summary JSON path
    ap.add_argument("--out_model", type=str, required=True, help="Where to write the surrogate bundle joblib file.")                                                  # Add argument for output model file path
    ap.add_argument("--min_rows_for_fit", type=int, default=20, help="Minimum labeled rows required before fitting learned surrogate models.")                          # Add argument for minimum row threshold
    ap.add_argument("--seed", type=int, default=42, help="Random seed for sklearn models.")                                                                           # Add argument to set random seed for reproducibility
    return ap.parse_args()                                                                                                                                            # Parse command-line inputs and return them



def main() -> None:
    # Read the surrogate dataset and the coverage summary that says whether there are enough structural labels to fit a real model.
    args = parse_args()                                         # Execute argument parser and store parsed arguments
    dataset_df = pd.read_csv(args.dataset_csv)                  # Read the specified CSV file into a pandas DataFrame
    summary = read_json(args.summary_json)                      # Read the specified JSON summary into a dictionary

    # Keep only the numeric feature columns declared by Stage 09b and drop rows missing all structural targets.
    feature_columns = list(summary.get("feature_columns", []))                                                                                                          # Extract the list of feature column names from the summary
    if not feature_columns:                                                                                                                                             # Check if the feature_columns list is empty
        raise ValueError("The Stage 09 surrogate summary JSON does not contain any feature columns.")                                                                   # Throw an error if no feature columns are found
    X = dataset_df.reindex(columns=feature_columns, fill_value=0.0).to_numpy(dtype=np.float32)                                                                          # Extract features, fill missing with 0.0, and convert to numpy array
    y_pass = dataset_df["stage08_pass"].astype(int).to_numpy(dtype=int) if "stage08_pass" in dataset_df.columns else np.zeros(len(dataset_df), dtype=int)               # Extract target labels as ints, defaulting to zeros if missing
    y_plddt = dataset_df["esmfold_mean_plddt"].to_numpy(dtype=np.float32) if "esmfold_mean_plddt" in dataset_df.columns else np.full(len(dataset_df), np.nan, dtype=np.float32) # Extract pLDDT scores, defaulting to NaNs if missing
    y_rmsd = dataset_df["rmsd_to_selected_seed"].to_numpy(dtype=np.float32) if "rmsd_to_selected_seed" in dataset_df.columns else np.full(len(dataset_df), np.nan, dtype=np.float32) # Extract RMSD scores, defaulting to NaNs if missing

    # Fall back to a rule-only bundle if the labeled dataset is too small or entirely one-sided.
    enough_rows = int(len(dataset_df)) >= int(args.min_rows_for_fit) # Check if the dataset size meets the minimum row threshold
    enough_class_variation = len(np.unique(y_pass)) > 1              # Check if both positive and negative labels exist in the dataset
    use_rule_bundle = not (enough_rows and enough_class_variation)   # Determine if we must fallback to rule-based logic

    bundle = {                                                       # Initialize the surrogate bundle dictionary
        "feature_columns": feature_columns,                          # Store the feature column names in the bundle
        "mode": "rule_based" if use_rule_bundle else "fitted",       # Set the bundle mode depending on data viability
        "pass_model": None,                                          # Initialize pass_model key as None
        "plddt_model": None,                                         # Initialize plddt_model key as None
        "rmsd_model": None,                                          # Initialize rmsd_model key as None
        "training_summary": {                                        # Initialize nested dictionary for training statistics
            "n_rows": int(len(dataset_df)),                          # Record total number of rows evaluated
            "n_positive": int(y_pass.sum()),                         # Record the count of positive instances (pass)
            "n_negative": int(len(y_pass) - y_pass.sum()),           # Record the count of negative instances (fail)
        },                                                           # Close the training_summary sub-dictionary
    }                                                                # Close the bundle dictionary

    # Fit lightweight tree-based models only when the Stage 08 label set is large enough to justify it.
    if not use_rule_bundle:                                                                                                         # Execute if data is sufficient for model fitting
        pass_model = RandomForestClassifier(n_estimators=300, max_depth=5, min_samples_leaf=2, random_state=args.seed)              # Initialize Random Forest classifier for pass/fail
        plddt_model = RandomForestRegressor(n_estimators=300, max_depth=5, min_samples_leaf=2, random_state=args.seed)              # Initialize Random Forest regressor for pLDDT
        rmsd_model = RandomForestRegressor(n_estimators=300, max_depth=5, min_samples_leaf=2, random_state=args.seed)               # Initialize Random Forest regressor for RMSD

        pass_model.fit(X, y_pass)                                                                                                   # Train the classification model on all features and labels
        valid_plddt = np.isfinite(y_plddt)                                                                                          # Create a boolean mask identifying non-NaN pLDDT scores
        valid_rmsd = np.isfinite(y_rmsd)                                                                                            # Create a boolean mask identifying non-NaN RMSD scores
        
        if valid_plddt.any():                                                                                                       # Check if there are any valid pLDDT values to train on
            plddt_model.fit(X[valid_plddt], y_plddt[valid_plddt])                                                                   # Train the pLDDT regressor on valid rows only
        else:                                                                                                                       # Handle case where no valid pLDDT data exists
            plddt_model = None                                                                                                      # Reset plddt_model to None
            
        if valid_rmsd.any():                                                                                                        # Check if there are any valid RMSD values to train on
            rmsd_model.fit(X[valid_rmsd], y_rmsd[valid_rmsd])                                                                       # Train the RMSD regressor on valid rows only
        else:                                                                                                                       # Handle case where no valid RMSD data exists
            rmsd_model = None                                                                                                       # Reset rmsd_model to None

        pred_pass = pass_model.predict_proba(X)[:, 1]                                                                               # Predict the probability of the positive class for the dataset
        auc = float(roc_auc_score(y_pass, pred_pass)) if len(np.unique(y_pass)) > 1 else float("nan")                               # Calculate the AUC metric, or assign NaN if impossible
        bundle.update(                                                                                                              # Update the bundle dictionary with fitted objects
            {                                                                                                                       # Open the update payload dictionary
                "mode": "fitted",                                                                                                   # Explicitly set the mode string to 'fitted'
                "pass_model": pass_model,                                                                                           # Add the trained pass_model object to the bundle
                "plddt_model": plddt_model,                                                                                         # Add the trained plddt_model object to the bundle
                "rmsd_model": rmsd_model,                                                                                           # Add the trained rmsd_model object to the bundle
                "training_summary": {                                                                                               # Update the nested training_summary dictionary
                    **bundle["training_summary"],                                                                                   # Unpack existing key-value pairs from the summary
                    "train_auc": auc,                                                                                               # Append the calculated AUC score
                },                                                                                                                  # Close the updated training_summary dictionary
            }                                                                                                                       # Close the update payload dictionary
        )                                                                                                                           # Complete the update method call

    # Persist the bundle in one joblib file so later Stage 09 scripts can load it with a single path argument.
    out_path = Path(args.out_model)                                      # Create a Path object for the target save destination
    out_path.parent.mkdir(parents=True, exist_ok=True)                   # Ensure the parent directory tree exists
    joblib.dump(bundle, out_path)                                        # Serialize and save the final bundle to disk
    print(f"Wrote: {out_path}")                                          # Print confirmation of the saved file path
    print(f"mode: {bundle['mode']}")                                     # Print the mode applied to the bundle
    print(f"training_summary: {bundle['training_summary']}")             # Print the final summary statistics


if __name__ == "__main__":                                               # Check if this script is being executed directly
    main()                                                               # Invoke the main execution function