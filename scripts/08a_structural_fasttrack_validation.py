#!/usr/bin/env python
"""
Robust Stage 08 structural validation for top Stage 07 candidates.

This script evaluates AI-engineered protein sequences by predicting their 3D 
structures using the ESMFold model. It compares these predictions against 
the original "seed" scaffold, calculating structural confidence (pLDDT), 
structural drift (RMSD), and various sequence health metrics (e.g., entropy, 
hydrophobicity). Ultimately, it scores and filters the candidates to identify 
which designs remain structurally plausible and biologically viable for lab testing.
"""

from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path
from typing import Iterable
import numpy as np
import pandas as pd


_MUT_RE = re.compile(r"(\d+):([A-Z])→([A-Z])")                                # Compile a regex to extract mutation positions from strings


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()                                        # Initialize the command-line argument parser
    parser.add_argument("--validated_csv", type=Path, required=True)          # Define the input path for the validated CSV file
    parser.add_argument("--ranked_csv", type=Path, required=True)             # Define the input path for the ranked CSV file
    parser.add_argument("--context_json", type=Path, required=True)           # Define the input path for the context JSON file
    parser.add_argument("--out_dir", type=Path, required=True)                # Define the output directory path for results
    parser.add_argument("--top_k", type=int, default=3)                       # Define the number of top candidates to process
    parser.add_argument("--device", type=str, default="auto")                 # Define the compute device to use (e.g., cuda, cpu)
    parser.add_argument("--chunk_size", type=int, default=128)                # Define the chunk size for ESMFold inference
    parser.add_argument("--num_recycles", type=int, default=1)                # Define the number of recycles for the folding model
    parser.add_argument("--resume", action="store_true", help="Reuse...")     # Define a flag to reuse existing PDBs to save time
    parser.add_argument("--max_candidates", type=int, default=None)           # Define a limit for the maximum number of candidates
    return parser.parse_args()                                                # Parse and return the supplied command-line arguments

def parse_mutation_positions(text: str) -> list[int]:
    """
    Extracts the raw integer positions of mutations from a formatted mutation string.

    It uses a pre-compiled regular expression (`_MUT_RE`) to find patterns matching 
    the specific "Position:OldAminoAcid→NewAminoAcid" format. It captures only the 
    numerical position group and ignores the amino acid characters.

    Example:
        parse_mutation_positions("12:A→T; 45:G→C")
        -> [12, 45]
    """
    if not isinstance(text, str) or not text.strip():                         # Check if the input is a valid, non-empty string
        return []                                                             # Return an empty list if the input is invalid
    return [int(m.group(1)) for m in _MUT_RE.finditer(text)]                  # Extract and return integer positions using the regex


def hydrophobic_fraction(seq: str) -> float:
    """
    Measures the proportion of the sequence made up of hydrophobic amino acids.

    It compares every character in the sequence against a hardcoded set of hydrophobic 
    residues (A, I, L, M, F, W, V, Y) and divides the count by the total sequence length. 
    High hydrophobicity can sometimes indicate sequences that will aggregate or misfold.

    Example:
        hydrophobic_fraction("AILMQ")
        -> 0.60 (Because A, I, L are hydrophobic; 3 out of 5)
    """
    hydrophobic = set("AILMFWVY")                                             # Define a set of characters representing hydrophobic amino acids
    return sum(aa in hydrophobic for aa in seq) / max(len(seq), 1)            # Calculate and return the ratio of hydrophobic residues


def charged_fraction(seq: str) -> float:
    """
    Measures the proportion of the sequence made up of electrically charged amino acids.

    It sums the occurrences of positive (K, R, H) and negative (D, E) residues. 
    This metric helps evaluate the solubility and surface characteristics of the protein.

    Example:
        charged_fraction("KRDEHQQ")
        -> 0.714 (Because K, R, D, E, H are charged; 5 out of 7)
    """
    charged = set("KRDEH")                                                    # Define a set of characters representing charged amino acids
    return sum(aa in charged for aa in seq) / max(len(seq), 1)                # Calculate and return the ratio of charged residues


def gly_fraction(seq: str) -> float:
    """
    Measures the exact proportion of Glycine (G) in the sequence.

    Glycine lacks a side chain, making it highly flexible. Generative models sometimes 
    "cheat" by using Glycine to force unstructured loops that artificially satisfy 
    mathematical constraints but fail to form rigid 3D binding structures in reality.

    Example:
        gly_fraction("AGGGTT")
        -> 0.50 (Because G appears 3 times out of 6)
    """
    return seq.count("G") / max(len(seq), 1)                                  # Calculate and return the fraction of Glycine (G) residues


def aromatic_fraction(seq: str) -> float:
    """
    Measures the proportion of bulky, aromatic amino acids in the sequence.

    It counts residues with aromatic rings (F, W, Y, H). Aromatic residues are 
    frequently found in protein-protein interfaces and binding pockets.

    Example:
        aromatic_fraction("FWYHA")
        -> 0.80 (4 out of 5 residues are aromatic)
    """
    aromatic = set("FWYH")                                                    # Define a set of characters representing aromatic amino acids
    return sum(aa in aromatic for aa in seq) / max(len(seq), 1)               # Calculate and return the ratio of aromatic residues


def shannon_entropy(seq: str) -> float:
    """
    Calculates the raw Shannon entropy of the amino acid composition.

    It computes the probability (frequency) of each unique amino acid present in the 
    sequence and applies the standard entropy formula: -sum(p * log2(p)). This flags 
    sequences that over-rely on a tiny subset of amino acids.

    Example:
        shannon_entropy("AAAAA") -> 0.0 (Zero uncertainty/diversity)
        shannon_entropy("ACDEF") -> ~2.32 (High diversity for 5 residues)
    """
    if not seq:                                                               # Check if the sequence is empty
        return 0.0                                                            # Return zero entropy for an empty sequence
    vals = np.array([seq.count(aa) for aa in sorted(set(seq))], dtype=float)  # Count occurrences of each unique amino acid
    probs = vals / vals.sum()                                                 # Convert the counts into probabilities
    return float(-(probs * np.log2(probs + 1e-12)).sum())                     # Calculate and return the Shannon entropy formula


def max_homopolymer_run(seq: str) -> int:
    """
    Finds the length of the longest contiguous stretch of the exact same amino acid.

    It slides through the sequence comparing adjacent characters. If they match, 
    the run counter increases; if they differ, it resets. This catches "stuttering" 
    generative artifacts.

    Example:
        max_homopolymer_run("MKGGGGGTY")
        -> 5 (Because 'G' repeats 5 times consecutively)
    """
    if not seq:                                                               # Check if the sequence is empty
        return 0                                                              # Return zero length for an empty sequence
    best = run = 1                                                            # Initialize the best and current run lengths to 1
    for a, b in zip(seq, seq[1:]):                                            # Iterate through adjacent pairs of characters
        run = run + 1 if a == b else 1                                        # Increment run if characters match, else reset to 1
        best = max(best, run)                                                 # Update the best run length found so far
    return best                                                               # Return the maximum homopolymer run length


def mutation_span(positions: Iterable[int]) -> int:
    """
    Calculates the total linear distance between the first and last mutated residues.

    It deduplicates and sorts the integer positions, then subtracts the minimum 
    from the maximum. This metric guarantees mutations aren't impossibly isolated 
    to a single point.

    Example:
        mutation_span([10, 15, 20])
        -> 10 (Because 20 - 10 = 10)
    """
    pos = sorted(set(int(x) for x in positions))                              # Deduplicate, convert to int, and sort the positions
    return 0 if not pos else pos[-1] - pos[0]                                 # Return the distance between first and last mutations


def choose_device(torch, wanted: str) -> str:
    """
    Dynamically assigns the best available hardware accelerator for neural network inference.

    It checks if the user forced a specific device. If set to "auto", it checks for 
    NVIDIA CUDA GPUs, then Apple Silicon MPS, and finally falls back to CPU to ensure 
    cross-platform compatibility without crashing.

    Example:
        choose_device(torch, "auto")
        -> "cuda" (If running on an AWS A10G or Google Colab T4)
    """
    if wanted != "auto":                                                      # Check if a specific device was requested
        return wanted                                                         # Return the explicitly requested device
    if torch.cuda.is_available():                                             # Check if a CUDA GPU is available
        return "cuda"                                                         # Return 'cuda' for NVIDIA GPU acceleration
    if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available(): # Check if Apple Silicon MPS is available
        return "mps"                                                          # Return 'mps' for Apple Silicon acceleration
    return "cpu"                                                              # Fallback and return 'cpu' if no accelerators exist


def import_esmfold():
    """
    Locally scopes the importation of heavy PyTorch and HuggingFace libraries.

    By deferring these imports into a function rather than putting them at the top 
    of the file, the script prevents massive memory consumption and slow startup times 
    until the folding models are strictly required.
    """
    import torch                                                              # Import the PyTorch library internally
    from transformers import AutoTokenizer, EsmForProteinFolding              # Import ESMFold classes from HuggingFace transformers
    return torch, AutoTokenizer, EsmForProteinFolding                         # Return the imported modules and classes


def load_esmfold(device: str):
    """
    Instantiates the ESMFold structural prediction model and prepares it for inference.

    It downloads or loads the `facebook/esmfold_v1` checkpoint. If a CUDA GPU is 
    detected, it aggressively optimizes the model by casting the trunk to Half-Precision 
    (FP16), massively speeding up inference and reducing VRAM usage.

    Example:
        torch, tokenizer, model = load_esmfold("cuda")
        -> (Loads model to GPU in FP16 mode)
    """
    torch, AutoTokenizer, EsmForProteinFolding = import_esmfold()             # Load required dependencies via the helper function
    model_name = "facebook/esmfold_v1"                                        # Specify the HuggingFace model identifier
    tokenizer = AutoTokenizer.from_pretrained(model_name)                     # Instantiate the tokenizer for the specified model
    model = EsmForProteinFolding.from_pretrained(model_name, low_cpu_mem_usage=True) # Instantiate the ESMFold model efficiently
    model.eval()                                                              # Set the model to evaluation mode (disable dropout)
    model.to(device)                                                          # Move the model weights to the selected hardware device
    # If the device is a CUDA GPU
    if device == "cuda":                                                      
    # Convert the trunk to FP16 to save memory and speed up
        model.esm = model.esm.half()                                          # Convert the trunk to FP16 to save memory and speed up
    return torch, tokenizer, model                                            # Return the torch module, tokenizer, and loaded model


def _scalar_tensor_to_float(x) -> float:
    """
    Safely unwraps 0-dimensional or 1-dimensional PyTorch tensors back into native Python floats.

    This prevents downstream Pandas operations from crashing due to unexpected 
    tensor shapes returned by the ESMFold `ptm` (Predicted Template Modeling) output.
    """
    arr = x.detach().cpu().numpy()                                            # Detach from graph, move to CPU, and convert to numpy
    # If the numpy array is a 0-dimensional scalar
    if np.ndim(arr) == 0:                                                     
        return float(arr.item())                                              # Extract and return the standard Python float value
    return float(np.asarray(arr).reshape(-1)[0])                              # Flatten array and return the first element as float


def _normalize_b_factors(plddt: np.ndarray, n_atoms: int = 37) -> np.ndarray:
    """
    Formats the 1D/2D pLDDT confidence array into the rigid 3D shape required by the PDB format.

    ESMFold generates a single confidence score per residue. However, standard `.pdb` 
    files require a B-factor (or temperature factor) score for the residue as a whole: 
    one value for every single atom of the residue (up to 37 atoms per residue). 
    This function intelligently repeats the residue's score across all its constituent atoms.
    
    For example:
        For a protein with 2 residues, and ESMFold giving them confidence scores of [90.0, 70.0].
        The PDB format requires a score for all 37 atoms per residue.
        _normalize_b_factors transforms that (2,) array into a (2, 37) matrix where the first row is thirty-seven 90.0s and the second row is thirty-seven 70.0s.
    """
    plddt = np.asarray(plddt)                                                 # Ensure the input pLDDT is a numpy array
    if plddt.ndim == 1:                                                       # Check if the pLDDT array is 1-dimensional
        return np.repeat(plddt[:, None], n_atoms, axis=1)                     # Expand dims and repeat for all atoms in the residue
    if plddt.ndim == 2:                                                       # Check if the pLDDT array is 2-dimensional
        if plddt.shape[1] == 1:                                               # Check if there is only 1 value per residue
            return np.repeat(plddt, n_atoms, axis=1)                          # Repeat that single value for all atoms
        if plddt.shape[1] == n_atoms:                                         # Check if there are already values for all atoms
            return plddt                                                      # Return the array unmodified
        return np.repeat(plddt.mean(axis=1, keepdims=True), n_atoms, axis=1)  # Average existing values and repeat for all atoms
    flat = plddt.reshape(plddt.shape[0], -1).mean(axis=1)                     # For higher dims, flatten and compute the mean per residue
    return np.repeat(flat[:, None], n_atoms, axis=1)                          # Expand dims and repeat the mean for all atoms


def tensor_to_pdb_lines(outputs) -> list[str]:
    """
    Translates raw 3D tensor spatial coordinates and masks generated by ESMFold 
    into the standardized, human-readable Protein Data Bank (.pdb) text format (strings).

    It uses atom14_to_atom37 to convert the AI’s compressed coordinate system into 
    the high-resolution system required by the structural software.
    Then it loops over the batch dimension, constructs a `Protein` data object using the 
    `openfold_utils` library, maps the pLDDT scores into its B-factor column, and 
    serializes the atomic coordinates.

    Example:
        The AI outputs a coordinate tensor for an Atom: [12.231, 45.321, 1.234].
        This function assembles the text line:
        ATOM      1  N   MET A   1      12.231  45.321   1.234  1.00 89.23

        1.	ATOM: Record type. Declares this line contains atomic coordinates.
        2.	1: Atom Serial Number. This is the 1st atom in the file.
        3.	N: Atom Name. This is the Nitrogen atom of the amino acid's backbone. (Other generated lines will have CA for Alpha-Carbon, C for Carboxyl Carbon, O for Oxygen, etc.)
        4.	MET: Residue Name. The amino acid type (Methionine).
        5.	A: Chain Identifier. If you had a multi-protein complex, they would be chains A, B, C.
        6.	1: Residue Sequence Number. This is the 1st amino acid in the protein chain.
        7.	12.231: X orthogonal coordinate (in Ångströms).
        8.	45.321: Y orthogonal coordinate (in Ångströms).
        9.	1.234: Z orthogonal coordinate (in Ångströms).
        10.	1.00: Occupancy. In crystal structures, an atom might alternate between two spots. In AI predictions, this is always 1.00 (100% occupying this spot).
        11.	89.23: The B-Factor column. ESMFold places its pLDDT confidence score here.
    """
    from transformers.models.esm.openfold_utils import protein as protein_module   # Import openfold protein utilities
    from transformers.models.esm.openfold_utils.feats import atom14_to_atom37      # Import coordinate conversion utility

    final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)     # Convert atom14 format coordinates to atom37
    batch = outputs["aatype"].shape[0]                                             # Determine the batch size from the output shape
    pdbs: list[str] = []                                                           # Initialize an empty list to store PDB strings
    for i in range(batch):                                                         # Iterate over each item in the batch
        aa = outputs["aatype"][i].cpu().numpy()                                    # Extract amino acid types to numpy array
        pos = final_atom_positions[i].detach().cpu().numpy()                       # Extract 3D positions to numpy array
        mask = outputs["atom37_atom_exists"][i].cpu().numpy()                      # Extract atom existence masks to numpy array
        residue_index = outputs["residue_index"][i].cpu().numpy() + 1              # Extract residue indices and make them 1-based
        plddt = outputs["plddt"][i].detach().cpu().numpy()                         # Extract pLDDT confidence scores to numpy array
        b_factors = _normalize_b_factors(plddt, n_atoms=37)                        # Compute formatted B-factors from the pLDDT scores
        protein = protein_module.Protein(                                          # Initialize a Protein data structure object
            aatype=aa,                                                             # Pass the amino acid types
            atom_positions=pos,                                                    # Pass the 3D atomic coordinates
            atom_mask=mask,                                                        # Pass the atom existence mask
            residue_index=residue_index,                                           # Pass the sequence indices
            b_factors=b_factors,                                                   # Pass the computed B-factors (pLDDT)
            chain_index=np.zeros_like(residue_index),                              # Assign all residues to a single chain (chain 0)
        )                                                                          # Finish initializing Protein object
        pdbs.append(protein_module.to_pdb(protein))                                # Convert the Protein object to a PDB string and append
    return pdbs                                                                    # Return the list of generated PDB strings


def fold_sequence(seq: str, tokenizer, model, torch, device: str, chunk_size: int, num_recycles: int):
    """
    Executes the end-to-end ESMFold neural network inference to predict the 3D structure of a protein.

    It tokenizes the input sequence, dynamically sets chunking limits to prevent VRAM 
    exhaustion (= model processing long sequences in smaller "bites" to save memory), 
    and runs the forward pass. It extracts the raw PDB text, the mean overall confidence (pLDDT), 
    and the optional template score (pTM).

    Example:
        pdb_text, mean_plddt, ptm = fold_sequence("MKAAGGG...", ...)
        -> ("ATOM 1 N...", 85.4, 0.72)
    """
    # Set the chunk size to optimize memory during attention
    model.trunk.set_chunk_size(chunk_size)                                         
    tokens = tokenizer([seq], return_tensors="pt", add_special_tokens=False)["input_ids"].to(device) # Tokenize sequence and send to device
    with torch.no_grad():                                                          # Disable gradient tracking for inference
        outputs = model(tokens, num_recycles=num_recycles)                         # Run the model forward pass to fold the protein
    # Convert the output tensors into a PDB formatted string
    pdb_text = tensor_to_pdb_lines(outputs)[0]                                     
    # Extract the per-residue pLDDT array
    plddt = outputs["plddt"][0].detach().cpu().numpy()                             
    # Calculate the average pLDDT across the whole protein
    mean_plddt = float(np.asarray(plddt).mean())                                   
    # Extract pTM score if present, else assign NaN
    ptm = _scalar_tensor_to_float(outputs["ptm"]) if "ptm" in outputs else math.nan 
    return pdb_text, mean_plddt, ptm            # Return the PDB text, mean confidence, and pTM


def parse_ca_coords_and_plddt(pdb_text: str):
    """
    Extracts the geometric and thermodynamic data from a raw PDB string.

    It scans every line of the PDB text looking only for standard "ATOM" declarations. 
    It isolates only the Alpha-Carbon ("CA") atoms (central "anchors" of every amino acid) 
    to form the geometric backbone of the protein. It returns their [X, Y, Z] spatial coordinates, 
    the pLDDT confidence score, and the integer residue index for every amino acid.
    
    Example:
        If a protein has 500 amino acids, the raw PDB file will contain roughly ~4,000 atoms (lines of text).
        By filtering the text file for only the CA lines, the script strips away all the messy side-chains.
        Returns a NumPy array of shape (500, 3)—representing the exact 3D locations of the 500 amino acids that make up that fiber's spine.
    """
    coords, plddts, residues = [], [], []                                     # Initialize lists to store coordinates, pLDDTs, and indices
    seen = set()                                                              # Initialize a set to track already processed residues
    for line in pdb_text.splitlines():                                        # Iterate line by line through the PDB text
        if not line.startswith("ATOM"):                                       # Skip lines that do not describe atomic coordinates
            continue                                                          # Move to the next line
        atom = line[12:16].strip()                                            # Extract and strip the atom name from the PDB format
        if atom != "CA":                                                      # Check if the atom is NOT an alpha-carbon (CA)
            continue                                                          # Skip non-alpha-carbon atoms
        chain = line[21].strip() or "A"                                       # Extract the chain ID, default to 'A' if blank
        resid = int(line[22:26].strip())                                      # Extract and parse the residue sequence integer
        key = (chain, resid)                                                  # Create a unique tuple key for this specific residue
        if key in seen:                                                       # Check if we have already recorded this residue's CA
            continue                                                          # Skip duplicates
        seen.add(key)                                                         # Mark this chain/residue combination as processed
        x = float(line[30:38].strip())                                        # Parse the X coordinate
        y = float(line[38:46].strip())                                        # Parse the Y coordinate
        z = float(line[46:54].strip())                                        # Parse the Z coordinate
        b = float(line[60:66].strip())                                        # Parse the B-factor column (which holds pLDDT here)
        coords.append([x, y, z])                                              # Append the 3D coordinate list
        plddts.append(b)                                                      # Append the pLDDT score
        residues.append(resid)                                                # Append the residue index
    return np.asarray(coords, dtype=float), np.asarray(plddts, dtype=float), residues # Convert to arrays and return


def parse_pdb_file(path: Path):
    """
    A simple file-reading wrapper that passes the contents of a .pdb file 
    directly to the `parse_ca_coords_and_plddt` function.
    """
    return parse_ca_coords_and_plddt(path.read_text())                        # Read the file contents and pass to the parsing function


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """
    Calculates the precise Root Mean Square Deviation (RMSD) between two sets of 3D protein coordinates.

    Since simply subtracting the coordinate matrices of two folded proteins that might 
    be facing different directions is misleading, this function implements the Kabsch algorithm. 
    It: 1) truncates the arrays to equal length, 
        2) centers them on their origins (both proteins to the coordinate (0,0,0)),
        3) applies Singular Value Decomposition (SVD) to find the mathematically 
           optimal rotation matrix that perfectly superimposes the candidate onto the seed, and 
        4) calculates the RMSD—the average distance (in Ångströms) between the atoms after they have been aligned.

    Example:
        The Input State:
            Seed: A triangle at coordinates (1,1), (2,1), (1,2).
            Candidate: An identical triangle, but generated shifted to the right and rotated 90°: (10,5), (10,6), (9,5).

        Truncation (The min() step):
            If the AI added 5 extra amino acids to the candidate's tail, we chop them off the array so both arrays are exactly the same length.

        Translation (Centering):
            The script calculates the geometric center (mean) of the Seed and subtracts it from all Seed points. The Seed moves to (0,0,0).
            It does the same for the Candidate. The Candidate moves from (10,5) down to (0,0,0). They are now on top of each other, but the Candidate is still facing the wrong way.

        Covariance Matrix (C = Pc.T @ Qc):
            It multiplies the centered coordinate matrices together. This creates a matrix representing how the X,Y,Z axes of the Seed correlate with the X,Y,Z axes of the Candidate.

        SVD & Optimal Rotation:
            It feeds the Covariance matrix into Singular Value Decomposition (np.linalg.svd). SVD mathematically isolates the exact 3D rotational angles required to untwist the Candidate so it perfectly faces the Seed.

        Calculate RMSD:
            The script applies that rotation matrix to the Candidate. Now, both proteins are at (0,0,0) and facing the exact same direction.
            Now it subtracts the coordinates. Since they are identical triangles, the distance between every point is 0.
            RMSD = 0.00 Ångströms.
            If the AI had mutated the candidate and accidentally bent the top of the triangle slightly, they would not superimpose perfectly (e.g. an average RMSD of 2.34 Ångströms).
    """
    if len(P) == 0 or len(Q) == 0:                                            # Guard against empty coordinate arrays
        return math.nan                                                       # Return NaN if inputs are empty
    n = min(len(P), len(Q))                                                   # Find the length of the shorter array to avoid mismatches
    P = P[:n]                                                                 # Truncate P to the minimum length
    Q = Q[:n]                                                                 # Truncate Q to the minimum length
    # Centers them around their origins
    Pc = P - P.mean(axis=0)                                                   # Center P around its origin
    Qc = Q - Q.mean(axis=0)                                                   # Center Q around its origin
    # Compute the covariance matrix between P and Q
    C = Pc.T @ Qc                                                             
    # Perform Singular Value Decomposition (SVD) on covariance
    V, _, Wt = np.linalg.svd(C)                                               
    d = np.sign(np.linalg.det(V @ Wt))                                        # Ensure right-handed coordinate system for rotation
    # Compute the optimal rotation matrix
    U = V @ np.diag([1.0, 1.0, d]) @ Wt                                       
    # Apply the rotation to the centered P coordinates
    aligned = Pc @ U                                                          
    return float(np.sqrt(np.mean(np.sum((aligned - Qc) ** 2, axis=1))))       # Calculate and return the Root Mean Square Deviation (RMSD)


def residue_confidence_fraction(residues: list[int], plddts: np.ndarray, selected: set[int], cutoff: float = 70.0) -> float:
    """
    Evaluates the localized structural stability of the specific mutations the AI generated.

    It scans the global pLDDT array, isolates only the indices that correspond to the 
    mutated residues, and returns the proportion of those specific edits that achieved 
    a safe confidence score (>= 70.0). This prevents globally stable proteins from passing 
    if their newly generated binding tip is an unstructured, unconfident mess.

    Example:
        residue_confidence_fraction(..., cutoff=70.0)
        -> 0.85 (85% of the mutated residues are structurally confident)
    """
    idx = [i for i, resid in enumerate(residues) if resid in selected]        # Get indices in the arrays matching the selected residues
    if not idx:                                                               # Check if no valid indices were found
        return math.nan                                                       # Return NaN if there's no data to compute
    vals = plddts[idx]                                                        # Extract the pLDDT scores for the identified indices
    return float(np.mean(vals >= cutoff))                                     # Return the proportion of scores equal or above cutoff


def residue_mean_confidence(residues: list[int], plddts: np.ndarray, selected: set[int]) -> float:
    """
    Calculates the strict average pLDDT confidence score localized solely to the 
    residues chosen for mutation by the generative model.
    """
    idx = [i for i, resid in enumerate(residues) if resid in selected]        # Get array indices for the selected residue numbers
    if not idx:                                                               # Check if no indices were matched
        return math.nan                                                       # Return NaN for empty selections
    return float(np.mean(plddts[idx]))                                        # Calculate and return the mean pLDDT for those indices


def stage08_decision_reason(row: pd.Series) -> str:
    """
    Acts as the final automated QA engine, evaluating a candidate's structural 
    and biochemical metrics against six strict physical safety rules.

    It checks the row data for violations (e.g., low global pLDDT, high RMSD drift, 
    low mutation site confidence, excess glycine). If any rule is broken, it appends 
    a specific error tag. It returns "pass" if the protein is flawless, or a 
    semicolon-separated string of every physical failure reason.

    Example:
        stage08_decision_reason(failed_candidate_row)
        -> "high_seed_drift;sequence_pathology_gly_rich"
    """
    reasons = []                                                              # Initialize an empty list to accumulate failure reasons
    # Check if the global confidence is below threshold
    if row["esmfold_mean_plddt"] < 70.0:                                      
        reasons.append("low_global_confidence")                               # Flag for low global confidence
    # Check if structural drift (RMSD) from seed exceeds threshold
    if row["rmsd_to_selected_seed"] > 3.5:                                    
        reasons.append("high_seed_drift")                                     # Flag for excessive structural deviation
    # Check if fewer than 50% of the newly mutated residues achieved a safe (>70.0) pLDDT score.
    if row["mutation_site_confidence_ge70_fraction"] < 0.5:                   
        reasons.append("low_mutation_site_confidence")                        # Flag for poor local prediction at mutation sites
    # Check if the model performed edits outside the pre-approved functional hotspots.
    if row["outside_editable_fraction"] > 0.0:                                
        reasons.append("outside_editable_region")                             # Flag for illegal editing
    # Check if the sequence contains a chain of >6 identical amino acids.
    if row["max_homopolymer_run"] > 6:                                        
        reasons.append("sequence_pathology_homopolymer")                      # Flag for homopolymer sequence pathology
    # Check if the sequence is >18% Glycine (too flexible/floppy to form a rigid structural binding tip).
    if row["gly_fraction"] > 0.18:                                            
        reasons.append("sequence_pathology_gly_rich")                         # Flag for high glycine sequence pathology
    return "pass" if not reasons else ";".join(reasons)                       # Return "pass" if no flags, else join flags with semicolons


def restore_candidate_sequence(validated: pd.DataFrame, ranked: pd.DataFrame) -> pd.DataFrame:
    """
    A vital safety mechanism that re-attaches the physical amino acid sequence strings 
    to the candidate metadata if they were dropped during intermediate CSV merging steps.

    It attempts a highly precise inner join using a massive key composed of 7 different 
    columns (scores, regimes, IDs). If that fails, it falls back to a simpler `sample_id` 
    merge. It guarantees the structural pipeline has actual sequence data to fold.
    """
    # Define a list of columns intended for the merge key
    merge_cols = [                                                            
        "sample_id",                                                          # Candidate identifier
        "generation_regime",                                                  # Mechanism of generation
        "final_multimodal_rank_score",                                        # Upstream ranking score
        "target_score",                                                       # Metric score for the target
        "strict_manifold_score",                                              # Plausibility score
        "structure_score",                                                    # Initial structural estimate
        "mutation_positions",                                                 # Annotated mutation list
    ]                                                                         # Close merge columns list
    # Intersect with actually present columns
    available_merge_cols = [c for c in merge_cols if c in ranked.columns and c in validated.columns] 
    # Log the columns being used for the merge operation
    print(f"[INFO] Attempting merge on columns: {available_merge_cols}")      

    # Drop any existing sequence column safely
    validated_base = validated.drop(columns=["candidate_sequence"], errors="ignore").copy() 
    # Perform a left join on the defined merge columns
    merged = validated_base.merge(                                            
        ranked[available_merge_cols + ["candidate_sequence"]],                # Right side: Ranked dataframe with the sequence added
        on=available_merge_cols,                                              # Join keys
        how="left",                                                           # Keep all rows from validated_base
    )                                                                         # Close merge statement

    # Check if merge failed to get sequences
    if "candidate_sequence" not in merged.columns or merged["candidate_sequence"].isna().any(): 
        print("[WARNING] Primary merge missed some sequences, falling back to sample_id merge") # Log the fallback operation
        # Create mapping using just sample_id
        fallback = ranked[["sample_id", "candidate_sequence"]].drop_duplicates("sample_id")     
        merged = merged.drop(columns=["candidate_sequence"], errors="ignore").merge(            # Strip old sequence data and merge again
            fallback,                                                         # Right side: Simple ID to sequence map
            on="sample_id",                                                   # Use only sample_id as the join key
            how="left",                                                       # Keep all rows from the primary table
        )                                                                     # Close merge statement

    # Final check if sequence column exists at all
    if "candidate_sequence" not in merged.columns:                            
        raise KeyError("candidate_sequence missing after merge")              # Abort execution with a missing key error
    # Check if literally every sequence is null and raise an error
    if merged["candidate_sequence"].isna().all():                             
        raise ValueError("All candidate_sequence values are missing after fallback merge")      # Abort indicating total failure
    # Check if partial missing sequences persist, gather their IDs, and raise an error
    if merged["candidate_sequence"].isna().any():                             
        missing_ids = merged.loc[merged["candidate_sequence"].isna(), "sample_id"].tolist()     # Gather the IDs of the missing rows
        raise ValueError(f"Some candidate_sequence values are still missing after merge. sample_id(s): {missing_ids}") # Abort and list IDs
    return merged                                                             # Return the successfully unified dataframe                                                        # Return the successfully unified dataframe


def main() -> None:
    # Parse CLI arguments, create output directories including a dedicated /pdbs folder exists to catch the raw structural outputs.
    args = parse_args()                                                                    # Parse CLI arguments into namespace
    args.out_dir.mkdir(parents=True, exist_ok=True)                                        # Create the root output directory
    pdb_dir = args.out_dir / "pdbs"                                                        # Define the nested path for PDB files
    pdb_dir.mkdir(parents=True, exist_ok=True)                                             # Create the nested PDB directory
    # Load validated_csv (Stage 07h QA check), and truncate it to user-defined top_k limit.
    validated = pd.read_csv(args.validated_csv).copy()                                     # Load the validated candidates CSV into memory
    if args.max_candidates is not None:                                                    # Check if a hard limit on candidates was set
        validated = validated.head(args.max_candidates)                                    # Truncate dataframe to the maximum limit
    validated = validated.head(args.top_k).copy()                                          # Keep only the top K candidates
    # Also load the broader ranked_csv and context_json to recover necessary background metadata.
    ranked = pd.read_csv(args.ranked_csv)                                                  # Load the ranked candidate data source
    context = json.loads(args.context_json.read_text())                                    # Load and parse the project context JSON
    # Run the safety protocol to ensure the physical amino acid sequence strings are perfectly attached to the metadata.
    merged = restore_candidate_sequence(validated, ranked)                                 # Recover candidate sequences missing from CSV
    # Load the original Canonical seed sequence and the allowed mutation positions from the context
    seed_seq = str(context["selected_seed"]["seed_sequence"])                              # Extract the original seed amino acid sequence
    editable = set(int(x) for x in context["editable_region"].get("hotspot_positions", []))# Extract allowed mutation positions into a set
    # Detect the optimal hardware and load the massive ESMFold neural network into VRAM.
    torch_mod, _, _ = import_esmfold()                                                     # Import just PyTorch to determine the device
    device = choose_device(torch_mod, args.device)                                         # Determine best compute device automatically
    torch, tokenizer, model = load_esmfold(device)                                         # Load ESMFold pipeline components to device
    # If the seed PDB already exists and we are resuming (args.resume), load the seed coordinates and pLDDTs, compute the mean pLDDT and assign NaN to pTM.
    seed_pdb_path = pdb_dir / "seed_selected_seed.pdb"                                     # Define path where seed PDB will be stored
    if args.resume and seed_pdb_path.exists():                                             # If resuming and the seed PDB exists
        seed_coords, seed_plddts, _ = parse_pdb_file(seed_pdb_path)                        # Load the cached seed structure
        seed_plddt = float(np.mean(seed_plddts)) if len(seed_plddts) else math.nan         # Compute cached mean confidence
        seed_ptm = math.nan                                                                # Assign NaN to pTM as it isn't in standard PDBs
    # If the seed PDB does not exist or we are not resuming
    else:                                                                                  # If not resuming or seed missing
    # Run the heavy ESMFold inference on the seed, save its PDB, and parse out the Alpha-Carbon coordinates and pLDDTs.
        seed_pdb, seed_plddt, seed_ptm = fold_sequence(                                    # Run the heavy ESMFold inference on seed
            seed_seq, tokenizer, model, torch, device=device,                              # Pass parameters
            chunk_size=args.chunk_size, num_recycles=args.num_recycles,                    # Pass runtime configurations
        )                                                                                  # Close function call
        seed_pdb_path.write_text(seed_pdb)                                                 # Save the newly generated PDB to disk
        seed_coords, seed_plddts, _ = parse_ca_coords_and_plddt(seed_pdb)                  # Extract coordinates and metrics from output

    # Iterate through the merged dataframe evaluating every candidate in the top-tier shortlist.
    rows = []                                                                              # Initialize list to hold metrics per candidate
    # For every candidate in the top-tier shortlist
    for row in merged.to_dict(orient="records"):                                           # Iterate through the candidate dataframe
        # Extract the sample ID, sequence, mutation positions, and PDB path.
        sample_id = int(row["sample_id"])                                                  # Extract sample ID
        seq = str(row.get("candidate_sequence", "") or "")                                 # Extract the sequence string
        muts = parse_mutation_positions(row.get("mutation_positions", ""))                 # Parse mutation positions into a list
        pdb_path = pdb_dir / f"candidate_{sample_id}.pdb"                                  # Define candidate-specific PDB path
        # If the candidate PDB already exists and we are resuming
        if args.resume and pdb_path.exists():                                              # Check if candidate PDB is cached and resuming
            # Load the cached PDB, compute the mean pLDDT, and assign NaN to pTM
            coords, plddts, residues = parse_pdb_file(pdb_path)                            # Load data from the cached file
            mean_plddt = float(np.mean(plddts)) if len(plddts) else math.nan               # Compute cached mean pLDDT
            ptm = math.nan                                                                 # Assign NaN to cached pTM
        # If the candidate PDB does not exist or we are not resuming
        else:                                                                              # If no cache or not resuming
            # Run ESMFold inference on the candidate, save its PDB, and parse out the Alpha-Carbon coordinates and pLDDTs
            pdb_text, mean_plddt, ptm = fold_sequence(                                     # Run ESMFold inference for the candidate
                seq, tokenizer, model, torch, device=device,                               # Pass parameters
                chunk_size=args.chunk_size, num_recycles=args.num_recycles,                # Pass configs
            )                                                                              # Close function call
            pdb_path.write_text(pdb_text)                                                  # Save predicted structure to disk
            coords, plddts, residues = parse_ca_coords_and_plddt(pdb_text)                 # Parse new structure data

        # Compute structural drift from seed (Kabsch RMSD), fraction of confident mutations (% of mutated residues with pLDDT > 70), 
        # mean pLDDT of mutations, and the fraction of illegal edits (outside the editable region).
        rmsd_to_seed = kabsch_rmsd(coords, seed_coords)                                         # Calculate structural drift vs the original seed
        mut_conf_frac = residue_confidence_fraction(residues, plddts, set(muts), cutoff=70.0)   # Check fraction of confident mutation sites
        mut_conf_mean = residue_mean_confidence(residues, plddts, set(muts))                    # Calculate average confidence across mutations
        outside_frac = float(np.mean([pos not in editable for pos in muts])) if muts else 0.0   # Measure illegal edit occurrences
        # Sequence-health checks: calculating exactly how much of the sequence is hydrophobic, charged, aromatic, or made of highly flexible Glycine.
        seq_metrics = {                                                                    # Create dictionary of sequence properties
            "sequence_entropy": shannon_entropy(seq),                                      # Record shannon entropy
            "gly_fraction": gly_fraction(seq),                                             # Record glycine ratio
            "hydrophobic_fraction": hydrophobic_fraction(seq),                             # Record hydrophobicity ratio
            "charged_fraction": charged_fraction(seq),                                     # Record charged residue ratio
            "aromatic_fraction": aromatic_fraction(seq),                                   # Record aromatic residue ratio
            "max_homopolymer_run": max_homopolymer_run(seq),                               # Record longest repetitive stretch
            "mutation_count": len(muts),                                                   # Record number of mutations
            "mutation_span": mutation_span(muts),                                          # Record the distance between extreme mutations
            "outside_editable_fraction": outside_frac,                                     # Record the fraction of unauthorized edits
        }                                                                                  # Close dictionary
        # Merge the original generative metadata, the new sequence health metrics, 
        # and the heavy 3D structural metrics (RMSD, global pLDDT, local mutation pLDDT) 
        # into a master dictionary and appends it to a rows list.
        rows.append({                                                                      # Build row data and append to master list
            **row,                                                                         # Merge original candidate metadata
            **seq_metrics,                                                                 # Merge calculated sequence metrics
            "esmfold_mean_plddt": mean_plddt,                                              # Store calculated global pLDDT
            "esmfold_ptm": ptm,                                                            # Store calculated pTM score
            "rmsd_to_selected_seed": rmsd_to_seed,                                         # Store structural drift metric
            "mutation_site_confidence_ge70_fraction": mut_conf_frac,                       # Store mutation-site confidence ratio
            "mutation_site_mean_plddt": mut_conf_mean,                                     # Store average mutation confidence
            "seed_esmfold_mean_plddt": seed_plddt,                                         # Record reference seed pLDDT for context
            "seed_esmfold_ptm": seed_ptm,                                                  # Record reference seed pTM for context
            "candidate_pdb": str(pdb_path),                                                # Store file path to the generated structure
        })                                                                                 # Close append statement
    # Convert the rows list (all candidates) to a DataFrame and sort it hierarchically, primarily by the original Stage 07 score,
    # but tie-broken by structural confidence (highest pLDDT) and minimal physical deformation (lowest RMSD). 
    out_df = pd.DataFrame(rows).sort_values(                                               # Convert list to DataFrame and sort it
        ["final_multimodal_rank_score", "esmfold_mean_plddt", "rmsd_to_selected_seed"],    # Define sort keys (rank, confidence, drift)
        ascending=[False, False, True],                                                    # Define sort directions (descending, asc)
    ).reset_index(drop=True)                                                               # Reset dataframe index after sorting
    # Create a formal structural rank column
    out_df["stage08_structural_rank"] = np.arange(1, len(out_df) + 1)                      
    # Enforces the Six Rules of Physical Safety in a new column ("stage08_pass").
    # If a protein passes every single rule, it earns a stage08_pass = True. If it fails, it explicitly logs the exact physical law it broke.
    out_df["stage08_pass"] = (                                                             # Define the boolean pass/fail column
        (out_df["esmfold_mean_plddt"] >= 70.0)                                             # Rule 1: Good overall structural confidence
        & (out_df["rmsd_to_selected_seed"] <= 3.5)                                         # Rule 2: Scaffold retains shape
        & (out_df["mutation_site_confidence_ge70_fraction"].fillna(0.0) >= 0.50)           # Rule 3: Mutated regions are confident
        & (out_df["gly_fraction"] <= 0.18)                                                 # Rule 4: Not pathologically floppy
        & (out_df["max_homopolymer_run"] <= 6)                                             # Rule 5: Avoids repetitive artifacts
        & (out_df["outside_editable_fraction"] <= 0.0)                                     # Rule 6: No mutations outside target zone
    )                                                                                      # Close boolean assignment
    # Populate reason column for failed designs.
    out_df["stage08_decision_reason"] = out_df.apply(stage08_decision_reason, axis=1)      

    # Export the dataFrame to CSV and JSON formats for downstream programmatic use (Stage 08b).
    csv_path = args.out_dir / "stage08_structural_fasttrack_summary.csv"                   # Define output CSV path
    json_path = args.out_dir / "stage08_structural_fasttrack_summary.json"                 # Define output JSON path
    md_path = args.out_dir / "stage08_structural_fasttrack_report.md"                      # Define output Markdown report path
    out_df.to_csv(csv_path, index=False)                                                   # Write dataframe to CSV
    json_path.write_text(out_df.to_json(orient="records", indent=2))                       # Write dataframe to formatted JSON
    
    # Create a detailed Markdown report building a dynamic summary table 
    md = [                                                                                 # Begin construction of the Markdown document list
        "# Stage 08 structural fast-track report",                                         # MD Title
        "",                                                                                # Blank line
        f"- candidates analyzed: **{len(out_df)}**",                                       # Log candidate count
        f"- seed ESMFold mean pLDDT: **{seed_plddt:.2f}**" if not math.isnan(seed_plddt) else "- seed ESMFold mean pLDDT: **nan**", # Log seed pLDDT
        f"- seed ESMFold pTM: **{seed_ptm:.4f}**" if not math.isnan(seed_ptm) else "- seed ESMFold pTM: **nan**",                   # Log seed pTM
        f"- resume mode: **{args.resume}**",                                               # Log cache usage
        "",                                                                                # Blank line
        "## Candidate summary",                                                            # Section header
        "",                                                                                # Blank line
        "| rank | sample_id | regime | final_score | mean_pLDDT | mut_mean_pLDDT | RMSD_to_seed | mut_conf>=70 | pass | reason |", # Table header
        "|---:|---:|---|---:|---:|---:|---:|---:|:---:|---|",                                                                        # Table alignment
    ]                                                                                      # Close list construction
    for _, r in out_df.iterrows():                                                         # Loop over dataframe rows to build MD table
        md.append(                                                                         # Add row to Markdown table
            f"| {int(r['stage08_structural_rank'])} | {int(r['sample_id'])} | {r['generation_regime']} | "                            # Add ID and regime columns
            f"{float(r['final_multimodal_rank_score']):.6f} | {float(r['esmfold_mean_plddt']):.2f} | "                              # Add score and global confidence columns
            f"{float(r['mutation_site_mean_plddt']):.2f} | {float(r['rmsd_to_selected_seed']):.3f} | "                               # Add local confidence and RMSD columns
            f"{float(r['mutation_site_confidence_ge70_fraction']):.3f} | {bool(r['stage08_pass'])} | {r['stage08_decision_reason']} |" # Add pass/fail and reason columns
        )                                                                                  # Close string formatting
    md.extend([                                                                            # Add interpretation footers
        "",                                                                                # Blank line
        "## Interpretation",                                                               # Section Header
        "- Prefer candidates with high multimodal score, high global pLDDT, high mutation-site confidence, and low RMSD to the selected seed.", # Advice line
        "- `stage08_pass` is intentionally conservative; near-pass candidates can still be useful backups.",                                    # Advice line
        "- If at least two candidates pass, those should become the immediate closeout handoff set.",                                           # Advice line
    ])                                                                                     # Close list extension
    md_path.write_text("\n".join(md))                                                      # Render and save the full Markdown report

    # Print the path to the output files and the dataframe to the terminal.
    print(f"Wrote: {csv_path}")                                                            # Log the writing of CSV
    print(f"Wrote: {json_path}")                                                           # Log the writing of JSON
    print(f"Wrote: {md_path}")                                                             # Log the writing of Markdown
    print(f"Wrote seed PDB: {seed_pdb_path}")                                              # Log the saving of the seed PDB
    for p in sorted(pdb_dir.glob("candidate_*.pdb")):                                      # Iterate through all generated candidate PDBs
        print(f"Wrote candidate PDB: {p}")                                                 # Log the saving of candidate PDBs

    cols = [                                                                               # Define subset of columns for terminal output
        "stage08_structural_rank",                                                         # Rank column
        "sample_id",                                                                       # Identifier column
        "generation_regime",                                                               # Generation mechanism column
        "final_multimodal_rank_score",                                                     # Model score column
        "esmfold_mean_plddt",                                                              # Confidence column
        "mutation_site_mean_plddt",                                                        # Mutated local confidence column
        "rmsd_to_selected_seed",                                                           # Drift metric column
        "mutation_site_confidence_ge70_fraction",                                          # Threshold ratio column
        "stage08_pass",                                                                    # Decision column
        "stage08_decision_reason",                                                         # Justification column
    ]                                                                                      # Close subset definition
    print("\nTop structural summary:")                                                     # Print spacer and header to terminal
    print(out_df[cols].to_string(index=False))                                             # Print dataframe table to terminal


if __name__ == "__main__":
    main()                                                                                 # Execute the primary application flow