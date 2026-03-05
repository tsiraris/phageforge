# PhageForge

Stage 1: Build an ESKAPEE phage RBP dataset from Virus-Host DB + GenBank.
Stage 2: Embed RBPs with ESM-2 and train a baseline host genus classifier.

## Quickstart
- Stage 1 (dataset): `python scripts/01_build_dataset.py`
- Stage 2 (embeddings): `python scripts/02_embed_rbps.py`
- Stage 2 (train): `python scripts/03_train_phi_mlp.py`