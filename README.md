# PhageForge

**PhageForge** is a research project exploring **host-conditioned phage receptor-binding protein (RBP) design** using protein language models.

The repository implements a full pipeline that:

1. builds curated phage–host datasets
2. extracts receptor-binding proteins (RBPs)
3. embeds RBPs using **ESM-2 protein language models**
4. trains host prediction models
5. benchmarks **broad vs strict datasets**
6. performs **host-conditioned RBP optimization**

The goal is to explore whether protein language models can guide the design of RBPs that preferentially target specific bacterial genera.

---

# Project Structure

```
phageforge/
│
├── data/
│   ├── raw/                     # raw downloaded datasets
│   └── processed/               # processed datasets and embeddings
│
├── phageforge/                  # reusable Python modules (future expansion)
│
├── results/
│   ├── benchmark/               # benchmark outputs
│   └── design_runs/             # RBP optimization outputs
│
├── scripts/
│   ├── 01_build_dataset.py          # build broad RBP dataset
│   ├── 01b_filter_strict_rbps.py    # create strict RBP dataset
│   ├── 02_embed_rbps.py             # embed RBPs with ESM-2
│   ├── 03_train_phi_mlp.py          # baseline MLP classifier
│   ├── 03b_linear_probe.py          # logistic regression probe
│   ├── 03c_run_benchmark.py         # benchmark runner
│   └── 04_optimize_rbp_for_host.py  # host-conditioned RBP optimization
│
├── pyproject.toml
└── README.md
```

---

# Pipeline Overview

## Stage 1 — Dataset construction

```
scripts/01_build_dataset.py
```

Builds the **broad RBP dataset** by:

* parsing Virus-Host DB
* filtering for ESKAPEE bacterial genera
* extracting RBPs from GenBank records using keyword matching

Output:

```
data/processed/rbp_dataset_eskapee_stage1.csv
```

---

## Stage 1b — Strict dataset

```
scripts/01b_filter_strict_rbps.py
```

Filters the broad dataset to create a **high-confidence RBP dataset**.

Strict RBPs include:

* tail fiber proteins
* receptor binding proteins
* tail spike proteins

Output:

```
data/processed/rbp_dataset_eskapee_strict.csv
```

---

# Stage 2 — Protein embeddings

```
scripts/02_embed_rbps.py
```

Embeds RBP amino-acid sequences using **ESM-2**.

Example:

```bash
python scripts/02_embed_rbps.py \
  --csv data/processed/rbp_dataset_eskapee_stage1.csv \
  --model facebook/esm2_t33_650M_UR50D \
  --batch_size 2
```

Outputs:

```
data/processed/esm2_embeddings.pt
data/processed/esm2_embeddings_index.csv
```

---

# Stage 3 — Host prediction models

Two baselines are implemented:

### MLP classifier

```
scripts/03_train_phi_mlp.py
```

A neural network classifier trained on ESM embeddings.

---

### Linear probe

```
scripts/03b_linear_probe.py
```

A logistic regression probe trained on frozen ESM embeddings.

This script also saves an **inference-ready model**:

```
model.joblib
label_classes.json
```

These artifacts are later used for **RBP design scoring**.

---

# Stage 3b — Benchmark

```
scripts/03c_run_benchmark.py
```

Runs a full benchmark comparing:

| Dataset | Model        |
| ------- | ------------ |
| Broad   | Linear probe |
| Broad   | MLP          |
| Strict  | Linear probe |
| Strict  | MLP          |

Across **multiple random seeds**.

Outputs:

```
results/summary.csv
results/summary_aggregate.csv
results/benchmark_report.md
```

---

# Stage 4 — Host-conditioned RBP optimization

```
scripts/04_optimize_rbp_for_host.py
```

This stage performs **guided sequence optimization**.

Given:

* a seed RBP sequence
* a target host genus

The script:

1. generates mutated candidate sequences
2. embeds them with ESM-2
3. scores them using the trained host classifier
4. ranks candidates by predicted host probability

Example:

```bash
python scripts/04_optimize_rbp_for_host.py \
  --seed_protein_id YP_004300729.1 \
  --target_host Klebsiella \
  --model_path results/broad/linear_probe/seed_42/model.joblib \
  --label_classes_path results/broad/linear_probe/seed_42/label_classes.json
```

Outputs:

```
results/design_runs/
    round_1_candidates.csv
    round_2_candidates.csv
    top_candidates.csv
```

These represent candidate RBPs predicted to target the specified host genus.

---

# Benchmark Summary (Current)

| Dataset | Model        | Accuracy |
| ------- | ------------ | -------- |
| Broad   | Linear probe | ~0.67    |
| Broad   | MLP          | ~0.36    |
| Strict  | Linear probe | ~0.74    |
| Strict  | MLP          | ~0.65    |

The **linear probe consistently performs best**, indicating that **ESM embeddings already encode host-specific information**.

---

# Requirements

Python ≥ 3.10

Main libraries:

```
torch
transformers
scikit-learn
pandas
numpy
biopython
joblib
```

Install locally:

```bash
pip install -e .
```

---

# Future Work

* structural RBP filtering
* diffusion-based sequence generation
* experimental host-range validation
* integration with phage genome design pipelines

---

# Author

Developed as a machine learning research project exploring **protein language models for phage host targeting**.
