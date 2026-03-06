# PhageForge

PhageForge is a research-oriented project for bacteriophage host prediction and early-stage host-conditioned design.

The project currently focuses on:

- building phage receptor-binding protein (RBP) datasets from public resources
- embedding RBPs with protein foundation models (ESM-2)
- benchmarking host genus prediction under different biological curation strategies
- comparing broad adsorption-module datasets against strict high-confidence RBP datasets

---

## Project motivation

Bacteriophage host specificity is strongly influenced by adsorption machinery, especially receptor-binding proteins such as tail fibers and tail spikes.

A key question in this project is:

> Does stricter biological curation improve host prediction more than simply scaling dataset size?

To test this, the project compares:

- **Broad dataset**: larger, noisier adsorption-module protein collection
- **Strict dataset**: smaller, high-confidence host-determining RBP set

and evaluates both with:

- **MLP classifier**
- **Linear probe (logistic regression)**

---

## Repository structure

```text
phageforge/
├── data/
│   ├── raw/
│   ├── processed/
│   └── cache/
├── phageforge/
│   └── data/
├── scripts/
│   ├── 01_build_dataset.py
│   ├── 01b_filter_strict_rbps.py
│   ├── 02_embed_rbps.py
│   ├── 03_train_phi_mlp.py
│   ├── 03b_linear_probe.py
│   ├── 03c_run_benchmark.py
│   └── 03d_make_benchmark_report.py
├── results/
├── pyproject.toml
└── README.md