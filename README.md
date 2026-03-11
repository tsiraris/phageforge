# PhageForge

**PhageForge** is a research project exploring **host-conditioned phage receptor-binding protein (RBP) design** using protein language models.

The repository implements a pipeline that:

1. builds curated phage-host datasets
2. extracts receptor-binding proteins (RBPs)
3. refines the dataset from **broad -> strict -> structural**
4. embeds RBPs with **ESM-2**
5. trains host prediction models
6. benchmarks model performance across dataset variants
7. performs **host-conditioned RBP optimization**
8. evaluates whether optimized candidates improve over the original seed

The goal is to test whether protein language models can support the design of RBPs that shift scoring toward a desired bacterial host genus.

---

# Project structure

```text
phageforge/
│
├── data/
│   ├── raw/                         # downloaded source data
│   └── processed/                   # processed datasets and embeddings
│
├── phageforge/                      # reusable Python package
│
├── results/
│   ├── broad/                       # trained broad-dataset models
│   ├── strict/                      # trained strict-dataset models
│   ├── structural/                  # optional structural-dataset models
│   ├── design_runs/                 # optimization run outputs
│   ├── summary.csv                  # benchmark results by run
│   ├── summary_aggregate.csv        # benchmark aggregates
│   └── benchmark_report.md          # generated benchmark report
│
├── scripts/
│   ├── 01_build_dataset.py
│   ├── 01b_filter_strict_rbps.py
│   ├── 01c_structural_filter_rbps.py
│   ├── 02_embed_rbps.py
│   ├── 03_train_phi_mlp.py
│   ├── 03b_linear_probe.py
│   ├── 03c_run_benchmark.py
│   ├── 03d_make_benchmark_report.py
│   ├── 04_optimize_rbp_for_host.py
│   ├── 04b_summarize_design_run.py
│   └── 04c_compare_seed_to_candidates.py
│
├── pyproject.toml
└── README.md