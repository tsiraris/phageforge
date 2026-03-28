# PhageForge

**PhageForge** is a research project exploring **validity-aware, host-conditioned phage receptor-binding protein (RBP) design** using protein language models, embedding-based host predictors, and constrained optimization.

The project focuses on a realistic design question:

> Can we retarget a validated phage RBP scaffold from one bacterial host genus toward another, while preserving RBP-like structure and family identity?

The current pipeline does **not** treat host transfer as a simple score-maximization problem. Instead, it combines:

* curated phage-host RBP datasets
* ESM-2 sequence embeddings
* host prediction probes
* validity-aware evaluation
* scaffold-family constraints
* host-ladder optimization
* controlled relaxation for harder host shifts

The main design case implemented in the repo is a **host ladder**:

```text
Klebsiella -> Enterobacter -> Acinetobacter
```

---

# Core idea

Initial optimization experiments showed that improving host prediction score alone is **not enough** to guarantee biologically plausible RBPs.

To address this, the project evolved into a **two-stage design framework**:

## Stage 05 — validity-aware evaluation

Evaluate optimized candidates with additional biological plausibility checks:

* RBP manifold similarity
* family-retention consistency
* mutation-aware plausibility / naturalness proxies
* corrected mutation accounting
* deduplicated ranking of unique candidates

## Stage 06 — family-conditioned host-ladder design

Use corrected Stage 05 outputs to:

* select a canonical seed scaffold
* build a tight seed-local family manifold
* define target-host anchor references
* optimize inside a hotspot-derived mutation window
* retarget stepwise across hosts with explicit structural constraints

This turns the project from:

* **single-objective host score optimization**

into:

* **validity-aware, scaffold-constrained protein retargeting**

---

# Main findings

## 1. Enterobacter transfer works under strict family constraints

A canonical Klebsiella seed can be shifted toward Enterobacter with:

* very high family similarity
* very small mutation count
* strong target-score improvement

## 2. Acinetobacter transfer stalls under overly strict constraints

The second ladder step initially fails to improve meaningfully when:

* mutation count is kept too low
* scaffold preservation is too strongly enforced

## 3. Controlled relaxation enables larger host shifts

Relaxing the optimization slightly by allowing:

* more mutations
* stronger target-anchor pull
* entropy-guided exploration
* weaker family penalty

improves Acinetobacter transfer while still preserving high scaffold similarity.

## 4. The relaxed behavior is reproducible across seeds

Running the relaxed Acinetobacter step from more than one Enterobacter-derived intermediate shows the same qualitative behavior.

### Main project insight

There is a real trade-off between:

* **structural/scaffold preservation**
* **functional movement toward a more distant host manifold**

This is the central design result of the repository.

---

# Repository structure

```text
phageforge/
│
├── data/
│   ├── raw/
│   └── processed/
│
├── phageforge/
│
├── results/
│   ├── broad/
│   ├── strict/
│   ├── structural/
│   ├── design_runs/
│   ├── analysis/
│   │   ├── validity/
│   │   └── validity_corrected/
│   ├── phaseA/
│   │   ├── phaseA_plan.json
│   │   ├── step1_enterobacter/
│   │   ├── step2_seed/
│   │   └── step2_acinetobacter/
│   ├── phaseB/
│   │   ├── relaxed_acinetobacter/
│   │   └── relaxed_acinetobacter_seed2/
│   ├── summary.csv
│   ├── summary_aggregate.csv
│   └── benchmark_report.md
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
│   ├── 04c_compare_seed_to_candidates.py
│   ├── 05_compute_validity_metrics.py
│   ├── 05b_rank_validated_candidates.py
│   ├── 05c_make_validity_report.py
│   ├── 06a_select_phaseA_family.py
│   ├── 06b_optimize_family_constrained.py
│   └── 06c_pick_phaseA_followup_seed.py
│
├── pyproject.toml
└── README.md
```

---

# Pipeline overview

## 1. Build and refine datasets

Construct and refine datasets into:

* broad
* strict
* structural
* structural_plus

## 2. Embed RBPs

Use **ESM-2** to generate reusable embeddings.

## 3. Train host predictors

Train classifiers to estimate host genus from embeddings.

## 4. Evaluate validity (Stage 05)

Re-rank candidates using biological plausibility signals beyond host score.

## 5. Run host-ladder design (Stage 06)

### 06a — Build design space

* choose seed
* define family
* derive hotspots and priors

### 06b — Optimize

* constrained mutation
* multi-objective scoring

### 06c — Select next seed

* pick best candidate
* export FASTA

---

# Phase interpretation

## Phase A — strict

* strong structural constraint
* limited mutation
* reliable short-distance transfer

## Phase B — relaxed

* increased mutation freedom
* improved long-distance transfer
* still structurally grounded

---

# Key outputs

### Phase A

* `results/phaseA/phaseA_plan.json`
* `results/phaseA/step1_enterobacter/`
* `results/phaseA/step2_seed/`

### Phase B

* `results/phaseB/relaxed_acinetobacter/`
* `results/phaseB/relaxed_acinetobacter_seed2/`

---

# Installation

```bash
pip install -e .
```

Dependencies:

```bash
pip install pandas numpy torch transformers scikit-learn joblib tabulate
```

---

# Minimal usage

### Build plan

```bash
python scripts/06a_select_phaseA_family.py ...
```

### Step 1

```bash
python scripts/06b_optimize_family_constrained.py ...
```

### Pick seed

```bash
python scripts/06c_pick_phaseA_followup_seed.py ...
```

### Step 2

```bash
python scripts/06b_optimize_family_constrained.py ...
```

---

# Summary

PhageForge now supports:

* validity-aware evaluation
* scaffold-constrained optimization
* host-ladder retargeting
* strict vs relaxed transfer analysis

The key result:

> Cross-host RBP transfer is possible — but larger host shifts require balancing structural preservation with functional adaptation.
