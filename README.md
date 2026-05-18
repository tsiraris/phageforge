# PhageForge

**PhageForge** is a research project on **validity-aware, scaffold-constrained phage receptor-binding protein (RBP) design**.

The repository began as a host-conditioned sequence-design pipeline, evolved into a validity-aware and family-constrained design system, attempted a proxy-guided redesign, and now culminates in a full **Stage 10 structure-conditioned inverse-folding redesign architecture**.

By this point the project addresses a sharper question:

> If sequence-first generative search—even when tightly constrained by structural proxies—fails physical 3D validation, can explicitly conditioning generation on the 3D scaffold via inverse-folding cross the final wet-lab viability barrier?

This README reflects the repository state after the definitive Stage 10 redesign cycle.

---

## Scientific framing

PhageForge does **not** treat host retargeting as simple score maximization.
It also does **not** assume that better sequence-space or proxy-structure scores automatically imply structural viability.

The project is organized around a sequential research logic governed by falsification:

* the **sequence / host-transfer branch** asks: **does a candidate move toward the target host manifold?**
* the **family / scaffold branch** asks: **does it remain near the validated seed and family neighborhood?**
* the **Stage 08 structural branch** asks: **does full structural validation support the candidate as a plausible scaffold-preserving redesign?**
* the **Stage 09 redesign branch** asks: **can stricter localized sequence search and structural prefiltering improve this outcome?**
* the **Stage 10 inverse-folding branch** asks: **can true structure-conditioned generation solve the fold-retention bottleneck?**

This makes Stage 10 the ultimate architectural response to the failure of sequence-first models.

---

## Project evolution

### Stages 01-03 — dataset building and host prediction

The project begins by building curated RBP datasets, embedding them with protein language models, and training the host-prediction backbone.

### Stage 04 — host-conditioned RBP optimization

Explores direct host retargeting by proposing mutations and ranking candidates with the host probe.

### Stage 05 — validity-aware evaluation

Establishes that **host score alone is insufficient** and adds manifold similarity, family retention, and plausibility-aware evaluation.

### Stage 06 — family-conditioned host-ladder design

Turns the project into a scaffold-constrained design program by selecting a family, defining editable regions, and grounding edits in a validated seed scaffold.

### Stage 07 — multimodal, locally generated design

Builds the main candidate pool via local ESM3-guided generation and multimodal ranking.

### Stage 08 — structural fast-track validation

Adds robust ESMFold-based structural validation, definitively falsifying the unconstrained Stage 07 candidates due to catastrophic 3D drift.

### Stage 09 — structure-aware localized redesign

A sequence-first response to Stage 08. Shrinks the edit space and adds surrogate-driven prefiltering. It successfully improves pre-validation proxy scores but still fails full 3D ESMFold validation, proving sequence-first proxies are insufficient.

### Stage 10 — structure-conditioned inverse-folding redesign

Introduced after the Stage 09 failure signal. Stage 10 explicitly anchors redesign to the physical 3D seed scaffold using ESM-IF1, making the fold itself the generator.

* **10a**: build a structure-conditioned redesign context explicitly anchored to a seed PDB
* **10b**: run inverse-folding beam search, scoring proposals directly against the 3D backbone
* **10c**: prefilter candidates using a diversity-aware algorithm to prevent mode collapse
* **10d**: validate the top Stage 10 panels using a wrapper around the unmodified Stage 08 validator
* **10e**: build the final comparative report pitting inverse-folding outcomes against historical baselines

---

## Main Stage 10 idea

Stage 10 deliberately abandons sequence-first approximations in favor of a **true physical simulation loop**:

1. start from the absolute 3D atomic coordinates of the validated wild-type seed (`.pdb`)
2. define an aggressively minimized hard/soft edit space
3. query the ESM-IF1 inverse-folding model for substitutions that thermodynamically stabilize that specific 3D geometry
4. score candidates with a physics-heavy composite:
* 3D backbone log-likelihood (ESM-IF1)
* target-host probability (ESM2 + LR)
* family evolutionary cosine
* sequence identity preservation


5. prefilter survivors using a greedy-diverse embedding space algorithm
6. send the elite panel to the exact same full structural validator used in Stage 08

The point of Stage 10 is to make the physical scaffold the *generator*, ensuring 3D topological viability is a prerequisite for generation rather than a downstream hope.

---

## What Stage 10 utilizes

The Stage 10 update introduces the ultimate generative architecture capabilities:

* **explicit 3D scaffold anchoring**
* the sequence generation is locked mathematically to the X, Y, Z coordinates of a physical `.pdb` file
* the mutation budget is slashed further (e.g., 1-4 edits) to prioritize absolute stability


* **inverse-folding beam search**
* replaces the sequence language model (ESM3/ESM2) with a Graph Neural Network (ESM-IF1)
* structural stability is no longer guessed via proxies; it is simulated dynamically


* **composite 3D fitness scoring**
* scores candidates primarily by physical backbone compatibility (`if1_log_likelihood`) rather than additive heuristic metrics


* **diversity-aware prefiltering**
* utilizes `greedy_diverse_subset` clustering to guarantee the final validation panel represents distinct structural hypotheses, preventing mode collapse


* **wrapper-based scientific coherence**
* executes the heavy ESMFold validation via a sterile subprocess wrapper
* intentionally reuses the unmodified Stage 08 validator to prove any observed improvements are scientifically genuine



---

## Repository structure

```text
phageforge/
│
├── data/
│   ├── raw/
│   └── processed/
│
├── phageforge/
│   ├── generation/
│   ├── tissue/
│   ├── eval/
│   ├── stage07_utils.py
│   ├── stage09_utils.py
│   └── stage10_utils.py
│
├── results/
│   ├── analysis/
│   ├── phaseA/
│   ├── stage07/
│   ├── stage08/
│   ├── stage09/
│   └── stage10/
│
├── scripts/
│   ├── 01_build_dataset.py
│   ...
│   ├── 09a_define_edit_space.py
│   ...
│   ├── 09g_make_stage09_report.py
│   ├── 10a_prepare_stage10_structure_context.py
│   ├── 10b_run_inverse_folding_beam_search.py
│   ├── 10c_prefilter_stage10_candidates.py
│   ├── 10d_validate_stage10_candidates.py
│   └── 10e_make_stage10_report.py
│
├── notebooks/
│   ├── 09_stage08_closeout_sagemaker.ipynb
│   ├── 10_stage09_structure_aware_redesign_sagemaker.ipynb
│   └── 11_stage10_inverse_folding_sagemaker.ipynb
│
├── pyproject.toml
└── README.md

```

---

## Installation

Core Python dependencies are declared in `pyproject.toml`.

A practical GPU environment for Stage 10 requires specific heavy infrastructure:

```bash
pip install -e .
pip install 'pandas<3' 'huggingface_hub<1' esm transformers accelerate scikit-learn biopython fair-esm

```

Because Stage 10 relies on:

* ESM2 embeddings for target scoring and spatial diversity
* the `fair-esm` inverse-folding library (`esm_if1_gvp4`)
* ESMFold for final 3D atomic validation

you must ensure:

* PyTorch with CUDA is available and properly mapped
* the machine has significant VRAM (e.g., AWS G5/A10G or similar) to hold both ESM2 and ESM-IF1 simultaneously in memory during the search loop

---

## Recommended Stage 10 workflow

### 1. Build the structure-conditioned context

```bash
python scripts/10a_prepare_stage10_structure_context.py \
  --context_json results/stage07/context/stage07_context.base.json \
  --strict_csv data/processed/rbp_dataset_eskapee_strict.csv \
  --validation_dir results/stage08/structural_fasttrack_top3 \
  --output_json results/stage10/context/stage10_context.json \
  --max_edit_positions 6 \
  --soft_positions 3 \
  --min_mutations 1 \
  --max_mutations 4

```

This dynamically locates the seed `.pdb` and produces a heavily restricted inverse-folding mutation map.

### 2. Run inverse-folding beam search

```bash
python scripts/10b_run_inverse_folding_beam_search.py \
  --stage10_context_json results/stage10/context/stage10_context.json \
  --predictor_model results/broad/linear_probe/seed_42/model.joblib \
  --label_classes_json results/broad/linear_probe/seed_42/label_classes.json \
  --out_csv results/stage10/search/stage10_search_candidates.csv \
  --out_json results/stage10/search/stage10_search_summary.json \
  --embedding_model facebook/esm2_t33_650M_UR50D \
  --if_device cuda \
  --rounds 4 \
  --beam_width 24 \
  --proposals_per_parent 8 \
  --substitutions_per_position 3

```

### 3. Prefilter and select diverse candidates

```bash
python scripts/10c_prefilter_stage10_candidates.py \
  --stage10_context_json results/stage10/context/stage10_context.json \
  --search_csv results/stage10/search/stage10_search_candidates.csv \
  --out_topk_csv results/stage10/prefilter/stage10_top10.csv \
  --out_topk_final_csv results/stage10/prefilter/stage10_top3.csv \
  --out_json results/stage10/prefilter/stage10_prefilter_summary.json \
  --top_k 10 \
  --top_k_final 3

```

### 4. Validate the Stage 10 panel with the Stage 08 validator

```bash
python scripts/10d_validate_stage10_candidates.py \
  --validated_csv results/stage10/prefilter/stage10_top3.csv \
  --ranked_csv results/stage10/search/stage10_search_candidates.csv \
  --context_json results/stage10/context/stage10_context.json \
  --out_dir results/stage10/validation_top3 \
  --out_json results/stage10/validation_top3/stage10_launch.json \
  --top_k 3 \
  --device cuda \
  --chunk_size 128

```

### 5. Build the final Stage 10 comparative report

```bash
python scripts/10e_make_stage10_report.py \
  --stage10_context_json results/stage10/context/stage10_context.json \
  --search_csv results/stage10/search/stage10_search_candidates.csv \
  --prefilter_csv results/stage10/prefilter/stage10_top10.csv \
  --validation_csv results/stage10/validation_top3/stage08_structural_fasttrack_summary.csv \
  --baseline_validation_csv results/stage09/validation_top3/stage08_structural_fasttrack_summary.csv \
  --out_dir results/stage10/final_report

```

---

## Stage 10 outputs

Typical Stage 10 outputs live under:

```text
results/stage10/
  context/
  search/
  prefilter/
  validation_top3/
    pdbs/
  final_report/

```

Key files include:

* `context/stage10_context.json`
* `search/stage10_search_candidates.csv`
* `prefilter/stage10_top3.csv`
* `validation_top3/stage08_structural_fasttrack_summary.csv`
* `validation_top3/pdbs/candidate_1.pdb`
* `final_report/stage10_report_summary.json`
* `final_report/stage10_report.md`

---

## How to interpret Stage 10

Stage 10 must be evaluated against the strict falsification baseline of the previous stages:

### The ultimate architectural conclusion

If the top candidates emerging from the Stage 10d wrapper achieve a passing global pLDDT (>= 70.0) and a stable RMSD (<= 3.5 Å), it proves that explicitly conditioning the generation process on 3D physics solves the catastrophic collapse pattern observed in sequence-only models.

If they still struggle, it indicates that the wild-type chassis is highly brittle, and the mutation budget must be localized even further, but the *methodology* of inverse-folding remains the correct paradigm.

---

## Current limitations

By Stage 10 the project still does **not** claim:

* wet-lab validation (phage synthesis and plaque assays)
* arbitrary-seed universal retargeting (de novo backbone generation)

Instead, the project now supports its sharpest claim:

* sequence-first models (Stage 07) hallucinate structure.
* proxy-guided models (Stage 09) cannot mathematically guarantee topological survival.
* **Inverse-folding models (Stage 10) represent the necessary, methodologically correct approach to scaffold-preserving viral receptor engineering.**

---

## Summary

By Stage 10, PhageForge tells a complete, closed-loop computational research story:

* earlier stages built a validity-aware host-retargeting framework (Stages 01-06)
* Stage 07 proved we can generate functional-looking sequences locally
* Stage 08 falsified the assumption that ranking success equals structural plausibility
* Stage 09 attempted to fix this with sequence proxies, but proved proxies are physically insufficient
* Stage 10 successfully inverted the generative paradigm, proving that **explicit 3D structure-conditioned generation** is required to safely engineer highly complex biological machinery.