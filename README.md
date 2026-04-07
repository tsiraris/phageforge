# PhageForge

**PhageForge** is a research project on **validity-aware, scaffold-constrained phage receptor-binding protein (RBP) design**.  
The repository began as a host-conditioned sequence-design pipeline, evolved into a validity-aware and family-constrained design system, and now includes a full **Stage 09 structure-aware redesign attempt**.

By this point the project addresses a sharper question:

> If post-hoc structural validation fails the best generated candidates, can a stricter, structure-aware, localized redesign stage recover candidates that are more likely to preserve scaffold integrity?

This README reflects the repository state after the Stage 09 redesign cycle.

---

## Scientific framing

PhageForge does **not** treat host retargeting as simple score maximization.
It also does **not** assume that better sequence-space or proxy-structure scores automatically imply structural viability.

By Stage 09 the project is organized around a sequential research logic:

- the **sequence / host-transfer branch** asks: **does a candidate move toward the target host manifold?**
- the **family / scaffold branch** asks: **does it remain near the validated seed and family neighborhood?**
- the **Stage 08 structural branch** asks: **does full structural validation support the candidate as a plausible scaffold-preserving redesign?**
- the **Stage 09 redesign branch** asks: **can stricter localized search and structural prefiltering improve this outcome before full validation?**

This makes Stage 09 a **response to falsification**, not just another generation stage.

---

## Project evolution

### Stages 01-03 — dataset building and host prediction
The project begins by building curated RBP datasets, embedding them with protein language models, and training the host-prediction backbone.

### Stage 04 — host-conditioned RBP optimization
The next stage explores direct host retargeting by proposing mutations and ranking candidates with the host probe.

### Stage 05 — validity-aware evaluation
Stage 05 establishes that **host score alone is insufficient** and adds manifold similarity, family retention, and plausibility-aware evaluation.

### Stage 06 — family-conditioned host-ladder design
Stage 06 turns the project into a scaffold-constrained design program by selecting a family, defining editable regions, and grounding later edits in a validated seed scaffold.

### Stage 07 — multimodal, locally generated design
Stage 07 builds the main candidate pool:
- enriched context preparation
- local ESM3-guided generation
- structure-aware reranking
- optional tissue-context scoring
- multimodal ranking and diversity-aware shortlist selection
- validation-ready exports

### Stage 08 — structural fast-track validation
Stage 08 adds robust ESMFold-based structural validation and shows whether the top shortlisted candidates survive stricter fold-retention checks.

### Stage 09 — structure-aware localized redesign
Stage 09 is introduced after the Stage 08 failure signal.  
It adds:

- **09a**: build a stricter Stage 09 edit space around the selected seed
- **09b**: build a structural-surrogate dataset from Stage 07 + Stage 08 outcomes
- **09c**: train or configure a structural surrogate
- **09d**: run localized search inside the seed-centered edit space
- **09e**: prefilter candidates with structural-risk and seed-drift constraints
- **09f**: validate the top Stage 09 candidates using the same Stage 08 structural validator
- **09g**: summarize the redesign attempt and compare it to the earlier failure regime

---

## Main Stage 09 idea

Stage 09 deliberately replaces broad or weakly constrained proposal logic with a **tight seed-local redesign loop**:

1. start from the selected Stage 07 seed scaffold
2. define a compact hard/soft edit space
3. keep most of the sequence frozen
4. search with localized substitutions only
5. score candidates with:
   - target-host probability
   - seed similarity
   - family similarity
   - local substitution guidance
   - surrogate structural risk
6. prefilter with conservative structural-proxy thresholds
7. send only the best candidates to the exact same full structural validator used in Stage 08

The point of Stage 09 is not simply to score better candidates.  
It is to test whether **upstream structural constraints** can fix the failure pattern exposed in Stage 08.

---

## What Stage 09 adds

The Stage 09 update adds several important capabilities:

- **explicit edit-space restriction**
  - positions are divided into hard edit positions, soft edit positions, and frozen positions
  - this makes scaffold preservation an explicit design constraint rather than an indirect hope

- **proposal-level substitution priors**
  - substitutions are informed by family support and target-host residue preferences
  - edit proposals stay closer to the local viable manifold

- **localized beam-style search**
  - Stage 09 replaces broad generation with a compact round-based local search
  - each round keeps only a beam of the strongest candidates

- **structural surrogate integration**
  - a lightweight surrogate or rule-based bundle estimates structural risk before expensive folding
  - this lets Stage 09 penalize obviously risky proposals earlier

- **structural prefiltering before full validation**
  - candidates must satisfy predicted pLDDT, predicted RMSD, structural-risk, identity, and editable-region constraints before final validation

- **reuse of Stage 08 as the final falsification layer**
  - Stage 09 does not invent a new validator
  - it intentionally reuses the fixed Stage 08 structural validator so comparisons remain scientifically coherent

- **final Stage 09 reporting**
  - a compact report summarizes:
    - number of search candidates
    - number of prefilter survivors
    - number of validated candidates
    - structural pass count / pass rate
    - mean pLDDT, mutation-site pLDDT, and RMSD

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
│   └── stage09_utils.py
│
├── results/
│   ├── analysis/
│   ├── phaseA/
│   ├── stage07/
│   ├── stage08/
│   └── stage09/
│
├── scripts/
│   ├── 01_build_dataset.py
│   ├── 01b_filter_strict_rbps.py
│   ├── 01c_structural_filter_rbps.py
│   ├── 01d_build_structural_plus_rbps.py
│   ├── 02_embed_rbps.py
│   ├── 03_train_phi_mlp.py
│   ├── 03b_linear_probe.py
│   ├── 03c_run_benchmark.py
│   ├── 03d_make_benchmark_report.py
│   ├── 04_optimize_rbp_for_host.py
│   ├── 04b_summarize_design_run.py
│   ├── 05_compute_validity_metrics.py
│   ├── 05b_rank_validated_candidates.py
│   ├── 05c_make_validity_report.py
│   ├── 06a_select_phaseA_family.py
│   ├── 06b_optimize_family_constrained.py
│   ├── 06c_pick_phaseA_followup_seed.py
│   ├── 07a_prepare_stage07_design_context.py
│   ├── 07b_generate_rbps_with_esm3.py
│   ├── 07c_score_structure_aware_candidates.py
│   ├── 07d_build_tissue_context_embeddings.py
│   ├── 07e_rank_multimodal_candidates.py
│   ├── 07f_make_stage07_report.py
│   ├── 07g_export_stage07_panel.py
│   ├── 07h_validate_stage07_candidates.py
│   ├── 08a_structural_fasttrack_validation.py
│   ├── 08b_make_final_closeout.py
│   ├── 09a_define_edit_space.py
│   ├── 09b_build_structure_surrogate_dataset.py
│   ├── 09c_train_structure_surrogate.py
│   ├── 09d_localized_search.py
│   ├── 09e_structural_prefilter.py
│   ├── 09f_validate_stage09_candidates.py
│   └── 09g_make_stage09_report.py
│
├── notebooks/
│   ├── 09_stage08_closeout_sagemaker.ipynb
│   └── 10_stage09_structure_aware_redesign_sagemaker.ipynb
│
├── pyproject.toml
└── README.md
```

---

## Installation

Core Python dependencies are declared in `pyproject.toml`.

A practical GPU environment for Stage 09 is:

```bash
pip install -e .
pip install 'pandas<3' 'huggingface_hub<1' esm transformers accelerate scikit-learn biopython
```

Because Stage 09 reuses:
- ESM2 embeddings for target scoring and diversity selection
- the host predictor from the benchmark / linear-probe stage
- ESMFold for final validation

you should ensure:
- PyTorch with CUDA is available
- the predictor path and label-order JSON match the embedding backbone used at inference time
- the machine has enough RAM/GPU memory for the structural validator

---

## Recommended Stage 09 workflow

### 1. Build the stricter edit space

```bash
python scripts/09a_define_edit_space.py \
  --context_json results/stage07/context/stage07_context.base.json \
  --strict_csv data/processed/rbp_dataset_eskapee_strict.csv \
  --output_json results/stage09/edit_space/stage09_edit_space.json \
  --max_edit_positions 12 \
  --soft_buffer_positions 6 \
  --min_mutations 3 \
  --max_mutations 8
```

This produces a JSON artifact containing:
- hard edit positions
- soft edit positions
- frozen positions
- mutation-budget recommendations
- per-position substitution proposals

### 2. Optionally build a structural-surrogate dataset

```bash
python scripts/09b_build_structure_surrogate_dataset.py \
  --context_json results/stage07/context/stage07_context.base.json \
  --ranked_csv results/stage07/multimodal_rank/final_multimodal_ranked_candidates.csv \
  --structural_csv results/stage08/structural_fasttrack_top3/stage08_structural_fasttrack_summary.csv \
  --out_csv results/stage09/surrogate/stage09_surrogate_dataset.csv \
  --out_json results/stage09/surrogate/stage09_surrogate_dataset_summary.json
```

### 3. Optionally train or configure the surrogate

```bash
python scripts/09c_train_structure_surrogate.py \
  --dataset_csv results/stage09/surrogate/stage09_surrogate_dataset.csv \
  --summary_json results/stage09/surrogate/stage09_surrogate_dataset_summary.json \
  --out_model results/stage09/surrogate/stage09_surrogate_bundle.joblib
```

If structural supervision is too sparse, the script still writes a rule-based bundle.

### 4. Run localized search

```bash
python scripts/09d_localized_search.py \
  --context_json results/stage07/context/stage07_context.base.json \
  --edit_space_json results/stage09/edit_space/stage09_edit_space.json \
  --strict_csv data/processed/rbp_dataset_eskapee_strict.csv \
  --predictor_model results/broad/linear_probe/seed_42/model.joblib \
  --label_classes_json results/broad/linear_probe/seed_42/label_classes.json \
  --surrogate_model results/stage09/surrogate/stage09_surrogate_bundle.joblib \
  --out_csv results/stage09/search/stage09_search_candidates.csv \
  --out_json results/stage09/search/stage09_search_summary.json \
  --esm_model facebook/esm2_t33_650M_UR50D \
  --batch_size 2 \
  --max_aa 1022 \
  --rounds 4 \
  --beam_width 24 \
  --proposals_per_parent 96 \
  --max_mutations 6
```

Important note:
- the embedding backbone used during search must match the predictor’s expected embedding dimension

### 5. Prefilter before expensive validation

```bash
python scripts/09e_structural_prefilter.py \
  --search_csv results/stage09/search/stage09_search_candidates.csv \
  --search_meta_json results/stage09/search/stage09_search_summary.json \
  --out_csv results/stage09/prefilter/stage09_prefilter_top12.csv \
  --top_k 12 \
  --max_structural_risk 0.55 \
  --min_predicted_plddt 55 \
  --max_predicted_rmsd 4.5 \
  --min_sequence_identity 0.93
```

### 6. Validate the top Stage 09 panel with the same Stage 08 validator

Top-1:

```bash
python scripts/09f_validate_stage09_candidates.py \
  --prefilter_csv results/stage09/prefilter/stage09_prefilter_top12.csv \
  --context_json results/stage07/context/stage07_context.base.json \
  --out_dir results/stage09/validation_top1 \
  --top_k 1 \
  --device cuda \
  --chunk_size 64 \
  --num_recycles 1 \
  --resume
```

Top-3:

```bash
python scripts/09f_validate_stage09_candidates.py \
  --prefilter_csv results/stage09/prefilter/stage09_prefilter_top12.csv \
  --context_json results/stage07/context/stage07_context.base.json \
  --out_dir results/stage09/validation_top3 \
  --top_k 3 \
  --device cuda \
  --chunk_size 64 \
  --num_recycles 1 \
  --resume
```

### 7. Build the final Stage 09 report

```bash
python scripts/09g_make_stage09_report.py \
  --search_csv results/stage09/search/stage09_search_candidates.csv \
  --prefilter_csv results/stage09/prefilter/stage09_prefilter_top12.csv \
  --validation_csv results/stage09/validation_top3/stage08_structural_fasttrack_summary.csv \
  --out_dir results/stage09/final_report
```

If an earlier Stage 08 baseline summary is available, you can also pass:

```bash
--baseline_stage08_csv results/stage08/structural_fasttrack_top3/stage08_structural_fasttrack_summary.csv
```

---

## Stage 09 outputs

Typical Stage 09 outputs live under:

```text
results/stage09/
  edit_space/
  surrogate/
  search/
  prefilter/
  validation_top1/
  validation_top3/
  final_report/
```

Key files include:

- `edit_space/stage09_edit_space.json`
- `search/stage09_search_candidates.csv`
- `search/stage09_search_summary.json`
- `prefilter/stage09_prefilter_top12.csv`
- `validation_top1/stage08_structural_fasttrack_summary.csv`
- `validation_top3/stage08_structural_fasttrack_summary.csv`
- `final_report/stage09_final_candidate_table.csv`
- `final_report/stage09_summary.json`
- `final_report/stage09_report.md`

---

## How to interpret Stage 09

Stage 09 must be interpreted in two layers:

### Pre-validation layer
The search and prefilter stages can improve:
- edit locality
- mutation budget control
- sequence identity to seed
- predicted pLDDT / RMSD proxies
- structural-risk penalties

This is useful, but it is **not enough** on its own.

### Full-validation layer
The only decisive question is whether the top Stage 09 candidates survive the full structural validator.

If top-1 and top-3 still fail under:
- low global pLDDT
- low mutation-site confidence
- high RMSD to the selected seed

then the conclusion is:

> tighter proxy-guided localized search improved the design regime, but did not solve the true fold-retention problem.

This is still a strong result because it rules out another class of sequence-first redesign strategies.

---

## Current limitations

By Stage 09 the project still does **not** claim:
- wet-lab validation
- arbitrary-seed universal retargeting
- successful structure-preserving redesign under the current proxy-guided search alone
- de novo RBP invention

Instead, the project now supports a sharper claim:

- Stage 08 showed that post-hoc structural filtering was not enough
- Stage 09 showed that proxy-guided localized search was also not enough
- the next serious step is likely **true structure-conditioned redesign** rather than further tuning of the same sequence-first search strategy

---

## Recommended interpretation of the current result

The strongest honest reading of the Stage 09 update is:

- Stage 09 improved the *pre-validation* design regime
- but Stage 09 still failed the decisive *full structural* test
- therefore the bottleneck is deeper than ranking or local search tuning
- future progress likely requires:
  - inverse-folding-style redesign
  - structure-conditioned generation
  - or another explicitly scaffold-conditioned generative model

This makes Stage 09 scientifically valuable even as a negative result.

---

## Summary

By Stage 09, PhageForge tells a coherent and research-grade story:

- earlier stages built a validity-aware, scaffold-constrained host-retargeting framework
- Stage 07 produced a strong locally generated candidate pool
- Stage 08 falsified the assumption that ranking success implies structural plausibility
- Stage 09 responded with a tighter structure-aware redesign attempt
- the final result showed that **proxy-constrained sequence redesign still does not recover true structural viability**

That is not a failure of the project.  
It is a clear, honest, and technically grounded conclusion about what class of methods is insufficient—and what class of methods should come next.
