# PhageForge

**PhageForge** is a research project on **validity-aware, scaffold-constrained phage receptor-binding protein (RBP) design**. The repository began as a host-conditioned sequence-design pipeline and now culminates in an upgraded **Stage 07** workflow that combines:

- **validated scaffold-family selection** from earlier stages
- **functionally prioritized, constrained RBP editing**
- **local gated ESM3-guided generation**
- **structure-aware reranking**
- **optional tissue-context scoring** from histopathology / omics / metadata
- **final multimodal candidate ranking with diversity-aware shortlist selection**

The central design question is:

> Can we retarget a validated phage RBP scaffold toward a new bacterial host while preserving family identity and biological plausibility, and optionally prioritize candidates for a specific tissue or infection context?

---

## Scientific framing

PhageForge does **not** treat host transfer as a simple score-maximization problem.
It also does **not** treat histopathology as a direct proxy for bacterial host range.

Instead, the project uses a strict division of roles:

- the **sequence / host-transfer branch** answers: **can this RBP plausibly move toward the target host manifold?**
- the **structure-aware branch** answers: **does the design remain plausible relative to the scaffold family and reference anchors?**
- the **optional tissue-context branch** answers: **is this candidate compatible with a tissue / infection context?**
- the **final Stage 07 reranker** combines these into one shortlist.

This is what makes the multimodal extension scientifically coherent rather than conceptually mixed.

---

## Project evolution

### Stages 01-03 — dataset building and host prediction
The early stages build curated RBP datasets and embed them with protein language models. These stages establish the host-prediction backbone that later design stages rely on.

### Stage 04 — host-conditioned RBP optimization
The repository then explores host retargeting by proposing mutations and ranking candidates by the host probe.

### Stage 05 — validity-aware evaluation
Stage 05 adds the crucial realization that **host score alone is not enough**. Candidates are evaluated with additional terms such as manifold similarity, family retention, novelty / naturalness, and corrected mutation accounting.

### Stage 06 — family-conditioned host-ladder design
Stage 06 turns the project into a scaffold-constrained design framework. It selects a canonical family, defines hotspot windows, uses target-anchor references, and demonstrates the trade-off between strict family preservation and functional movement toward a more distant host.

### Stage 07 — multimodal, locally generated, frontier RBP design
Stage 07 is the final extension of the project. It keeps the validated Stage 05/06 core fixed, then adds:

- **07a**: preparation of one self-contained Stage 07 design context
- **07b**: constrained candidate generation with **local gated ESM3**
- **07c**: structure-aware candidate reranking
- **07d**: optional tissue-context embedding construction
- **07e**: final multimodal ranking with diversity-aware panel selection
- **07f**: compact report generation
- **07g**: validation-panel export to CSV / FASTA / JSON

---

## Main Stage 07 idea

Stage 07 is intentionally **not** whole-genome generation and **not** unconstrained protein generation.

It performs **scaffold-constrained editing**:

1. start from a validated seed scaffold from Stage 06
2. score editable positions by functional weight rather than using naive uniform windows
3. define hotspot or structured-block editing regimes
4. generate multiple candidate RBPs with **local gated ESM3**
5. rerank them with structure-aware and optional tissue-aware terms
6. apply diversity-aware selection to produce a compact validation panel

This is the shortest realistic path to a scientifically defensible and impressive final portfolio piece.

---

## What changed in the upgraded Stage 07 workflow

The current repository version moves beyond the earlier Stage 07 implementation in several important ways:

- **local gated `esm3-open` is now the default practical backend**
  - this avoids Forge API daily-credit failures
  - it also keeps provenance clean and reproducible

- **ESM2 fallback can be disabled for provenance-complete runs**
  - final shortlisted candidates can now come from fully local ESM3 generation only

- **functional hotspot weighting is built into context preparation**
  - editable positions are ranked using family mutability and target-oriented priors
  - this replaces weaker “uniform window” logic for primary runs

- **structured-block regimes replace the old failing window modes**
  - contiguous high-value regions are now explored explicitly
  - this improves regime diversity without collapsing generation quality

- **multi-attempt local generation improves stability**
  - Stage 07 generation now retries guided local ESM3 sampling per sample
  - this removes the earlier regime failure mode where some runs returned no usable sequence output

- **diversity-aware ranking is now part of the final shortlist**
  - the final panel is not just top-scoring clones from one regime
  - it is selected from a stronger candidate pool while preserving regime diversity

- **validation-ready exports are now standard outputs**
  - top candidates are written directly to CSV / FASTA / JSON
  - a single Stage 07 bundle can be archived and downloaded cleanly

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
│   └── stage07_utils.py
│
├── results/
│   ├── analysis/
│   │   ├── validity/
│   │   └── validity_corrected/
│   ├── phaseA/
│   ├── phaseB/
│   └── stage07/
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
│   ├── 07b_train_conditional_generator.py
│   ├── 07b_generate_rbps_with_esm3.py
│   ├── 07c_score_structure_aware_candidates.py
│   ├── 07d_build_tissue_context_embeddings.py
│   ├── 07e_rank_multimodal_candidates.py
│   ├── 07f_make_stage07_report.py
│   └── 07g_export_stage07_panel.py
│
├── notebooks/
│   ├── run_stage07_with_upgraded_esm3_colab_v2.ipynb
│   └── run_stage07_with_stage071_upgrades.ipynb
│
├── pyproject.toml
└── README.md
```

---

## Installation

Core Python dependencies are declared in `pyproject.toml`.
For the upgraded Stage 07 Colab workflow, install the repo and the required model dependencies:

```bash
pip install -e .
pip install 'pandas<3' 'huggingface_hub<1' esm transformers accelerate scikit-learn biopython
```

### Hugging Face access for local gated ESM3

The recommended Stage 07 path uses the gated **`EvolutionaryScale/esm3-sm-open-v1`** release through local inference.
You therefore need:

1. a Hugging Face account
2. access approval for the gated model
3. a read token available in the environment

Example:

```bash
export HF_TOKEN=your_real_token_here
export HF_HUB_DISABLE_XET=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
```

---

## Recommended Stage 07 workflow

### 1. Build the enriched Stage 07 context
This stage packages the validated seed, editable region, family context, ranked functional positions, and structured editing windows into one reusable JSON.

```bash
python scripts/07a_prepare_stage07_design_context.py \
  --phaseA_plan_json results/phaseA/phaseA_plan.json \
  --phase06c_followup_summary_json results/phaseA/step2_seed/phaseA_followup_seed_summary.json \
  --strict_csv data/processed/rbp_dataset_eskapee_strict.csv \
  --target_host Acinetobacter \
  --output_json results/stage07/context/stage07_context.base.json
```

### 2. Generate candidates with local gated ESM3
This is the preferred Stage 07 generation mode. It performs constrained editing over the validated seed using hotspot-priority or structured-block regimes.

```bash
python scripts/07b_generate_rbps_with_esm3.py \
  --context_json results/stage07/context/stage07_context.base.json \
  --out_csv results/stage07/generation/all_generated_candidates.csv \
  --n_samples 32 \
  --temperature 0.7 \
  --top_k 5 \
  --esm3_backend local \
  --esm3_model esm3-open \
  --sampling_seed 42 \
  --esm3_num_steps 8 \
  --max_attempts_per_sample 3 \
  --esm3_error_fallback none
```

Notes:
- **local gated ESM3** is now the recommended default
- **ESM2 fallback is best left disabled** for provenance-complete final runs
- the strongest results currently come from **functional hotspot** and **structured-block** regimes rather than naive contiguous windows

### 3. Score candidates with structure-aware proxies

```bash
python scripts/07c_score_structure_aware_candidates.py \
  --context_json results/stage07/context/stage07_context.base.json \
  --generated_csv results/stage07/generation/all_generated_candidates.csv \
  --scored_csv results/stage07/structure_rerank/structure_scored_candidates.csv
```

### 4. Optionally build tissue-context embeddings
If you have pathology, omics, or metadata context, build the optional tissue branch.

```bash
python scripts/07d_build_tissue_context_embeddings.py \
  --input_csv data/processed/tissue_context_inputs.csv \
  --out_pt results/stage07/tissue_context/tissue_embeddings.pt \
  --out_csv results/stage07/tissue_context/tissue_metadata.csv
```

### 5. Final multimodal ranking

Without tissue context:

```bash
python scripts/07e_rank_multimodal_candidates.py \
  --generated_csv results/stage07/generation/all_generated_candidates.csv \
  --structure_scored_csv results/stage07/structure_rerank/structure_scored_candidates.csv \
  --diverse_top_k 5 \
  --per_regime_pool 3 \
  --out_csv results/stage07/multimodal_rank/final_multimodal_ranked_candidates.csv
```

With tissue context:

```bash
python scripts/07e_rank_multimodal_candidates.py \
  --generated_csv results/stage07/generation/all_generated_candidates.csv \
  --structure_scored_csv results/stage07/structure_rerank/structure_scored_candidates.csv \
  --tissue_embeddings_pt results/stage07/tissue_context/tissue_embeddings.pt \
  --tissue_metadata_csv results/stage07/tissue_context/tissue_metadata.csv \
  --diverse_top_k 5 \
  --per_regime_pool 3 \
  --out_csv results/stage07/multimodal_rank/final_multimodal_ranked_candidates_with_tissue.csv
```

### 6. Write the compact report

```bash
python scripts/07f_make_stage07_report.py \
  --context_json results/stage07/context/stage07_context.base.json \
  --ranked_csv results/stage07/multimodal_rank/final_multimodal_ranked_candidates.csv \
  --report_md results/stage07/report/stage07_report.md
```

### 7. Export the validation panel

```bash
python scripts/07g_export_stage07_panel.py \
  --ranked_csv results/stage07/multimodal_rank/final_multimodal_ranked_candidates.csv \
  --output_dir results/stage07/exports \
  --top_k 5
```

---

## Colab workflow

The current recommended notebook is:

```text
notebooks/run_stage07_with_stage071_upgrades.ipynb
```

It is designed to:

1. upload your zipped repo to Colab
2. unzip it inside `/content`
3. install the local package plus model dependencies
4. authenticate gated local ESM3 through Hugging Face
5. build the enriched Stage 07 context
6. run multiple hotspot and structured-block regimes
7. merge, deduplicate, rerank, and export a validation panel
8. write a single downloadable Stage 07 bundle

The older Colab notebook is still useful historically, but the upgraded notebook is the preferred path for final Stage 07 runs.

---

## Outputs produced by the upgraded Stage 07 workflow

Typical Stage 07 outputs live under:

```text
results/stage07/
  context/
  generation/
  structure_rerank/
  tissue_context/
  multimodal_rank/
  report/
  exports/
```

Key files include:

- `context/stage07_context.base.json`
- `generation/all_generated_candidates.csv`
- `generation/generation_regime_summary.csv`
- `structure_rerank/structure_scored_candidates.csv`
- `multimodal_rank/final_multimodal_ranked_candidates.csv`
- `multimodal_rank/top_validation_panel.csv`
- `report/stage07_report.md`
- `exports/top5_validation_panel.csv`
- `exports/top5_validation_panel.fasta`
- `exports/top5_validation_panel.json`
- `results/stage07_bundle.zip`

---

## Recommended final demo

The strongest Stage 07 demo remains a **single polished case study**, for example:

- source scaffold family rooted in **Klebsiella**
- stepwise transfer through **Enterobacter**
- hard target **Acinetobacter**
- optional comparison of rankings with and without tissue context

For each final candidate, report:

- host-transfer score
- family similarity
- structure-aware score
- tissue score (if used)
- final rank
- mutation summary / provenance

For the current upgraded Stage 07 workflow, the most useful practical handoff is:

- a **primary top-3 validation subset** for immediate downstream screening
- an **expanded top-5 panel** for broader structural / biological follow-up

---

## Current limitations

This repository does **not** yet claim:

- whole-genome phage design
- universal arbitrary-seed retargeting
- fully trained pathology foundation models from scratch
- experimentally validated wet-lab success

Stage 07 is intended as a **frontier but finishable** multimodal extension of a validated design core.

---

## Practical note on Stage 07 generation backends

`07b_generate_rbps_with_esm3.py` supports multiple practical modes, but the recommended final mode is now clear:

1. **local gated ESM3 (`esm3-open`)** — preferred
2. **local conditional generator** — optional experimental fallback
3. **ESM2 masked-LM fallback** — available, but typically disabled for final provenance-clean runs

This means the pipeline is still robust, but final reported candidates can now be generated under a cleaner and more reproducible local setup.

---

## Summary

PhageForge now tells one coherent story:

- earlier stages proved that **validity-aware, scaffold-constrained host retargeting** is possible
- the upgraded Stage 07 workflow extends that core into a **multimodal, locally generated, ESM3-enabled RBP design system**
- functional hotspot weighting, structured-block regimes, multi-attempt generation, and diversity-aware ranking now make the final shortlist more stable and more scientifically defensible
- the final output is not just a mutated sequence, but a **ranked shortlist of biologically plausible candidates** ready for downstream structural and biological screening
