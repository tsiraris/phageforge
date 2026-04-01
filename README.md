# PhageForge

**PhageForge** is a research project on **validity-aware, scaffold-constrained phage receptor-binding protein (RBP) design**. The repository began as a host-conditioned sequence-design pipeline and now culminates in **Stage 07**, a multimodal extension that combines:

- **validated scaffold-family selection** from earlier stages
- **ESM-guided / ESM3-guided constrained RBP editing**
- **structure-aware reranking**
- **optional tissue-context scoring** from histopathology / omics / metadata
- **final multimodal candidate ranking**

The central design question is:

> Can we retarget a validated phage RBP scaffold toward a new bacterial host while preserving family identity and biological plausibility, and optionally prioritize candidates for a specific tissue or infection context?

---

## Scientific framing

PhageForge does **not** treat host transfer as a simple score maximization problem.
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

### Stage 07 — multimodal, tissue-aware, frontier RBP design
Stage 07 is the final extension of the project. It keeps the validated Stage 05/06 core fixed, then adds:

- **07a**: preparation of one self-contained Stage 07 design context
- **07b**: constrained candidate generation with **ESM3** (or a local / ESM2 fallback)
- **07c**: structure-aware candidate reranking
- **07d**: optional tissue-context embedding construction
- **07e**: final multimodal ranking
- **07f**: compact report generation

---

## Main Stage 07 idea

Stage 07 is intentionally **not** whole-genome generation and **not** unconstrained protein generation.

It performs **scaffold-constrained editing**:

1. start from a validated seed scaffold from Stage 06
2. restrict edits to a hotspot-derived or window-defined region
3. generate multiple candidate RBPs with ESM3 or a fallback generator
4. rerank them with structure-aware and optional tissue-aware terms
5. produce final shortlists for one target host case study

This is the shortest realistic path to a scientifically defensible and impressive final portfolio piece.

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
│   └── eval/
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
│   └── 07f_make_stage07_report.py
│
├── notebooks/
│   └── run_stage07_with_esm3_colab.ipynb
│
├── pyproject.toml
└── README.md
```

---

## Requirements

Core Python dependencies are declared in `pyproject.toml`.
For Stage 07 with ESM3 Forge, install the repo and the ESM client:

```bash
pip install -e .
pip install esm tabulate
```

If you want to use the ESM3 Forge API, set your token in the environment:

**Windows CMD**
```bat
set ESM_API_KEY=your_real_key_here
```

**PowerShell**
```powershell
$env:ESM_API_KEY="your_real_key_here"
```

---

## Recommended Stage 07 workflow

### 1. Build the Stage 07 context
This stage packages the validated seed, editable window, family context, and optional tissue metadata into one reusable JSON.

```bash
python scripts/07a_prepare_stage07_design_context.py   --phaseA_plan_json results/phaseA/phaseA_plan.json   --phase06c_followup_summary_json results/phaseA/step2_seed/phaseA_followup_seed_summary.json   --strict_csv data/processed/rbp_dataset_eskapee_strict.csv   --target_host Acinetobacter   --output_json results/stage07/context/stage07_context.json
```

### 2. Generate candidates with ESM3
This is the preferred Stage 07 generation mode. It performs constrained editing over the validated seed window.

```bash
python scripts/07b_generate_rbps_with_esm3.py   --context_json results/stage07/context/stage07_context.json   --out_csv results/stage07/generation/all_generated_candidates.csv   --n_samples 64   --temperature 0.7   --top_k 5   --use_esm3_api   --esm3_model esm3-medium-2024-08   --esm3_num_steps 8   --sampling_seed 42
```

Fallback modes are also supported:
- **local conditional generator checkpoint** from `07b_train_conditional_generator.py`
- **ESM2 masked-LM fallback** if neither ESM3 nor a local checkpoint is used

### 3. Score candidates with structure-aware proxies

```bash
python scripts/07c_score_structure_aware_candidates.py   --context_json results/stage07/context/stage07_context.json   --generated_csv results/stage07/generation/all_generated_candidates.csv   --reference_embeddings_pt data/processed/strict/esm2_embeddings.pt   --reference_index_csv data/processed/strict/esm2_embeddings_index.csv   --scored_csv results/stage07/structure_rerank/structure_scored_candidates.csv
```

### 4. Optionally build tissue-context embeddings
If you have pathology, omics, or metadata context, build the optional tissue branch.

```bash
python scripts/07d_build_tissue_context_embeddings.py   --input_csv data/processed/tissue_context_inputs.csv   --out_pt results/stage07/tissue_context/tissue_embeddings.pt   --out_csv results/stage07/tissue_context/tissue_metadata.csv
```

### 5. Final multimodal ranking

Without tissue context:

```bash
python scripts/07e_rank_multimodal_candidates.py   --generated_csv results/stage07/generation/all_generated_candidates.csv   --structure_scored_csv results/stage07/structure_rerank/structure_scored_candidates.csv   --out_csv results/stage07/multimodal_rank/final_multimodal_ranked_candidates.csv
```

With tissue context:

```bash
python scripts/07e_rank_multimodal_candidates.py   --generated_csv results/stage07/generation/all_generated_candidates.csv   --structure_scored_csv results/stage07/structure_rerank/structure_scored_candidates.csv   --tissue_embeddings_pt results/stage07/tissue_context/tissue_embeddings.pt   --tissue_metadata_csv results/stage07/tissue_context/tissue_metadata.csv   --out_csv results/stage07/multimodal_rank/final_multimodal_ranked_candidates_with_tissue.csv
```

### 6. Write the compact report

```bash
python scripts/07f_make_stage07_report.py   --context_json results/stage07/context/stage07_context.json   --ranked_csv results/stage07/multimodal_rank/final_multimodal_ranked_candidates.csv   --report_md results/stage07/report/stage07_report.md
```

---

## Colab workflow

A ready-to-run notebook is provided at:

```text
notebooks/run_stage07_with_esm3_colab.ipynb
```

It is designed for the exact workflow:

1. upload your zipped repo to Colab
2. unzip it inside `/content`
3. install the local package plus `esm`
4. set the Forge API key
5. run Stage 07 end to end with default paths from the current repository layout

---

## Outputs produced by Stage 07

Typical Stage 07 outputs live under:

```text
results/stage07/
  context/
  generation/
  structure_rerank/
  tissue_context/
  multimodal_rank/
  report/
```

Key files include:
- `context/stage07_context.json`
- `generation/all_generated_candidates.csv`
- `structure_rerank/structure_scored_candidates.csv`
- `multimodal_rank/final_multimodal_ranked_candidates.csv`
- `report/stage07_report.md`

---

## Recommended final demo

The strongest Stage 07 demo is still a **single polished case study**, for example:

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

---

## Current limitations

This repository does **not** yet claim:
- whole-genome phage design
- universal arbitrary-seed retargeting
- fully trained pathology foundation models from scratch
- experimentally validated wet-lab success

Stage 07 is intended as a **frontier but finishable** multimodal extension of a validated design core.

---

## Practical note on Stage 07b

`07b_generate_rbps_with_esm3.py` now supports three generation backends:

1. **ESM3 Forge API** (preferred)
2. **local conditional generator** trained by `07b_train_conditional_generator.py`
3. **ESM2 masked-LM fallback**

This makes the pipeline robust: it can run with your API key, but it also remains executable without remote access.

---

## Summary

PhageForge now tells one coherent story:

- earlier stages proved that **validity-aware, scaffold-constrained host retargeting** is possible
- Stage 07 upgrades that core into a **multimodal, tissue-aware, ESM3-enabled RBP design system**
- the final output is not just a mutated protein sequence, but a **ranked shortlist of biologically plausible candidates** for a specific host and, optionally, a specific tissue / infection context
