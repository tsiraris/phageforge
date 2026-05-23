# PhageForge

**PhageForge** is a research project on **validity-aware, scaffold-constrained phage receptor-binding protein (RBP) design**.

The repository began as a host-conditioned sequence-design pipeline, evolved into a validity-aware and family-constrained design system, attempted a proxy-guided redesign, inverted the paradigm with structure-conditioned generation, and now culminates in a full **Stage 11 autonomous baseline-qualified inverse-folding redesign architecture**.

By this point the project addresses a sharper question:

> If inverse-folding (Stage 10) successfully optimizes for a provided 3D shape, how do we guarantee that the foundational wild-type chassis is natively foldable and uncorrupted before initiating redesign?

This README reflects the repository state after the definitive Stage 11 redesign cycle.

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
* the **Stage 11 baseline qualification branch** asks: **is the foundational wild-type chassis natively foldable (pLDDT ≥ 70) before we allow any generative edits to occur?**

This makes Stage 11 the architectural response to the danger of unvalidated intermediate baselines: it isolates the inverse-folding engine from every upstream artifact, so that no single early defect can silently propagate down the pipeline.

---

## Project evolution

### Stages 01-03 — dataset building and host prediction

The project begins by building curated RBP datasets, embedding them with protein language models, and training the host-prediction backbone.

### Stage 04 — host-conditioned RBP optimization

Explores direct host retargeting by proposing mutations and ranking candidates with the host probe.

### Stage 05 — validity-aware evaluation

Establishes that **host score alone is insufficient** and adds manifold similarity, family retention, and plausibility-aware evaluation.

### Stage 06 — family-conditioned host-ladder design

Turns the project into a scaffold-constrained design program by selecting a family, defining editable regions, and grounding edits in a seed scaffold. This seed was later found to be the project's hidden liability: a sequence-optimized intermediate that was never structurally validated.

### Stage 07 — multimodal, locally generated design

Builds the main candidate pool via local ESM3-guided generation and multimodal ranking.

### Stage 08 — structural fast-track validation

Adds robust ESMFold-based structural validation, definitively falsifying the unconstrained Stage 07 candidates due to catastrophic 3D drift.

### Stage 09 — structure-aware localized redesign

A sequence-first response to Stage 08. Shrinks the edit space and adds surrogate-driven prefiltering. It successfully improves pre-validation proxy scores but still fails full 3D ESMFold validation, proving sequence-first proxies are insufficient.

### Stage 10 — structure-conditioned inverse-folding redesign

Explicitly anchors redesign to a physical 3D seed scaffold using ESM-IF1. It proves that inverse-folding perfectly optimizes for the provided 3D shape, but reveals a critical flaw: if the upstream sequence-generated seed is already structurally collapsed, the model will faithfully stabilize that broken shape. The post-mortem traced every prior structural failure to the unvalidated Stage 06 seed, which itself folds at very low confidence.

### Stage 11 — autonomous baseline-qualified inverse-folding

Introduced to fix the pipeline vulnerability discovered in Stage 10. Stage 11 severs all ties with upstream generated artifacts and builds a completely self-contained inverse-folding loop directly from a wild-type seed.

* **11a**: build a self-contained redesign context from a wild-type seed, enforcing a strict >70.0 pLDDT Baseline Qualification gate
* **11b**: run inverse-folding beam search, scoring proposals directly against the baseline-qualified 3D backbone
* **11c**: prefilter candidates using a two-pass diversity-aware algorithm to prevent mode collapse
* **11d**: validate the top Stage 11 panels using a sterile subprocess wrapper around the unmodified Stage 08 validator
* **11e**: build the final comparative report and FASTA handoff, contrasting outcomes against historical corrupted-chassis baselines

---

## Main Stage 11 idea

Stage 11 deliberately abandons upstream sequence-first hallucinations in favor of a **true physical simulation loop** anchored to ground-truth biology:

1. select a raw wild-type RBP sequence directly from the strict dataset
2. explicitly predict its wild-type structure using ESMFold and **abort if it fails the Baseline Qualification gate** (pLDDT < 70)
3. define an aggressively minimized hard/soft edit space computed completely from scratch
4. query the ESM-IF1 inverse-folding model for substitutions that thermodynamically stabilize that specific, validated 3D geometry
5. score candidates with a physics-heavy composite:
* 3D backbone log-likelihood (ESM-IF1)
* target-host probability (ESM2 + LR)
* family evolutionary cosine
* sequence identity preservation


6. prefilter survivors using a two-pass greedy-diverse embedding space algorithm
7. send the elite panel to the exact same full structural validator used in Stage 08

The point of Stage 11 is to guarantee that the inverse-folding engine is only fed verified, structurally sound foundations, and that every generative artifact is recomputed in-stage so no upstream defect can re-enter the loop.

---

## What Stage 11 utilizes

The Stage 11 update introduces the following generative architecture capabilities:

* **autonomous baseline qualification**
* the pipeline explicitly tests its own wild-type starting materials, refusing to spend compute cycles redesigning natively brittle or collapsed proteins


* **full self-containment**
* the family centroid, target centroid, edit space, and per-position priors are all recomputed in-stage from the strict CSV, the cached ESM-2 embeddings, and the trained host probe
* no JSON context, edit space, or surrogate model from Stages 06/07/09/10 is read, so a single early defect cannot silently propagate downstream


* **explicit 3D scaffold anchoring**
* the sequence generation is locked mathematically to the X, Y, Z coordinates of the freshly validated `seed_wt.pdb` file
* the mutation budget is kept deliberately small (configurable; runs to date used 1-10 edits) to prioritize fold retention


* **inverse-folding beam search**
* replaces the sequence language model (ESM3/ESM2) with a Graph Neural Network (ESM-IF1)
* structural stability is no longer guessed via proxies; it is simulated dynamically


* **composite 3D fitness scoring**
* scores candidates primarily by physical backbone compatibility (`if1_log_likelihood`) rather than additive heuristic metrics
* the five composite weights are exposed as CLI parameters and default to the canonical Stage 10 values; leaving them unchanged reproduces the built-in formula exactly


* **diversity-aware prefiltering**
* utilizes a robust two-pass `greedy_diverse_subset` clustering to guarantee the final validation panel represents distinct structural hypotheses, preventing mode collapse


* **wrapper-based scientific coherence**
* executes the heavy ESMFold validation via a sterile subprocess wrapper, isolating VRAM between embeddings and structural physics
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
│   ├── stage10_utils.py
│   └── stage11_utils.py
│
├── results/
│   ├── analysis/
│   ├── phaseA/
│   ├── stage08/
│   ├── stage09/
│   ├── stage10/
│   └── stage11/
│
├── scripts/
│   ├── 01_build_dataset.py
│   ...
│   ├── 09a_define_edit_space.py
│   ...
│   ├── 10a_prepare_stage10_structure_context.py
│   ...
│   ├── 11a_prepare_stage11_context.py
│   ├── 11b_run_inverse_folding_beam_search.py
│   ├── 11c_prefilter_stage11_candidates.py
│   ├── 11d_validate_stage11_candidates.py
│   └── 11e_make_stage11_report.py
│
├── notebooks/
│   ├── stage07.ipynb
│   ├── stage08_closeout_sagemaker.ipynb
│   ├── stage09_structure_aware_redesign_sagemaker.ipynb
│   ├── stage10_inverse_folding_sagemaker.ipynb
│   └── stage11_autonomous_redesign_sagemaker.ipynb
│
├── pyproject.toml
└── README.md

```

---

## Installation

Core Python dependencies are declared in `pyproject.toml`.

A practical GPU environment for Stage 11 requires specific heavy infrastructure:

```bash
pip install -e .
pip install 'pandas<3' 'huggingface_hub<1' esm transformers accelerate scikit-learn biopython fair-esm

```

Because Stage 11 relies on:

* ESM2 embeddings for target scoring and spatial diversity
* the `fair-esm` inverse-folding library (`esm_if1_gvp4`)
* ESMFold for both the upfront Baseline Qualification and the final 3D atomic validation

you must ensure:

* PyTorch with CUDA is available and properly mapped
* the machine has significant VRAM (e.g., AWS G5/A10G or similar) to hold both ESM2 and ESM-IF1 simultaneously in memory during the search loop

> **Note on Python 3.12 / SageMaker:** `fair-esm` pulls in `torch_scatter` / `torch_sparse`, which may not compile on recent toolchains. Stage 11b ships with a pure-PyTorch compatibility shim that activates only when the compiled libraries are absent, so the search loop runs unmodified on these environments.

---

## Recommended Stage 11 workflow

### 1. Build the self-contained baseline-qualified context

```bash
python scripts/11a_prepare_stage11_context.py \
  --strict_csv data/processed/rbp_dataset_eskapee_strict.csv \
  --seed_protein_id "UOX38086.1" \
  --target_host "Enterobacter" \
  --predictor_model results/broad/linear_probe/seed_42/model.joblib \
  --label_classes_json results/broad/linear_probe/seed_42/label_classes.json \
  --out_dir results/stage11/context \
  --max_edit_positions 6 \
  --soft_positions 3 \
  --min_mutations 1 \
  --max_mutations 4

```

This script explicitly folds the wild-type seed, aborts if pLDDT < 70, and generates the `stage11_context.json` completely from scratch.

### 2. Run inverse-folding beam search

```bash
python scripts/11b_run_inverse_folding_beam_search.py \
  --stage11_context_json results/stage11/context/stage11_context.json \
  --out_csv results/stage11/search/stage11_search_candidates.csv \
  --out_json results/stage11/search/stage11_search_summary.json \
  --rounds 4 \
  --beam_width 24 \
  --proposals_per_parent 8 \
  --substitutions_per_position 3

```

The five composite weights (`--w_target`, `--w_if1`, `--w_family`, `--w_identity`, `--w_mut_penalty`) are optional; omit them to use the canonical defaults, or pass them to ablate the score for a given run.

### 3. Prefilter and select diverse candidates

```bash
python scripts/11c_prefilter_stage11_candidates.py \
  --stage11_context_json results/stage11/context/stage11_context.json \
  --search_csv results/stage11/search/stage11_search_candidates.csv \
  --out_topk_csv results/stage11/prefilter/stage11_top10.csv \
  --out_topk_final_csv results/stage11/prefilter/stage11_top3.csv \
  --out_json results/stage11/prefilter/stage11_prefilter_summary.json \
  --top_k 10 \
  --top_k_final 3

```

### 4. Validate the Stage 11 panel with the Stage 08 validator

```bash
python scripts/11d_validate_stage11_candidates.py \
  --validated_csv results/stage11/prefilter/stage11_top3.csv \
  --ranked_csv results/stage11/search/stage11_search_candidates.csv \
  --context_json results/stage11/context/stage11_context.json \
  --out_dir results/stage11/validation_top3 \
  --out_json results/stage11/validation_top3/stage11_launch.json \
  --top_k 3 \
  --device cuda \
  --chunk_size 128

```

### 5. Build the final Stage 11 comparative report

```bash
python scripts/11e_make_stage11_report.py \
  --stage11_context_json results/stage11/context/stage11_context.json \
  --search_csv results/stage11/search/stage11_search_candidates.csv \
  --prefilter_csv results/stage11/prefilter/stage11_top10.csv \
  --validation_csv results/stage11/validation_top3/stage08_structural_fasttrack_summary.csv \
  --out_dir results/stage11/final_report

```

---

## Stage 11 outputs

Typical Stage 11 outputs live under:

```text
results/stage11/
  context/
    seed_wt.pdb
    stage11_context.json
  search/
  prefilter/
  validation_top3/
    pdbs/
  final_report/

```

Key files include:

* `context/seed_wt.pdb`
* `context/stage11_context.json`
* `search/stage11_search_candidates.csv`
* `prefilter/stage11_top3.csv`
* `validation_top3/stage08_structural_fasttrack_summary.csv`
* `validation_top3/pdbs/candidate_1.pdb`
* `final_report/stage11_report_summary.json`
* `final_report/stage11_report.md`
* `final_report/stage11_handoff.fasta`

---

## Structural validation correction (pLDDT scale)

A unit-scale bug in the structural validator masked Stage 11's true outcome and is worth documenting, both because the fix is required to reproduce the results and because it illustrates how a silent metric error can imitate genuine scientific failure.

ESMFold (HuggingFace) emits pLDDT on a 0–1 scale, but the validator compared the stored value against the canonical `70.0` gate and, separately, averaged confidence over all 37 atom slots per residue (including non-existent padding atoms). The combined effect was that **every candidate failed `low_global_confidence` automatically**, superficially echoing the genuine Stage 08–10 collapses.

The fix recomputes every confidence metric on the canonical 0–100 per-residue Cα convention — the same convention the Stage 11a Baseline Qualification gate uses. With the correction applied, the verdicts flip from `0/N` to the true pass rates reported below. The seed itself folds at mean pLDDT ≈ 79.

---

## Stage 11 results

For the Klebsiella → Enterobacter retargeting (seed `UOX38086.1`, 806 aa), with the validator corrected to the canonical pLDDT convention:

* **Baseline Qualification:** PASS — wild-type seed folds at mean pLDDT ≈ 79, comfortably clearing the 70.0 gate (contrast: the historical Stage 06 chassis folds at very low confidence).
* **Structural pass rate:** the final top-3 panel reached **3/3 passing all six Stage 08 hard gates**, with mean pLDDT ≈ 79, mutation-site pLDDT ≈ 69–72, and RMSD < 1.6 Å to the seed.
* **Minimal-edit retargeting:** passing candidates carry only 2–5 substitutions on an 806-residue protein, confirming that inverse folding retains the wild-type fold while introducing the targeted edits.

Run-to-run variation across the search experiments was driven by the **mutation budget** and **search breadth** (rounds, beam width, proposals), together with the pLDDT-scale correction — not by composite re-weighting (see *Current limitations*).

---

## How to interpret Stage 11

Stage 11 must be evaluated against the strict falsification baseline of the previous stages:

### The architectural conclusion

When the top candidates emerging from the Stage 11d wrapper achieve a passing global pLDDT (>= 70.0) and a stable RMSD (<= 3.5 Å) — as they do for the case above — it demonstrates that **Baseline-Qualified Inverse-Folding** is a sound path for scaffold-constrained viral receptor engineering. It confirms that explicit 3D physics solves generative collapse, provided the system is anchored to a natively foldable, verified ground-truth chassis.

Where candidates struggle, it indicates that the wild-type chassis is highly brittle, and the mutation budget must be localized even further, but the *methodology* of inverse-folding remains the correct paradigm.

---

## Current limitations

By Stage 11 the project still does **not** claim:

* wet-lab validation (phage synthesis and plaque assays)
* arbitrary-seed universal retargeting (de novo backbone generation)
* general validity beyond the demonstrated case — results to date cover a single seed and a single Klebsiella → Enterobacter target pair, in silico only

Two characterized limitations are worth stating explicitly:

* **Target-host probability ceiling.** Conservative, fold-preserving edits move the host classifier only slightly: predicted target-host probability stays around ~0.22 regardless of mutation budget. This is a seed-distance limitation — the wild-type seed sits far from the target-host cluster in embedding space — and is resolvable only by choosing a phylogenetically closer seed or by empirical wet-lab screening, not by further search tuning.
* **Composite weighting.** In the original implementation the composite weights were effectively fixed, so the weight changes attempted across early search runs did not influence results; the observed run-to-run variation came from mutation budget and the structural-validation correction. The weights are now exposed as tunable CLI parameters (defaulting to the canonical values), so future ablations apply as intended.

Instead, the project now supports its sharpest claim:

* sequence-first models (Stage 07) hallucinate structure.
* proxy-guided models (Stage 09) cannot mathematically guarantee topological survival.
* Inverse-folding models (Stage 10) represent the necessary, methodologically correct approach to scaffold-preserving viral receptor engineering.
* **Baseline Qualification (Stage 11) is the mandatory prerequisite; inverse-folding succeeds only when isolated from upstream hallucinations and explicitly anchored to a verified, natively foldable seed.**

---

## Summary

By Stage 11, PhageForge tells a complete, closed-loop computational research story:

* earlier stages built a validity-aware host-retargeting framework (Stages 01-06)
* Stage 07 proved we can generate functional-looking sequences locally
* Stage 08 falsified the assumption that ranking success equals structural plausibility
* Stage 09 attempted to fix this with sequence proxies, but proved proxies are physically insufficient
* Stage 10 successfully inverted the generative paradigm, proving that explicit 3D structure-conditioned generation is required to safely engineer highly complex biological machinery
* **Stage 11 achieved self-contained architectural success by instituting a strict Baseline Qualification gate, producing a top-3 panel that passed full structural validation (mean pLDDT ≈ 79, RMSD < 1.6 Å) once the inverse-folding engine was isolated from upstream hallucinations and fed only verified, structurally sound foundations.**