# Stage 11 results — pLDDT scale correction

## What was wrong
Every `stage08_pass` was `False` with reason `low_global_confidence` (and, in runs 1–2,
`low_mutation_site_confidence`). This was **not** a science failure — it was a unit bug in the
Stage 08a validator:

1. **Global pLDDT stored on the wrong scale.** `esmfold_mean_plddt` was written on a 0–1 scale
   (e.g. 0.71) but compared against the gate threshold of `70.0`, so it always failed.
2. **Global pLDDT averaged over non-existent atom slots.** The validator used
   `outputs["plddt"].mean()` (all 37 atom slots per residue, including padding), which is ~8–10
   points lower than the canonical per-residue Cα pLDDT.
3. **Per-residue confidence fraction** (`mutation_site_confidence_ge70_fraction`) was 0.0 in runs
   1–2 for the same 0–1-vs-70 reason. (Run 3 had a partial live patch, so its fraction was already
   correct, but its global pLDDT was not.)

## How it was corrected
All pLDDT-derived metrics were recomputed **directly from the PDB Cα B-factors × 100** — the same
0–100 Cα convention the Stage 11a baseline gate uses (which reported the seed at ~81). RMSD,
glycine fraction, homopolymer run, and outside-editable fraction are scale-independent and were left
unchanged. The six `stage08_pass` rules and `stage08_decision_reason` were then re-evaluated, and the
CSV / JSON / Markdown reports (Stage 08a and Stage 11e) regenerated.

## Corrected outcome (seed folds at mean pLDDT ≈ 79)

| Run | top-10 pass | top-3 pass | verdict |
|---|---:|---:|---|
| 1 (`T101007Z`, 1–3 mut) | 7/10 | 1/3 | partially_supported |
| 2 (`2nd_run`, 3–6 mut)  | 9/10 | 2/3 | supported |
| 3 (`3rd_run`, 4–10 mut) | 7/10 | **3/3** | **supported** |

Headline: structure-conditioned inverse-folding on a structurally-qualified wild-type seed yields
RBP redesigns that retain the fold (mean pLDDT ≈ 79, RMSD < 1.6 Å) while carrying the targeted edits.
Predicted target-host probability stays ~0.22–0.23 (a known seed-distance limitation, unchanged by
this correction).

_Note: stale nested archives (`*.tar.gz`) and `.ipynb_checkpoints` were removed from the tree._
