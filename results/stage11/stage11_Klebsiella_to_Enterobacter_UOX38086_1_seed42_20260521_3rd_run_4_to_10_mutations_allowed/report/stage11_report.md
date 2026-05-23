# Stage 11 — Universal Structure-Conditioned RBP Redesign Report

## Why Stage 11 exists

Stages 07–10 all failed ESMFold validation. The Stage 10 post-mortem traced this to a corrupted Stage 06 seed chassis (`round8_cand473` itself folds at very low pLDDT). Stage 11 mandates a Baseline Qualification gate (default pLDDT ≥ 70.0 on the wild-type seed) before any redesign compute is spent, and recomputes every family/target/edit-space artifact in-stage so no failed downstream artifact can re-enter the loop.

## Headline outcome

- **Verdict:** `supported`
- **Interpretation:** 3/3 Stage 11 candidates passed all six 08a hard thresholds. The Baseline Qualification hypothesis is supported: anchoring inverse-folding on a structurally healthy wild-type seed yields candidates that retain the fold (high pLDDT, low RMSD) while carrying the targeted edits.
- **Pass rate:** 100%
- **Best Stage 11 mean pLDDT:** 79.04
- **Best Stage 11 RMSD to seed:** 1.079 Å

## Baseline Qualification

- **Gate result:** PASS
- **Seed mean pLDDT:** 81.057
- **Threshold:** 70.000

## Core redesign setup

- **Source host:** Klebsiella
- **Target host:** Enterobacter
- **Seed protein_id:** `UOX38086.1`
- **Seed length:** 806 AA
- **Hard editable positions:** [300, 328, 354, 422, 561, 615, 656, 768]
- **Soft editable positions:** [360, 574, 642]
- **Mutation budget:** 4 … 10

## Search statistics

- **Unique candidates evaluated:** 662
- **Rounds completed:** 6
- **Best composite score:** 0.932
- **Best target probability:** 0.230
- **Best seed identity:** 0.999

## Stage 11 structural validation (top-3, pLDDT corrected to 0–100 Cα convention)

- **Validated rows:** 3
- **Pass count:** 3 (100% pass rate)
- **Best mean pLDDT:** 79.04
- **Mean mean pLDDT:** 78.90
- **Best mutation-site mean pLDDT:** 71.75
- **Best RMSD to seed:** 1.079 Å
- **Mean RMSD to seed:** 1.183 Å

## Passing candidates

| sample_id | mutation_count | mutation_text | mean_pLDDT | mut_site_pLDDT | RMSD (Å) | target_prob | pass |
|---:|---:|---|---:|---:|---:|---:|:---:|
| 1 | 4 | 354:D→T;561:A→V;615:K→Q;656:Y→T | 79.03 | 71.75 | 1.257 | 0.225 | True |
| 2 | 4 | 354:D→N;561:A→V;615:K→Q;656:Y→T | 79.04 | 71.25 | 1.212 | 0.225 | True |
| 3 | 5 | 300:N→R;354:D→T;561:A→V;615:K→Q;656:Y→T | 78.62 | 68.60 | 1.079 | 0.227 | True |

## Conclusion

3/3 Stage 11 candidates passed all six 08a hard thresholds. The Baseline Qualification hypothesis is supported: anchoring inverse-folding on a structurally healthy wild-type seed yields candidates that retain the fold (high pLDDT, low RMSD) while carrying the targeted edits.

_pLDDT values corrected to the canonical 0–100 Cα convention used by the Stage 11a baseline gate. The original report applied a 0–1 global scale (and averaged over non-existent atom slots), which spuriously failed every candidate despite folds with mean pLDDT ≈ 79 and RMSD < 1.6 Å._