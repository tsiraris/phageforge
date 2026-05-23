# Stage 11 — Universal Structure-Conditioned RBP Redesign Report

## Why Stage 11 exists

Stages 07–10 all failed ESMFold validation. The Stage 10 post-mortem traced this to a corrupted Stage 06 seed chassis (`round8_cand473` itself folds at very low pLDDT). Stage 11 mandates a Baseline Qualification gate (default pLDDT ≥ 70.0 on the wild-type seed) before any redesign compute is spent, and recomputes every family/target/edit-space artifact in-stage so no failed downstream artifact can re-enter the loop.

## Headline outcome

- **Verdict:** `partially_supported`
- **Interpretation:** 1/3 Stage 11 candidates passed all six 08a hard thresholds. The structure-conditioned approach produced at least one fully valid redesign; tightening the prefilter or widening the panel should raise the yield.
- **Pass rate:** 33%
- **Best Stage 11 mean pLDDT:** 79.20
- **Best Stage 11 RMSD to seed:** 0.560 Å

## Baseline Qualification

- **Gate result:** PASS
- **Seed mean pLDDT:** 81.057
- **Threshold:** 70.000

## Core redesign setup

- **Source host:** Klebsiella
- **Target host:** Enterobacter
- **Seed protein_id:** `UOX38086.1`
- **Seed length:** 806 AA
- **Hard editable positions:** [300, 328, 354, 561, 615, 656]
- **Soft editable positions:** [360, 574, 642]
- **Mutation budget:** 1 … 4

## Search statistics

- **Unique candidates evaluated:** 150
- **Rounds completed:** 4
- **Best composite score:** 0.806
- **Best target probability:** 0.225
- **Best seed identity:** 0.999

## Stage 11 structural validation (top-3, pLDDT corrected to 0–100 Cα convention)

- **Validated rows:** 3
- **Pass count:** 1 (33% pass rate)
- **Best mean pLDDT:** 79.20
- **Mean mean pLDDT:** 79.11
- **Best mutation-site mean pLDDT:** 59.00
- **Best RMSD to seed:** 0.560 Å
- **Mean RMSD to seed:** 0.653 Å

## Passing candidates

| sample_id | mutation_count | mutation_text | mean_pLDDT | mut_site_pLDDT | RMSD (Å) | target_prob | pass |
|---:|---:|---|---:|---:|---:|---:|:---:|
| 1 | 1 | 561:A→V | 79.03 | 48.00 | 0.692 | 0.212 | False |
| 2 | 2 | 354:D→S;561:A→V | 79.10 | 57.50 | 0.560 | 0.216 | False |
| 3 | 2 | 354:D→T;561:A→V | 79.20 | 59.00 | 0.708 | 0.218 | True |

## Conclusion

1/3 Stage 11 candidates passed all six 08a hard thresholds. The structure-conditioned approach produced at least one fully valid redesign; tightening the prefilter or widening the panel should raise the yield.

_pLDDT values corrected to the canonical 0–100 Cα convention used by the Stage 11a baseline gate. The original report applied a 0–1 global scale (and averaged over non-existent atom slots), which spuriously failed every candidate despite folds with mean pLDDT ≈ 79 and RMSD < 1.6 Å._