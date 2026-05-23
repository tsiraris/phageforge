# Stage 08 structural fast-track report

- candidates analyzed: **3**
- seed ESMFold mean pLDDT: **79.39**
- seed ESMFold pTM: **0.7257**
- candidates passing all six gates: **1/3**
- resume mode: **False**

## Candidate summary

| rank | sample_id | regime | final_score | mean_pLDDT | mut_mean_pLDDT | RMSD_to_seed | mut_conf>=70 | pass | reason |
|---:|---:|---|---:|---:|---:|---:|---:|:---:|---|
| 1 | 1 | stage11_inverse_folding | 0.805972 | 79.03 | 48.00 | 0.692 | 0.000 | False | low_mutation_site_confidence |
| 2 | 2 | stage11_inverse_folding | 0.797846 | 79.10 | 57.50 | 0.560 | 0.000 | False | low_mutation_site_confidence |
| 3 | 3 | stage11_inverse_folding | 0.785109 | 79.20 | 59.00 | 0.708 | 0.500 | True | pass |

## Interpretation
- Prefer candidates with high multimodal score, high global pLDDT, high mutation-site confidence, and low RMSD to the selected seed.
- `stage08_pass` is intentionally conservative; near-pass candidates can still be useful backups.
- If at least two candidates pass, those should become the immediate closeout handoff set.

_Note: pLDDT corrected to the canonical 0–100 Cα convention (matching the Stage 11a baseline gate). Prior values were on a 0–1 scale (global) or averaged over non-existent atom slots, which spuriously failed every candidate._