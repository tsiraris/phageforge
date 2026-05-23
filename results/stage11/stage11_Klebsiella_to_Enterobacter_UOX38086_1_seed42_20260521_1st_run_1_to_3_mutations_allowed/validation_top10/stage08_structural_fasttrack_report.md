# Stage 08 structural fast-track report

- candidates analyzed: **10**
- seed ESMFold mean pLDDT: **79.39**
- seed ESMFold pTM: **0.7257**
- candidates passing all six gates: **7/10**
- resume mode: **False**

## Candidate summary

| rank | sample_id | regime | final_score | mean_pLDDT | mut_mean_pLDDT | RMSD_to_seed | mut_conf>=70 | pass | reason |
|---:|---:|---|---:|---:|---:|---:|---:|:---:|---|
| 1 | 1 | stage11_inverse_folding | 0.805972 | 79.03 | 48.00 | 0.692 | 0.000 | False | low_mutation_site_confidence |
| 2 | 2 | stage11_inverse_folding | 0.797846 | 79.10 | 57.50 | 0.560 | 0.000 | False | low_mutation_site_confidence |
| 3 | 3 | stage11_inverse_folding | 0.785109 | 79.20 | 59.00 | 0.708 | 0.500 | True | pass |
| 4 | 4 | stage11_inverse_folding | 0.778649 | 78.67 | 68.50 | 1.322 | 0.500 | True | pass |
| 5 | 5 | stage11_inverse_folding | 0.760309 | 78.75 | 69.25 | 1.348 | 0.750 | True | pass |
| 6 | 6 | stage11_inverse_folding | 0.705663 | 79.32 | 68.67 | 1.432 | 0.333 | False | low_mutation_site_confidence |
| 7 | 7 | stage11_inverse_folding | 0.699097 | 79.60 | 71.67 | 1.369 | 0.667 | True | pass |
| 8 | 8 | stage11_inverse_folding | 0.694714 | 79.25 | 77.50 | 0.462 | 0.500 | True | pass |
| 9 | 9 | stage11_inverse_folding | 0.682938 | 79.25 | 78.50 | 0.453 | 1.000 | True | pass |
| 10 | 10 | stage11_inverse_folding | 0.672115 | 78.90 | 68.75 | 1.559 | 0.500 | True | pass |

## Interpretation
- Prefer candidates with high multimodal score, high global pLDDT, high mutation-site confidence, and low RMSD to the selected seed.
- `stage08_pass` is intentionally conservative; near-pass candidates can still be useful backups.
- If at least two candidates pass, those should become the immediate closeout handoff set.

_Note: pLDDT corrected to the canonical 0–100 Cα convention (matching the Stage 11a baseline gate). Prior values were on a 0–1 scale (global) or averaged over non-existent atom slots, which spuriously failed every candidate._