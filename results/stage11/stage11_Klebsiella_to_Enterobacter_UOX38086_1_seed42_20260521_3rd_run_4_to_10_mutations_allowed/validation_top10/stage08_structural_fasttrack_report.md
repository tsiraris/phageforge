# Stage 08 structural fast-track report

- candidates analyzed: **10**
- seed ESMFold mean pLDDT: **79.39**
- seed ESMFold pTM: **0.7257**
- candidates passing all six gates: **7/10**
- resume mode: **False**

## Candidate summary

| rank | sample_id | regime | final_score | mean_pLDDT | mut_mean_pLDDT | RMSD_to_seed | mut_conf>=70 | pass | reason |
|---:|---:|---|---:|---:|---:|---:|---:|:---:|---|
| 1 | 1 | stage11_inverse_folding | 0.854635 | 79.03 | 71.75 | 1.257 | 0.750 | True | pass |
| 2 | 2 | stage11_inverse_folding | 0.834216 | 79.04 | 71.25 | 1.212 | 0.500 | True | pass |
| 3 | 3 | stage11_inverse_folding | 0.822953 | 78.62 | 68.60 | 1.079 | 0.600 | True | pass |
| 4 | 4 | stage11_inverse_folding | 0.820713 | 78.63 | 68.33 | 1.509 | 0.500 | True | pass |
| 5 | 5 | stage11_inverse_folding | 0.816063 | 78.67 | 68.60 | 1.051 | 0.400 | False | low_mutation_site_confidence |
| 6 | 6 | stage11_inverse_folding | 0.812728 | 78.99 | 70.25 | 1.308 | 0.500 | True | pass |
| 7 | 7 | stage11_inverse_folding | 0.811726 | 78.90 | 70.50 | 1.183 | 0.500 | True | pass |
| 8 | 8 | stage11_inverse_folding | 0.807555 | 78.70 | 68.50 | 1.469 | 0.333 | False | low_mutation_site_confidence |
| 9 | 9 | stage11_inverse_folding | 0.788686 | 78.79 | 68.60 | 1.226 | 0.400 | False | low_mutation_site_confidence |
| 10 | 10 | stage11_inverse_folding | 0.786424 | 78.85 | 69.60 | 1.583 | 0.600 | True | pass |

## Interpretation
- Prefer candidates with high multimodal score, high global pLDDT, high mutation-site confidence, and low RMSD to the selected seed.
- `stage08_pass` is intentionally conservative; near-pass candidates can still be useful backups.
- If at least two candidates pass, those should become the immediate closeout handoff set.

_Note: pLDDT corrected to the canonical 0–100 Cα convention (matching the Stage 11a baseline gate). Prior values were on a 0–1 scale (global) or averaged over non-existent atom slots, which spuriously failed every candidate._