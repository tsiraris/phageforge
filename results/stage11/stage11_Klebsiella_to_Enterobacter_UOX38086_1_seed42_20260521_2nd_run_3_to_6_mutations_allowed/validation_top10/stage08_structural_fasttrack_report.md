# Stage 08 structural fast-track report

- candidates analyzed: **10**
- seed ESMFold mean pLDDT: **79.39**
- seed ESMFold pTM: **0.7257**
- candidates passing all six gates: **9/10**
- resume mode: **False**

## Candidate summary

| rank | sample_id | regime | final_score | mean_pLDDT | mut_mean_pLDDT | RMSD_to_seed | mut_conf>=70 | pass | reason |
|---:|---:|---|---:|---:|---:|---:|---:|:---:|---|
| 1 | 1 | stage11_inverse_folding | 0.910658 | 78.73 | 65.33 | 0.934 | 0.667 | True | pass |
| 2 | 2 | stage11_inverse_folding | 0.906756 | 78.75 | 69.25 | 1.348 | 0.750 | True | pass |
| 3 | 3 | stage11_inverse_folding | 0.903672 | 78.60 | 63.67 | 0.797 | 0.333 | False | low_mutation_site_confidence |
| 4 | 4 | stage11_inverse_folding | 0.890361 | 79.00 | 76.00 | 1.123 | 1.000 | True | pass |
| 5 | 5 | stage11_inverse_folding | 0.887567 | 78.67 | 68.50 | 1.322 | 0.500 | True | pass |
| 6 | 6 | stage11_inverse_folding | 0.884209 | 79.00 | 75.33 | 1.098 | 0.667 | True | pass |
| 7 | 7 | stage11_inverse_folding | 0.878848 | 79.03 | 71.75 | 1.257 | 0.750 | True | pass |
| 8 | 8 | stage11_inverse_folding | 0.865334 | 78.68 | 69.33 | 1.271 | 0.667 | True | pass |
| 9 | 9 | stage11_inverse_folding | 0.861079 | 78.90 | 70.50 | 1.183 | 0.500 | True | pass |
| 10 | 10 | stage11_inverse_folding | 0.848521 | 78.46 | 67.80 | 1.280 | 0.600 | True | pass |

## Interpretation
- Prefer candidates with high multimodal score, high global pLDDT, high mutation-site confidence, and low RMSD to the selected seed.
- `stage08_pass` is intentionally conservative; near-pass candidates can still be useful backups.
- If at least two candidates pass, those should become the immediate closeout handoff set.

_Note: pLDDT corrected to the canonical 0–100 Cα convention (matching the Stage 11a baseline gate). Prior values were on a 0–1 scale (global) or averaged over non-existent atom slots, which spuriously failed every candidate._