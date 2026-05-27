# Stage 08 structural fast-track report

- candidates analyzed: **3**
- seed ESMFold mean pLDDT: **73.85**
- seed ESMFold pTM: **0.6279**
- resume mode: **False**

## Candidate summary

| rank | sample_id | regime | final_score | mean_pLDDT | mut_mean_pLDDT | RMSD_to_seed | mut_conf>=70 | pass | reason |
|---:|---:|---|---:|---:|---:|---:|---:|:---:|---|
| 1 | 1 | stage11_inverse_folding | 0.779082 | 72.72 | 75.83 | 1.814 | 0.833 | True | pass |
| 2 | 2 | stage11_inverse_folding | 0.772378 | 72.94 | 76.50 | 1.809 | 0.833 | True | pass |
| 3 | 3 | stage11_inverse_folding | 0.768013 | 72.81 | 75.33 | 1.278 | 0.833 | True | pass |

## Interpretation
- Prefer candidates with high multimodal score, high global pLDDT, high mutation-site confidence, and low RMSD to the selected seed.
- `stage08_pass` is intentionally conservative; near-pass candidates can still be useful backups.
- If at least two candidates pass, those should become the immediate closeout handoff set.