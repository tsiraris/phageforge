# Stage 08 structural fast-track report

- candidates analyzed: **3**
- seed ESMFold mean pLDDT: **88.57**
- seed ESMFold pTM: **0.9135**
- resume mode: **False**

## Candidate summary

| rank | sample_id | regime | final_score | mean_pLDDT | mut_mean_pLDDT | RMSD_to_seed | mut_conf>=70 | pass | reason |
|---:|---:|---|---:|---:|---:|---:|---:|:---:|---|
| 1 | 1 | stage11_inverse_folding | 0.679925 | 88.59 | 96.00 | 0.214 | 1.000 | True | pass |
| 2 | 2 | stage11_inverse_folding | 0.671870 | 88.71 | 96.00 | 0.188 | 1.000 | True | pass |
| 3 | 3 | stage11_inverse_folding | 0.660072 | 88.49 | 96.00 | 0.243 | 1.000 | True | pass |

## Interpretation
- Prefer candidates with high multimodal score, high global pLDDT, high mutation-site confidence, and low RMSD to the selected seed.
- `stage08_pass` is intentionally conservative; near-pass candidates can still be useful backups.
- If at least two candidates pass, those should become the immediate closeout handoff set.