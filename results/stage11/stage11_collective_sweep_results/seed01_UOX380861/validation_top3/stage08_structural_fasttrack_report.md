# Stage 08 structural fast-track report

- candidates analyzed: **3**
- seed ESMFold mean pLDDT: **79.39**
- seed ESMFold pTM: **0.7257**
- resume mode: **False**

## Candidate summary

| rank | sample_id | regime | final_score | mean_pLDDT | mut_mean_pLDDT | RMSD_to_seed | mut_conf>=70 | pass | reason |
|---:|---:|---|---:|---:|---:|---:|---:|:---:|---|
| 1 | 1 | stage11_inverse_folding | 0.760494 | 78.67 | 68.50 | 1.322 | 0.500 | True | pass |
| 2 | 2 | stage11_inverse_folding | 0.748007 | 78.75 | 69.25 | 1.348 | 0.750 | True | pass |
| 3 | 3 | stage11_inverse_folding | 0.722062 | 78.68 | 69.33 | 1.271 | 0.667 | True | pass |

## Interpretation
- Prefer candidates with high multimodal score, high global pLDDT, high mutation-site confidence, and low RMSD to the selected seed.
- `stage08_pass` is intentionally conservative; near-pass candidates can still be useful backups.
- If at least two candidates pass, those should become the immediate closeout handoff set.