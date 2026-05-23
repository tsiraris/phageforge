# Stage 08 structural fast-track report

- candidates analyzed: **3**
- seed ESMFold mean pLDDT: **0.22**
- seed ESMFold pTM: **0.1475**
- resume mode: **False**

## Candidate summary

| rank | sample_id | regime | final_score | mean_pLDDT | mut_mean_pLDDT | RMSD_to_seed | mut_conf>=70 | pass | reason |
|---:|---:|---|---:|---:|---:|---:|---:|:---:|---|
| 1 | 3 | InverseFolding | 0.000000 | 0.22 | nan | 11.792 | nan | False | low_global_confidence;high_seed_drift |
| 2 | 2 | InverseFolding | 0.000000 | 0.22 | nan | 15.812 | nan | False | low_global_confidence;high_seed_drift |
| 3 | 1 | InverseFolding | 0.000000 | 0.22 | nan | 15.902 | nan | False | low_global_confidence;high_seed_drift |

## Interpretation
- Prefer candidates with high multimodal score, high global pLDDT, high mutation-site confidence, and low RMSD to the selected seed.
- `stage08_pass` is intentionally conservative; near-pass candidates can still be useful backups.
- If at least two candidates pass, those should become the immediate closeout handoff set.