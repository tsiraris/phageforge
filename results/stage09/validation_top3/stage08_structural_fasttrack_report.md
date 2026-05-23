# Stage 08 structural fast-track report

- candidates analyzed: **3**
- seed ESMFold mean pLDDT: **0.22**
- seed ESMFold pTM: **0.1475**
- resume mode: **True**

## Candidate summary

| rank | sample_id | regime | final_score | mean_pLDDT | mut_mean_pLDDT | RMSD_to_seed | mut_conf>=70 | pass | reason |
|---:|---:|---|---:|---:|---:|---:|---:|:---:|---|
| 1 | 1 | stage09_localized_search | 0.335713 | 0.22 | 0.22 | 13.341 | 0.000 | False | low_global_confidence;high_seed_drift;low_mutation_site_confidence |
| 2 | 2 | stage09_localized_search | 0.335432 | 0.23 | 0.23 | 14.073 | 0.000 | False | low_global_confidence;high_seed_drift;low_mutation_site_confidence |
| 3 | 3 | stage09_localized_search | 0.335353 | 0.23 | 0.22 | 13.129 | 0.000 | False | low_global_confidence;high_seed_drift;low_mutation_site_confidence |

## Interpretation
- Prefer candidates with high multimodal score, high global pLDDT, high mutation-site confidence, and low RMSD to the selected seed.
- `stage08_pass` is intentionally conservative; near-pass candidates can still be useful backups.
- If at least two candidates pass, those should become the immediate closeout handoff set.