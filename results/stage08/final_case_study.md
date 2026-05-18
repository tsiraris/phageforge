# Final closeout case study

## Project framing
A validity-aware, scaffold-constrained phage RBP retargeting workflow was closed out with local ESM3 generation, Stage 07 multimodal reranking, and Stage 08 structural fast-track validation.

## Selected seed
- protein_id: **round8_cand473**
- source host: **Klebsiella**
- virus accession: **round8_cand473**
- sequence length: **658**

## Structural closeout summary
- final candidate count: **3**
- Stage 08 passes: **0**
- primary candidates: **none**

## Final ranking

| final_rank | sample_id | decision | pass | final_score | mean_pLDDT | mut_mean_pLDDT | RMSD_to_seed | reason |
|---:|---:|---|:---:|---:|---:|---:|---:|---|
| 1 | 11 | near_pass_or_fail | False | 0.784310 | 0.24 | 0.24 | 17.628 | low_global_confidence;high_seed_drift;low_mutation_site_confidence |
| 2 | 21 | near_pass_or_fail | False | 0.523553 | 0.24 | 0.26 | 17.059 | low_global_confidence;high_seed_drift;low_mutation_site_confidence |
| 3 | 28 | near_pass_or_fail | False | 0.519171 | 0.26 | 0.49 | 14.749 | low_global_confidence;high_seed_drift;low_mutation_site_confidence |

## Deliverables
- final table: `/home/sagemaker-user/results/final_closeout/final_candidate_table.csv`
- final FASTA: `/home/sagemaker-user/results/final_closeout/final_top_candidates.fasta`
- summary JSON: `/home/sagemaker-user/results/final_closeout/final_closeout_summary.json`