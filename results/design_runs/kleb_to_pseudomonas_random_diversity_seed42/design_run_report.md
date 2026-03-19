# Design Run Summary

## Run metadata

- Seed protein ID: `QFR57578.1`
- Virus accession: `MN434096`
- Source host: `Klebsiella`
- Target host: `Pseudomonas`
- ESM model: `facebook/esm2_t33_650M_UR50D`
- Rounds: `6`
- Candidates per round: `64`
- Min mutations: `1`
- Max mutations: `4`
- Top-K kept per round: `10`
- Proposal top-K: `8`
- Seed: `42`
- Total candidates evaluated: `3263`
- Top candidates saved: `50`
- Diversity enabled: `True`
- Diversity min distance: `8`
- Position selection strategy: `random`
- Position pool size: `N/A`

## Target score summary

| metric   |   target_score |
|:---------|---------------:|
| min      |       0.343598 |
| mean     |       0.382832 |
| median   |       0.372433 |
| max      |       0.464707 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Pseudomonas  |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       1 |
|             2 |       1 |
|             3 |       4 |
|             4 |      44 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  26.0327 |
| unique_mutation_site_count_top_candidates | 172      |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000797749 |
| mean_novelty_distance                     | 0.00156089  |
| median_novelty_distance                   | 0.00157472  |
| max_novelty_distance                      | 0.00281113  |
| min_nearest_neighbor_cosine_similarity    | 0.997189    |
| mean_nearest_neighbor_cosine_similarity   | 0.998439    |
| median_nearest_neighbor_cosine_similarity | 0.998425    |
| max_nearest_neighbor_cosine_similarity    | 0.999202    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round6_cand150 | QFR57578.1                    | Klebsiella                    |                             0.997189 |         0.00281113 |
| round6_cand528 | QFR57578.1                    | Klebsiella                    |                             0.997824 |         0.00217593 |
| round6_cand199 | QFR57578.1                    | Klebsiella                    |                             0.997745 |         0.00225496 |
| round6_cand117 | QFR57578.1                    | Klebsiella                    |                             0.998351 |         0.0016489  |
| round6_cand564 | QFR57578.1                    | Klebsiella                    |                             0.99811  |         0.00188965 |
| round6_cand100 | QFR57578.1                    | Klebsiella                    |                             0.9985   |         0.00150037 |
| round6_cand78  | QFR57578.1                    | Klebsiella                    |                             0.998434 |         0.00156647 |
| round6_cand585 | QFR57578.1                    | Klebsiella                    |                             0.99866  |         0.00133955 |
| round6_cand331 | QFR57578.1                    | Klebsiella                    |                             0.998733 |         0.00126719 |
| round6_cand392 | QFR57578.1                    | Klebsiella                    |                             0.998222 |         0.00177824 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score | mutations               |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|:------------------------|--------------:|
| round6_cand150 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.464707 | G150Y                   |             1 |
| round6_cand528 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.445993 | Y77G;R175A;I192G;E526V  |             4 |
| round6_cand199 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.437639 | N25V;D67A;I250F;F577Q   |             4 |
| round6_cand117 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.434377 | G122A;N367D;T413G;V433A |             4 |
| round6_cand564 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.433784 | Y200A;H228N;Q560A;C587V |             4 |
| round6_cand100 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.429396 | V168A;I381A;Y383L;S558V |             4 |
| round6_cand78  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.41917  | Y336A;V398R;F401L;W469I |             4 |
| round6_cand585 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.415962 | V141I;G167R;Q339V;V571N |             4 |
| round6_cand331 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.415439 | Q71R;G97D               |             2 |
| round6_cand392 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.412138 | F94R;R477Y;D549A;M596G  |             4 |
