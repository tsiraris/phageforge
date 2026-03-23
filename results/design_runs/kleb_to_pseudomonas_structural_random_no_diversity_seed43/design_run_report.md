# Design Run Summary

## Run metadata

- Seed protein ID: `QFR57578.1`
- Virus accession: `MN434096`
- Source host: `Klebsiella`
- Target host: `Pseudomonas`
- ESM model: `facebook/esm2_t33_650M_UR50D`
- Rounds: `4`
- Candidates per round: `64`
- Min mutations: `1`
- Max mutations: `3`
- Top-K kept per round: `10`
- Proposal top-K: `8`
- Seed: `43`
- Total candidates evaluated: `1983`
- Top candidates saved: `50`
- Diversity enabled: `False`
- Diversity min distance: `20`
- Position selection strategy: `random`
- Position pool size: `32`
- Proposal batch size: `16`
- Novelty penalty lambda: `0.0`
- Novelty penalty schedule: `constant`

## Target score summary

| metric   |   target_score |   selection_score |
|:---------|---------------:|------------------:|
| min      |       0.2432   |          0.2432   |
| mean     |       0.251892 |          0.251892 |
| median   |       0.250289 |          0.250289 |
| max      |       0.274391 |          0.274391 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       4 |
|             2 |      17 |
|             3 |      29 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  16.1731 |
| unique_mutation_site_count_top_candidates | 104      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998914    |
| mean_seed_cosine_similarity   | 0.999367    |
| median_seed_cosine_similarity | 0.999379    |
| max_seed_cosine_similarity    | 0.999757    |
| min_seed_novelty_distance     | 0.000242591 |
| mean_seed_novelty_distance    | 0.000633045 |
| median_seed_novelty_distance  | 0.000621408 |
| max_seed_novelty_distance     | 0.00108635  |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000242412 |
| mean_novelty_distance                     | 0.000633138 |
| median_novelty_distance                   | 0.000621676 |
| max_novelty_distance                      | 0.00108618  |
| min_nearest_neighbor_cosine_similarity    | 0.998914    |
| mean_nearest_neighbor_cosine_similarity   | 0.999367    |
| median_nearest_neighbor_cosine_similarity | 0.999378    |
| max_nearest_neighbor_cosine_similarity    | 0.999758    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand572 | QFR57578.1                    | Klebsiella                    |                             0.999698 |        0.000301957 |
| round4_cand115 | QFR57578.1                    | Klebsiella                    |                             0.999531 |        0.000468612 |
| round4_cand544 | QFR57578.1                    | Klebsiella                    |                             0.999659 |        0.000340879 |
| round4_cand15  | QFR57578.1                    | Klebsiella                    |                             0.999335 |        0.000665188 |
| round4_cand336 | QFR57578.1                    | Klebsiella                    |                             0.999436 |        0.000563741 |
| round4_cand122 | QFR57578.1                    | Klebsiella                    |                             0.999532 |        0.000468433 |
| round4_cand446 | QFR57578.1                    | Klebsiella                    |                             0.999698 |        0.000302196 |
| round4_cand267 | QFR57578.1                    | Klebsiella                    |                             0.999589 |        0.000410676 |
| round4_cand284 | QFR57578.1                    | Klebsiella                    |                             0.999733 |        0.00026679  |
| round4_cand448 | QFR57578.1                    | Klebsiella                    |                             0.999214 |        0.000785589 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand572 | Klebsiella    | Pseudomonas   | Enterobacter |       0.274391 |          0.274391 |                 0.999698 |             0.0003016   | S546A;D626P       |             2 |
| round4_cand115 | Klebsiella    | Pseudomonas   | Enterobacter |       0.271929 |          0.271929 |                 0.999531 |             0.000468791 | S38P;Y144A        |             2 |
| round4_cand544 | Klebsiella    | Pseudomonas   | Enterobacter |       0.268598 |          0.268598 |                 0.999659 |             0.000340939 | S213P;R289A;N424L |             3 |
| round4_cand15  | Klebsiella    | Pseudomonas   | Enterobacter |       0.266001 |          0.266001 |                 0.999335 |             0.000664949 | E86R;G176P;M598A  |             3 |
| round4_cand336 | Klebsiella    | Pseudomonas   | Enterobacter |       0.265745 |          0.265745 |                 0.999437 |             0.000563323 | H49T;F321A;E353P  |             3 |
| round4_cand122 | Klebsiella    | Pseudomonas   | Enterobacter |       0.262979 |          0.262979 |                 0.999532 |             0.000468075 | Y251L;E636A       |             2 |
| round4_cand446 | Klebsiella    | Pseudomonas   | Enterobacter |       0.262494 |          0.262494 |                 0.999698 |             0.000302255 | Y77P;Y200T;S262G  |             3 |
| round4_cand267 | Klebsiella    | Pseudomonas   | Enterobacter |       0.257839 |          0.257839 |                 0.99959  |             0.000410438 | D67G              |             1 |
| round4_cand284 | Klebsiella    | Pseudomonas   | Enterobacter |       0.257118 |          0.257118 |                 0.999733 |             0.00026679  | N367R;T473A       |             2 |
| round4_cand448 | Klebsiella    | Pseudomonas   | Enterobacter |       0.256833 |          0.256833 |                 0.999214 |             0.000785649 | K260P;Y336I;T583S |             3 |
