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
- Seed: `44`
- Total candidates evaluated: `1982`
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
| min      |       0.231999 |          0.231999 |
| mean     |       0.242567 |          0.242567 |
| median   |       0.239781 |          0.239781 |
| max      |       0.277427 |          0.277427 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      49 |
| Klebsiella   |       1 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       6 |
|             2 |      15 |
|             3 |      29 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  16.3584 |
| unique_mutation_site_count_top_candidates | 100      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998679    |
| mean_seed_cosine_similarity   | 0.999559    |
| median_seed_cosine_similarity | 0.999617    |
| max_seed_cosine_similarity    | 0.999766    |
| min_seed_novelty_distance     | 0.000234306 |
| mean_seed_novelty_distance    | 0.000441414 |
| median_seed_novelty_distance  | 0.000382751 |
| max_seed_novelty_distance     | 0.0013209   |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000234723 |
| mean_novelty_distance                     | 0.000441496 |
| median_novelty_distance                   | 0.000382751 |
| max_novelty_distance                      | 0.00132096  |
| min_nearest_neighbor_cosine_similarity    | 0.998679    |
| mean_nearest_neighbor_cosine_similarity   | 0.999559    |
| median_nearest_neighbor_cosine_similarity | 0.999617    |
| max_nearest_neighbor_cosine_similarity    | 0.999765    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand254 | QFR57578.1                    | Klebsiella                    |                             0.999285 |        0.00071466  |
| round4_cand321 | QFR57578.1                    | Klebsiella                    |                             0.999515 |        0.000485182 |
| round4_cand230 | QFR57578.1                    | Klebsiella                    |                             0.999512 |        0.000487924 |
| round4_cand288 | QFR57578.1                    | Klebsiella                    |                             0.999567 |        0.00043273  |
| round4_cand243 | QFR57578.1                    | Klebsiella                    |                             0.999422 |        0.000577927 |
| round4_cand350 | QFR57578.1                    | Klebsiella                    |                             0.9995   |        0.000499845 |
| round4_cand263 | QFR57578.1                    | Klebsiella                    |                             0.999441 |        0.000559151 |
| round4_cand134 | QFR57578.1                    | Klebsiella                    |                             0.999529 |        0.000471354 |
| round4_cand126 | QFR57578.1                    | Klebsiella                    |                             0.999621 |        0.000378609 |
| round4_cand51  | QFR57578.1                    | Klebsiella                    |                             0.999701 |        0.000299096 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand254 | Klebsiella    | Pseudomonas   | Enterobacter |       0.277427 |          0.277427 |                 0.999285 |             0.0007146   | Y200A;R232E;E526P |             3 |
| round4_cand321 | Klebsiella    | Pseudomonas   | Enterobacter |       0.27413  |          0.27413  |                 0.999515 |             0.000485003 | Y77P;E526N;G582A  |             3 |
| round4_cand230 | Klebsiella    | Pseudomonas   | Enterobacter |       0.263513 |          0.263513 |                 0.999512 |             0.000487804 | G101P;D566E;W602A |             3 |
| round4_cand288 | Klebsiella    | Pseudomonas   | Enterobacter |       0.25737  |          0.25737  |                 0.999567 |             0.00043273  | S384P;G588I;M596S |             3 |
| round4_cand243 | Klebsiella    | Pseudomonas   | Enterobacter |       0.256491 |          0.256491 |                 0.999422 |             0.000577569 | N313A;N402S;T413L |             3 |
| round4_cand350 | Klebsiella    | Pseudomonas   | Enterobacter |       0.253227 |          0.253227 |                 0.9995   |             0.000499785 | K374L;G411R;F559L |             3 |
| round4_cand263 | Klebsiella    | Pseudomonas   | Enterobacter |       0.252042 |          0.252042 |                 0.999441 |             0.000559092 | E50R;K374L;H564G  |             3 |
| round4_cand134 | Klebsiella    | Pseudomonas   | Enterobacter |       0.251067 |          0.251067 |                 0.999528 |             0.000471711 | S357G;T473R;C574N |             3 |
| round4_cand126 | Klebsiella    | Pseudomonas   | Enterobacter |       0.24979  |          0.24979  |                 0.999622 |             0.000378489 | D423G;A529P;Q586V |             3 |
| round4_cand51  | Klebsiella    | Pseudomonas   | Enterobacter |       0.249175 |          0.249175 |                 0.999701 |             0.000299335 | G6K;Y144A         |             2 |
