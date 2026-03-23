# Design Run Summary

## Run metadata

- Seed protein ID: `QFR57578.1`
- Virus accession: `MN434096`
- Source host: `Klebsiella`
- Target host: `Enterobacter`
- ESM model: `facebook/esm2_t33_650M_UR50D`
- Rounds: `4`
- Candidates per round: `64`
- Min mutations: `1`
- Max mutations: `3`
- Top-K kept per round: `10`
- Proposal top-K: `8`
- Seed: `43`
- Total candidates evaluated: `1964`
- Top candidates saved: `50`
- Diversity enabled: `True`
- Diversity min distance: `20`
- Position selection strategy: `entropy`
- Position pool size: `32`
- Proposal batch size: `16`
- Novelty penalty lambda: `5.0`
- Novelty penalty schedule: `constant`

## Target score summary

| metric   |   target_score |   selection_score |
|:---------|---------------:|------------------:|
| min      |       0.604255 |          -4.39048 |
| mean     |       0.614356 |          -4.38088 |
| median   |       0.612152 |          -4.38325 |
| max      |       0.645079 |          -4.34995 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       7 |
|             2 |      16 |
|             3 |      27 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  7.63184 |
| unique_mutation_site_count_top_candidates | 31       |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998566    |
| mean_seed_cosine_similarity   | 0.999047    |
| median_seed_cosine_similarity | 0.999065    |
| max_seed_cosine_similarity    | 0.999308    |
| min_seed_novelty_distance     | 0.000692129 |
| mean_seed_novelty_distance    | 0.000953016 |
| median_seed_novelty_distance  | 0.000934869 |
| max_seed_novelty_distance     | 0.00143355  |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.00069195  |
| mean_novelty_distance                     | 0.000949882 |
| median_novelty_distance                   | 0.000934839 |
| max_novelty_distance                      | 0.00143349  |
| min_nearest_neighbor_cosine_similarity    | 0.998567    |
| mean_nearest_neighbor_cosine_similarity   | 0.99905     |
| median_nearest_neighbor_cosine_similarity | 0.999065    |
| max_nearest_neighbor_cosine_similarity    | 0.999308    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand26  | QFR57578.1                    | Klebsiella                    |                             0.999007 |        0.000993013 |
| round4_cand92  | WWD14686.1                    | Klebsiella                    |                             0.999058 |        0.000941634 |
| round4_cand3   | WWD14686.1                    | Klebsiella                    |                             0.998976 |        0.00102425  |
| round4_cand34  | QFR57578.1                    | Klebsiella                    |                             0.998567 |        0.00143349  |
| round4_cand74  | QFR57578.1                    | Klebsiella                    |                             0.999031 |        0.00096941  |
| round4_cand54  | QFR57578.1                    | Klebsiella                    |                             0.998924 |        0.00107551  |
| round4_cand88  | QFR57578.1                    | Klebsiella                    |                             0.99911  |        0.000890017 |
| round4_cand94  | QFR57578.1                    | Klebsiella                    |                             0.999308 |        0.00069195  |
| round4_cand309 | QFR57578.1                    | Klebsiella                    |                             0.998809 |        0.00119066  |
| round4_cand328 | QFR57578.1                    | Klebsiella                    |                             0.999066 |        0.000934482 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand26  | Klebsiella    | Enterobacter  | Enterobacter |       0.645079 |          -4.34995 |                 0.999007 |             0.000993133 | C132S;K523G;W527I |             3 |
| round4_cand92  | Klebsiella    | Enterobacter  | Enterobacter |       0.632285 |          -4.363   |                 0.999057 |             0.000942707 | K22P;C516T;K523T  |             3 |
| round4_cand3   | Klebsiella    | Enterobacter  | Enterobacter |       0.630249 |          -4.36459 |                 0.998967 |             0.00103271  | H27A;K523N;C568S  |             3 |
| round4_cand34  | Klebsiella    | Enterobacter  | Enterobacter |       0.626794 |          -4.36604 |                 0.998566 |             0.00143355  | C162R;W527T;C587T |             3 |
| round4_cand74  | Klebsiella    | Enterobacter  | Enterobacter |       0.627344 |          -4.36781 |                 0.99903  |             0.000969708 | C132V;K335A;C516D |             3 |
| round4_cand54  | Klebsiella    | Enterobacter  | Enterobacter |       0.626678 |          -4.36794 |                 0.998925 |             0.00107545  | E427V;K523T;K599R |             3 |
| round4_cand88  | Klebsiella    | Enterobacter  | Enterobacter |       0.625687 |          -4.36986 |                 0.99911  |             0.000890076 | C132S;C516G;K599Q |             3 |
| round4_cand94  | Klebsiella    | Enterobacter  | Enterobacter |       0.62352  |          -4.37302 |                 0.999308 |             0.000692129 | C516D             |             1 |
| round4_cand309 | Klebsiella    | Enterobacter  | Enterobacter |       0.620361 |          -4.37369 |                 0.998809 |             0.00119054  | K22T;C162R;W505T  |             3 |
| round4_cand328 | Klebsiella    | Enterobacter  | Enterobacter |       0.620236 |          -4.37509 |                 0.999066 |             0.000934482 | M226Y;W505D;K523N |             3 |
