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
- Seed: `42`
- Total candidates evaluated: `1952`
- Top candidates saved: `50`
- Diversity enabled: `True`
- Diversity min distance: `20`
- Position selection strategy: `entropy`
- Position pool size: `32`
- Proposal batch size: `16`
- Novelty penalty lambda: `0.0`
- Novelty penalty schedule: `constant`

## Target score summary

| metric   |   target_score |   selection_score |
|:---------|---------------:|------------------:|
| min      |       0.596256 |          0.596256 |
| mean     |       0.603935 |          0.603935 |
| median   |       0.602803 |          0.602803 |
| max      |       0.626677 |          0.626677 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       9 |
|             2 |      20 |
|             3 |      21 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 11.2359 |
| unique_mutation_site_count_top_candidates | 34      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998582    |
| mean_seed_cosine_similarity   | 0.999092    |
| median_seed_cosine_similarity | 0.999107    |
| max_seed_cosine_similarity    | 0.999469    |
| min_seed_novelty_distance     | 0.000531375 |
| mean_seed_novelty_distance    | 0.000908493 |
| median_seed_novelty_distance  | 0.000892878 |
| max_seed_novelty_distance     | 0.00141823  |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000523865 |
| mean_novelty_distance                     | 0.000896707 |
| median_novelty_distance                   | 0.000892937 |
| max_novelty_distance                      | 0.00141835  |
| min_nearest_neighbor_cosine_similarity    | 0.998582    |
| mean_nearest_neighbor_cosine_similarity   | 0.999103    |
| median_nearest_neighbor_cosine_similarity | 0.999107    |
| max_nearest_neighbor_cosine_similarity    | 0.999476    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand102 | WWD14686.1                    | Klebsiella                    |                             0.999211 |        0.000788689 |
| round4_cand97  | WWD14686.1                    | Klebsiella                    |                             0.999366 |        0.000633717 |
| round4_cand46  | QFR57578.1                    | Klebsiella                    |                             0.998582 |        0.00141835  |
| round4_cand149 | QFR57578.1                    | Klebsiella                    |                             0.999268 |        0.000731647 |
| round4_cand134 | QFR57578.1                    | Klebsiella                    |                             0.999218 |        0.000781894 |
| round4_cand22  | QFR57578.1                    | Klebsiella                    |                             0.998613 |        0.00138682  |
| round4_cand205 | QFR57578.1                    | Klebsiella                    |                             0.999048 |        0.000951648 |
| round4_cand37  | QFR57578.1                    | Klebsiella                    |                             0.998745 |        0.00125486  |
| round4_cand595 | QFR57578.1                    | Klebsiella                    |                             0.998913 |        0.00108743  |
| round4_cand315 | WWD14686.1                    | Klebsiella                    |                             0.999037 |        0.000962973 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand102 | Klebsiella    | Enterobacter  | Enterobacter |       0.626677 |          0.626677 |                 0.999164 |             0.000836492 | K22V;W329N;K335A  |             3 |
| round4_cand97  | Klebsiella    | Enterobacter  | Enterobacter |       0.61574  |          0.61574  |                 0.999344 |             0.000655651 | C240N;W505N;H593L |             3 |
| round4_cand46  | Klebsiella    | Enterobacter  | Enterobacter |       0.613975 |          0.613975 |                 0.998582 |             0.00141823  | K335T;W505T       |             2 |
| round4_cand149 | Klebsiella    | Enterobacter  | Enterobacter |       0.611886 |          0.611886 |                 0.999268 |             0.000731885 | K128R;W329D;K374L |             3 |
| round4_cand134 | Klebsiella    | Enterobacter  | Enterobacter |       0.611271 |          0.611271 |                 0.999218 |             0.000781894 | K335S;W505N;C574S |             3 |
| round4_cand22  | Klebsiella    | Enterobacter  | Enterobacter |       0.610218 |          0.610218 |                 0.998613 |             0.00138652  | K22G;H271L;C516N  |             3 |
| round4_cand205 | Klebsiella    | Enterobacter  | Enterobacter |       0.609873 |          0.609873 |                 0.999048 |             0.000951648 | W505R;K523T;M637S |             3 |
| round4_cand37  | Klebsiella    | Enterobacter  | Enterobacter |       0.60973  |          0.60973  |                 0.998745 |             0.00125462  | K22G;M643R        |             2 |
| round4_cand595 | Klebsiella    | Enterobacter  | Enterobacter |       0.609675 |          0.609675 |                 0.998913 |             0.00108713  | K260P;W505G       |             2 |
| round4_cand315 | Klebsiella    | Enterobacter  | Enterobacter |       0.608736 |          0.608736 |                 0.999027 |             0.000973463 | K128R;P291A;W329N |             3 |
