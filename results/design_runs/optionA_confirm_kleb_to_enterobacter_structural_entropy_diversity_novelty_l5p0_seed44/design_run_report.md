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
- Seed: `44`
- Total candidates evaluated: `1966`
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
| min      |       0.60056  |          -4.3949  |
| mean     |       0.609993 |          -4.38523 |
| median   |       0.607647 |          -4.38765 |
| max      |       0.644352 |          -4.3497  |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       7 |
|             2 |      17 |
|             3 |      26 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 13.0147 |
| unique_mutation_site_count_top_candidates | 35      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998452    |
| mean_seed_cosine_similarity   | 0.999044    |
| median_seed_cosine_similarity | 0.999065    |
| max_seed_cosine_similarity    | 0.999439    |
| min_seed_novelty_distance     | 0.000561476 |
| mean_seed_novelty_distance    | 0.000956339 |
| median_seed_novelty_distance  | 0.000935435 |
| max_seed_novelty_distance     | 0.00154787  |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000561357 |
| mean_novelty_distance                     | 0.00094448  |
| median_novelty_distance                   | 0.000926167 |
| max_novelty_distance                      | 0.00154299  |
| min_nearest_neighbor_cosine_similarity    | 0.998457    |
| mean_nearest_neighbor_cosine_similarity   | 0.999056    |
| median_nearest_neighbor_cosine_similarity | 0.999074    |
| max_nearest_neighbor_cosine_similarity    | 0.999439    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand58  | WWD14686.1                    | Klebsiella                    |                             0.998879 |        0.00112134  |
| round4_cand60  | QFR57578.1                    | Klebsiella                    |                             0.998808 |        0.00119221  |
| round4_cand0   | QFR57578.1                    | Klebsiella                    |                             0.99906  |        0.000939727 |
| round4_cand498 | QFR57578.1                    | Klebsiella                    |                             0.998979 |        0.00102103  |
| round4_cand243 | QFR57578.1                    | Klebsiella                    |                             0.999069 |        0.000931144 |
| round4_cand4   | WWD14686.1                    | Klebsiella                    |                             0.999079 |        0.000920534 |
| round4_cand45  | WWD14686.1                    | Klebsiella                    |                             0.99922  |        0.000779986 |
| round4_cand375 | QFR57578.1                    | Klebsiella                    |                             0.999292 |        0.000707984 |
| round4_cand12  | QFR57578.1                    | Klebsiella                    |                             0.999023 |        0.000977159 |
| round4_cand111 | QFR57578.1                    | Klebsiella                    |                             0.99928  |        0.000720322 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand58  | Klebsiella    | Enterobacter  | Enterobacter |       0.644352 |          -4.3497  |                 0.99881  |             0.00118965  | K128A;W329N;K599R |             3 |
| round4_cand60  | Klebsiella    | Enterobacter  | Enterobacter |       0.638388 |          -4.35565 |                 0.998808 |             0.00119215  | K325T;W329D;C568R |             3 |
| round4_cand0   | Klebsiella    | Enterobacter  | Enterobacter |       0.626493 |          -4.36881 |                 0.99906  |             0.000939786 | H27A;M226N;W329N  |             3 |
| round4_cand498 | Klebsiella    | Enterobacter  | Enterobacter |       0.624352 |          -4.37054 |                 0.998979 |             0.00102109  | W197Q;K455Q;W505G |             3 |
| round4_cand243 | Klebsiella    | Enterobacter  | Enterobacter |       0.621526 |          -4.37382 |                 0.999069 |             0.000931084 | K22A;K523D        |             2 |
| round4_cand4   | Klebsiella    | Enterobacter  | Enterobacter |       0.617284 |          -4.37789 |                 0.999035 |             0.000965297 | K523R;K599V       |             2 |
| round4_cand45  | Klebsiella    | Enterobacter  | Enterobacter |       0.617689 |          -4.37801 |                 0.99914  |             0.000859618 | K523T;W620T;M637P |             3 |
| round4_cand375 | Klebsiella    | Enterobacter  | Enterobacter |       0.617313 |          -4.37915 |                 0.999292 |             0.000707746 | M334T;W373N       |             2 |
| round4_cand12  | Klebsiella    | Enterobacter  | Enterobacter |       0.615324 |          -4.37979 |                 0.999023 |             0.00097692  | K335P;K523N;E540L |             3 |
| round4_cand111 | Klebsiella    | Enterobacter  | Enterobacter |       0.61502  |          -4.38138 |                 0.99928  |             0.000720322 | K325G;W620R       |             2 |
