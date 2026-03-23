# Design Run Summary

## Run metadata

- Seed protein ID: `QFR57578.1`
- Virus accession: `MN434096`
- Source host: `Klebsiella`
- Target host: `Acinetobacter`
- ESM model: `facebook/esm2_t33_650M_UR50D`
- Rounds: `4`
- Candidates per round: `64`
- Min mutations: `1`
- Max mutations: `3`
- Top-K kept per round: `10`
- Proposal top-K: `8`
- Seed: `43`
- Total candidates evaluated: `1962`
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
| min      |       0.086115 |          0.086115 |
| mean     |       0.101538 |          0.101538 |
| median   |       0.100379 |          0.100379 |
| max      |       0.115573 |          0.115573 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |      12 |
|             2 |      14 |
|             3 |      24 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      |  11.249 |
| unique_mutation_site_count_top_candidates |  41     |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998707    |
| mean_seed_cosine_similarity   | 0.999314    |
| median_seed_cosine_similarity | 0.999363    |
| max_seed_cosine_similarity    | 0.999574    |
| min_seed_novelty_distance     | 0.000426412 |
| mean_seed_novelty_distance    | 0.000686257 |
| median_seed_novelty_distance  | 0.000637412 |
| max_seed_novelty_distance     | 0.00129294  |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000426292 |
| mean_novelty_distance                     | 0.000686272 |
| median_novelty_distance                   | 0.000637323 |
| max_novelty_distance                      | 0.00129282  |
| min_nearest_neighbor_cosine_similarity    | 0.998707    |
| mean_nearest_neighbor_cosine_similarity   | 0.999314    |
| median_nearest_neighbor_cosine_similarity | 0.999363    |
| max_nearest_neighbor_cosine_similarity    | 0.999574    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand72  | QFR57578.1                    | Klebsiella                    |                             0.999059 |        0.000940621 |
| round4_cand563 | QFR57578.1                    | Klebsiella                    |                             0.999077 |        0.000923038 |
| round4_cand52  | QFR57578.1                    | Klebsiella                    |                             0.99933  |        0.000670433 |
| round4_cand49  | QFR57578.1                    | Klebsiella                    |                             0.999275 |        0.000725269 |
| round4_cand71  | QFR57578.1                    | Klebsiella                    |                             0.998844 |        0.00115585  |
| round4_cand38  | QFR57578.1                    | Klebsiella                    |                             0.999324 |        0.000675619 |
| round4_cand53  | QFR57578.1                    | Klebsiella                    |                             0.999282 |        0.000718117 |
| round4_cand104 | QFR57578.1                    | Klebsiella                    |                             0.999156 |        0.000844121 |
| round4_cand365 | QFR57578.1                    | Klebsiella                    |                             0.999112 |        0.000887752 |
| round4_cand403 | QFR57578.1                    | Klebsiella                    |                             0.999462 |        0.000538349 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand72  | Klebsiella    | Acinetobacter | Enterobacter |       0.115573 |          0.115573 |                 0.99906  |             0.000940382 | H271L;C516K;M643T |             3 |
| round4_cand563 | Klebsiella    | Acinetobacter | Enterobacter |       0.086115 |          0.086115 |                 0.999077 |             0.000923097 | M369F;C574A;I580Q |             3 |
| round4_cand52  | Klebsiella    | Acinetobacter | Enterobacter |       0.113117 |          0.113117 |                 0.99933  |             0.000670075 | R317I             |             1 |
| round4_cand49  | Klebsiella    | Acinetobacter | Enterobacter |       0.111992 |          0.111992 |                 0.999275 |             0.000725389 | E154A;R317I       |             2 |
| round4_cand71  | Klebsiella    | Acinetobacter | Enterobacter |       0.1096   |          0.1096   |                 0.998844 |             0.00115615  | P256L;C574V;W602F |             3 |
| round4_cand38  | Klebsiella    | Acinetobacter | Enterobacter |       0.109469 |          0.109469 |                 0.999324 |             0.000675559 | C62T;H271A;R317G  |             3 |
| round4_cand53  | Klebsiella    | Acinetobacter | Enterobacter |       0.108204 |          0.108204 |                 0.999282 |             0.000717938 | M226V;R317I;M637R |             3 |
| round4_cand104 | Klebsiella    | Acinetobacter | Enterobacter |       0.107159 |          0.107159 |                 0.999156 |             0.000843704 | M139D;M334E;E427I |             3 |
| round4_cand365 | Klebsiella    | Acinetobacter | Enterobacter |       0.106071 |          0.106071 |                 0.999112 |             0.000887692 | K325E;M334E;H499T |             3 |
| round4_cand403 | Klebsiella    | Acinetobacter | Enterobacter |       0.104935 |          0.104935 |                 0.999462 |             0.000538468 | R317G;C568E       |             2 |
