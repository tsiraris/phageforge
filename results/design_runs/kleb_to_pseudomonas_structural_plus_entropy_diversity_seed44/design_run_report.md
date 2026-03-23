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
- Total candidates evaluated: `1964`
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
| min      |       0.227252 |          0.227252 |
| mean     |       0.281352 |          0.281352 |
| median   |       0.279489 |          0.279489 |
| max      |       0.324938 |          0.324938 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       5 |
|             2 |      20 |
|             3 |      25 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 11.2367 |
| unique_mutation_site_count_top_candidates | 34      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998577    |
| mean_seed_cosine_similarity   | 0.999092    |
| median_seed_cosine_similarity | 0.999103    |
| max_seed_cosine_similarity    | 0.999601    |
| min_seed_novelty_distance     | 0.000399351 |
| mean_seed_novelty_distance    | 0.000907792 |
| median_seed_novelty_distance  | 0.000896901 |
| max_seed_novelty_distance     | 0.00142288  |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000399411 |
| mean_novelty_distance                     | 0.000907829 |
| median_novelty_distance                   | 0.00089699  |
| max_novelty_distance                      | 0.00142306  |
| min_nearest_neighbor_cosine_similarity    | 0.998577    |
| mean_nearest_neighbor_cosine_similarity   | 0.999092    |
| median_nearest_neighbor_cosine_similarity | 0.999103    |
| max_nearest_neighbor_cosine_similarity    | 0.999601    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand180 | QFR57578.1                    | Klebsiella                    |                             0.998864 |        0.00113606  |
| round3_cand249 | QFR57578.1                    | Klebsiella                    |                             0.999601 |        0.000399411 |
| round4_cand6   | QFR57578.1                    | Klebsiella                    |                             0.999154 |        0.000846267 |
| round4_cand219 | QFR57578.1                    | Klebsiella                    |                             0.998749 |        0.00125074  |
| round4_cand7   | QFR57578.1                    | Klebsiella                    |                             0.998982 |        0.00101781  |
| round4_cand90  | QFR57578.1                    | Klebsiella                    |                             0.999105 |        0.000895143 |
| round4_cand63  | QFR57578.1                    | Klebsiella                    |                             0.999431 |        0.000568628 |
| round4_cand135 | QFR57578.1                    | Klebsiella                    |                             0.998853 |        0.00114679  |
| round4_cand65  | QFR57578.1                    | Klebsiella                    |                             0.999126 |        0.000873864 |
| round4_cand197 | QFR57578.1                    | Klebsiella                    |                             0.998957 |        0.00104272  |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand180 | Klebsiella    | Pseudomonas   | Enterobacter |       0.324938 |          0.324938 |                 0.998864 |             0.00113589  | M40R;K105P;C574P  |             3 |
| round3_cand249 | Klebsiella    | Pseudomonas   | Enterobacter |       0.227252 |          0.227252 |                 0.999601 |             0.000399351 | M65R;E427R;M637V  |             3 |
| round4_cand6   | Klebsiella    | Pseudomonas   | Enterobacter |       0.303732 |          0.303732 |                 0.999154 |             0.000846446 | K128R;K374V       |             2 |
| round4_cand219 | Klebsiella    | Pseudomonas   | Enterobacter |       0.295185 |          0.295185 |                 0.998749 |             0.00125068  | H271L;K374L       |             2 |
| round4_cand7   | Klebsiella    | Pseudomonas   | Enterobacter |       0.290666 |          0.290666 |                 0.998982 |             0.00101757  | M65L;K260A;K325D  |             3 |
| round4_cand90  | Klebsiella    | Pseudomonas   | Enterobacter |       0.290484 |          0.290484 |                 0.999105 |             0.000895083 | K374I             |             1 |
| round4_cand63  | Klebsiella    | Pseudomonas   | Enterobacter |       0.289322 |          0.289322 |                 0.999432 |             0.00056833  | W197A;K284P       |             2 |
| round4_cand135 | Klebsiella    | Pseudomonas   | Enterobacter |       0.288943 |          0.288943 |                 0.998853 |             0.00114691  | K105P;C574N;M643T |             3 |
| round4_cand65  | Klebsiella    | Pseudomonas   | Enterobacter |       0.288866 |          0.288866 |                 0.999126 |             0.000873864 | M40A;K374G        |             2 |
| round4_cand197 | Klebsiella    | Pseudomonas   | Enterobacter |       0.288444 |          0.288444 |                 0.998958 |             0.00104237  | E154V;M334T;C574P |             3 |
