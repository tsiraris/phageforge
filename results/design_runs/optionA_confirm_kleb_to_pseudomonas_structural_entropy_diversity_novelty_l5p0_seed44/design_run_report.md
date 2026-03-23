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
- Novelty penalty lambda: `5.0`
- Novelty penalty schedule: `constant`

## Target score summary

| metric   |   target_score |   selection_score |
|:---------|---------------:|------------------:|
| min      |       0.227252 |          -4.77075 |
| mean     |       0.282654 |          -4.71278 |
| median   |       0.281375 |          -4.71433 |
| max      |       0.324938 |          -4.66938 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       2 |
|             2 |      21 |
|             3 |      27 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 11.2669 |
| unique_mutation_site_count_top_candidates | 38      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998577    |
| mean_seed_cosine_similarity   | 0.999087    |
| median_seed_cosine_similarity | 0.999074    |
| max_seed_cosine_similarity    | 0.999601    |
| min_seed_novelty_distance     | 0.000399351 |
| mean_seed_novelty_distance    | 0.000913043 |
| median_seed_novelty_distance  | 0.000926286 |
| max_seed_novelty_distance     | 0.00142288  |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000399411 |
| mean_novelty_distance                     | 0.000913054 |
| median_novelty_distance                   | 0.000926405 |
| max_novelty_distance                      | 0.00142306  |
| min_nearest_neighbor_cosine_similarity    | 0.998577    |
| mean_nearest_neighbor_cosine_similarity   | 0.999087    |
| median_nearest_neighbor_cosine_similarity | 0.999074    |
| max_nearest_neighbor_cosine_similarity    | 0.999601    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand180 | QFR57578.1                    | Klebsiella                    |                             0.998864 |        0.00113606  |
| round3_cand249 | QFR57578.1                    | Klebsiella                    |                             0.999601 |        0.000399411 |
| round4_cand33  | QFR57578.1                    | Klebsiella                    |                             0.999084 |        0.000915825 |
| round4_cand59  | QFR57578.1                    | Klebsiella                    |                             0.998657 |        0.00134319  |
| round4_cand219 | QFR57578.1                    | Klebsiella                    |                             0.998749 |        0.00125074  |
| round4_cand125 | QFR57578.1                    | Klebsiella                    |                             0.998976 |        0.00102443  |
| round4_cand80  | QFR57578.1                    | Klebsiella                    |                             0.999373 |        0.00062716  |
| round4_cand35  | QFR57578.1                    | Klebsiella                    |                             0.999065 |        0.000935316 |
| round4_cand135 | QFR57578.1                    | Klebsiella                    |                             0.998853 |        0.00114679  |
| round4_cand189 | QFR57578.1                    | Klebsiella                    |                             0.998577 |        0.00142306  |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand180 | Klebsiella    | Pseudomonas   | Enterobacter |       0.324938 |          -4.66938 |                 0.998864 |             0.00113589  | M40R;K105P;C574P  |             3 |
| round3_cand249 | Klebsiella    | Pseudomonas   | Enterobacter |       0.227252 |          -4.77075 |                 0.999601 |             0.000399351 | M65R;E427R;M637V  |             3 |
| round4_cand33  | Klebsiella    | Pseudomonas   | Enterobacter |       0.3099   |          -4.68552 |                 0.999084 |             0.000915825 | M40R;W197A;M226G  |             3 |
| round4_cand59  | Klebsiella    | Pseudomonas   | Enterobacter |       0.298053 |          -4.69523 |                 0.998657 |             0.00134301  | C162N;K374L;C587P |             3 |
| round4_cand219 | Klebsiella    | Pseudomonas   | Enterobacter |       0.295185 |          -4.69856 |                 0.998749 |             0.00125068  | H271L;K374L       |             2 |
| round4_cand125 | Klebsiella    | Pseudomonas   | Enterobacter |       0.295121 |          -4.69976 |                 0.998976 |             0.00102437  | K284T;K325T;K374G |             3 |
| round4_cand80  | Klebsiella    | Pseudomonas   | Enterobacter |       0.294009 |          -4.70285 |                 0.999373 |             0.000627279 | E427R;K599A       |             2 |
| round4_cand35  | Klebsiella    | Pseudomonas   | Enterobacter |       0.290554 |          -4.70477 |                 0.999065 |             0.000935316 | M65G;K374L        |             2 |
| round4_cand135 | Klebsiella    | Pseudomonas   | Enterobacter |       0.288943 |          -4.70532 |                 0.998853 |             0.00114691  | K105P;C574N;M643T |             3 |
| round4_cand189 | Klebsiella    | Pseudomonas   | Enterobacter |       0.287047 |          -4.70584 |                 0.998577 |             0.00142288  | K155R;C240A;K455R |             3 |
