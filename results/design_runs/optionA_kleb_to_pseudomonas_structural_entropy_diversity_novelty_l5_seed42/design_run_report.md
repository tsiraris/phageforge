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
- Seed: `42`
- Total candidates evaluated: `1951`
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
| min      |       0.261745 |          -4.73295 |
| mean     |       0.272832 |          -4.7228  |
| median   |       0.270746 |          -4.72565 |
| max      |       0.297622 |          -4.69693 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       4 |
|             2 |      30 |
|             3 |      16 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 12.3216 |
| unique_mutation_site_count_top_candidates | 36      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998509    |
| mean_seed_cosine_similarity   | 0.999127    |
| median_seed_cosine_similarity | 0.999165    |
| max_seed_cosine_similarity    | 0.999432    |
| min_seed_novelty_distance     | 0.000568211 |
| mean_seed_novelty_distance    | 0.000872772 |
| median_seed_novelty_distance  | 0.000835001 |
| max_seed_novelty_distance     | 0.00149137  |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000567794 |
| mean_novelty_distance                     | 0.000872778 |
| median_novelty_distance                   | 0.000835001 |
| max_novelty_distance                      | 0.00149131  |
| min_nearest_neighbor_cosine_similarity    | 0.998509    |
| mean_nearest_neighbor_cosine_similarity   | 0.999127    |
| median_nearest_neighbor_cosine_similarity | 0.999165    |
| max_nearest_neighbor_cosine_similarity    | 0.999432    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand148 | QFR57578.1                    | Klebsiella                    |                             0.998715 |        0.00128472  |
| round4_cand227 | QFR57578.1                    | Klebsiella                    |                             0.998912 |        0.00108778  |
| round4_cand48  | QFR57578.1                    | Klebsiella                    |                             0.999051 |        0.000949144 |
| round4_cand116 | QFR57578.1                    | Klebsiella                    |                             0.999095 |        0.00090462  |
| round4_cand286 | QFR57578.1                    | Klebsiella                    |                             0.999195 |        0.000804722 |
| round4_cand71  | QFR57578.1                    | Klebsiella                    |                             0.999012 |        0.000987649 |
| round4_cand178 | QFR57578.1                    | Klebsiella                    |                             0.998509 |        0.00149131  |
| round4_cand173 | QFR57578.1                    | Klebsiella                    |                             0.998942 |        0.00105751  |
| round4_cand574 | QFR57578.1                    | Klebsiella                    |                             0.999358 |        0.000642478 |
| round4_cand22  | QFR57578.1                    | Klebsiella                    |                             0.999215 |        0.000785351 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand148 | Klebsiella    | Pseudomonas   | Enterobacter |       0.296648 |          -4.69693 |                 0.998715 |             0.00128478  | K325T;K374I;M643R |             3 |
| round4_cand227 | Klebsiella    | Pseudomonas   | Enterobacter |       0.297142 |          -4.69742 |                 0.998913 |             0.00108743  | K22P;K374A        |             2 |
| round4_cand48  | Klebsiella    | Pseudomonas   | Enterobacter |       0.297622 |          -4.69763 |                 0.999051 |             0.000949204 | E154F;K284A;K455G |             3 |
| round4_cand116 | Klebsiella    | Pseudomonas   | Enterobacter |       0.293147 |          -4.70233 |                 0.999095 |             0.00090462  | M226G;K374I;C568S |             3 |
| round4_cand286 | Klebsiella    | Pseudomonas   | Enterobacter |       0.291787 |          -4.70419 |                 0.999195 |             0.000804901 | M226A;K374G;C587R |             3 |
| round4_cand71  | Klebsiella    | Pseudomonas   | Enterobacter |       0.287611 |          -4.70745 |                 0.999012 |             0.00098753  | K374I;K599I       |             2 |
| round4_cand178 | Klebsiella    | Pseudomonas   | Enterobacter |       0.279786 |          -4.71276 |                 0.998509 |             0.00149137  | H253V;M643P       |             2 |
| round4_cand173 | Klebsiella    | Pseudomonas   | Enterobacter |       0.281055 |          -4.71366 |                 0.998942 |             0.00105774  | K325A;M643R       |             2 |
| round4_cand574 | Klebsiella    | Pseudomonas   | Enterobacter |       0.282502 |          -4.71429 |                 0.999358 |             0.000642478 | H27P;K374A        |             2 |
| round4_cand22  | Klebsiella    | Pseudomonas   | Enterobacter |       0.28165  |          -4.71442 |                 0.999215 |             0.000785351 | H27P;M139G;M641P  |             3 |
