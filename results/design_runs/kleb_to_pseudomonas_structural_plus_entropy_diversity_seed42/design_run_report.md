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
- Novelty penalty lambda: `0.0`
- Novelty penalty schedule: `constant`

## Target score summary

| metric   |   target_score |   selection_score |
|:---------|---------------:|------------------:|
| min      |       0.262199 |          0.262199 |
| mean     |       0.272894 |          0.272894 |
| median   |       0.271641 |          0.271641 |
| max      |       0.297622 |          0.297622 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       4 |
|             2 |      25 |
|             3 |      21 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 12.4531 |
| unique_mutation_site_count_top_candidates | 39      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998509    |
| mean_seed_cosine_similarity   | 0.99916     |
| median_seed_cosine_similarity | 0.999193    |
| max_seed_cosine_similarity    | 0.999492    |
| min_seed_novelty_distance     | 0.000508487 |
| mean_seed_novelty_distance    | 0.000839578 |
| median_seed_novelty_distance  | 0.000806689 |
| max_seed_novelty_distance     | 0.00149137  |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000508428 |
| mean_novelty_distance                     | 0.000839548 |
| median_novelty_distance                   | 0.000806749 |
| max_novelty_distance                      | 0.00149131  |
| min_nearest_neighbor_cosine_similarity    | 0.998509    |
| mean_nearest_neighbor_cosine_similarity   | 0.99916     |
| median_nearest_neighbor_cosine_similarity | 0.999193    |
| max_nearest_neighbor_cosine_similarity    | 0.999492    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand48  | QFR57578.1                    | Klebsiella                    |                             0.999051 |        0.000949144 |
| round4_cand148 | QFR57578.1                    | Klebsiella                    |                             0.998715 |        0.00128472  |
| round4_cand116 | QFR57578.1                    | Klebsiella                    |                             0.999095 |        0.00090462  |
| round4_cand229 | QFR57578.1                    | Klebsiella                    |                             0.999129 |        0.000870943 |
| round4_cand71  | QFR57578.1                    | Klebsiella                    |                             0.999012 |        0.000987649 |
| round4_cand230 | QFR57578.1                    | Klebsiella                    |                             0.999132 |        0.00086832  |
| round4_cand22  | QFR57578.1                    | Klebsiella                    |                             0.999215 |        0.000785351 |
| round4_cand173 | QFR57578.1                    | Klebsiella                    |                             0.998942 |        0.00105751  |
| round4_cand178 | QFR57578.1                    | Klebsiella                    |                             0.998509 |        0.00149131  |
| round4_cand41  | QFR57578.1                    | Klebsiella                    |                             0.999261 |        0.000739455 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand48  | Klebsiella    | Pseudomonas   | Enterobacter |       0.297622 |          0.297622 |                 0.999051 |             0.000949204 | E154F;K284A;K455G |             3 |
| round4_cand148 | Klebsiella    | Pseudomonas   | Enterobacter |       0.296648 |          0.296648 |                 0.998715 |             0.00128478  | K325T;K374I;M643R |             3 |
| round4_cand116 | Klebsiella    | Pseudomonas   | Enterobacter |       0.293147 |          0.293147 |                 0.999095 |             0.00090462  | M226G;K374I;C568S |             3 |
| round4_cand229 | Klebsiella    | Pseudomonas   | Enterobacter |       0.288437 |          0.288437 |                 0.999129 |             0.000871003 | M226A;R350V;M643P |             3 |
| round4_cand71  | Klebsiella    | Pseudomonas   | Enterobacter |       0.287611 |          0.287611 |                 0.999012 |             0.00098753  | K374I;K599I       |             2 |
| round4_cand230 | Klebsiella    | Pseudomonas   | Enterobacter |       0.284085 |          0.284085 |                 0.999132 |             0.000868499 | K455R;C568P;M641A |             3 |
| round4_cand22  | Klebsiella    | Pseudomonas   | Enterobacter |       0.28165  |          0.28165  |                 0.999215 |             0.000785351 | H27P;M139G;M641P  |             3 |
| round4_cand173 | Klebsiella    | Pseudomonas   | Enterobacter |       0.281055 |          0.281055 |                 0.998942 |             0.00105774  | K325A;M643R       |             2 |
| round4_cand178 | Klebsiella    | Pseudomonas   | Enterobacter |       0.279786 |          0.279786 |                 0.998509 |             0.00149137  | H253V;M643P       |             2 |
| round4_cand41  | Klebsiella    | Pseudomonas   | Enterobacter |       0.279763 |          0.279763 |                 0.999261 |             0.000739455 | W602I;M637R;M641Q |             3 |
