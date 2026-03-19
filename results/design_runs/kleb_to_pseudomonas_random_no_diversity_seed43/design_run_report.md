# Design Run Summary

## Run metadata

- Seed protein ID: `QFR57578.1`
- Virus accession: `MN434096`
- Source host: `Klebsiella`
- Target host: `Pseudomonas`
- ESM model: `facebook/esm2_t33_650M_UR50D`
- Rounds: `6`
- Candidates per round: `64`
- Min mutations: `1`
- Max mutations: `4`
- Top-K kept per round: `10`
- Proposal top-K: `8`
- Seed: `43`
- Total candidates evaluated: `3263`
- Top candidates saved: `50`
- Diversity enabled: `False`
- Diversity min distance: `8`
- Position selection strategy: `random`
- Position pool size: `N/A`

## Target score summary

| metric   |   target_score |
|:---------|---------------:|
| min      |       0.375615 |
| mean     |       0.385217 |
| median   |       0.383528 |
| max      |       0.408484 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Pseudomonas  |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       4 |
|             2 |      12 |
|             3 |      19 |
|             4 |      15 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  23.4571 |
| unique_mutation_site_count_top_candidates | 123      |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000528276 |
| mean_novelty_distance                     | 0.000798683 |
| median_novelty_distance                   | 0.000741422 |
| max_novelty_distance                      | 0.00192565  |
| min_nearest_neighbor_cosine_similarity    | 0.998074    |
| mean_nearest_neighbor_cosine_similarity   | 0.999201    |
| median_nearest_neighbor_cosine_similarity | 0.999259    |
| max_nearest_neighbor_cosine_similarity    | 0.999472    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round6_cand583 | QFR57578.1                    | Klebsiella                    |                             0.999263 |        0.000737309 |
| round6_cand385 | QFR57578.1                    | Klebsiella                    |                             0.999026 |        0.000974    |
| round6_cand118 | QFR57578.1                    | Klebsiella                    |                             0.998488 |        0.00151169  |
| round6_cand28  | QFR57578.1                    | Klebsiella                    |                             0.999387 |        0.000613332 |
| round6_cand495 | QFR57578.1                    | Klebsiella                    |                             0.999132 |        0.000868082 |
| round6_cand64  | QFR57578.1                    | Klebsiella                    |                             0.998949 |        0.00105131  |
| round6_cand622 | QFR57578.1                    | Klebsiella                    |                             0.999219 |        0.00078094  |
| round6_cand128 | QFR57578.1                    | Klebsiella                    |                             0.998862 |        0.00113833  |
| round6_cand239 | QFR57578.1                    | Klebsiella                    |                             0.999124 |        0.000875652 |
| round6_cand21  | QFR57578.1                    | Klebsiella                    |                             0.99937  |        0.00063014  |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score | mutations               |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|:------------------------|--------------:|
| round6_cand583 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.408484 | G129V;N367A;N416R       |             3 |
| round6_cand385 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.404177 | L199E;F415I;N463I;P576A |             4 |
| round6_cand118 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.403989 | S18P;N153L;I192G;H499I  |             4 |
| round6_cand28  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.401737 | S112A;G173D;G243L;G309F |             4 |
| round6_cand495 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.399442 | E50R;G90P;I381P;H593V   |             4 |
| round6_cand64  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.394871 | F94R;S205T;G480P        |             3 |
| round6_cand622 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.39436  | C162R;K335G;Y412A       |             3 |
| round6_cand128 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.393986 | D52F;G311S;N432G;C574R  |             4 |
| round6_cand239 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.393384 | F143L;N313I;G588Y       |             3 |
| round6_cand21  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.393306 | E360G;A426V             |             2 |
