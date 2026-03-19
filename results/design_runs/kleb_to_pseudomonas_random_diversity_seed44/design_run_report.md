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
- Seed: `44`
- Total candidates evaluated: `3263`
- Top candidates saved: `50`
- Diversity enabled: `True`
- Diversity min distance: `8`
- Position selection strategy: `random`
- Position pool size: `N/A`

## Target score summary

| metric   |   target_score |
|:---------|---------------:|
| min      |       0.370167 |
| mean     |       0.40289  |
| median   |       0.399884 |
| max      |       0.471159 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Pseudomonas  |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             3 |       2 |
|             4 |      48 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  25.7527 |
| unique_mutation_site_count_top_candidates | 172      |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000803471 |
| mean_novelty_distance                     | 0.00111987  |
| median_novelty_distance                   | 0.00104854  |
| max_novelty_distance                      | 0.00161457  |
| min_nearest_neighbor_cosine_similarity    | 0.998385    |
| mean_nearest_neighbor_cosine_similarity   | 0.99888     |
| median_nearest_neighbor_cosine_similarity | 0.998951    |
| max_nearest_neighbor_cosine_similarity    | 0.999197    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round6_cand28  | QFR57578.1                    | Klebsiella                    |                             0.998385 |        0.00161457  |
| round6_cand110 | QFR57578.1                    | Klebsiella                    |                             0.998401 |        0.00159919  |
| round6_cand81  | QFR57578.1                    | Klebsiella                    |                             0.998848 |        0.0011518   |
| round6_cand261 | QFR57578.1                    | Klebsiella                    |                             0.998562 |        0.00143778  |
| round6_cand101 | QFR57578.1                    | Klebsiella                    |                             0.998979 |        0.00102133  |
| round6_cand152 | QFR57578.1                    | Klebsiella                    |                             0.999083 |        0.000916839 |
| round6_cand188 | QFR57578.1                    | Klebsiella                    |                             0.998989 |        0.00101089  |
| round6_cand198 | QFR57578.1                    | Klebsiella                    |                             0.998824 |        0.00117576  |
| round6_cand161 | QFR57578.1                    | Klebsiella                    |                             0.998958 |        0.00104243  |
| round6_cand171 | QFR57578.1                    | Klebsiella                    |                             0.998949 |        0.00105101  |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score | mutations               |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|:------------------------|--------------:|
| round6_cand28  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.471159 | G122I;N409T;V468T       |             3 |
| round6_cand110 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.455021 | P256T;I366A;N382G;Q560G |             4 |
| round6_cand81  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.453388 | D140A;P229G;K374G;E547V |             4 |
| round6_cand261 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.435258 | R68G;Y77A;D449A;N557L   |             4 |
| round6_cand101 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.432993 | S15R;G97S;Y269V;W505T   |             4 |
| round6_cand152 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.432898 | G159V;N524P;I554F;K561I |             4 |
| round6_cand188 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.431662 | N20V;R400V;N409P;Y453G  |             4 |
| round6_cand198 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.43166  | I35P;G243L;G462R;D521Q  |             4 |
| round6_cand161 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.431365 | G2P;V72P;G307V;I407V    |             4 |
| round6_cand171 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.429374 | N272L;N445L;G480P;T503N |             4 |
