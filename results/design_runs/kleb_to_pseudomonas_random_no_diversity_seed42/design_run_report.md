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
- Seed: `42`
- Total candidates evaluated: `3263`
- Top candidates saved: `50`
- Diversity enabled: `False`
- Diversity min distance: `8`
- Position selection strategy: `random`
- Position pool size: `N/A`

## Target score summary

| metric   |   target_score |
|:---------|---------------:|
| min      |       0.421203 |
| mean     |       0.436321 |
| median   |       0.430975 |
| max      |       0.487475 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Pseudomonas  |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       8 |
|             2 |      12 |
|             3 |      18 |
|             4 |      12 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  15.1029 |
| unique_mutation_site_count_top_candidates | 115      |

## Nearest-neighbor novelty summary

| metric                                    |      value |
|:------------------------------------------|-----------:|
| min_novelty_distance                      | 0.00121677 |
| mean_novelty_distance                     | 0.00182781 |
| median_novelty_distance                   | 0.00186813 |
| max_novelty_distance                      | 0.00234461 |
| min_nearest_neighbor_cosine_similarity    | 0.997655   |
| mean_nearest_neighbor_cosine_similarity   | 0.998172   |
| median_nearest_neighbor_cosine_similarity | 0.998132   |
| max_nearest_neighbor_cosine_similarity    | 0.998783   |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round6_cand20  | QFR57578.1                    | Klebsiella                    |                             0.997913 |         0.00208724 |
| round6_cand8   | QFR57578.1                    | Klebsiella                    |                             0.99803  |         0.00197029 |
| round6_cand4   | QFR57578.1                    | Klebsiella                    |                             0.998077 |         0.00192279 |
| round6_cand33  | QFR57578.1                    | Klebsiella                    |                             0.998053 |         0.00194657 |
| round6_cand15  | QFR57578.1                    | Klebsiella                    |                             0.997854 |         0.00214571 |
| round6_cand160 | QFR57578.1                    | Klebsiella                    |                             0.997869 |         0.00213099 |
| round6_cand171 | QFR57578.1                    | Klebsiella                    |                             0.99864  |         0.0013597  |
| round6_cand454 | QFR57578.1                    | Klebsiella                    |                             0.997835 |         0.0021646  |
| round6_cand250 | QFR57578.1                    | Klebsiella                    |                             0.998124 |         0.00187635 |
| round6_cand30  | QFR57578.1                    | Klebsiella                    |                             0.998272 |         0.00172782 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score | mutations               |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|:------------------------|--------------:|
| round6_cand20  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.487475 | G167L;T234P;A405T       |             3 |
| round6_cand8   | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.469624 | D67R;D158V;Q586A        |             3 |
| round6_cand4   | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.459942 | D136T;N351A;K355G;Q417K |             4 |
| round6_cand33  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.455421 | T17A;T98R               |             2 |
| round6_cand15  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.454737 | D52Y;N216P              |             2 |
| round6_cand160 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.453339 | A75V;N169R;E526Q;I580Q  |             4 |
| round6_cand171 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.450593 | M296A;E526P;S609A       |             3 |
| round6_cand454 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.450361 | F94P;N608T              |             2 |
| round6_cand250 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.448593 | Y96A;W632H              |             2 |
| round6_cand30  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.445752 | E146P;A395G             |             2 |
