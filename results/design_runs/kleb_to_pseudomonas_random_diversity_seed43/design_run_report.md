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
- Diversity enabled: `True`
- Diversity min distance: `8`
- Position selection strategy: `random`
- Position pool size: `N/A`

## Target score summary

| metric   |   target_score |
|:---------|---------------:|
| min      |       0.305131 |
| mean     |       0.349559 |
| median   |       0.33798  |
| max      |       0.492084 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      26 |
| Pseudomonas  |      24 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       1 |
|             2 |       2 |
|             3 |      12 |
|             4 |      35 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      |  29.898 |
| unique_mutation_site_count_top_candidates | 164     |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000342429 |
| mean_novelty_distance                     | 0.000974201 |
| median_novelty_distance                   | 0.000716805 |
| max_novelty_distance                      | 0.00237817  |
| min_nearest_neighbor_cosine_similarity    | 0.997622    |
| mean_nearest_neighbor_cosine_similarity   | 0.999026    |
| median_nearest_neighbor_cosine_similarity | 0.999283    |
| max_nearest_neighbor_cosine_similarity    | 0.999658    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round6_cand78  | QFR57578.1                    | Klebsiella                    |                             0.997825 |         0.00217462 |
| round6_cand21  | QFR57578.1                    | Klebsiella                    |                             0.997922 |         0.00207794 |
| round6_cand43  | QFR57578.1                    | Klebsiella                    |                             0.997948 |         0.00205171 |
| round6_cand142 | QFR57578.1                    | Klebsiella                    |                             0.998877 |         0.00112307 |
| round6_cand50  | QFR57578.1                    | Klebsiella                    |                             0.998006 |         0.00199425 |
| round6_cand629 | QFR57578.1                    | Klebsiella                    |                             0.998835 |         0.00116527 |
| round6_cand547 | QFR57578.1                    | Klebsiella                    |                             0.998524 |         0.00147605 |
| round6_cand1   | QFR57578.1                    | Klebsiella                    |                             0.998161 |         0.0018388  |
| round6_cand11  | QFR57578.1                    | Klebsiella                    |                             0.997622 |         0.00237817 |
| round6_cand336 | QFR57578.1                    | Klebsiella                    |                             0.998851 |         0.00114864 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score | mutations               |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|:------------------------|--------------:|
| round6_cand78  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.492084 | G150L;G167L;L569P       |             3 |
| round6_cand21  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.449294 | S112T;G173A;G243L;G309A |             4 |
| round6_cand43  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.418009 | G31L;C162A;A497T;Y615P  |             4 |
| round6_cand142 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.409782 | Y269L;Q352P;L501R       |             3 |
| round6_cand50  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.404431 | H49T;N263G;K335P;V575T  |             4 |
| round6_cand629 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.401324 | N21L;N294I;D525S        |             3 |
| round6_cand547 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.392015 | N130A;I185S;R326E;N445L |             4 |
| round6_cand1   | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.390717 | I315A;M334S;A530S;M601F |             4 |
| round6_cand11  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.388218 | G9A;R64D;P475V;S558T    |             4 |
| round6_cand336 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.382645 | N216G;N382L             |             2 |
