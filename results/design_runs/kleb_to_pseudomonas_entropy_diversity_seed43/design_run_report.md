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
- Total candidates evaluated: `3241`
- Top candidates saved: `50`
- Diversity enabled: `True`
- Diversity min distance: `8`
- Position selection strategy: `entropy`
- Position pool size: `N/A`

## Target score summary

| metric   |   target_score |
|:---------|---------------:|
| min      |       0.307266 |
| mean     |       0.359079 |
| median   |       0.349567 |
| max      |       0.456488 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      33 |
| Pseudomonas  |      17 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       1 |
|             2 |       1 |
|             3 |       6 |
|             4 |      42 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 23.3478 |
| unique_mutation_site_count_top_candidates | 56      |

## Nearest-neighbor novelty summary

| metric                                    |      value |
|:------------------------------------------|-----------:|
| min_novelty_distance                      | 0.00121993 |
| mean_novelty_distance                     | 0.00231878 |
| median_novelty_distance                   | 0.00234833 |
| max_novelty_distance                      | 0.00337815 |
| min_nearest_neighbor_cosine_similarity    | 0.996622   |
| mean_nearest_neighbor_cosine_similarity   | 0.997681   |
| median_nearest_neighbor_cosine_similarity | 0.997652   |
| max_nearest_neighbor_cosine_similarity    | 0.99878    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round6_cand97  | QFR57578.1                    | Klebsiella                    |                             0.997344 |         0.00265598 |
| round6_cand13  | QFR57578.1                    | Klebsiella                    |                             0.996879 |         0.00312054 |
| round6_cand260 | QFR57578.1                    | Klebsiella                    |                             0.996823 |         0.00317717 |
| round6_cand65  | QFR57578.1                    | Klebsiella                    |                             0.997154 |         0.00284564 |
| round6_cand161 | QFR57578.1                    | Klebsiella                    |                             0.997415 |         0.00258505 |
| round6_cand369 | QFR57578.1                    | Klebsiella                    |                             0.997669 |         0.00233078 |
| round6_cand157 | QFR57578.1                    | Klebsiella                    |                             0.99765  |         0.00234973 |
| round6_cand254 | QFR57578.1                    | Klebsiella                    |                             0.99701  |         0.00299001 |
| round6_cand511 | QFR57578.1                    | Klebsiella                    |                             0.997392 |         0.00260764 |
| round6_cand272 | QFR57578.1                    | Klebsiella                    |                             0.997049 |         0.00295103 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score | mutations               |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|:------------------------|--------------:|
| round6_cand97  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.456488 | H27P;C132P;C240A;M641A  |             4 |
| round6_cand13  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.433889 | K155R;P256Y;K284V       |             3 |
| round6_cand260 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.425109 | K260A;H271A;K455G;K599Q |             4 |
| round6_cand65  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.424971 | K128V;C162R;H253I;C587V |             4 |
| round6_cand161 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.420565 | K22N;C132R;H253I;F577A  |             4 |
| round6_cand369 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.415803 | H253I;H271A;K455A       |             3 |
| round6_cand157 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.414651 | K128R;H271N;K325A;M573V |             4 |
| round6_cand254 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.410737 | K455R;C535V;M637R       |             3 |
| round6_cand511 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.404666 | M334P;K455A;C587I       |             3 |
| round6_cand272 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.402733 | M40R;M567A;C568T;H578T  |             4 |
