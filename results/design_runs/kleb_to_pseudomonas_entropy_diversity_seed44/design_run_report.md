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
- Total candidates evaluated: `3236`
- Top candidates saved: `50`
- Diversity enabled: `True`
- Diversity min distance: `8`
- Position selection strategy: `entropy`
- Position pool size: `N/A`

## Target score summary

| metric   |   target_score |
|:---------|---------------:|
| min      |       0.28787  |
| mean     |       0.340128 |
| median   |       0.335305 |
| max      |       0.409178 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      45 |
| Pseudomonas  |       5 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             2 |       2 |
|             3 |       6 |
|             4 |      42 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 26.5722 |
| unique_mutation_site_count_top_candidates | 60      |

## Nearest-neighbor novelty summary

| metric                                    |      value |
|:------------------------------------------|-----------:|
| min_novelty_distance                      | 0.00134051 |
| mean_novelty_distance                     | 0.00256844 |
| median_novelty_distance                   | 0.00253624 |
| max_novelty_distance                      | 0.00425559 |
| min_nearest_neighbor_cosine_similarity    | 0.995744   |
| mean_nearest_neighbor_cosine_similarity   | 0.997432   |
| median_nearest_neighbor_cosine_similarity | 0.997464   |
| max_nearest_neighbor_cosine_similarity    | 0.998659   |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round6_cand26  | QFR57578.1                    | Klebsiella                    |                             0.997237 |         0.00276309 |
| round6_cand256 | QFR57578.1                    | Klebsiella                    |                             0.997031 |         0.00296926 |
| round6_cand218 | QFR57578.1                    | Klebsiella                    |                             0.997079 |         0.0029211  |
| round6_cand447 | QFR57578.1                    | Klebsiella                    |                             0.997113 |         0.00288731 |
| round6_cand379 | QFR57578.1                    | Klebsiella                    |                             0.997543 |         0.00245702 |
| round6_cand141 | QFR57578.1                    | Klebsiella                    |                             0.996123 |         0.00387675 |
| round6_cand488 | QFR57578.1                    | Klebsiella                    |                             0.997228 |         0.00277215 |
| round6_cand92  | QFR57578.1                    | Klebsiella                    |                             0.99808  |         0.00191981 |
| round6_cand529 | QFR57578.1                    | Klebsiella                    |                             0.996332 |         0.00366843 |
| round6_cand162 | QFR57578.1                    | Klebsiella                    |                             0.996122 |         0.00387764 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score | mutations               |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|:------------------------|--------------:|
| round6_cand26  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.409178 | M40R;K523A;K561A        |             3 |
| round6_cand256 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.408348 | K105R;K128A;K561A;C587R |             4 |
| round6_cand218 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.404613 | H253I;K455A;M567T;M573V |             4 |
| round6_cand447 | Klebsiella    | Pseudomonas   | Enterobacter |       0.402663 | H253F;R350P;E427R;M641A |             4 |
| round6_cand379 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.398427 | M40R;M65V;K105N;H593W   |             4 |
| round6_cand141 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.39829  | H271G;K284R;M334A;M369A |             4 |
| round6_cand488 | Klebsiella    | Pseudomonas   | Enterobacter |       0.393506 | H27A;M139P;C240V;C568P  |             4 |
| round6_cand92  | Klebsiella    | Pseudomonas   | Enterobacter |       0.393077 | Y444L;H499G;W602P       |             3 |
| round6_cand529 | Klebsiella    | Pseudomonas   | Enterobacter |       0.391936 | K335P;W527G;C535R;K599V |             4 |
| round6_cand162 | Klebsiella    | Pseudomonas   | Enterobacter |       0.38512  | H27R;M65G;K128R;M226G   |             4 |
