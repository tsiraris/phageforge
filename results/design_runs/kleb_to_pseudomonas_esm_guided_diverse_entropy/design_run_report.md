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
- Total candidates evaluated: `1962`
- Top candidates saved: `50`
- Diversity enabled: `True`
- Diversity min distance: `8`
- Position selection strategy: `entropy`
- Position pool size: `N/A`

## Target score summary

| metric   |   target_score |
|:---------|---------------:|
| min      |       0.151137 |
| mean     |       0.244319 |
| median   |       0.27177  |
| max      |       0.307651 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       9 |
|             2 |       7 |
|             3 |      34 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 12.7918 |
| unique_mutation_site_count_top_candidates | 36      |

## Nearest-neighbor novelty summary

| metric                                    |      value |
|:------------------------------------------|-----------:|
| min_novelty_distance                      | 0.00011915 |
| mean_novelty_distance                     | 0.00109571 |
| median_novelty_distance                   | 0.00106439 |
| max_novelty_distance                      | 0.0028069  |
| min_nearest_neighbor_cosine_similarity    | 0.997193   |
| mean_nearest_neighbor_cosine_similarity   | 0.998904   |
| median_nearest_neighbor_cosine_similarity | 0.998936   |
| max_nearest_neighbor_cosine_similarity    | 0.999881   |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand26  | QFR57578.1                    | Klebsiella                    |                             0.99807  |        0.00192958  |
| round4_cand296 | QFR57578.1                    | Klebsiella                    |                             0.998695 |        0.00130463  |
| round4_cand160 | QFR57578.1                    | Klebsiella                    |                             0.998969 |        0.00103146  |
| round4_cand230 | QFR57578.1                    | Klebsiella                    |                             0.999015 |        0.000985026 |
| round4_cand335 | QFR57578.1                    | Klebsiella                    |                             0.998899 |        0.00110126  |
| round4_cand82  | QFR57578.1                    | Klebsiella                    |                             0.998525 |        0.0014751   |
| round4_cand502 | QFR57578.1                    | Klebsiella                    |                             0.999366 |        0.000634134 |
| round4_cand422 | QFR57578.1                    | Klebsiella                    |                             0.998377 |        0.00162292  |
| round4_cand560 | QFR57578.1                    | Klebsiella                    |                             0.999286 |        0.000714362 |
| round4_cand593 | QFR57578.1                    | Klebsiella                    |                             0.998854 |        0.0011459   |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|:------------------|--------------:|
| round4_cand26  | Klebsiella    | Pseudomonas   | Enterobacter |       0.307651 | C132P;K455T;H593F |             3 |
| round4_cand296 | Klebsiella    | Pseudomonas   | Enterobacter |       0.291649 | K284A;K335A;H593V |             3 |
| round4_cand160 | Klebsiella    | Pseudomonas   | Enterobacter |       0.289649 | M40R;E154V;C568P  |             3 |
| round4_cand230 | Klebsiella    | Pseudomonas   | Enterobacter |       0.277863 | M334G;K561P       |             2 |
| round4_cand335 | Klebsiella    | Pseudomonas   | Enterobacter |       0.274878 | K335S;C574P;M637N |             3 |
| round4_cand82  | Klebsiella    | Pseudomonas   | Enterobacter |       0.270464 | K455G;K599A;M641Q |             3 |
| round4_cand502 | Klebsiella    | Pseudomonas   | Enterobacter |       0.263795 | M334G;K374A;K455S |             3 |
| round4_cand422 | Klebsiella    | Pseudomonas   | Enterobacter |       0.263588 | K284R;K325A;C568P |             3 |
| round4_cand560 | Klebsiella    | Pseudomonas   | Enterobacter |       0.257474 | C132P;K260N;K561P |             3 |
| round4_cand593 | Klebsiella    | Pseudomonas   | Enterobacter |       0.24468  | K325A;M637P       |             2 |
