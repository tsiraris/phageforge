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
- Total candidates evaluated: `3235`
- Top candidates saved: `50`
- Diversity enabled: `True`
- Diversity min distance: `8`
- Position selection strategy: `entropy`
- Position pool size: `N/A`

## Target score summary

| metric   |   target_score |
|:---------|---------------:|
| min      |       0.289236 |
| mean     |       0.360793 |
| median   |       0.360796 |
| max      |       0.432541 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      37 |
| Pseudomonas  |      13 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             2 |       2 |
|             3 |       9 |
|             4 |      39 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 25.2441 |
| unique_mutation_site_count_top_candidates | 58      |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000785828 |
| mean_novelty_distance                     | 0.00233917  |
| median_novelty_distance                   | 0.0023309   |
| max_novelty_distance                      | 0.00342369  |
| min_nearest_neighbor_cosine_similarity    | 0.996576    |
| mean_nearest_neighbor_cosine_similarity   | 0.997661    |
| median_nearest_neighbor_cosine_similarity | 0.997669    |
| max_nearest_neighbor_cosine_similarity    | 0.999214    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round6_cand22  | QFR57578.1                    | Klebsiella                    |                             0.997706 |         0.00229383 |
| round6_cand188 | QFR57578.1                    | Klebsiella                    |                             0.99724  |         0.00276035 |
| round6_cand87  | QFR57578.1                    | Klebsiella                    |                             0.996693 |         0.00330681 |
| round6_cand31  | QFR57578.1                    | Klebsiella                    |                             0.997018 |         0.00298154 |
| round6_cand202 | QFR57578.1                    | Klebsiella                    |                             0.997378 |         0.00262213 |
| round6_cand324 | QFR57578.1                    | Klebsiella                    |                             0.996952 |         0.00304806 |
| round6_cand570 | QFR57578.1                    | Klebsiella                    |                             0.996788 |         0.00321221 |
| round6_cand290 | QFR57578.1                    | Klebsiella                    |                             0.997309 |         0.00269079 |
| round6_cand63  | QFR57578.1                    | Klebsiella                    |                             0.996775 |         0.00322509 |
| round6_cand373 | WWD14686.1                    | Klebsiella                    |                             0.997257 |         0.00274301 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score | mutations               |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|:------------------------|--------------:|
| round6_cand22  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.432541 | H253G;K523R;M567A;H578P |             4 |
| round6_cand188 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.432061 | M163A;H253T;H499A;M641P |             4 |
| round6_cand87  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.42767  | M124P;K128R;H253I;C587R |             4 |
| round6_cand31  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.415805 | C240A;W329A;K335N;M637R |             4 |
| round6_cand202 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.413948 | K335Q;C574P             |             2 |
| round6_cand324 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.410447 | H271G;K455T;C535V;M573A |             4 |
| round6_cand570 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.404582 | K128S;K455T;E540N;K561V |             4 |
| round6_cand290 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.404357 | M65V;M226A;C568T        |             3 |
| round6_cand63  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.401551 | C62R;C162V;C535T;C568P  |             4 |
| round6_cand373 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.400569 | H27A;K155L;W197A;K284S  |             4 |
