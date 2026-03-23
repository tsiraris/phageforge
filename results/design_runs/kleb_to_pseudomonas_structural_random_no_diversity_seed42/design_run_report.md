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
- Total candidates evaluated: `1983`
- Top candidates saved: `50`
- Diversity enabled: `False`
- Diversity min distance: `20`
- Position selection strategy: `random`
- Position pool size: `32`
- Proposal batch size: `16`
- Novelty penalty lambda: `0.0`
- Novelty penalty schedule: `constant`

## Target score summary

| metric   |   target_score |   selection_score |
|:---------|---------------:|------------------:|
| min      |       0.252561 |          0.252561 |
| mean     |       0.261295 |          0.261295 |
| median   |       0.258071 |          0.258071 |
| max      |       0.285625 |          0.285625 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       5 |
|             2 |      19 |
|             3 |      26 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  12.8359 |
| unique_mutation_site_count_top_candidates | 108      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.999131    |
| mean_seed_cosine_similarity   | 0.999465    |
| median_seed_cosine_similarity | 0.999478    |
| max_seed_cosine_similarity    | 0.999772    |
| min_seed_novelty_distance     | 0.000227988 |
| mean_seed_novelty_distance    | 0.000534766 |
| median_seed_novelty_distance  | 0.000521839 |
| max_seed_novelty_distance     | 0.000869334 |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000227928 |
| mean_novelty_distance                     | 0.000534844 |
| median_novelty_distance                   | 0.000522107 |
| max_novelty_distance                      | 0.000869691 |
| min_nearest_neighbor_cosine_similarity    | 0.99913     |
| mean_nearest_neighbor_cosine_similarity   | 0.999465    |
| median_nearest_neighbor_cosine_similarity | 0.999478    |
| max_nearest_neighbor_cosine_similarity    | 0.999772    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand200 | QFR57578.1                    | Klebsiella                    |                             0.999189 |        0.00081116  |
| round4_cand27  | QFR57578.1                    | Klebsiella                    |                             0.999311 |        0.00068891  |
| round4_cand492 | QFR57578.1                    | Klebsiella                    |                             0.999759 |        0.000240922 |
| round4_cand464 | QFR57578.1                    | Klebsiella                    |                             0.999569 |        0.000430584 |
| round4_cand319 | QFR57578.1                    | Klebsiella                    |                             0.999457 |        0.000543416 |
| round4_cand103 | QFR57578.1                    | Klebsiella                    |                             0.999312 |        0.000687718 |
| round4_cand8   | QFR57578.1                    | Klebsiella                    |                             0.999568 |        0.000432372 |
| round4_cand249 | QFR57578.1                    | Klebsiella                    |                             0.999374 |        0.000625551 |
| round4_cand467 | QFR57578.1                    | Klebsiella                    |                             0.999747 |        0.00025332  |
| round4_cand537 | QFR57578.1                    | Klebsiella                    |                             0.999416 |        0.000584126 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand200 | Klebsiella    | Pseudomonas   | Enterobacter |       0.285625 |          0.285625 |                 0.999189 |             0.0008111   | S126V;T482S;E540T |             3 |
| round4_cand27  | Klebsiella    | Pseudomonas   | Enterobacter |       0.283021 |          0.283021 |                 0.999311 |             0.00068903  | D52W;S134T;K455Q  |             3 |
| round4_cand492 | Klebsiella    | Pseudomonas   | Enterobacter |       0.279974 |          0.279974 |                 0.999759 |             0.000240862 | N20P;T63A;D541P   |             3 |
| round4_cand464 | Klebsiella    | Pseudomonas   | Enterobacter |       0.276576 |          0.276576 |                 0.99957  |             0.000430346 | G6R;G19A;A405G    |             3 |
| round4_cand319 | Klebsiella    | Pseudomonas   | Enterobacter |       0.275701 |          0.275701 |                 0.999456 |             0.000543535 | N108P;T406L;D449S |             3 |
| round4_cand103 | Klebsiella    | Pseudomonas   | Enterobacter |       0.275609 |          0.275609 |                 0.999312 |             0.000687599 | G101P;S262N       |             2 |
| round4_cand8   | Klebsiella    | Pseudomonas   | Enterobacter |       0.272227 |          0.272227 |                 0.999568 |             0.000432193 | S38G;D136P;S211V  |             3 |
| round4_cand249 | Klebsiella    | Pseudomonas   | Enterobacter |       0.267774 |          0.267774 |                 0.999375 |             0.000625372 | Y200N;E540V;M598A |             3 |
| round4_cand467 | Klebsiella    | Pseudomonas   | Enterobacter |       0.267766 |          0.267766 |                 0.999747 |             0.000253201 | F415A;I617A       |             2 |
| round4_cand537 | Klebsiella    | Pseudomonas   | Enterobacter |       0.266069 |          0.266069 |                 0.999416 |             0.000584006 | N108V;M254V;Y336A |             3 |
