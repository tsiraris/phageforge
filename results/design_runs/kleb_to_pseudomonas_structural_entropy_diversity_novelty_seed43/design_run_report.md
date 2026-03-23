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
- Seed: `43`
- Total candidates evaluated: `1962`
- Top candidates saved: `50`
- Diversity enabled: `True`
- Diversity min distance: `20`
- Position selection strategy: `entropy`
- Position pool size: `32`
- Proposal batch size: `16`
- Novelty penalty lambda: `0.03`
- Novelty penalty schedule: `constant`

## Target score summary

| metric   |   target_score |   selection_score |
|:---------|---------------:|------------------:|
| min      |       0.259944 |          0.229969 |
| mean     |       0.27338  |          0.243407 |
| median   |       0.268381 |          0.238405 |
| max      |       0.305155 |          0.275177 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       3 |
|             2 |      19 |
|             3 |      28 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 13.8016 |
| unique_mutation_site_count_top_candidates | 41      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998701    |
| mean_seed_cosine_similarity   | 0.999131    |
| median_seed_cosine_similarity | 0.999164    |
| max_seed_cosine_similarity    | 0.999408    |
| min_seed_novelty_distance     | 0.000591993 |
| mean_seed_novelty_distance    | 0.000869195 |
| median_seed_novelty_distance  | 0.000835955 |
| max_seed_novelty_distance     | 0.00129861  |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000592113 |
| mean_novelty_distance                     | 0.000869222 |
| median_novelty_distance                   | 0.000835806 |
| max_novelty_distance                      | 0.00129837  |
| min_nearest_neighbor_cosine_similarity    | 0.998702    |
| mean_nearest_neighbor_cosine_similarity   | 0.999131    |
| median_nearest_neighbor_cosine_similarity | 0.999164    |
| max_nearest_neighbor_cosine_similarity    | 0.999408    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand258 | QFR57578.1                    | Klebsiella                    |                             0.999273 |        0.000726521 |
| round4_cand620 | QFR57578.1                    | Klebsiella                    |                             0.998935 |        0.00106454  |
| round4_cand218 | QFR57578.1                    | Klebsiella                    |                             0.998743 |        0.00125706  |
| round4_cand181 | QFR57578.1                    | Klebsiella                    |                             0.998754 |        0.00124621  |
| round4_cand24  | QFR57578.1                    | Klebsiella                    |                             0.999013 |        0.000987172 |
| round4_cand35  | QFR57578.1                    | Klebsiella                    |                             0.999175 |        0.000824571 |
| round4_cand34  | QFR57578.1                    | Klebsiella                    |                             0.999125 |        0.000874758 |
| round4_cand157 | QFR57578.1                    | Klebsiella                    |                             0.999233 |        0.000767052 |
| round4_cand130 | QFR57578.1                    | Klebsiella                    |                             0.998702 |        0.00129837  |
| round4_cand349 | QFR57578.1                    | Klebsiella                    |                             0.999258 |        0.000741839 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand258 | Klebsiella    | Pseudomonas   | Enterobacter |       0.305155 |          0.275177 |                 0.999274 |             0.000726342 | W197A;K335P;M573A |             3 |
| round4_cand620 | Klebsiella    | Pseudomonas   | Enterobacter |       0.274655 |          0.244687 |                 0.998936 |             0.00106442  | K284S;K374L;M641A |             3 |
| round4_cand218 | Klebsiella    | Pseudomonas   | Enterobacter |       0.300256 |          0.270293 |                 0.998743 |             0.00125694  | W197A;H253I;K374V |             3 |
| round4_cand181 | Klebsiella    | Pseudomonas   | Enterobacter |       0.296757 |          0.266794 |                 0.998754 |             0.0012458   | M226N;K335P;C568A |             3 |
| round4_cand24  | Klebsiella    | Pseudomonas   | Enterobacter |       0.295739 |          0.265768 |                 0.999013 |             0.000987411 | K128A;C132R;K374L |             3 |
| round4_cand35  | Klebsiella    | Pseudomonas   | Enterobacter |       0.290578 |          0.260603 |                 0.999176 |             0.000824451 | K374A;H593M       |             2 |
| round4_cand34  | Klebsiella    | Pseudomonas   | Enterobacter |       0.290303 |          0.26033  |                 0.999126 |             0.000874341 | C516P;C535V;C568P |             3 |
| round4_cand157 | Klebsiella    | Pseudomonas   | Enterobacter |       0.28808  |          0.258103 |                 0.999233 |             0.000766814 | M163R;C574P       |             2 |
| round4_cand130 | Klebsiella    | Pseudomonas   | Enterobacter |       0.286451 |          0.25649  |                 0.998701 |             0.00129861  | C132T;K335A;M369A |             3 |
| round4_cand349 | Klebsiella    | Pseudomonas   | Enterobacter |       0.28445  |          0.254473 |                 0.999258 |             0.00074178  | W197G;C574A;H593A |             3 |
