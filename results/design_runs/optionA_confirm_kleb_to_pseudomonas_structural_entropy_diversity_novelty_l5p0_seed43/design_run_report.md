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
- Novelty penalty lambda: `5.0`
- Novelty penalty schedule: `constant`

## Target score summary

| metric   |   target_score |   selection_score |
|:---------|---------------:|------------------:|
| min      |       0.24993  |          -4.74555 |
| mean     |       0.265689 |          -4.72978 |
| median   |       0.262451 |          -4.73268 |
| max      |       0.295739 |          -4.69932 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       5 |
|             2 |      16 |
|             3 |      29 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 11.9037 |
| unique_mutation_site_count_top_candidates | 34      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998457    |
| mean_seed_cosine_similarity   | 0.999095    |
| median_seed_cosine_similarity | 0.999131    |
| max_seed_cosine_similarity    | 0.999432    |
| min_seed_novelty_distance     | 0.000567675 |
| mean_seed_novelty_distance    | 0.000905473 |
| median_seed_novelty_distance  | 0.000869483 |
| max_seed_novelty_distance     | 0.00154287  |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000567675 |
| mean_novelty_distance                     | 0.000905521 |
| median_novelty_distance                   | 0.000869721 |
| max_novelty_distance                      | 0.00154293  |
| min_nearest_neighbor_cosine_similarity    | 0.998457    |
| mean_nearest_neighbor_cosine_similarity   | 0.999094    |
| median_nearest_neighbor_cosine_similarity | 0.99913     |
| max_nearest_neighbor_cosine_similarity    | 0.999432    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round4_cand24  | QFR57578.1                    | Klebsiella                    |                             0.999013 |        0.000987172 |
| round4_cand335 | QFR57578.1                    | Klebsiella                    |                             0.999096 |        0.000903845 |
| round4_cand35  | QFR57578.1                    | Klebsiella                    |                             0.999175 |        0.000824571 |
| round4_cand34  | QFR57578.1                    | Klebsiella                    |                             0.999125 |        0.000874758 |
| round4_cand37  | QFR57578.1                    | Klebsiella                    |                             0.999112 |        0.00088799  |
| round4_cand22  | QFR57578.1                    | Klebsiella                    |                             0.999162 |        0.000837922 |
| round4_cand11  | QFR57578.1                    | Klebsiella                    |                             0.999363 |        0.000636518 |
| round4_cand287 | QFR57578.1                    | Klebsiella                    |                             0.998935 |        0.00106496  |
| round4_cand417 | QFR57578.1                    | Klebsiella                    |                             0.999304 |        0.000695825 |
| round4_cand183 | QFR57578.1                    | Klebsiella                    |                             0.998567 |        0.00143266  |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand24  | Klebsiella    | Pseudomonas   | Enterobacter |       0.295739 |          -4.69932 |                 0.999013 |             0.000987411 | K128A;C132R;K374L |             3 |
| round4_cand335 | Klebsiella    | Pseudomonas   | Enterobacter |       0.24993  |          -4.74555 |                 0.999096 |             0.000903726 | K155R;C162V;C516A |             3 |
| round4_cand35  | Klebsiella    | Pseudomonas   | Enterobacter |       0.290578 |          -4.7053  |                 0.999176 |             0.000824451 | K374A;H593M       |             2 |
| round4_cand34  | Klebsiella    | Pseudomonas   | Enterobacter |       0.290303 |          -4.70533 |                 0.999126 |             0.000874341 | C516P;C535V;C568P |             3 |
| round4_cand37  | Klebsiella    | Pseudomonas   | Enterobacter |       0.283076 |          -4.71248 |                 0.999112 |             0.000888288 | K325A;C516A;M641V |             3 |
| round4_cand22  | Klebsiella    | Pseudomonas   | Enterobacter |       0.281781 |          -4.71403 |                 0.999162 |             0.00083828  | C62T;K105R;E427R  |             3 |
| round4_cand11  | Klebsiella    | Pseudomonas   | Enterobacter |       0.28236  |          -4.71446 |                 0.999363 |             0.000636637 | C574A             |             1 |
| round4_cand287 | Klebsiella    | Pseudomonas   | Enterobacter |       0.278988 |          -4.71569 |                 0.998935 |             0.00106466  | K335A;H593L       |             2 |
| round4_cand417 | Klebsiella    | Pseudomonas   | Enterobacter |       0.276463 |          -4.72006 |                 0.999304 |             0.000695884 | K374G;C574P       |             2 |
| round4_cand183 | Klebsiella    | Pseudomonas   | Enterobacter |       0.272072 |          -4.72077 |                 0.998567 |             0.0014326   | C132P;K155R;H271I |             3 |
