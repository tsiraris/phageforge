# Design Run Summary

## Run metadata

- Seed protein ID: `QFR57578.1`
- Virus accession: `MN434096`
- Source host: `Klebsiella`
- Target host: `Acinetobacter`
- ESM model: `facebook/esm2_t33_650M_UR50D`
- Rounds: `4`
- Candidates per round: `64`
- Min mutations: `1`
- Max mutations: `3`
- Top-K kept per round: `10`
- Proposal top-K: `8`
- Seed: `44`
- Total candidates evaluated: `1984`
- Top candidates saved: `50`
- Diversity enabled: `True`
- Diversity min distance: `20`
- Position selection strategy: `random`
- Position pool size: `32`
- Proposal batch size: `16`
- Novelty penalty lambda: `0.0`
- Novelty penalty schedule: `constant`

## Target score summary

| metric   |   target_score |   selection_score |
|:---------|---------------:|------------------:|
| min      |       0.135731 |          0.135731 |
| mean     |       0.14532  |          0.14532  |
| median   |       0.140882 |          0.140882 |
| max      |       0.219928 |          0.219928 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Klebsiella   |      34 |
| Enterobacter |      16 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       6 |
|             2 |      17 |
|             3 |      27 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  11.3633 |
| unique_mutation_site_count_top_candidates | 100      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.997574    |
| mean_seed_cosine_similarity   | 0.999251    |
| median_seed_cosine_similarity | 0.999396    |
| max_seed_cosine_similarity    | 0.999634    |
| min_seed_novelty_distance     | 0.000365734 |
| mean_seed_novelty_distance    | 0.000749059 |
| median_seed_novelty_distance  | 0.000604153 |
| max_seed_novelty_distance     | 0.00242609  |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand475 | Klebsiella    | Acinetobacter | Klebsiella   |       0.219928 |          0.219928 |                 0.997783 |             0.00221699  | T233N;Q417I;G495L |             3 |
| round4_cand523 | Klebsiella    | Acinetobacter | Enterobacter |       0.14234  |          0.14234  |                 0.999401 |             0.000598729 | W59A;L520K        |             2 |
| round4_cand9   | Klebsiella    | Acinetobacter | Klebsiella   |       0.186496 |          0.186496 |                 0.998011 |             0.00198877  | L476I;G495V;M567S |             3 |
| round4_cand274 | Klebsiella    | Acinetobacter | Klebsiella   |       0.165414 |          0.165414 |                 0.999343 |             0.000656843 | E86D;R492Y;G514K  |             3 |
| round4_cand218 | Klebsiella    | Acinetobacter | Klebsiella   |       0.161178 |          0.161178 |                 0.998944 |             0.00105584  | D55N;R242D;N416L  |             3 |
| round4_cand7   | Klebsiella    | Acinetobacter | Klebsiella   |       0.160197 |          0.160197 |                 0.999119 |             0.000880659 | I69A;Q417V;N594Y  |             3 |
| round4_cand124 | Klebsiella    | Acinetobacter | Klebsiella   |       0.152265 |          0.152265 |                 0.9994   |             0.000599682 | R232D;M254A;V560G |             3 |
| round4_cand54  | Klebsiella    | Acinetobacter | Klebsiella   |       0.152163 |          0.152163 |                 0.99931  |             0.000689507 | P358E;N594I;W632K |             3 |
| round4_cand17  | Klebsiella    | Acinetobacter | Klebsiella   |       0.151465 |          0.151465 |                 0.999391 |             0.000608623 | G90A;C240E;R651V  |             3 |
| round4_cand30  | Klebsiella    | Acinetobacter | Klebsiella   |       0.150422 |          0.150422 |                 0.998975 |             0.00102532  | F401N;R553Q       |             2 |
