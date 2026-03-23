# Design Run Summary

## Run metadata

- Seed protein ID: `QFR57578.1`
- Virus accession: `MN434096`
- Source host: `Klebsiella`
- Target host: `Enterobacter`
- ESM model: `facebook/esm2_t33_650M_UR50D`
- Rounds: `4`
- Candidates per round: `64`
- Min mutations: `1`
- Max mutations: `3`
- Top-K kept per round: `10`
- Proposal top-K: `8`
- Seed: `43`
- Total candidates evaluated: `1964`
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
| min      |       0.604799 |          0.574821 |
| mean     |       0.614387 |          0.584415 |
| median   |       0.612152 |          0.58218  |
| max      |       0.645079 |          0.615109 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       9 |
|             2 |      15 |
|             3 |      26 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  7.49796 |
| unique_mutation_site_count_top_candidates | 31       |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998566    |
| mean_seed_cosine_similarity   | 0.999064    |
| median_seed_cosine_similarity | 0.999071    |
| max_seed_cosine_similarity    | 0.99934     |
| min_seed_novelty_distance     | 0.000659764 |
| mean_seed_novelty_distance    | 0.000936322 |
| median_seed_novelty_distance  | 0.000929117 |
| max_seed_novelty_distance     | 0.00143355  |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand26  | Klebsiella    | Enterobacter  | Enterobacter |       0.645079 |          0.615109 |                 0.999007 |             0.000993133 | C132S;K523G;W527I |             3 |
| round4_cand92  | Klebsiella    | Enterobacter  | Enterobacter |       0.632285 |          0.602314 |                 0.999057 |             0.000942707 | K22P;C516T;K523T  |             3 |
| round4_cand3   | Klebsiella    | Enterobacter  | Enterobacter |       0.630249 |          0.60028  |                 0.998967 |             0.00103271  | H27A;K523N;C568S  |             3 |
| round4_cand74  | Klebsiella    | Enterobacter  | Enterobacter |       0.627344 |          0.597373 |                 0.99903  |             0.000969708 | C132V;K335A;C516D |             3 |
| round4_cand34  | Klebsiella    | Enterobacter  | Enterobacter |       0.626794 |          0.596837 |                 0.998566 |             0.00143355  | C162R;W527T;C587T |             3 |
| round4_cand54  | Klebsiella    | Enterobacter  | Enterobacter |       0.626678 |          0.596711 |                 0.998925 |             0.00107545  | E427V;K523T;K599R |             3 |
| round4_cand88  | Klebsiella    | Enterobacter  | Enterobacter |       0.625687 |          0.595713 |                 0.99911  |             0.000890076 | C132S;C516G;K599Q |             3 |
| round4_cand94  | Klebsiella    | Enterobacter  | Enterobacter |       0.62352  |          0.59354  |                 0.999308 |             0.000692129 | C516D             |             1 |
| round4_cand309 | Klebsiella    | Enterobacter  | Enterobacter |       0.620361 |          0.590397 |                 0.998809 |             0.00119054  | K22T;C162R;W505T  |             3 |
| round4_cand328 | Klebsiella    | Enterobacter  | Enterobacter |       0.620236 |          0.590264 |                 0.999066 |             0.000934482 | M226Y;W505D;K523N |             3 |
