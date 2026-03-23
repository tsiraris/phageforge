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
- Seed: `42`
- Total candidates evaluated: `1982`
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
| min      |       0.522196 |          0.522196 |
| mean     |       0.533759 |          0.533759 |
| median   |       0.53115  |          0.53115  |
| max      |       0.566061 |          0.566061 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       9 |
|             2 |      11 |
|             3 |      30 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  15.5453 |
| unique_mutation_site_count_top_candidates | 103      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.999311    |
| mean_seed_cosine_similarity   | 0.999625    |
| median_seed_cosine_similarity | 0.999634    |
| max_seed_cosine_similarity    | 0.999793    |
| min_seed_novelty_distance     | 0.000207126 |
| mean_seed_novelty_distance    | 0.000374879 |
| median_seed_novelty_distance  | 0.000365704 |
| max_seed_novelty_distance     | 0.000689209 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand9   | Klebsiella    | Enterobacter  | Enterobacter |       0.566061 |          0.566061 |                 0.99963  |             0.000369728 | K105D;W527D;C574R |             3 |
| round4_cand28  | Klebsiella    | Enterobacter  | Enterobacter |       0.560412 |          0.560412 |                 0.999744 |             0.000255704 | W632Q             |             1 |
| round4_cand574 | Klebsiella    | Enterobacter  | Enterobacter |       0.558763 |          0.558763 |                 0.999472 |             0.000527501 | W329P;Y412R;V622R |             3 |
| round4_cand158 | Klebsiella    | Enterobacter  | Enterobacter |       0.556632 |          0.556632 |                 0.999493 |             0.000506639 | A343F;W632R       |             2 |
| round4_cand141 | Klebsiella    | Enterobacter  | Enterobacter |       0.551719 |          0.551719 |                 0.999533 |             0.000467181 | Y278N;G462V;Y536R |             3 |
| round4_cand21  | Klebsiella    | Enterobacter  | Enterobacter |       0.547763 |          0.547763 |                 0.999664 |             0.000335693 | G365N;N457R;W632S |             3 |
| round4_cand106 | Klebsiella    | Enterobacter  | Enterobacter |       0.544377 |          0.544377 |                 0.999511 |             0.000488758 | W620Q             |             1 |
| round4_cand311 | Klebsiella    | Enterobacter  | Enterobacter |       0.543213 |          0.543213 |                 0.999603 |             0.000396788 | W373N             |             1 |
| round4_cand159 | Klebsiella    | Enterobacter  | Enterobacter |       0.540219 |          0.540219 |                 0.999638 |             0.000361681 | V45G;Y107R        |             2 |
| round4_cand89  | Klebsiella    | Enterobacter  | Enterobacter |       0.538103 |          0.538103 |                 0.999441 |             0.000559032 | V43I;L100D;W545V  |             3 |
