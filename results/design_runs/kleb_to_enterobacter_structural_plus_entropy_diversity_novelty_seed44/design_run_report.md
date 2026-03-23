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
- Seed: `44`
- Total candidates evaluated: `1966`
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
| min      |       0.600079 |          0.570107 |
| mean     |       0.609565 |          0.579592 |
| median   |       0.607037 |          0.577064 |
| max      |       0.644352 |          0.614387 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       9 |
|             2 |      16 |
|             3 |      25 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 12.6416 |
| unique_mutation_site_count_top_candidates | 35      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998457    |
| mean_seed_cosine_similarity   | 0.999094    |
| median_seed_cosine_similarity | 0.999096    |
| max_seed_cosine_similarity    | 0.999439    |
| min_seed_novelty_distance     | 0.000561476 |
| mean_seed_novelty_distance    | 0.000905672 |
| median_seed_novelty_distance  | 0.000904024 |
| max_seed_novelty_distance     | 0.00154322  |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand58  | Klebsiella    | Enterobacter  | Enterobacter |       0.644352 |          0.614387 |                 0.99881  |             0.00118965  | K128A;W329N;K599R |             3 |
| round4_cand60  | Klebsiella    | Enterobacter  | Enterobacter |       0.638388 |          0.608423 |                 0.998808 |             0.00119215  | K325T;W329D;C568R |             3 |
| round4_cand631 | Klebsiella    | Enterobacter  | Enterobacter |       0.631685 |          0.601714 |                 0.999027 |             0.000973463 | W505G;W527R       |             2 |
| round4_cand0   | Klebsiella    | Enterobacter  | Enterobacter |       0.626493 |          0.596521 |                 0.99906  |             0.000939786 | H27A;M226N;W329N  |             3 |
| round4_cand243 | Klebsiella    | Enterobacter  | Enterobacter |       0.621526 |          0.591554 |                 0.999069 |             0.000931084 | K22A;K523D        |             2 |
| round4_cand45  | Klebsiella    | Enterobacter  | Enterobacter |       0.617689 |          0.587714 |                 0.99914  |             0.000859618 | K523T;W620T;M637P |             3 |
| round4_cand375 | Klebsiella    | Enterobacter  | Enterobacter |       0.617313 |          0.587335 |                 0.999292 |             0.000707746 | M334T;W373N       |             2 |
| round4_cand4   | Klebsiella    | Enterobacter  | Enterobacter |       0.617284 |          0.587313 |                 0.999035 |             0.000965297 | K523R;K599V       |             2 |
| round4_cand12  | Klebsiella    | Enterobacter  | Enterobacter |       0.615324 |          0.585353 |                 0.999023 |             0.00097692  | K335P;K523N;E540L |             3 |
| round4_cand111 | Klebsiella    | Enterobacter  | Enterobacter |       0.61502  |          0.585041 |                 0.99928  |             0.000720322 | K325G;W620R       |             2 |
