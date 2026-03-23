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
| min      |       0.533782 |          0.533782 |
| mean     |       0.541086 |          0.541086 |
| median   |       0.538504 |          0.538504 |
| max      |       0.563627 |          0.563627 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       8 |
|             2 |      18 |
|             3 |      24 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 10.9306 |
| unique_mutation_site_count_top_candidates | 97      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.999275    |
| mean_seed_cosine_similarity   | 0.999668    |
| median_seed_cosine_similarity | 0.999688    |
| max_seed_cosine_similarity    | 0.999776    |
| min_seed_novelty_distance     | 0.000223577 |
| mean_seed_novelty_distance    | 0.000332361 |
| median_seed_novelty_distance  | 0.000311732 |
| max_seed_novelty_distance     | 0.000724912 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand23  | Klebsiella    | Enterobacter  | Enterobacter |       0.563627 |          0.563627 |                 0.999572 |             0.000427783 | E86R;W310S;E630R  |             3 |
| round4_cand428 | Klebsiella    | Enterobacter  | Enterobacter |       0.556532 |          0.556532 |                 0.999722 |             0.000278413 | T482G;W527S       |             2 |
| round4_cand526 | Klebsiella    | Enterobacter  | Enterobacter |       0.555727 |          0.555727 |                 0.999638 |             0.000362098 | W59V;L520N        |             2 |
| round4_cand42  | Klebsiella    | Enterobacter  | Enterobacter |       0.553985 |          0.553985 |                 0.999598 |             0.000402153 | D78G;G202N;W505G  |             3 |
| round4_cand68  | Klebsiella    | Enterobacter  | Enterobacter |       0.552265 |          0.552265 |                 0.999715 |             0.000285447 | H253N;G429R       |             2 |
| round4_cand114 | Klebsiella    | Enterobacter  | Enterobacter |       0.550555 |          0.550555 |                 0.999718 |             0.000281513 | M598L;W620Q       |             2 |
| round4_cand92  | Klebsiella    | Enterobacter  | Enterobacter |       0.549915 |          0.549915 |                 0.999703 |             0.000297427 | W59R;N345E;N377S  |             3 |
| round4_cand383 | Klebsiella    | Enterobacter  | Enterobacter |       0.549526 |          0.549526 |                 0.999677 |             0.000322878 | A109L;S117T;K155R |             3 |
| round4_cand278 | Klebsiella    | Enterobacter  | Enterobacter |       0.549087 |          0.549087 |                 0.999594 |             0.00040561  | W59A;I407F        |             2 |
| round4_cand632 | Klebsiella    | Enterobacter  | Enterobacter |       0.546156 |          0.546156 |                 0.999605 |             0.000394583 | P66R;I554N        |             2 |
