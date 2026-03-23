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
- Seed: `43`
- Total candidates evaluated: `1980`
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
| min      |       0.130775 |          0.130775 |
| mean     |       0.13716  |          0.13716  |
| median   |       0.136041 |          0.136041 |
| max      |       0.159028 |          0.159028 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      38 |
| Klebsiella   |      12 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       6 |
|             2 |      19 |
|             3 |      25 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 16.9967 |
| unique_mutation_site_count_top_candidates | 97      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998785    |
| mean_seed_cosine_similarity   | 0.999367    |
| median_seed_cosine_similarity | 0.999445    |
| max_seed_cosine_similarity    | 0.999681    |
| min_seed_novelty_distance     | 0.000318706 |
| mean_seed_novelty_distance    | 0.000633256 |
| median_seed_novelty_distance  | 0.000554532 |
| max_seed_novelty_distance     | 0.0012148   |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand228 | Klebsiella    | Acinetobacter | Klebsiella   |       0.159028 |          0.159028 |                 0.999305 |             0.000694811 | P184V;R350T;Q586G |             3 |
| round4_cand22  | Klebsiella    | Acinetobacter | Klebsiella   |       0.151078 |          0.151078 |                 0.999192 |             0.000808179 | G57V;I354G;A605L  |             3 |
| round4_cand544 | Klebsiella    | Acinetobacter | Enterobacter |       0.141648 |          0.141648 |                 0.999633 |             0.000367105 | S213A;R289A;N424K |             3 |
| round4_cand157 | Klebsiella    | Acinetobacter | Enterobacter |       0.139264 |          0.139264 |                 0.999343 |             0.000656962 | E146T;G340Y;W632Q |             3 |
| round4_cand530 | Klebsiella    | Acinetobacter | Enterobacter |       0.154173 |          0.154173 |                 0.999262 |             0.000737607 | V393D;Q565G       |             2 |
| round4_cand492 | Klebsiella    | Acinetobacter | Enterobacter |       0.149347 |          0.149347 |                 0.999554 |             0.000445962 | R392N;W522N       |             2 |
| round4_cand245 | Klebsiella    | Acinetobacter | Enterobacter |       0.144848 |          0.144848 |                 0.999413 |             0.000587046 | R392G             |             1 |
| round4_cand323 | Klebsiella    | Acinetobacter | Enterobacter |       0.144528 |          0.144528 |                 0.999442 |             0.00055778  | V13L;R289A        |             2 |
| round4_cand50  | Klebsiella    | Acinetobacter | Enterobacter |       0.144383 |          0.144383 |                 0.999176 |             0.000824094 | S238D;N351V;W452A |             3 |
| round4_cand42  | Klebsiella    | Acinetobacter | Enterobacter |       0.143737 |          0.143737 |                 0.998867 |             0.00113291  | W32A;N249D;P292Y  |             3 |
