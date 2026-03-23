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
- Seed: `42`
- Total candidates evaluated: `1950`
- Top candidates saved: `50`
- Diversity enabled: `True`
- Diversity min distance: `20`
- Position selection strategy: `entropy`
- Position pool size: `32`
- Proposal batch size: `16`
- Novelty penalty lambda: `0.0`
- Novelty penalty schedule: `constant`

## Target score summary

| metric   |   target_score |   selection_score |
|:---------|---------------:|------------------:|
| min      |      0.0998439 |         0.0998439 |
| mean     |      0.105893  |         0.105893  |
| median   |      0.103188  |         0.103188  |
| max      |      0.130846  |         0.130846  |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |      15 |
|             2 |      19 |
|             3 |      16 |

## Diversity summary

| metric                                    |   value |
|:------------------------------------------|--------:|
| avg_pairwise_distance_top_candidates      | 10.6988 |
| unique_mutation_site_count_top_candidates | 36      |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998037    |
| mean_seed_cosine_similarity   | 0.998999    |
| median_seed_cosine_similarity | 0.999111    |
| max_seed_cosine_similarity    | 0.999573    |
| min_seed_novelty_distance     | 0.00042665  |
| mean_seed_novelty_distance    | 0.00100092  |
| median_seed_novelty_distance  | 0.000888735 |
| max_seed_novelty_distance     | 0.00196284  |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand54  | Klebsiella    | Acinetobacter | Enterobacter |       0.130846 |          0.130846 |                 0.999124 |             0.000876427 | M65L;R350D        |             2 |
| round4_cand49  | Klebsiella    | Acinetobacter | Enterobacter |       0.127258 |          0.127258 |                 0.999337 |             0.00066328  | R350V;W527A       |             2 |
| round4_cand38  | Klebsiella    | Acinetobacter | Enterobacter |       0.126821 |          0.126821 |                 0.999114 |             0.00088644  | H271Y;R350T       |             2 |
| round4_cand29  | Klebsiella    | Acinetobacter | Enterobacter |       0.120373 |          0.120373 |                 0.999192 |             0.000808358 | R350A             |             1 |
| round4_cand19  | Klebsiella    | Acinetobacter | Enterobacter |       0.11602  |          0.11602  |                 0.999172 |             0.000828326 | R350P;M641P       |             2 |
| round4_cand8   | Klebsiella    | Acinetobacter | Enterobacter |       0.114328 |          0.114328 |                 0.999207 |             0.000792682 | R350S;F577T       |             2 |
| round4_cand58  | Klebsiella    | Acinetobacter | Enterobacter |       0.113467 |          0.113467 |                 0.998993 |             0.00100696  | C240A;H271A       |             2 |
| round4_cand241 | Klebsiella    | Acinetobacter | Enterobacter |       0.108515 |          0.108515 |                 0.998704 |             0.00129575  | H253L;E427V       |             2 |
| round4_cand9   | Klebsiella    | Acinetobacter | Enterobacter |       0.108373 |          0.108373 |                 0.998863 |             0.00113714  | P256V;F577P       |             2 |
| round4_cand140 | Klebsiella    | Acinetobacter | Enterobacter |       0.108337 |          0.108337 |                 0.998409 |             0.00159138  | P179N;M226D;M643N |             3 |
