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
- Total candidates evaluated: `1965`
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
| min      |       0.10314  |         0.0731598 |
| mean     |       0.107609 |         0.0776353 |
| median   |       0.106491 |         0.0765142 |
| max      |       0.125518 |         0.0955404 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |      15 |
|             2 |      14 |
|             3 |      21 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  9.54694 |
| unique_mutation_site_count_top_candidates | 37       |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.998402    |
| mean_seed_cosine_similarity   | 0.999123    |
| median_seed_cosine_similarity | 0.999173    |
| max_seed_cosine_similarity    | 0.999581    |
| min_seed_novelty_distance     | 0.000418544 |
| mean_seed_novelty_distance    | 0.000877436 |
| median_seed_novelty_distance  | 0.000826508 |
| max_seed_novelty_distance     | 0.00159794  |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand97  | Klebsiella    | Acinetobacter | Enterobacter |       0.125518 |         0.0955404 |                 0.999267 |             0.00073278  | P179V;E427N;C568P |             3 |
| round4_cand108 | Klebsiella    | Acinetobacter | Enterobacter |       0.122024 |         0.0920574 |                 0.998884 |             0.00111598  | K128G;P179L;P256V |             3 |
| round4_cand67  | Klebsiella    | Acinetobacter | Enterobacter |       0.11562  |         0.0856465 |                 0.9991   |             0.000899673 | P179L;E427S;F577Q |             3 |
| round4_cand117 | Klebsiella    | Acinetobacter | Enterobacter |       0.113258 |         0.0832877 |                 0.998999 |             0.00100064  | P179V;V441T;K561A |             3 |
| round4_cand94  | Klebsiella    | Acinetobacter | Enterobacter |       0.113044 |         0.083076  |                 0.99894  |             0.00105983  | H271F             |             1 |
| round4_cand101 | Klebsiella    | Acinetobacter | Enterobacter |       0.111682 |         0.0817054 |                 0.999229 |             0.000771344 | H27V;W522G        |             2 |
| round4_cand95  | Klebsiella    | Acinetobacter | Enterobacter |       0.111518 |         0.0815417 |                 0.999212 |             0.000787914 | P179S             |             1 |
| round4_cand58  | Klebsiella    | Acinetobacter | Enterobacter |       0.111063 |         0.0810917 |                 0.999039 |             0.000960886 | W505N;M637S;M643T |             3 |
| round4_cand271 | Klebsiella    | Acinetobacter | Enterobacter |       0.110897 |         0.0809448 |                 0.998402 |             0.00159794  | M226D;H253I;M567T |             3 |
| round4_cand10  | Klebsiella    | Acinetobacter | Enterobacter |       0.110495 |         0.0805264 |                 0.998959 |             0.00104088  | W197E;P268V;M643G |             3 |
