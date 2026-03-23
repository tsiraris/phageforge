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
- Total candidates evaluated: `1981`
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
| min      |       0.558301 |          0.558301 |
| mean     |       0.567089 |          0.567089 |
| median   |       0.564796 |          0.564796 |
| max      |       0.592952 |          0.592952 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |      19 |
|             2 |      12 |
|             3 |      19 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  5.81469 |
| unique_mutation_site_count_top_candidates | 92       |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.999298    |
| mean_seed_cosine_similarity   | 0.999497    |
| median_seed_cosine_similarity | 0.999504    |
| max_seed_cosine_similarity    | 0.999641    |
| min_seed_novelty_distance     | 0.00035876  |
| mean_seed_novelty_distance    | 0.000502689 |
| median_seed_novelty_distance  | 0.000495613 |
| max_seed_novelty_distance     | 0.000701666 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand41  | Klebsiella    | Enterobacter  | Enterobacter |       0.592952 |          0.592952 |                 0.999365 |             0.000634789 | Y283S;Y412T;G483N |             3 |
| round4_cand10  | Klebsiella    | Enterobacter  | Enterobacter |       0.586316 |          0.586316 |                 0.99941  |             0.000589669 | S112R;G495A;A502N |             3 |
| round4_cand13  | Klebsiella    | Enterobacter  | Enterobacter |       0.584729 |          0.584729 |                 0.999349 |             0.000651419 | E460G;I533R;V650Y |             3 |
| round4_cand21  | Klebsiella    | Enterobacter  | Enterobacter |       0.584385 |          0.584385 |                 0.999298 |             0.000701666 | W197T;A450R;I554A |             3 |
| round4_cand30  | Klebsiella    | Enterobacter  | Enterobacter |       0.583828 |          0.583828 |                 0.999485 |             0.000515342 | W172P             |             1 |
| round4_cand36  | Klebsiella    | Enterobacter  | Enterobacter |       0.581284 |          0.581284 |                 0.9994   |             0.00060004  | W59G;K260S;D525S  |             3 |
| round4_cand68  | Klebsiella    | Enterobacter  | Enterobacter |       0.579334 |          0.579334 |                 0.999486 |             0.000513852 | A195D;L491M;W505R |             3 |
| round4_cand9   | Klebsiella    | Enterobacter  | Enterobacter |       0.575224 |          0.575224 |                 0.99943  |             0.000569582 | V327D;N524P;G534T |             3 |
| round4_cand37  | Klebsiella    | Enterobacter  | Enterobacter |       0.572928 |          0.572928 |                 0.999448 |             0.000552356 | E276D;L569A;S609A |             3 |
| round4_cand5   | Klebsiella    | Enterobacter  | Enterobacter |       0.572563 |          0.572563 |                 0.999486 |             0.000513554 | G297F;R421T;A531T |             3 |
