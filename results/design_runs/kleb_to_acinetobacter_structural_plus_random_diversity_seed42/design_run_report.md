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
- Total candidates evaluated: `1982`
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
| min      |       0.142265 |          0.142265 |
| mean     |       0.197128 |          0.197128 |
| median   |       0.198568 |          0.198568 |
| max      |       0.24584  |          0.24584  |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      49 |
| Klebsiella   |       1 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |      18 |
|             2 |      17 |
|             3 |      15 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  5.38694 |
| unique_mutation_site_count_top_candidates | 89       |

## Seed similarity summary

| metric                        |       value |
|:------------------------------|------------:|
| min_seed_cosine_similarity    | 0.99369     |
| mean_seed_cosine_similarity   | 0.994515    |
| median_seed_cosine_similarity | 0.994269    |
| max_seed_cosine_similarity    | 0.999488    |
| min_seed_novelty_distance     | 0.000512004 |
| mean_seed_novelty_distance    | 0.00548519  |
| median_seed_novelty_distance  | 0.00573102  |
| max_seed_novelty_distance     | 0.00631022  |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score |   selection_score |   seed_cosine_similarity |   seed_novelty_distance | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|------------------:|-------------------------:|------------------------:|:------------------|--------------:|
| round4_cand45  | Klebsiella    | Acinetobacter | Enterobacter |       0.24584  |          0.24584  |                 0.99411  |             0.00589001  | V81P;R477A;E547R  |             3 |
| round4_cand395 | Klebsiella    | Acinetobacter | Klebsiella   |       0.148443 |          0.148443 |                 0.999488 |             0.000512004 | R289N;Y336D;N608S |             3 |
| round4_cand335 | Klebsiella    | Acinetobacter | Enterobacter |       0.142265 |          0.142265 |                 0.999052 |             0.000947773 | L241A;Y251V;M254N |             3 |
| round4_cand324 | Klebsiella    | Acinetobacter | Enterobacter |       0.22652  |          0.22652  |                 0.996228 |             0.00377196  | F401T             |             1 |
| round4_cand50  | Klebsiella    | Acinetobacter | Enterobacter |       0.223769 |          0.223769 |                 0.99369  |             0.00631022  | W632N             |             1 |
| round4_cand48  | Klebsiella    | Acinetobacter | Enterobacter |       0.222588 |          0.222588 |                 0.994281 |             0.00571883  | R477S;T482N       |             2 |
| round4_cand37  | Klebsiella    | Acinetobacter | Enterobacter |       0.222301 |          0.222301 |                 0.993764 |             0.00623566  | D52A;S134V;K455A  |             3 |
| round4_cand32  | Klebsiella    | Acinetobacter | Enterobacter |       0.22115  |          0.22115  |                 0.993909 |             0.00609094  | G556A             |             1 |
| round4_cand21  | Klebsiella    | Acinetobacter | Enterobacter |       0.215573 |          0.215573 |                 0.994311 |             0.00568873  | L46S;S609T        |             2 |
| round4_cand63  | Klebsiella    | Acinetobacter | Enterobacter |       0.215144 |          0.215144 |                 0.994281 |             0.00571895  | R562S             |             1 |
