# PhageForge corrected validity closeout report

## Executive summary

- Total candidate rows analyzed: **4800**
- Total unique sequences after deduplication: **1587**
- Corrected validated rows: **3786**
- Corrected validated unique sequences: **978**
- Unique sequences in final shortlist: **20**
- Strongest target direction after correction: **Pseudomonas**

## What was corrected

1. Family-retention matching was broadened to correctly recognize both hyphenated and non-hyphenated RBP annotations, including receptor-binding/receptor binding, tail-fiber/tail fiber, and tail-spike/tail spike.
2. Mutation counts were recomputed directly from candidate sequences and seed sequences, rather than relying on the original logged mutation field alone.
3. Final ranking was performed on unique sequences per target so repeated copies of the same sequence across runs do not dominate the shortlist.

## Patch metadata

```json
{
  "source_master_csv": "results/analysis/validity/master_candidates_validated.csv",
  "source_validated_pool_csv": "results/analysis/validity/validated_pool_ranked.csv",
  "max_actual_mutations": 12,
  "strict_cosine_floor": 0.995,
  "strict_top5_floor": 0.992,
  "pll_delta_floor": 0.0,
  "preferred_targets": [
    "Enterobacter",
    "Acinetobacter"
  ],
  "mutation_disagreement_summary": {
    "n_rows": 4800,
    "n_disagreements": 4800,
    "fraction_disagreement": 1.0,
    "mean_reported_n_mutations": 2.4175,
    "mean_actual_n_mutations": 10.981666666666667,
    "mean_original_seq_dist_from_seed": 10.981666666666667,
    "mean_patched_seq_dist_from_seed": 10.981666666666667
  },
  "n_rows_full": 4800,
  "n_rows_unique": 1587,
  "n_rows_unique_preferred": 619,
  "n_validated_patched_full": 3786,
  "n_validated_patched_unique": 978
}
```

## Summary by target

| target_host   |   n_rows |   n_unique_sequences |   n_validated_original |   n_validated_patched |   mean_target_score |   best_target_score |   mean_actual_seq_dist |   median_actual_seq_dist |   mean_pll_delta |   best_corrected_rank_score |
|:--------------|---------:|---------------------:|-----------------------:|----------------------:|--------------------:|--------------------:|-----------------------:|-------------------------:|-----------------:|----------------------------:|
| Pseudomonas   |     2100 |                  968 |                      0 |                  1472 |            0.289858 |            0.492084 |                12.659  |                       11 |         0.993909 |                    1.18253  |
| Enterobacter  |     1500 |                  317 |                    242 |                  1344 |            0.583972 |            0.645079 |                10.0387 |                       10 |         1.00206  |                    0.879159 |
| Acinetobacter |     1200 |                  302 |                     60 |                   970 |            0.132538 |            0.24584  |                 9.225  |                        9 |         0.768132 |                    0.904089 |

## Final shortlist

| target_host   | validated_pass_patched   |   target_score |   corrected_rank_score |   strict_nn_cosine |   pll_delta_candidate_minus_seed |   actual_n_mutations | actual_mutation_list                                              | strict_nn_product             | method                                                           | run_name                                                         |
|:--------------|:-------------------------|---------------:|-----------------------:|-------------------:|---------------------------------:|---------------------:|:------------------------------------------------------------------|:------------------------------|:-----------------------------------------------------------------|:-----------------------------------------------------------------|
| Enterobacter  | True                     |       0.603959 |               0.879159 |           0.999468 |                         1.75241  |                    8 | P179N;W197R;K260P;W373T;W505R;W522P;K523D;W527I                   | receptor-binding protein      | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 |
| Enterobacter  | True                     |       0.61574  |               0.876392 |           0.999366 |                         1.93237  |                    9 | W197R;C240N;K260P;W373T;W505N;W522P;K523D;W527I;H593L             | receptor binding tail protein | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 |
| Enterobacter  | True                     |       0.601474 |               0.875199 |           0.999448 |                         1.93682  |                    8 | W197R;K260P;W329N;W373T;W522P;K523D;W527I;C568T                   | receptor binding tail protein | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 |
| Enterobacter  | True                     |       0.596277 |               0.779577 |           0.999436 |                         1.94463  |                    8 | W197R;K260P;W373T;W522P;K523D;W527I;K561S;C574R                   | receptor binding tail protein | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 |
| Enterobacter  | True                     |       0.597279 |               0.771271 |           0.999476 |                         1.22073  |                    7 | W197R;K260P;W373T;C516D;W522P;K523D;W527I                         | receptor binding tail protein | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 |
| Enterobacter  | True                     |       0.626677 |               0.766905 |           0.999212 |                         1.53883  |                    9 | K22V;W197R;K260P;W329N;K335A;W373T;W522P;K523D;W527I              | receptor binding tail protein | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 |
| Enterobacter  | True                     |       0.60811  |               0.745481 |           0.999399 |                         1.49269  |                    9 | C162N;W197R;K260P;W373T;W505S;C516T;W522P;K523D;W527I             | receptor binding tail protein | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 |
| Enterobacter  | True                     |       0.605467 |               0.728431 |           0.99934  |                         2.60774  |                   10 | K128V;W197D;K260A;W329N;W373T;K374F;W505S;C516N;W522T;C574R       | receptor-binding protein      | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed43 | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed43 |
| Enterobacter  | True                     |       0.603442 |               0.709978 |           0.999402 |                         1.69595  |                    9 | P179G;W197R;K260P;W329P;K335G;W373T;W522P;K523D;W527I             | receptor-binding protein      | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed42 |
| Enterobacter  | True                     |       0.62352  |               0.6839   |           0.999308 |                         0.994361 |                    9 | K260A;W329N;W373T;K374F;W505R;C516D;W522T;W527T;C574R             | receptor-binding protein      | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed43 | kleb_to_enterobacter_structural_entropy_diversity_novelty_seed43 |
| Acinetobacter | True                     |       0.149347 |               0.606315 |           0.999554 |                         1.33529  |                    9 | Y96G;E171N;R232D;G267V;R289A;R392N;W522N;R562A;S609L              | receptor-binding protein      | kleb_to_acinetobacter_structural_plus_random_diversity_seed43    | kleb_to_acinetobacter_structural_plus_random_diversity_seed43    |
| Acinetobacter | True                     |       0.152265 |               0.538463 |           0.9994   |                         0.87468  |                   10 | P177T;R232D;M254A;R322G;G338D;G372P;R400G;Q560G;Y615N;V633C       | receptor-binding protein      | kleb_to_acinetobacter_structural_plus_random_diversity_seed44    | kleb_to_acinetobacter_structural_plus_random_diversity_seed44    |
| Acinetobacter | True                     |       0.140755 |               0.512615 |           0.999549 |                         1.10639  |                   10 | P177T;R322G;G338D;G372P;R400G;Q560V;M598L;Y615N;W620S;V633C       | receptor-binding protein      | kleb_to_acinetobacter_structural_plus_random_diversity_seed44    | kleb_to_acinetobacter_structural_plus_random_diversity_seed44    |
| Acinetobacter | True                     |       0.14721  |               0.499911 |           0.99962  |                         0.201523 |                    9 | D58K;R157L;R201T;G338D;G372P;R400G;Q560V;Y615N;V633C              | receptor-binding protein      | kleb_to_acinetobacter_structural_plus_random_diversity_seed44    | kleb_to_acinetobacter_structural_plus_random_diversity_seed44    |
| Acinetobacter | True                     |       0.154173 |               0.472561 |           0.999263 |                         1.34662  |                    9 | R232D;T233V;V323D;G386T;D390L;V393D;D525N;Q565G;P591N             | receptor-binding protein      | kleb_to_acinetobacter_structural_plus_random_diversity_seed43    | kleb_to_acinetobacter_structural_plus_random_diversity_seed43    |
| Acinetobacter | True                     |       0.159028 |               0.469348 |           0.999305 |                         0.594917 |                   11 | P184V;G235P;R248S;Y336D;N337A;R350T;G386V;G408T;L419A;R553A;Q586G | receptor-binding protein      | kleb_to_acinetobacter_structural_plus_random_diversity_seed43    | kleb_to_acinetobacter_structural_plus_random_diversity_seed43    |
| Acinetobacter | True                     |       0.144848 |               0.417966 |           0.999413 |                         0.501558 |                    9 | G235P;R248S;Y336D;N337A;G386V;R392G;G408T;L419A;R553A             | receptor-binding protein      | kleb_to_acinetobacter_structural_plus_random_diversity_seed43    | kleb_to_acinetobacter_structural_plus_random_diversity_seed43    |
| Acinetobacter | True                     |       0.140289 |               0.395776 |           0.999584 |                         0.63743  |                   10 | P177T;R322G;G338D;G372P;R400G;E540N;Q560V;Y615N;V633C;R644S       | receptor-binding protein      | kleb_to_acinetobacter_structural_plus_random_diversity_seed44    | kleb_to_acinetobacter_structural_plus_random_diversity_seed44    |
| Acinetobacter | True                     |       0.138788 |               0.393693 |           0.999611 |                         0.159508 |                    9 | R232S;G338D;G372P;R400G;W469V;R553N;Q560V;Y615N;V633C             | receptor-binding protein      | kleb_to_acinetobacter_structural_plus_random_diversity_seed44    | kleb_to_acinetobacter_structural_plus_random_diversity_seed44    |
| Acinetobacter | True                     |       0.132469 |               0.373019 |           0.999403 |                         1.20484  |                    8 | P179V;R217V;R248A;M254N;P268G;V323D;R328G;Q464N                   | receptor-binding protein      | kleb_to_acinetobacter_structural_plus_random_diversity_seed43    | kleb_to_acinetobacter_structural_plus_random_diversity_seed43    |

## Top runs after correction

| run_name                                                               | target_host   | method_family     |   n_validated_patched |   best_target_score |   best_corrected_rank_score |   best_strict_nn_cosine |   best_pll_delta |   mean_actual_seq_dist |
|:-----------------------------------------------------------------------|:--------------|:------------------|----------------------:|--------------------:|----------------------------:|------------------------:|-----------------:|-----------------------:|
| kleb_to_acinetobacter_structural_plus_random_diversity_seed44          | Acinetobacter | diversity         |                    28 |            0.219928 |                    0.904089 |                0.999635 |         1.54663  |                  10.8  |
| kleb_to_acinetobacter_structural_plus_random_no_diversity_seed44       | Acinetobacter | diversity         |                    28 |            0.219928 |                    0.904089 |                0.999635 |         1.54663  |                  10.8  |
| kleb_to_acinetobacter_structural_random_diversity_seed44               | Acinetobacter | diversity         |                    28 |            0.219928 |                    0.904089 |                0.999635 |         1.54663  |                  10.8  |
| kleb_to_acinetobacter_structural_random_no_diversity_seed44            | Acinetobacter | diversity         |                    28 |            0.219928 |                    0.904089 |                0.999635 |         1.54663  |                  10.8  |
| kleb_to_acinetobacter_structural_plus_random_diversity_seed43          | Acinetobacter | diversity         |                    37 |            0.159028 |                    0.606315 |                0.999681 |         1.34662  |                   9.86 |
| kleb_to_acinetobacter_structural_plus_random_no_diversity_seed43       | Acinetobacter | diversity         |                    37 |            0.159028 |                    0.606315 |                0.999681 |         1.34662  |                   9.86 |
| kleb_to_acinetobacter_structural_random_diversity_seed43               | Acinetobacter | diversity         |                    37 |            0.159028 |                    0.606315 |                0.999681 |         1.34662  |                   9.86 |
| kleb_to_acinetobacter_structural_random_no_diversity_seed43            | Acinetobacter | diversity         |                    37 |            0.159028 |                    0.606315 |                0.999681 |         1.34662  |                   9.86 |
| kleb_to_acinetobacter_structural_plus_random_diversity_seed42          | Acinetobacter | diversity         |                    28 |            0.24584  |                    0.450824 |                0.999488 |         0.827961 |                   8.02 |
| kleb_to_acinetobacter_structural_plus_random_no_diversity_seed42       | Acinetobacter | diversity         |                    27 |            0.24584  |                    0.450824 |                0.997727 |         0.827961 |                   7.88 |
| kleb_to_acinetobacter_structural_random_diversity_seed42               | Acinetobacter | diversity         |                    28 |            0.24584  |                    0.450824 |                0.999488 |         0.827961 |                   8.02 |
| kleb_to_acinetobacter_structural_random_no_diversity_seed42            | Acinetobacter | diversity         |                    27 |            0.24584  |                    0.450824 |                0.997727 |         0.827961 |                   7.88 |
| kleb_to_acinetobacter_structural_entropy_diversity_novelty_seed44      | Acinetobacter | diversity+novelty |                    50 |            0.125518 |                    0.274165 |                0.999581 |         2.37211  |                   8.22 |
| kleb_to_acinetobacter_structural_entropy_diversity_seed44              | Acinetobacter | diversity         |                    50 |            0.125518 |                    0.274165 |                0.999581 |         2.37211  |                   8.22 |
| kleb_to_acinetobacter_structural_plus_entropy_diversity_novelty_seed44 | Acinetobacter | diversity+novelty |                    50 |            0.125518 |                    0.274165 |                0.999581 |         2.37211  |                   8.22 |
| kleb_to_acinetobacter_structural_plus_entropy_diversity_seed44         | Acinetobacter | diversity         |                    50 |            0.125518 |                    0.274165 |                0.999581 |         2.37211  |                   8.22 |
| kleb_to_acinetobacter_structural_entropy_diversity_novelty_seed42      | Acinetobacter | diversity+novelty |                    50 |            0.130846 |                    0.241036 |                0.999573 |         2.45068  |                   8.7  |
| kleb_to_acinetobacter_structural_entropy_diversity_seed42              | Acinetobacter | diversity         |                    50 |            0.130846 |                    0.241036 |                0.999573 |         2.45068  |                   8.7  |
| kleb_to_acinetobacter_structural_plus_entropy_diversity_novelty_seed42 | Acinetobacter | diversity+novelty |                    50 |            0.130846 |                    0.241036 |                0.999573 |         2.45068  |                   8.7  |
| kleb_to_acinetobacter_structural_plus_entropy_diversity_seed42         | Acinetobacter | diversity         |                    50 |            0.130846 |                    0.241036 |                0.999573 |         2.45068  |                   8.7  |

## Interpretation

The corrected analysis should be used instead of the initial pass because the initial family-retention logic undercounted legitimate RBP-like annotations and the original mutation counts did not reliably match the actual sequence differences. The corrected closeout should emphasize unique-sequence evidence, not duplicated copies of the same sequence across multiple runs.

For the final project story, focus on the strongest corrected target direction first, keep Acinetobacter as secondary if it remains plausible after correction, and remove noisy targets from the headline if they do not hold up under the corrected validity layer.
