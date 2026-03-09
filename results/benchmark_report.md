# PhageForge Benchmark Report

This report compares:
- broad vs strict dataset definitions
- MLP vs linear probe baselines
- 3 random seeds (42, 43, 44)

## Aggregate benchmark table

| dataset_name   | model_name   | accuracy      | macro_f1      | weighted_f1   |   n_total |   n_train |    n_val |   n_test |
|:---------------|:-------------|:--------------|:--------------|:--------------|----------:|----------:|---------:|---------:|
| broad          | linear_probe | 0.667 ± 0.050 | 0.673 ± 0.011 | 0.671 ± 0.051 |      3396 |   2674.67 | 386      | 335.333  |
| broad          | mlp          | 0.356 ± 0.018 | 0.293 ± 0.033 | 0.321 ± 0.052 |      3396 |   2674.67 | 386      | 335.333  |
| strict         | linear_probe | 0.743 ± 0.087 | 0.596 ± 0.101 | 0.689 ± 0.083 |       154 |    125    |  14.6667 |  14.3333 |
| strict         | mlp          | 0.651 ± 0.014 | 0.378 ± 0.052 | 0.563 ± 0.019 |       154 |    125    |  14.6667 |  14.3333 |

## Best run per dataset/model

| dataset_name   | model_name   |   seed |   accuracy |   macro_f1 |   weighted_f1 |
|:---------------|:-------------|-------:|-----------:|-----------:|--------------:|
| broad          | linear_probe |     43 |   0.634921 |   0.683369 |      0.641401 |
| broad          | mlp          |     44 |   0.364261 |   0.329396 |      0.373693 |
| strict         | linear_probe |     42 |   0.8      |   0.654762 |      0.744444 |
| strict         | mlp          |     43 |   0.642857 |   0.438095 |      0.57619  |

## Broad dataset split counts (seed 42)

### Virus counts

| host_genus     |   test |   train |   val |
|:---------------|-------:|--------:|------:|
| Acinetobacter  |     21 |     164 |    20 |
| Enterobacter   |      8 |      64 |     8 |
| Enterococcus   |      4 |      38 |     5 |
| Escherichia    |     11 |      94 |    12 |
| Klebsiella     |     13 |     104 |    13 |
| Pseudomonas    |     11 |      90 |    11 |
| Staphylococcus |      8 |      65 |     8 |

### Protein counts

| host_genus     |   test |   train |   val |
|:---------------|-------:|--------:|------:|
| Acinetobacter  |     60 |     734 |    72 |
| Enterobacter   |     55 |     473 |    93 |
| Enterococcus   |      4 |      67 |    15 |
| Escherichia    |    104 |     517 |    85 |
| Klebsiella     |     53 |     518 |    27 |
| Pseudomonas    |     36 |     251 |    35 |
| Staphylococcus |     25 |     148 |    24 |

## Strict dataset split counts (seed 42)

### Virus counts

| host_genus     |   test |   train |   val |
|:---------------|-------:|--------:|------:|
| Acinetobacter  |      5 |      42 |     5 |
| Enterobacter   |      1 |      10 |     1 |
| Enterococcus   |      1 |       9 |     1 |
| Escherichia    |      1 |       7 |     1 |
| Klebsiella     |      4 |      30 |     4 |
| Pseudomonas    |      1 |       8 |     1 |
| Staphylococcus |      1 |       7 |     1 |

### Protein counts

| host_genus     |   test |   train |   val |
|:---------------|-------:|--------:|------:|
| Acinetobacter  |      5 |      44 |     5 |
| Enterobacter   |      1 |      13 |     1 |
| Enterococcus   |      1 |       9 |     1 |
| Escherichia    |      1 |       8 |     1 |
| Klebsiella     |      4 |      32 |     4 |
| Pseudomonas    |      1 |       9 |     1 |
| Staphylococcus |      2 |       9 |     2 |

## Interpretation

- The broad dataset tests whether larger but noisier adsorption-module collections support host prediction.
- The strict dataset tests whether high-confidence RBPs provide cleaner host-specific signal.
- Comparing MLP vs linear probe helps distinguish representation quality from classifier instability.
- Macro-F1 is the most important metric here because host genera are not equally represented.
