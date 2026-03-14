# Design Run Summary

## Run metadata

- Seed protein ID: `QFR57578.1`
- Virus accession: `MN434096`
- Source host: `Klebsiella`
- Target host: `Pseudomonas`
- ESM model: `facebook/esm2_t33_650M_UR50D`
- Rounds: `4`
- Candidates per round: `64`
- Min mutations: `1`
- Max mutations: `3`
- Top-K kept per round: `10`
- Proposal top-K: `8`
- Seed: `42`
- Total candidates evaluated: `1984`
- Top candidates saved: `50`

## Target score summary

| metric   |   target_score |
|:---------|---------------:|
| min      |       0.15581  |
| mean     |       0.234041 |
| median   |       0.254156 |
| max      |       0.286722 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       7 |
|             2 |      16 |
|             3 |      27 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|:------------------|--------------:|
| round4_cand12  | Klebsiella    | Pseudomonas   | Enterobacter |       0.286722 | W310V;N463V       |             2 |
| round4_cand100 | Klebsiella    | Pseudomonas   | Enterobacter |       0.268047 | P324A;T543L       |             2 |
| round4_cand140 | Klebsiella    | Pseudomonas   | Enterobacter |       0.258882 | N153L;S550L;N594T |             3 |
| round4_cand347 | Klebsiella    | Pseudomonas   | Enterobacter |       0.249429 | N20P;I185N;K561P  |             3 |
| round4_cand212 | Klebsiella    | Pseudomonas   | Enterobacter |       0.248504 | N300L;N463T;S546I |             3 |
| round4_cand299 | Klebsiella    | Pseudomonas   | Enterobacter |       0.244175 | N320A;Y615P       |             2 |
| round4_cand549 | Klebsiella    | Pseudomonas   | Enterobacter |       0.236383 | N153L;M334P;E460Y |             3 |
| round4_cand497 | Klebsiella    | Pseudomonas   | Enterobacter |       0.232044 | D52N;D449P;W469A  |             3 |
| round3_cand300 | Klebsiella    | Pseudomonas   | Enterobacter |       0.229886 | E88A;N313A;K561A  |             3 |
| round4_cand413 | Klebsiella    | Pseudomonas   | Enterobacter |       0.227163 | E146D;G267T;Y444A |             3 |
