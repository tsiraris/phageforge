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
- Total candidates evaluated: `1982`
- Top candidates saved: `50`

## Target score summary

| metric   |   target_score |
|:---------|---------------:|
| min      |       0.272466 |
| mean     |       0.281248 |
| median   |       0.279834 |
| max      |       0.302903 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Enterobacter |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       6 |
|             2 |      15 |
|             3 |      29 |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score | mutations         |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|:------------------|--------------:|
| round4_cand176 | Klebsiella    | Pseudomonas   | Enterobacter |       0.302903 | S456A;N594A;T607V |             3 |
| round4_cand353 | Klebsiella    | Pseudomonas   | Enterobacter |       0.300055 | N445L;S546L;E630R |             3 |
| round4_cand241 | Klebsiella    | Pseudomonas   | Enterobacter |       0.2959   | T151V;Y453T;F579A |             3 |
| round4_cand236 | Klebsiella    | Pseudomonas   | Enterobacter |       0.295787 | F36P;K325E;E353A  |             3 |
| round4_cand100 | Klebsiella    | Pseudomonas   | Enterobacter |       0.294331 | N409P;G488A       |             2 |
| round4_cand386 | Klebsiella    | Pseudomonas   | Enterobacter |       0.294223 | N29T;K335G;W465A  |             3 |
| round4_cand133 | Klebsiella    | Pseudomonas   | Enterobacter |       0.292662 | Y342E;L368A;D549S |             3 |
| round4_cand66  | Klebsiella    | Pseudomonas   | Enterobacter |       0.28859  | I391F;R392A;N416R |             3 |
| round4_cand69  | Klebsiella    | Pseudomonas   | Enterobacter |       0.287262 | K12L;T129F;H138G  |             3 |
| round4_cand16  | Klebsiella    | Pseudomonas   | Enterobacter |       0.286791 | W310V;N463V       |             2 |
