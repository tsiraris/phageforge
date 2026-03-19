# Design Run Summary

## Run metadata

- Seed protein ID: `QFR57578.1`
- Virus accession: `MN434096`
- Source host: `Klebsiella`
- Target host: `Pseudomonas`
- ESM model: `facebook/esm2_t33_650M_UR50D`
- Rounds: `6`
- Candidates per round: `64`
- Min mutations: `1`
- Max mutations: `4`
- Top-K kept per round: `10`
- Proposal top-K: `8`
- Seed: `44`
- Total candidates evaluated: `3263`
- Top candidates saved: `50`
- Diversity enabled: `False`
- Diversity min distance: `8`
- Position selection strategy: `random`
- Position pool size: `N/A`

## Target score summary

| metric   |   target_score |
|:---------|---------------:|
| min      |       0.410205 |
| mean     |       0.425513 |
| median   |       0.423241 |
| max      |       0.462639 |

## Predicted label distribution in top candidates

| pred_label   |   count |
|:-------------|--------:|
| Pseudomonas  |      50 |

## Mutation count distribution in top candidates

|   n_mutations |   count |
|--------------:|--------:|
|             1 |       6 |
|             2 |      11 |
|             3 |      14 |
|             4 |      19 |

## Diversity summary

| metric                                    |    value |
|:------------------------------------------|---------:|
| avg_pairwise_distance_top_candidates      |  23.1159 |
| unique_mutation_site_count_top_candidates | 125      |

## Nearest-neighbor novelty summary

| metric                                    |       value |
|:------------------------------------------|------------:|
| min_novelty_distance                      | 0.000788331 |
| mean_novelty_distance                     | 0.0010894   |
| median_novelty_distance                   | 0.00105718  |
| max_novelty_distance                      | 0.00147247  |
| min_nearest_neighbor_cosine_similarity    | 0.998528    |
| mean_nearest_neighbor_cosine_similarity   | 0.998911    |
| median_nearest_neighbor_cosine_similarity | 0.998943    |
| max_nearest_neighbor_cosine_similarity    | 0.999212    |

## Nearest-neighbor novelty preview (top 10)

| candidate_id   | nearest_neighbor_protein_id   | nearest_neighbor_host_genus   |   nearest_neighbor_cosine_similarity |   novelty_distance |
|:---------------|:------------------------------|:------------------------------|-------------------------------------:|-------------------:|
| round6_cand67  | QFR57578.1                    | Klebsiella                    |                             0.999212 |        0.000788331 |
| round6_cand2   | QFR57578.1                    | Klebsiella                    |                             0.998918 |        0.00108182  |
| round6_cand35  | QFR57578.1                    | Klebsiella                    |                             0.998698 |        0.00130236  |
| round6_cand25  | QFR57578.1                    | Klebsiella                    |                             0.998846 |        0.00115448  |
| round6_cand43  | QFR57578.1                    | Klebsiella                    |                             0.998528 |        0.00147247  |
| round6_cand268 | QFR57578.1                    | Klebsiella                    |                             0.998593 |        0.00140727  |
| round6_cand172 | QFR57578.1                    | Klebsiella                    |                             0.998783 |        0.00121725  |
| round6_cand40  | QFR57578.1                    | Klebsiella                    |                             0.99864  |        0.00136018  |
| round6_cand272 | QFR57578.1                    | Klebsiella                    |                             0.998938 |        0.00106174  |
| round6_cand297 | QFR57578.1                    | Klebsiella                    |                             0.999065 |        0.00093472  |

## Top 10 candidates

| candidate_id   | source_host   | target_host   | pred_label   |   target_score | mutations               |   n_mutations |
|:---------------|:--------------|:--------------|:-------------|---------------:|:------------------------|--------------:|
| round6_cand67  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.462639 | N108P;E344P;D423A;N524A |             4 |
| round6_cand2   | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.458929 | F94P;G159R;Q560P;N642R  |             4 |
| round6_cand35  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.451613 | F94R;G297F;G431A        |             3 |
| round6_cand25  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.44557  | F94P;N263T;T607P        |             3 |
| round6_cand43  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.444681 | S38G;M65D;F94R;R289A    |             4 |
| round6_cand268 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.444263 | H138A;R257V;S456L;M596G |             4 |
| round6_cand172 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.438691 | G2P;V72R;G307P;I407V    |             4 |
| round6_cand40  | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.43433  | R317I;A529R             |             2 |
| round6_cand272 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.434133 | R68G;Y77N;D449A;N557I   |             4 |
| round6_cand297 | Klebsiella    | Pseudomonas   | Pseudomonas  |       0.433131 | K155R;N367R             |             2 |
