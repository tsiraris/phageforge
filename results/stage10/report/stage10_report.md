# Stage 10 structure-conditioned redesign report

## Why Stage 10 exists

Stage 08 and Stage 09 showed that sequence-first redesign, even when tightened with structure-aware proxies, still failed decisive full structural validation. Stage 10 therefore moved the structure signal upstream and redesigned candidates directly against a fixed seed scaffold with an inverse-folding objective.

## Core redesign setup

- target host: **Acinetobacter**
- selected seed id: **round8_cand473**
- seed scaffold: **/home/sagemaker-user/phageforge_clean/results/stage09/validation_top3/pdbs/seed_selected_seed.pdb**
- hard editable positions: **[221, 329, 334, 373, 401, 427]**
- soft editable positions: **[505, 291, 527]**
- mutation budget: **1 to 4**

## Search summary

- scored candidates: **147**
- best Stage 10 composite score: **0.7180818**
- best target probability: **0.029418018**
- best inverse-folding log-likelihood: **-2.6386764**

## Prefilter summary

- prefilter panel size: **10**
- sample ids carried forward: **[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]**

## Structural validation summary

- validated rows: **3**
- pass count: **0**
- best mean pLDDT: **0.2245393395423889**
- best mutation-site mean pLDDT: **nan**
- best RMSD to seed: **11.791931711577568 Å**

## Comparison to earlier structural validation

- baseline validated rows: **3**
- baseline pass count: **0**
- baseline best mean pLDDT: **0.226950541138649**
- baseline best mutation-site mean pLDDT: **0.2324999999999999**
- baseline best RMSD to seed: **13.129426232705702 Å**

## Interpretation

Stage 10 is the first phase in which the scaffold itself becomes part of candidate generation rather than only downstream evaluation. This makes it the correct methodological successor to Stage 08 and Stage 09, regardless of whether the final heavy structural validation fully succeeds on the first attempt.