# P3 Logbook Entry — Data and Results

**Date:** 27 October 2026  
**Experiment:** Consolidated P3 optimization and retrospective validation audit  
**Purpose:** Preserve the defensible existing results in one reproducible evidence package.  
**Package:** `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/logbook/p3_final_defensible_package_20261027T174229`

## Design and objective

P3 addresses two distinct questions:

1. Does adaptive guidance improve the optimizer's internal objective relative to the matched uniform baseline?
2. Do the internal score, mutation-position signal, or structural-scoring method correspond to independently measured binding-affinity outcomes?

These questions are not interchangeable. An optimizer can improve its assigned score without improving real binding.

This package reuses the existing paired run, candidate records, metric files, logs, and completed validation results. **It does not launch an optimizer, add seeds, perform another seed-by-seed comparison, or rerun the validation scripts.** The figures are summaries of existing results, not new biological experiments.

## A. Existing optimizer and candidate inventory

The existing paired run contains 480 candidate records per arm.

| Measure | Uniform baseline | Adaptive guidance |
|---|---:|---:|
| Candidate records | 480 | 480 |
| Unique candidate sequences | 275 | 450 |
| Mean native score | 0.65224 | 0.74933 |
| Best native score | 0.71454 | 0.84553 |

The armwise score difference is approximately +0.0971 for adaptive minus baseline. The existing paired-run summary reports the same direction for mean score, best score, and final-generation mean across its three matched runs. The underlying run and original summary files are preserved in `supporting_data/` and `source_reports/`.

These are internal objective scores. They are not measured binding affinities.

## B. Mutation behavior audit

The existing parent-child mutation audit reported:

| Measure | Baseline | Adaptive |
|---|---:|---:|
| Parent-child records checked | 480 | 480 |
| Records with an observed sequence change | 480 | 450 |
| Substitution events counted | 480 | 510 |
| Unique parent-child transitions | 275 | 473 |
| Mutation-count consistency | 480/480 | 480/480 |

The metadata marked 432/480 adaptive records (90%) as hotspot-guided; the remaining 48 were not marked guided. All 480 baseline records used the `uniform_matched` mode.

The selected-position frequency differences are descriptive. The presence of a guidance flag does not demonstrate that a particular mutation was caused by a hotspot or that it improves binding.

## C. Broad native C-score versus experimental affinity

The retrospective benchmark included 693 matched single-substitution pairs, 621 unique peptides, peptide lengths of 9–20 amino acids, and 466 informative pairs.

- Spearman correlation between native score change and measured pKd change: **-0.1071**.
- Cluster-bootstrap 95% confidence interval: **[-0.2022, -0.0175]**.
- Direction accuracy on the reported informative comparison: **42.49%**, with a reported 95% interval of approximately **[36.55%, 47.77%]**.

The result does not support treating the native sequence-only score as a reliable predictor of experimental affinity changes. The native score has no receptor input; this benchmark tests the scoring function, not whether the adaptive optimizer produces experimentally better candidates.

## D. Exact-21-amino-acid native-score audit

The filtered audit retained 129 rows, 67 distinct peptide sequences, and 73 measurement contexts. Only seven contexts had enough distinct peptides and variability for the reported within-context correlation analysis.

- Mean within-context Spearman correlation: **-0.4872**.
- Group-bootstrap 95% interval: **[-0.7714, -0.0301]**.
- Informative single-substitution pairs: **9**.
- Direction accuracy on those pairs: **1/9 (11.11%)**.

The small, retrospective sample limits generalization. The source also reports a positive exploratory correlation on a different subset. That result is documented in `unresolved_evidence/EXACT_21AA_SUBSET_NOTE.md` and is not used as a headline finding until the differing subsets and sign conventions are reconciled.

## E. Mutation-position external validation

Across all 693 external mutation pairs, the positions 1–3 comparison reported a centered effect of **-0.0115**, with a target-cluster bootstrap interval of **[-0.1690, +0.1543]** and two-sided permutation p = **0.9218**. This does not show a general position 1–3 effect across all pairs.

Within the 466 informative pairs, the centered difference was **+0.1716**, interval **[+0.0148, +0.3419]**. The one-sided permutation p-value was **0.0200**, and the two-sided value was **0.0528**. This is a narrower exploratory signal and is not independent confirmation of a universal hotspot effect. The data are from a retrospective peptide benchmark, not direct HERV-K Env-specific experimental testing.

## F. Grouped sequence-only Ridge comparison

The recorded grouped validation used 693 pairs, 123 targets, 111 leakage components, and five folds. The audit reported zero target overlap, zero exact-pair overlap, zero peptide-sequence overlap, and one out-of-fold prediction per row.

| Metric | Training-fold mean | Ridge |
|---|---:|---:|
| MAE (lower is better) | 0.8182 | 0.8643 |
| RMSE (lower is better) | 1.3656 | 1.3957 |
| Spearman correlation | 0.0459 | 0.1087 |
| Direction accuracy | 62.66% | 58.96% |

Ridge ranked changes somewhat better by correlation, but the mean baseline had lower error and higher direction accuracy. This sequence-only benchmark is separate from the native C-score and the optimizer. The provenance path discrepancy is documented in `unresolved_evidence/GROUPED_VALIDATION_PROVENANCE_NOTE.md`.

## G. PEPBI strict-matched benchmark

The strict-matched audit reported 290 eligible rows, 69 distinct sequences, 18 matched contexts, and 382 pairwise comparisons.

Across 11 independent groups, mean concordance was:

- Native score: **0.3762**, 95% interval **[0.2110, 0.5315]**, reported one-sided sign-test p = **0.8867**.
- Rosetta delta G: **0.3927**, 95% interval **[0.2345, 0.5476]**, reported one-sided sign-test p = **0.7256**.

The target-reported experimental-delta-G subset had only one independent group, which is insufficient for a dependable general conclusion. These retrospective concordance values do not establish that either score reliably ranks the generated candidates.

The exact candidate-match count remains unresolved: one audit reported zero matches in both arms, while a separate Numbers worksheet was summarized as 125 adaptive matches and zero baseline matches. The latter count is not treated as confirmed evidence because a sequence-level match list and linked affinity records were not supplied in the audit. See `unresolved_evidence/PEPBI_MATCH_CONFLICT.md`.

## H. Structural-scoring benchmark

Of 693 considered experimental mutation pairs, only **18** were successfully mapped and scored.

- Overall Spearman correlation: **0.0588**, p = **0.8167**.
- Overall direction accuracy: **50%**.
- Informative subset: 14 pairs, correlation **0.0593**, p = **0.8403**.
- Same-receptor subset: 10 pairs across five targets, correlation **0.2849**, p = **0.425**, direction accuracy **60%**.

The pose sanity check reported a crystal-pose score of 0.96537 and a docked-pose score of 0.97751; under the recorded score orientation, the unfavorable docked pose received the preferred score. This indicates that the scoring setup needs more validation. The small mapped sample is not sufficient to claim that the structural score reliably predicts affinity.

## I. Receptor-context control

The recorded sequence-score receptor-context audit compared 24 matched receptor-variant pairs. The sequence-only native score was unchanged for zero pairs, while seven pairs had reported experimental energy differences greater than 0.5 kcal/mol.

This illustrates that a score without receptor input cannot directly encode receptor-specific differences. It is not a direct test of adaptive versus baseline binding.

## Files and figures

- `graphs/`: curated, fully labeled evidence plots.
- `poster_graphs/`: three suggested core poster figures.
- `supporting_data/`: source metrics, existing candidate tables, and mutation-audit outputs.
- `source_reports/`: existing report and execution artifacts.
- `source_archive/`: unchanged copy of the prior full audit.
- `unresolved_evidence/`: conflicts that must not be presented as confirmed findings.
- `evidence_manifest.csv`: checksums and paths for files in this package.
