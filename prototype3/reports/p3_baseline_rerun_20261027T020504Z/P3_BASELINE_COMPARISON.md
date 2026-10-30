# P3 Baseline Comparison Report

Generated: 2026-10-27T02:10:04.812201+00:00

This report consolidates the current metrics. Check the run-status section below: if a rerun failed, the corresponding result file may still contain an earlier result. This report does not claim that failed steps were reproduced.

## 1. ARISE mutation-position validation

**Primary comparison — positions 1–3 vs. random equal-width windows**

- All pairs: context-adjusted difference = -0.0115; 95% CI [-0.1690, 0.1543]; two-sided permutation p = 0.9218.
- Informative subset: difference = 0.1716; 95% CI [0.0148, 0.3419]; two-sided permutation p = 0.0528.
- Previous-window sensitivity check, positions 2–3: difference = -0.0184; two-sided p = 0.8454.

Interpretation: the all-pairs result is the primary result. The filtered-subset finding is exploratory.

## 2. Native C score vs. experimental mutation effects

- Matched mutation pairs: 693.
- Spearman correlation, score change vs. pKd change: -0.1071; 95% cluster-bootstrap CI [-0.2022, -0.0175].
- Informative-pair direction accuracy: 0.4249; 95% CI [0.3655, 0.4777]; nominal direction reference = 50%.

Interpretation: compare against the 50% direction reference; do not describe this as a target-aware score.

## 3. Exact 21-amino-acid benchmark

- Eligible measurement rows: 129.
- Distinct scored peptides: 67.
- Within-context mean Spearman rho: -0.4872; 95% group-bootstrap CI [-0.7714, -0.0301].
- Informative mutation pairs: 0.1111111111111111 direction accuracy (one observed correct-direction pair out of nine in the recorded analysis); reference = 50%.

Interpretation: the small informative-pair sample makes the direction estimate unstable.

## 4. Grouped sequence-model baseline

**training_fold_mean_baseline**
- MAE: 0.8182; RMSE: 1.3656; Spearman rho: 0.0459; direction accuracy: 0.6266.

**sequence_only_ridge**
- MAE: 0.8643; RMSE: 1.3957; Spearman rho: 0.1087; direction accuracy: 0.5896.

Source: `prototype3/results/sequence_only_grouped_target_validation/run_20261027T020510Z/metrics.json`.
Interpretation: Ridge is a separate comparator, not the production native C score.

## 5. PEPBI strict matched comparison

**Target: Calculated ΔG from K_D**
- native_score: mean group concordance 0.3762, 11 independent binding groups, 95% group-bootstrap CI [0.2110, 0.5315]; chance reference = 0.5000.
- rosetta_dg: mean group concordance 0.3927, 11 independent binding groups, 95% group-bootstrap CI [0.2345, 0.5475]; chance reference = 0.5000.

**Target: Reported experimental ΔG**
- native_score: mean group concordance 0.0000, 1 independent binding groups, 95% group-bootstrap CI [0.0000, 0.0000]; chance reference = 0.5000.
- rosetta_dg: mean group concordance 0.0000, 1 independent binding groups, 95% group-bootstrap CI [0.0000, 0.0000]; chance reference = 0.5000.

Interpretation: preserve the strict matched result and the Rosetta comparison. Do not infer successful affinity prediction from the pooled PEPBI correlations alone.

## 6. Receptor-context invariance control

- Matched peptide contexts: 4.
- Matched receptor-variant pairs: 24.
- Pairs with any native-score change: 0.
- Pairs with an experimental ΔG difference greater than 0.5 kcal/mol: 7.

Interpretation: a peptide-only score is invariant when the peptide sequence remains unchanged; this is a limitation for receptor-specific prediction, not evidence of binding accuracy.

## 7. Structural benchmark carried forward from existing results

- Scored pairs: 18; target groups: 7.
- Spearman rho: 0.0588; direction accuracy: 0.5000; chance reference = 50%.
- Source: `prototype3/results/experimental_mutation_benchmark/structural_scoring_20261026T223525Z/benchmark_metrics.json`.
- Status: existing result retained; it was not rerun because no corresponding standalone runner appears in the listed `prototype3/scripts` directory.

## Run status

```text
PASS	ppikb_preprocessing
PASS	grouped_sequence_validation
PASS	native_score_mutation_pairs
PASS	exact_21aa_benchmark
PASS	arise_position_validation
PASS	pepbi_corrected_benchmark
PASS	pepbi_strict_matched
PASS	pepbi_target_invariance
PASS	pepbi_verified_substitutions
PASS	pepbi_native_exploratory
```

## Overall interpretation

Baselines make the comparison more interpretable; they do not change the observed metrics. These retrospective tests assess score/affinity relationships and mutation-position associations. They do not establish prospective HERV-K Env binding or biological efficacy.

Preserve the backup and original logs. Do not replace a failed rerun with an older result without labeling it.
