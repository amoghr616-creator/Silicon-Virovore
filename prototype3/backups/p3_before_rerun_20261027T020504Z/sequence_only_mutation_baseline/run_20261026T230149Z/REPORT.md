# Sequence-only mutation-effect benchmark

- Run: `run_20261026T230149Z`
- Input: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/ppikb_kd_benchmark/filtered_kd_measurements.csv`
- Peptide column: `Peptide_Sequence`
- Target grouping column: `Protein_Sequence`
- Affinity column: `_p_affinity`
- Context columns: `_assay, _affinity_type`
- Split: five-fold GroupKFold grouped by target; no target overlap across train/test folds.
- Model: fixed Ridge regression (`alpha=10.0`); no hyperparameter selection on the held-out data.

## Data counts

- Input rows: 3439
- Valid rows after filtering: 3439
- Unique target/context/peptide measurements after median aggregation: 2496
- Matched single-substitution pairs: 693
- Distinct target groups: 123
- Pairs with absolute experimental change >= 0.3 pKd: 466

## Out-of-fold results

| Method | Pairs | MAE (pKd) | RMSE (pKd) | Spearman rho | Direction accuracy for |delta pKd| >= 0.3 |
|---|---:|---:|---:|---:|---:|
| Predict no change | 693 | 0.8254 | 1.3840 | None | 0.0 |
| Training-fold mean | 693 | 0.8204 | 1.3734 | -0.08411059247873977 | 0.6266094420600858 |
| Sequence-only Ridge | 693 | 0.8418 | 1.3634 | 0.24293544687990765 | 0.6673819742489271 |

## Cluster-bootstrap uncertainty

- Spearman rho 95% CI: `[0.14729049978348047, 0.33721265232300474]`
- Informative direction accuracy 95% CI: `[0.601446451050644, 0.7329791099000909]`

## Interpretation

This benchmark tests a simple sequence-only comparator on retrospective experimental mutation pairs. It does not use or validate the native `c_score`, and it does not validate the target-aware structural score. The target-held-out split reduces target leakage, but assay heterogeneity, data selection, and limited generalization remain relevant. Do not claim real-world binding or efficacy from this benchmark.
