# Paired sequence-model comparison

- Predictions: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/sequence_only_mutation_baseline/run_20261026T234823Z`
- Resampling unit: target group, not individual pair.
- Bootstrap: 2,000 resamples; intervals condition on fixed out-of-fold predictions.

## Observed performance differences

| Difference (Ridge minus fold-mean baseline) | Estimate | 95% target-bootstrap interval |
|---|---:|---:|
| mae_delta_ridge_minus_fold_mean | 0.0213 | [-0.0193, 0.0588] |
| rmse_delta_ridge_minus_fold_mean | -0.0100 | [-0.0511, 0.0298] |
| spearman_delta_ridge_minus_fold_mean | 0.3270 | [0.1870, 0.4567] |
| direction_delta_ridge_minus_fold_mean | 0.0408 | [-0.0141, 0.0951] |

## Reading the differences

- Negative MAE or RMSE differences favor Ridge.
- Positive Spearman or direction-score differences favor Ridge.
- If an interval contains zero, this comparison does not clearly establish a difference.
- This compares a sequence-only comparator with a fold-mean baseline. It does not validate `c_score` or the structural score.
