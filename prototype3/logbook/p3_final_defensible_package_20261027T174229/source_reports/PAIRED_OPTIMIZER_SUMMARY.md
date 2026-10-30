# P3 paired optimizer run analysis

Source run directory: `prototype3/results/p3_paired_arm_20261027T022810Z`

This report summarizes the actual output files; it does not interpret an internal sequence score as experimentally measured affinity.

## Run-level summary

| Arm | Seed | Records | Unique sequences | Primary score field | Score mean | Score max | Final-generation mean | Parent-child comparisons |
|---|---:|---:|---:|---|---:|---:|---:|---:|
| baseline | 2601 | 160 | 138 | native_fitness | 0.65109 | 0.70555 | 0.64102 | 160 |
| baseline | 2602 | 160 | 135 | native_fitness | 0.65268 | 0.71454 | 0.64842 | 160 |
| baseline | 2603 | 160 | 125 | native_fitness | 0.65295 | 0.7072 | 0.65695 | 160 |
| adaptive | 2601 | 160 | 151 | native_fitness | 0.75933 | 0.8433 | 0.80516 | 160 |
| adaptive | 2602 | 160 | 151 | native_fitness | 0.75315 | 0.84553 | 0.80448 | 160 |
| adaptive | 2603 | 160 | 151 | native_fitness | 0.73549 | 0.83327 | 0.792 | 160 |

## Paired adaptive-minus-baseline differences

| Seed | Mean score Δ | Max score Δ | Final-generation mean Δ | Unique-sequence Δ |
|---:|---:|---:|---:|---:|
| 2601 | 0.10824 | 0.13775 | 0.16415 | 13 |
| 2602 | 0.10047 | 0.131 | 0.15606 | 16 |
| 2603 | 0.082544 | 0.12607 | 0.13504 | 26 |

## Mutation-guidance signal

- Metadata fields matching hotspot, selected-position, entropy, or guidance terms: **12**.
- Distinct observed parent-child mutation positions: **42**.
- Position count CSV: `mutation_position_counts.csv`.
- Metadata signal CSV: `adaptive_metadata_signals.csv`.

Inspect the metadata signal CSV to verify that the adaptive arm actually recorded selected hotspot positions and that those positions were used. A requested mode name alone does not prove effective guidance.

## Grouped Ridge versus training-fold mean

Metrics source: `prototype3/results/p3_paired_arm_20261027T022810Z/grouped_sequence_validation/run_20261027T023029Z/metrics.json`
- `informative_direction_accuracy`: 0.589622641509434
- `informative_direction_accuracy_ties_half`: 0.5815450643776824
- `informative_direction_coverage`: 0.9098712446351931
- `mae`: 0.8643157903240979
- `rmse`: 1.395672869959655
- `spearman_rho`: 0.10874027339042665
- `informative_direction_accuracy`: 0.6266094420600858
- `informative_direction_accuracy_ties_half`: 0.6266094420600858
- `informative_direction_coverage`: 1.0
- `mae`: 0.818212554563056
- `rmse`: 1.3655814586084296
- `spearman_rho`: 0.045870601293509276
- `informative_direction_accuracy_ties_half`: 0.5
- `informative_direction_coverage`: 0.0
- `mae`: 0.825359312626174
- `rmse`: 1.3840433664872516
- `valid_replicates_direction_tie_half`: 2000
- `valid_replicates_mae`: 2000
- `valid_replicates_spearman`: 2000

## What this run does not establish

- These three matched seeds compare optimizer outputs; they are not the prior P2 100-seed experiment.
- Native C fitness is a sequence-level heuristic, not experimental affinity or a validated target-specific binding score.
- ARISE validation using the positions actually identified in these new runs remains a separate analysis.
- PEPBI adaptive-versus-baseline prediction accuracy requires distinct arm-specific predictions on the same measured cases; copying one fixed score into two labels is invalid.
- Structural scoring of generated candidates requires valid candidate structures/poses and identical scoring/eligibility rules for both arms.
- Three seeds are a small exploratory optimizer comparison; do not claim statistical superiority from descriptive differences alone.
