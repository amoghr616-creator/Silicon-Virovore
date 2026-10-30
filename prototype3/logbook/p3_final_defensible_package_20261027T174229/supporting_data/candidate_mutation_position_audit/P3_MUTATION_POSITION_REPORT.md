# Candidate Mutation-Position Audit

Source: six existing candidate_history.json files from `p3_paired_arm_20261027T022810Z`.
No optimizer runs, new seeds, or P2 seed comparisons were performed.

## Overall comparison

| Metric | Baseline | Adaptive |
|---|---:|---:|
| Candidate records | 480 | 480 |
| Unique child sequences | 275 | 450 |
| Valid parent-child comparisons | 480 | 480 |
| Changed comparisons | 480 | 450 |
| Unchanged comparisons | 0 | 30 |
| Missing parent sequence | 0 | 0 |
| Invalid/unusable comparisons | 0 | 0 |
| Observed substitution events | 480 | 510 |
| Unique parent-child transitions | 275 | 473 |

Mean substitutions per valid parent-child comparison: baseline 1; adaptive 1.0625.

## Mutation rate by seed

| Arm | Seed | Valid comparisons | Changed comparisons | Changed % | Substitution events | Mean substitutions/valid comparison |
|---|---|---:|---:|---:|---:|---:|
| adaptive | seed2601 | 160 | 150 | 93.75% | 170 | 1.0625 |
| baseline | seed2601 | 160 | 160 | 100.00% | 160 | 1 |
| adaptive | seed2602 | 160 | 150 | 93.75% | 170 | 1.0625 |
| baseline | seed2602 | 160 | 160 | 100.00% | 160 | 1 |
| adaptive | seed2603 | 160 | 150 | 93.75% | 170 | 1.0625 |
| baseline | seed2603 | 160 | 160 | 100.00% | 160 | 1 |

## Position-level mutation frequencies

Percentages are the fraction of valid parent-child record comparisons that changed at that position. Positions are one-based.

| Position | Baseline changed | Baseline % | Adaptive changed | Adaptive % | Adaptive − baseline (percentage points) |
|---:|---:|---:|---:|---:|---:|
| 1 | 20 | 4.17% | 20 | 4.17% | +0.00 |
| 2 | 30 | 6.25% | 16 | 3.33% | -2.92 |
| 3 | 20 | 4.17% | 33 | 6.88% | +2.71 |
| 4 | 20 | 4.17% | 22 | 4.58% | +0.42 |
| 5 | 33 | 6.88% | 17 | 3.54% | -3.33 |
| 6 | 14 | 2.92% | 21 | 4.38% | +1.46 |
| 7 | 28 | 5.83% | 32 | 6.67% | +0.83 |
| 8 | 17 | 3.54% | 33 | 6.88% | +3.33 |
| 9 | 24 | 5.00% | 21 | 4.38% | -0.62 |
| 10 | 16 | 3.33% | 28 | 5.83% | +2.50 |
| 11 | 12 | 2.50% | 15 | 3.12% | +0.62 |
| 12 | 19 | 3.96% | 19 | 3.96% | +0.00 |
| 13 | 33 | 6.88% | 30 | 6.25% | -0.62 |
| 14 | 21 | 4.38% | 34 | 7.08% | +2.71 |
| 15 | 23 | 4.79% | 21 | 4.38% | -0.42 |
| 16 | 31 | 6.46% | 28 | 5.83% | -0.62 |
| 17 | 26 | 5.42% | 20 | 4.17% | -1.25 |
| 18 | 23 | 4.79% | 30 | 6.25% | +1.46 |
| 19 | 27 | 5.62% | 24 | 5.00% | -0.62 |
| 20 | 23 | 4.79% | 23 | 4.79% | +0.00 |
| 21 | 20 | 4.17% | 23 | 4.79% | +0.62 |

## Largest average per-seed position differences

Sorted by the absolute mean adaptive-minus-baseline difference across matched seeds. Positive means adaptive changed that position more often.

| Position | Mean difference (pp) | Median (pp) | Seed differences (pp) |
|---:|---:|---:|---|
| 5 | -3.33 | -5.00 | +0.62, -5.62, -5.00 |
| 8 | +3.33 | +2.50 | +1.88, +5.62, +2.50 |
| 2 | -2.92 | -3.12 | -4.38, -3.12, -1.25 |
| 3 | +2.71 | +1.88 | -0.62, +6.88, +1.88 |
| 14 | +2.71 | +1.88 | +1.88, +6.25, +0.00 |
| 10 | +2.50 | +1.88 | +5.00, +1.88, +0.62 |
| 6 | +1.46 | +0.62 | +0.62, -0.62, +4.38 |
| 18 | +1.46 | +1.25 | +0.62, +1.25, +2.50 |
| 17 | -1.25 | +0.00 | -4.38, +0.00, +0.62 |
| 7 | +0.83 | +2.50 | +2.50, -3.75, +3.75 |

## Unique parent-child transition check

The file `position_frequency_unique_transitions.csv` repeats the position analysis after counting each distinct parent-to-child sequence transition once per arm. Compare it with record-level frequencies to see whether repeated records are driving apparent differences.

## Reported mutation-count consistency

| Arm | Check | Comparable rows | Matches | Mismatches | Match % |
|---|---|---:|---:|---:|---:|
| baseline | mutation_count_matches_observed | 480 | 480 | 0 | 100.00% |
| baseline | changed_positions_count_matches_observed | 480 | 480 | 0 | 100.00% |
| adaptive | mutation_count_matches_observed | 480 | 480 | 0 | 100.00% |
| adaptive | changed_positions_count_matches_observed | 480 | 480 | 0 | 100.00% |

## Guidance metadata

| Arm | hotspot_guided field | mutation_guidance_mode | Records |
|---|---|---|---:|
| adaptive | False | hotspot_entropy_substitution_epistasis_exploratory_v2 | 48 |
| adaptive | True | hotspot_entropy_substitution_epistasis_exploratory_v2 | 432 |
| baseline | False | uniform_matched | 480 |

## Interpretation limits

- These are mutation-pattern measurements, not experimental binding-affinity measurements.
- A changed position is not automatically an ARISE-selected hotspot. The guidance flags are included as metadata, but their presence alone does not validate hotspot selection.
- Records may repeat sequences or related parent-child transitions; inspect both record-level and unique-transition summaries.
- The three matched seeds are summarized descriptively. Position-level percentages and their differences are not independent observations, so this report does not treat them as proof of statistical significance.
- Mutation positions are one-based sequence indices. They are not automatically residue numbers in the HERV-K Env receptor structure.

## Input files
- `prototype3/results/p3_paired_arm_20261027T022810Z/runs/adaptive/seed2601/p3_adaptive_2601/candidate_history.json: 160 records`
- `prototype3/results/p3_paired_arm_20261027T022810Z/runs/adaptive/seed2602/p3_adaptive_2602/candidate_history.json: 160 records`
- `prototype3/results/p3_paired_arm_20261027T022810Z/runs/adaptive/seed2603/p3_adaptive_2603/candidate_history.json: 160 records`
- `prototype3/results/p3_paired_arm_20261027T022810Z/runs/baseline/seed2601/p3_baseline_2601/candidate_history.json: 160 records`
- `prototype3/results/p3_paired_arm_20261027T022810Z/runs/baseline/seed2602/p3_baseline_2602/candidate_history.json: 160 records`
- `prototype3/results/p3_paired_arm_20261027T022810Z/runs/baseline/seed2603/p3_baseline_2603/candidate_history.json: 160 records`
