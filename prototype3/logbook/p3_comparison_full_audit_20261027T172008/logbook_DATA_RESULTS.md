# Logbook Entry — P3 Comparison: Data and Results

**Date:** 27 October 2026  
**Audit generated:** 2026-10-27T17:20:18-04:00  
**Experiment name:** P3 paired adaptive-versus-baseline comparison, candidate diversity,
mutation-position audit, and external-validation reconciliation  
**Primary run:** `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/p3_paired_arm_20261027T022810Z`  
**Audit output:** `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/logbook/p3_comparison_full_audit_20261027T172008`

## Design and objective

The purpose of this audit was to consolidate the most recent paired adaptive-versus-baseline
optimizer outputs with the available independent affinity-validation benchmarks. The analysis
keeps two questions separate:

1. Did adaptive guidance improve the optimizer's internal objective and alter the candidate
   sequence set compared with the uniform baseline?
2. Do the scoring and candidate-selection results correspond to experimentally measured
   binding-affinity changes?

The paired optimizer run used three matched seeds (2601, 2602, and 2603). The external
benchmarks use different datasets and sample sizes; they are reported separately rather
than pooled into a single performance number.

## Source files and reproducibility

- Selected paired run: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/p3_paired_arm_20261027T022810Z`
- Six-area comparison report: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/logbook/p3_direct_armwise_20261027T204305Z/P3_DIRECT_ARMWISE_REPORT.md`
- Source files inventoried: **175**
- Source files copied into the snapshot: **175**
- Files by extension: `.csv` 61, `.json` 61, `.log` 8, `.md` 19, `.png` 25, `.tsv` 1
- Candidate sequence-bearing records parsed from the paired run: **2916**
- JSON scalar values extracted: **413597**
- Supporting evidence `native_mutation_benchmark`: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/native_score_all_mutation_pairs`
- Supporting evidence `21aa_affinity_audit`: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/p3_21aa_affinity_audit`
- Supporting evidence `hotspot_position_audit`: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/arise_experimental_position_audit`
- Supporting evidence `pepbi_strict_matched`: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/pepbi_strict_matched`
- Supporting evidence `structural_benchmark`: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/experimental_mutation_benchmark`
- Supporting evidence `grouped_sequence_validation`: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/sequence_only_grouped_target_validation/run_20261027T205129Z`
- Supporting evidence `candidate_mutation_position_audit`: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/logbook/p3_mutation_position_audit_20261027T205650Z`

The `source_file_inventory.csv` records file paths, file sizes, modification times, SHA-256
hashes, and whether each file was copied into `source_snapshot/`. JSON scalar values are
indexed in `json_scalar_inventory.csv`. Sequence-bearing records that the parser could
identify from the paired run are exported to `candidate_records_extracted.csv`; the
extraction summary is `candidate_extraction_summary.csv`. The raw snapshots remain the
authoritative record if a data structure was not recognized by the parser.

## Results A — Paired optimizer outputs

| Seed | Baseline mean | Adaptive mean | Baseline best | Adaptive best | Baseline final-gen mean | Adaptive final-gen mean | Baseline unique | Adaptive unique |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2601 | 0.65109 | 0.75933 | 0.70555 | 0.84330 | 0.64102 | 0.80516 | 138 | 151 |
| 2602 | 0.65268 | 0.75315 | 0.71454 | 0.84553 | 0.64842 | 0.80448 | 135 | 151 |
| 2603 | 0.65295 | 0.73549 | 0.70720 | 0.83327 | 0.65695 | 0.79200 | 125 | 151 |

Across the available 480 candidate records in each arm, the previously reported inventory
contained 275 unique baseline sequences and 450 unique adaptive sequences. The corresponding
reported mean native scores were 0.65224 for baseline and 0.74933 for adaptive; the best
reported scores were 0.71454 and 0.84553, respectively. These are optimizer-objective
results, not direct experimental binding measurements.

The per-seed summary shows the same direction of the internal-score difference in all three
matched seeds. The experiment contains three seed pairs, so this is a small, exploratory
optimizer comparison; the seed summaries are not a substitute for an independent affinity
test.

## Results B — Mutation behavior and sequence diversity

The candidate mutation-position audit reported:

| Measure | Baseline | Adaptive |
|---|---:|---:|
| Parent-child records checked | 480 | 480 |
| Records with a sequence change | 480 (100%) | 450 (93.75%) |
| Mutation events counted | 480 | 510 |
| Unique parent-child sequence transitions | 275 | 473 |
| Unique candidate sequences in the candidate inventory | 275 | 450 |
| Mutation-count consistency | 480/480 | 480/480 |

The reported metadata identified 432/480 adaptive records (90%) as hotspot-guided and
48/480 (10%) as not hotspot-guided. All 480 baseline records used `uniform_matched`.
The metadata flag by itself does not prove that each particular substitution was caused by
or located within a selected hotspot.

The largest reported adaptive-minus-baseline position-frequency differences were: position
5, -3.33 percentage points; position 8, +3.33; position 2, -2.92; positions 3 and 14,
+2.71 each; and position 10, +2.50. These position differences are descriptive across only
three seeds and should not be represented as statistically significant.

## Results C — New PEPBI exact-sequence audit (provisional)

The newly supplied Numbers worksheet reports:

| Arm | Unique candidate sequences | Reported exact PEPBI matches | Reported match rate |
|---|---:|---:|---:|
| Baseline | 275 | 0 | 0% |
| Adaptive | 450 | 125 | 27.78% |

**Status: provisional and in conflict with an earlier audit that reported zero exact
matches for both arms.** The worksheet's visible `match_sequences` field did not identify
the 125 matching sequences. This audit therefore preserves the new counts as reported but
does not treat them as confirmed candidate-level experimental validation. Reconciliation
requires the actual matched sequence list, the source PEPBI files, matching rules, duplicate
handling, and the affinity/receptor-context records associated with each match.

## Results D — Native sequence-score experimental benchmark

The broad mutation benchmark contained 693 matched single-substitution pairs, 621 unique
peptides, and 466 informative pairs. The reported Spearman correlation between score change
and experimental pKd change was -0.1071, with cluster-bootstrap 95% CI
[-0.2022, -0.0175]. Direction accuracy was 42.49%, with 95% CI approximately
[36.55%, 47.78%]. A shuffle test was reported as p=0.223.

The exact-21-amino-acid audit contained 129 rows, 67 distinct peptides, and 73 measurement
contexts, but only seven contexts had enough variation for within-context analysis. The
reported mean within-context Spearman correlation was -0.487. Twenty substitution pairs
were available, only nine of them informative, and direction accuracy was 1/9 (11.11%).
The same report also contained a positive correlation summary for an informative subset
(rho=+0.678); the discrepancy between these subset summaries should be audited before
either is quoted as a definitive result.

## Results E — Mutation-position external validation

In the external PPIKB mutation dataset, the all-pair position 1–3 comparison reported a
centered effect of -0.0115, bootstrap 95% CI [-0.1690, +0.1543], and two-sided permutation
p=0.9218. This does not show a general effect across all pairs.

For the informative subset, the reported centered difference in absolute affinity change
was +0.1716, CI [+0.0148, +0.3419], with two-sided permutation p=0.05279 and one-sided
p=0.020. This narrower subset is exploratory and should not be presented as independent
confirmation of a general mutation-position effect.

## Results F — Grouped sequence-only prediction benchmark

The grouped validation used 693 pairs from 123 targets and 111 components across five
folds. The recorded leakage checks reported zero target overlap, zero exact-pair overlap,
zero peptide-sequence overlap, and one out-of-fold prediction per row.

| Metric | Training-fold mean | Sequence-only Ridge |
|---|---:|---:|
| MAE | 0.81821 | 0.86432 |
| RMSE | 1.36558 | 1.39567 |
| Spearman rho | 0.04587 | 0.10874 |
| Direction accuracy | 62.66% | 58.96% |

The paired bootstrap interval for Ridge minus mean MAE was [+0.0173, +0.0782].
The Ridge model had a higher reported rank correlation, but worse MAE and direction
accuracy. No single model dominated every metric.

## Results G — PEPBI strict-matched benchmark

The strict-matched audit reported 290 eligible rows, 69 distinct sequences, 18 strict
contexts, and 382 pairwise comparisons. Across 11 independent groups, mean concordance
was 0.3762 for the native score (95% CI [0.2110, 0.5315]) and 0.3927 for Rosetta delta G
(95% CI [0.2345, 0.5476]). Reported one-sided sign-test p-values were 0.8867 and 0.7256,
respectively. The target-reported experimental-delta-G subset had only one independent
group, which is insufficient for a reliable conclusion.

## Results H — Structural scoring benchmark

The structural benchmark began with 693 candidate mutation pairs, but only 18 mapped pairs
were successfully scored, across seven targets. Overall correlation was 0.0588 (p=0.8167)
and directional accuracy was 50%. The informative subset of 14 pairs had correlation
0.0593 (p=0.8403). The same-receptor subset had ten pairs across five targets, correlation
0.2849 (p=0.425), and directional accuracy 60%; uncertainty is substantial at this sample
size.

The pose sanity check reported a crystal-pose score of 0.96537 and a docked-pose score of
0.97751. Under the score orientation used in that audit, the score preferred the unfavorable
docked pose, so the scoring setup requires further validation before strong claims are made.

## Figures generated

- `graphs/01_paired_native_score_by_seed.png`
- `graphs/02_candidate_diversity_and_provisional_pepbi_matches.png`
- `graphs/03_mutation_position_frequency_differences.png`
- `graphs/04_external_affinity_correlations.png`
- `graphs/05_external_direction_accuracy.png`
- `graphs/06_ridge_vs_training_fold_mean_mae.png`
- `graphs/07_pepbi_strict_matched_concordance.png`

## Data integrity notes

1. The 125-match adaptive PEPBI count remains provisional pending source-level reconciliation.
2. The 21-aa native-score summary contains subset-level results that need reconciliation.
3. The internal native score is not itself an independent receptor-aware binding measurement.
4. External mutation and structural benchmarks are retrospective method-validation studies;
   they do not automatically validate the actual candidates generated in the optimizer run.
5. Where the extraction parser could not recognize a file's schema, the preserved raw source
   snapshot and manifest should be used rather than inferring that the data were absent.
