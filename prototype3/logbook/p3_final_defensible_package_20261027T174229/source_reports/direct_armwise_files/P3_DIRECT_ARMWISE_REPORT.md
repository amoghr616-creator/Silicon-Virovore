# P3 Direct Armwise Audit

Generated: 2026-10-27T20:51:36.765768+00:00

## Scope and exclusions

This patch did not launch the optimizer, create new seeds, or rerun the P2 baseline-versus-adaptive seed comparison. It reuses existing candidate artifacts and reruns the listed P3 validation scripts.

## Direct arm-specific candidate inventory

The inventory reports only records whose arm can be identified from the file path or row fields. Records labeled `unknown` are not silently assigned to an arm.

| Arm | Candidate records | Unique sequences | Scored records | Mean native score | Best native score |
|---|---:|---:|---:|---:|---:|
| baseline | 480 | 275 | 480 | 0.6522404876005785 | 0.7145354444591977 |
| adaptive | 480 | 450 | 480 | 0.7493252275618555 | 0.8455310975436344 |
| unknown | 0 | 0 | 0 | n/a | n/a |

## Direct candidate-level PEPBI sequence audit

| Arm | Unique candidates | Exact PEPBI sequence matches | Interpretation |
|---|---:|---:|---|
| baseline | 275 | 0 | No exact matches located in scanned PEPBI files; accuracy is unavailable, not zero. |
| adaptive | 450 | 0 | No exact matches located in scanned PEPBI files; accuracy is unavailable, not zero. |

An exact match count is not an accuracy metric. If matching sequences exist, their measured outcomes and experimental contexts must be validated before calculating arm-specific prediction metrics.

## Status of every requested P3 area

| Area | Direct baseline-vs-adaptive evidence available? | What this run can defensibly report |
|---|---|---|
| Native C-score / broad mutation benchmark | No arm-specific measured affinity labels established by this audit | Shared score-validity benchmark only; cannot infer which arm is more accurate. |
| Exact 21-aa affinity audit | No arm-specific measured affinity labels established by this audit | Shared length-restricted score-validity test; not two arm accuracy values. |
| Ridge vs fold-mean | No, unless candidate-specific predictions and labels exist | Prediction-method benchmark; fold-mean is not the baseline optimizer. |
| Mutation positions | Partly | Arm-specific mutation frequencies can be compared only when candidate records expose mutation positions; external positional validation remains shared. |
| PEPBI | Exact sequence overlap audit performed by arm | Arm-specific accuracy is calculable only if exact matches have valid measurements. |
| Receptor-context control | No | Shared test of sequence-only score invariance; not arm-specific binding accuracy. |
| Structural scoring | Only if comparable candidate-run poses are found and verified | Potential PDBs discovered: 0. File presence alone does not establish a valid pose or score. |

## P3 validation script execution status
## broad_native_c_score
Script: `benchmark_native_score_all_mutation_pairs.py`
Status: PASS
Exit code: 0
Log: `prototype3/logbook/p3_direct_armwise_20261027T204305Z/broad_native_c_score.log`

## exact_21aa
Script: `benchmark_native_c_score_21aa.py`
Status: PASS
Exit code: 0
Log: `prototype3/logbook/p3_direct_armwise_20261027T204305Z/exact_21aa.log`

## mutation_position
Script: `validate_arise_hotspot_vs_experimental_positions.py`
Status: PASS
Exit code: 0
Log: `prototype3/logbook/p3_direct_armwise_20261027T204305Z/mutation_position.log`

## ridge_vs_fold_mean
Script: `run_grouped_sequence_pair_validation.py`
Status: PASS
Exit code: 0
Log: `prototype3/logbook/p3_direct_armwise_20261027T204305Z/ridge_vs_fold_mean.log`

## pepbi_strict_matched
Script: `pepbi_strict_matched_benchmark.py`
Status: PASS
Exit code: 0
Log: `prototype3/logbook/p3_direct_armwise_20261027T204305Z/pepbi_strict_matched.log`

## receptor_context_control
Script: `pepbi_target_mutation_invariance.py`
Status: PASS
Exit code: 0
Log: `prototype3/logbook/p3_direct_armwise_20261027T204305Z/receptor_context_control.log`


## Existing metric files discovered

### `prototype3/results/native_score_all_mutation_pairs/metrics.json`
```json
{
  "benchmark": "Native sequence fitness versus experimental mutation effects",
  "pair_file": "/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/sequence_only_mutation_baseline/run_20261026T234823Z/single_substitution_pairs.csv",
  "score_function": "native_score(sequence) -> c_evaluate_sequence -> fitness_score",
  "n_matched_single_substitution_pairs": 693,
  "n_unique_peptide_sequences_freshly_scored": 621,
  "peptide_length_min": 9,
  "peptide_length_max": 20,
  "n_exact_21aa_pairs": 0,
  "bootstrap_unit": "Connected components formed by shared target groups or individual peptide sequences",
  "n_independent_components": 111,
  "bootstrap_replicates_requested": 2000,
  "observed_metrics": {
    "n_pairs": 693,
    "spearman_delta_score_vs_delta_pKd": -0.10713922442596775,
    "n_informative_pairs": 466,
    "direction_accuracy_ties_half": 0.4248927038626609,
    "spearman_abs_delta_score_vs_abs_delta_pKd": -0.041827830390897214
  },
  "cluster_bootstrap_95ci": {
    "spearman_delta_score_vs_delta_pKd": [
      -0.20219641390743148,
      -0.017549489786270734
    ],
    "direction_accuracy_ties_half": [
      0.3654991147711051,
      0.4777496065745759
    ],
    "spearman_abs_delta_score_vs_abs_delta_pKd": [
      -0.12619585069276026,
      0.05186891862131836
    ]
  },
  "position_profile_test": {
    "eligible_length_position_cells": 71,
    "minimum_pairs_per_cell": 3,
    "position_profile_spearman": -0.13543037752301498,
    "position_shuffle_pvalue_two_sided": 0.22315536892621476,
    "permutations": 5000
  },
  "interpretation_limits": [
    "The C score is a sequence-property heuristic with no receptor input.",
    "This benchmark tests the native score, not the adaptive optimizer.",
    "The broad dataset includes peptide lengths 9-20; it is not an exact-21-aa-only result.",
    "Position profiling is exploratory and is not validation of ARISE/ALE hotspot selection.",
    "Pairs overlap in peptides and target groups; connected-component bootstrap partly addresses this dependence."
  ]
}
```

### `prototype3/results/p3_21aa_affinity_audit/native_c_score_21aa/metrics.json`
```json
{
  "benchmark": "PPIKB exact canonical 21-aa native sequence-score audit",
  "input_candidate_rows": 247,
  "rows_passing_exact21_kd_linear_no_explicit_modification_filters": 129,
  "unique_peptide_sequences_eligible": 67,
  "unique_peptides_scored_successfully": 67,
  "score_failures": {},
  "filter_rejection_counts": {
    "not_exact_canonical_21aa": 0,
    "not_numeric_uncensored_Kd": 81,
    "explicit_modification_or_not_linear": 37,
    "missing_target_sequence": 0
  },
  "measurement_contexts": 73,
  "contexts_with_at_least_3_distinct_scored_peptides_and_variable_values": 7,
  "within_context_spearman_mean": -0.48719876297629444,
  "within_context_spearman_median": -0.7999999999999999,
  "within_context_spearman_group_bootstrap_95ci_for_mean": [
    -0.7714285714285715,
    -0.030055905833437393
  ],
  "single_substitution_pairs_same_context": 20,
  "informative_single_substitution_pairs_abs_delta_pKd_ge_0_3": 9,
  "informative_mutation_direction_accuracy_ties_half": 0.1111111111111111,
  "informative_delta_score_vs_delta_pKd_spearman_exploratory": 0.6779661016949153,
  "interpretation_warning": "These are retrospective associations, not prospective prediction. The native score is sequence-only and has no receptor input. Pairs overlap in peptides and are not independent observations. The position summary describes this filtered dataset, not a validated biological hotspot selector. The grouped correlations are exploratory."
}
```

### `prototype3/results/arise_experimental_position_audit/metrics.json`
```json
{
  "benchmark": "Retrospective PPIKB position-window validation",
  "input": "/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/sequence_only_mutation_baseline/run_20261026T234823Z/single_substitution_pairs.csv",
  "n_single_substitution_pairs": 693,
  "n_informative_pairs_abs_delta_pKd_ge_0_3": 466,
  "hotspot_coordinate_mapping": "The logs specify zero-based start and end-exclusive indices; [0,3) maps to one-based positions 1,2,3; [1,3) maps to positions 2,3.",
  "primary_window_test": {
    "label": "latest_logged_hotspot_positions_1_to_3",
    "positions_1based": [
      1,
      2,
      3
    ],
    "informative_pairs_only": false,
    "pairs_total": 693,
    "pairs_in_window": 109,
    "pairs_outside_window": 584,
    "contexts_after_matching_filter": 86,
    "unique_target_sequences": 123,
    "raw_mean_abs_delta_pKd_in_window": 0.8582792040816429,
    "raw_mean_abs_delta_pKd_outside_window": 0.8192150178168482,
    "raw_median_abs_delta_pKd_in_window": 0.5201081171791433,
    "raw_median_abs_delta_pKd_outside_window": 0.5237776248695273,
    "centered_effect_difference_window_minus_elsewhere": -0.011501530384395456,
    "target_cluster_bootstrap_95ci": [
      -0.16901286911387303,
      0.15432571704219142
    ],
    "bootstrap_valid_replicates": 2000,
    "random_contiguous_window_permutation_p_greater": 0.5860827834433113,
    "random_contiguous_window_permutation_p_two_sided": 0.9218156368726255,
    "permutation_valid_replicates": 5000,
    "permutation_count_requested": 5000,
    "bootstrap_count_requested": 2000,
    "interpretation": "Positive centered difference means mutations in the selected window have larger absolute measured affinity effects than the within-context expectation. This is retrospective and exploratory."
  },
  "previous_generation_window_sensitivity": {
    "label": "previous_generation_hotspot_positions_2_to_3",
    "positions_1based": [
      2,
      3
    ],
    "informative_pairs_only": false,
    "pairs_total": 693,
    "pairs_in_window": 72,
    "pairs_outside_window": 621,
    "contexts_after_matching_filter": 86,
    "unique_target_sequences": 123,
    "raw_mean_abs_delta_pKd_in_window": 0.8558699062833472,
    "raw_mean_abs_delta_pKd_outside_window": 0.8218218524920089,
    "raw_median_abs_delta_pKd_in_window": 0.49121637042493926,
    "raw_median_abs_delta_pKd_outside_window": 0.5288841700741336,
    "centered_effect_difference_window_minus_elsewhere": -0.018426000846649134,
    "target_cluster_bootstrap_95ci": [
      -0.20258150000062086,
      0.1683673545532448
    ],
    "bootstrap_valid_replicates": 2000,
    "random_contiguous_window_permutation_p_greater": 0.6196760647870426,
    "random_contiguous_window_permutation_p_two_sided": 0.8454309138172366,
    "permutation_valid_replicates": 5000,
    "permutation_count_requested": 5000,
    "bootstrap_count_requested": 2000,
    "interpretation": "Positive centered difference means mutations in the selected window have larger absolute measured affinity effects than the within-context expectation. This is retrospective and exploratory."
  },
  "informative_subset_sensitivity": {
    "label": "latest_hotspot_positions_1_to_3_informative_pairs",
    "positions_1based": [
      1,
      2,
      3
    ],
    "informative_pairs_only": true,
    "pairs_total": 466,
    "pairs_in_window": 68,
    "pairs_outside_window": 398,
    "contexts_after_matching_filter": 71,
    "unique_target_sequences": 111,
    "raw_mean_abs_delta_pKd_in_window": 1.3042212755964766,
    "raw_mean_abs_delta_pKd_outside_window": 1.1423108474271246,
    "raw_median_abs_delta_pKd_in_window": 0.9463189699149308,
    "raw_median_abs_delta_pKd_outside_window": 0.7755395751650886,
    "centered_effect_difference_window_minus_elsewhere": 0.17159960568330337,
    "target_cluster_bootstrap_95ci": [
      0.01481250014849741,
      0.3418983806200788
    ],
    "bootstrap_valid_replicates": 2000,
    "random_contiguous_window_permutation_p_greater": 0.01999600079984003,
    "random_contiguous_window_permutation_p_two_sided": 0.052789442111577684,
    "permutation_valid_replicates": 5000,
    "permutation_count_requested": 5000,
    "bootstrap_count_requested": 2000,
    "interpretation": "Positive centered difference means mutations in the selected window have larger absolute measured affinity effects than the within-context expectation. This is retrospective and exploratory."
  },
  "limitations": [
    "The experimental dataset was explored in prior analyses; this is not an untouched prospective test.",
    "PPIKB peptides have different targets and sequences from the 21-aa search seed.",
    "This tests a general peptide-position window, not HERV-K Env-specific residue importance.",
    "It does not evaluate target-aware structural contacts or prove that ARISE guided the mutations in the logged smoke run.",
    "The permutation test compares the observed window with random contiguous windows of equal width, stratified by peptide length."
  ]
}
```

### `prototype3/results/sequence_only_grouped_target_validation/run_20261027T164512Z/metrics.json`
```json
{
  "benchmark": "experimental_single_substitution_pKd_difference_sequence_only_grouped_validation",
  "bootstrap_seed": 20261026,
  "claim_scope": "separate sequence-only comparator; not native c_score, not the adaptive optimizer, and not target-aware structural score",
  "collapsed_target_assay_peptide_records": 2496,
  "independent_leakage_components": 111,
  "input_csv": "/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/ppikb_kd_benchmark/filtered_kd_measurements.csv",
  "input_sha256": "48aafbf55007d258b07c6075820c7dd3a2b01491fa8a8e322789ab193792e418",
  "limitations": [
    "Retrospective benchmark on deposited affinity data; not prospective validation.",
    "Sequence-only comparator is separate from the production native c_score and is not integrated into candidate generation.",
    "Holding target sequences and exact peptide pairs together prevents those forms of split leakage, but does not establish generalization to HERV-K Env.",
    "Pair orientation is lexicographic; delta_pKd is pKd_b minus pKd_a and positive values mean peptide_b has stronger measured affinity.",
    "Experimental annotations may differ in construct, assay conditions, and uncertainty even after exact target-sequence grouping."
  ],
  "metrics": {
    "sequence_only_ridge": {
      "informative_direction_accuracy": 0.589622641509434,
      "informative_direction_accuracy_ties_half": 0.5815450643776824,
      "informative_direction_coverage": 0.9098712446351931,
      "mae": 0.8643157903240979,
      "n": 693,
      "n_informative": 466,
      "rmse": 1.395672869959655,
      "spearman_rho": 0.10874027339042665
    },
    "training_fold_mean_baseline": {
      "informative_direction_accuracy": 0.6266094420600858,
      "informative_direction_accuracy_ties_half": 0.6266094420600858,
      "informative_direction_coverage": 1.0,
      "mae": 0.818212554563056,
      "n": 693,
      "n_informative": 466,
      "rmse": 1.3655814586084296,
      "spearman_rho": 0.045870601293509276
    },
    "zero_change_baseline": {
      "informative_direction_accuracy": null,
      "informative_direction_accuracy_ties_half": 0.5,
      "informative_direction_coverage": 0.0,
      "mae": 0.825359312626174,
      "n": 693,
      "n_informative": 466,
      "rmse": 1.3840433664872516,
      "spearman_rho": null
    }
  },
  "model": {
    "alpha": 10.0,
    "feature_count": 453,
    "feature_definition": "directional amino-acid substitution, peptide length, normalized position, and flanking residues; sequence-only; no target/assay identifiers",
    "hyperparameter_tuning": "none; alpha fixed before CV",
    "name": "Ridge"
  },
  "n_folds": 5,
  "paired_bootstrap_ridge_minus_fold_mean": {
    "bootstrap_unit": "connected leakage component, resampled as a cluster",
    "cluster_count": 111,
    "direction_tie_half_difference_a_minus_b_ci95": {
      "ci95_high": 2.399232245668356e-05,
      "ci95_low": -0.09796700753735615
    },
    "mae_difference_a_minus_b_ci95": {
      "ci95_high": 0.07815702909143259,
      "ci95_low": 0.01731472836918097
    },
    "replicates_requested": 2000,
    "spearman_difference_a_minus_b_ci95": {
      "ci95_high": 0.1726406814692609,
      "ci95_low": -0.0636806757785614
    },
    "valid_replicates_direction_tie_half": 2000,
    "valid_replicates_mae": 2000,
    "valid_replicates_spearman": 2000
  },
  "result_directory": "/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/sequence_only_grouped_target_validation/run_20261027T164512Z",
  "single_substitution_pairs": 693,
  "split_leakage_checks": {
    "all_folds_have_zero_exact_pair_overlap": true,
    "all_folds_have_zero_peptide_sequence_overlap": true,
    "all_folds_have_zero_target_overlap": true,
    "all_rows_have_one_oof_prediction": true
  },
  "split_rule": "GroupKFold over connected components; every exact target sequence, individual peptide sequence, and unordered peptide pair is confined to one fold",
  "unique_exact_pair_signatures": 633,
  "unique_target_sequences": 123
}
```

## Interpretation boundary

The direct comparison is limited to properties measured for the generated candidates themselves, such as score distributions, sequence counts, and mutation behavior. External PPIKB, exact-21-aa, PEPBI, receptor-context, and structural benchmarks remain shared validation unless candidate-level matched measurements are available. Missing experimental labels or valid poses are reported as unavailable evidence, never as zero predictive performance.

Do not claim that native C-score predicts or optimizes binding to HERV-K Env 9SYA or ARG203/ARG205: the native score has no receptor input.
