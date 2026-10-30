# P3 baseline-versus-adaptive: six-area comparison

Generated: 2026-10-27T16:45:18.625886+00:00

## 1. What this report compares

The three-seed paired optimizer experiment is a direct arm comparison on internal
scores and sequence diversity. External experimental benchmarks are separate:
they validate scoring methods on measured cases, but do not automatically validate
the candidates created by either optimizer arm.

## 2. Direct optimizer comparison: matched three-seed run

| Seed | Baseline mean native score | Adaptive mean native score | Difference | Baseline best | Adaptive best | Baseline final-gen mean | Adaptive final-gen mean | Baseline unique sequences | Adaptive unique sequences |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2601 | 0.65109 | 0.75933 | +0.10824 | 0.70555 | 0.84330 | 0.64102 | 0.80516 | 138 | 151 |
| 2602 | 0.65268 | 0.75315 | +0.10047 | 0.71454 | 0.84553 | 0.64842 | 0.80448 | 135 | 151 |
| 2603 | 0.65295 | 0.73549 | +0.08254 | 0.70720 | 0.83327 | 0.65695 | 0.79200 | 125 | 151 |

Across the three seeds, the prior paired-run summary reported mean final-generation
native-score difference **+0.15175**, mean best-score difference **+0.1316**, and
mean unique-sequence difference **+18.33 per run**, all adaptive minus baseline.
This is preliminary evidence for higher internal scores/diversity in these three runs;
it is not evidence of higher experimental affinity.

## 3. Requested validation areas: arm-by-arm status

| Area | Baseline arm | Adaptive arm | Shared external evidence | Defensible conclusion |
|---|---|---|---|---|
| Ridge versus training-fold mean | No experimental labels for generated baseline candidates; arm-specific MAE/RMSE/correlation unavailable. | No experimental labels for generated adaptive candidates; arm-specific MAE/RMSE/correlation unavailable. | Ridge and training-fold mean tested on the same strict grouped held-out pairs; zero target, exact-pair, or peptide-sequence overlap. | This compares prediction methods on the external dataset, not optimizer arms. Do not duplicate the same test metrics as two arm scores. |
| PEPBI | Prior exact-sequence audit: 0 exact PEPBI matches among 275 unique baseline candidate sequences. | Prior exact-sequence audit: 0 exact PEPBI matches among 450 unique adaptive candidate sequences. | Strict matched PEPBI benchmark: native-score group concordance 0.3762 (11 groups); Rosetta 0.3927 (11 groups), target-calculated ΔG from Kd. Target-reported ΔG analysis has only 1 independent group. | Score-method comparison is possible on PEPBI. Candidate-arm experimental accuracy is not estimable without measured labels for those sequences. |
| Receptor-context control | The native sequence-only score has no receptor input; it cannot change with receptor context for a fixed baseline peptide. | The same limitation applies to every fixed adaptive peptide sequence. | Prior control: native score unchanged in 24 matched receptor-variant pairs; 7 pairs had experimental ΔG differences >0.5 kcal/mol. | This demonstrates receptor-context blindness of the native score, not successful receptor-aware prediction. Candidate-level structural comparison remains unavailable without valid complex poses. |
| Mutation-position validation | 12.08% of baseline mutation events fall at proxy positions [14, 3, 8]. | 19.61% of adaptive mutation events fall at the same proxy positions [14, 3, 8]. | Shared external context-adjusted effect for these positions: -0.08283; random-position permutation p=0.54948. Separate all-pairs window test for positions 1–3: centered effect -0.01150, two-sided p=0.92182. | Arm mutation frequencies differ, but the external test did not validate the proxy positions. The proxies must not be described as verified ARISE-selected hotspots. |
| Exact 21-aa native C-score audit | No exact candidate-to-experimental pKd matches in the previous paired-arm audit; per-arm measured-affinity accuracy unavailable. | No exact candidate-to-experimental pKd matches in the previous paired-arm audit; per-arm measured-affinity accuracy unavailable. | 129 eligible measurement rows; 67 distinct peptides; 73 measurement contexts; within-context mean Spearman -0.4872. Only 9 informative single-substitution pairs; direction accuracy 1/9 = 11.11%. | Shared retrospective score-validity test, not an adaptive-versus-baseline accuracy test. |
| Broad native C-score mutation benchmark | No exact candidate-to-experimental pKd matches in the previous paired-arm audit; per-arm measured-affinity accuracy unavailable. | No exact candidate-to-experimental pKd matches in the previous paired-arm audit; per-arm measured-affinity accuracy unavailable. | 693 matched pairs; Spearman Δscore versus ΔpKd -0.10714 (95% CI [-0.20220, -0.01755]); informative-pair direction accuracy 42.49% (95% CI [36.55%, 47.77%]). | The shared benchmark argues against treating native C-score as a reliable experimental-affinity predictor. It does not establish which optimizer arm is more accurate experimentally. |
| Structural scoring of generated candidates | No candidate-run PDBs available in the prior audit; candidate structural score comparison unavailable. | No candidate-run PDBs available in the prior audit; candidate structural score comparison unavailable. | Prior structural benchmark scored 18 verified mapped pairs; Δscore versus ΔpKd rho 0.05882, p=0.81665; direction accuracy 50%. | Do not report candidate structural accuracy until both arms have valid, comparable poses and the same scoring protocol. |

## 4. Fresh metrics extracted from rerun outputs

- Broad native C-score JSON: `prototype3/results/native_score_all_mutation_pairs/metrics.json`
  - Pairs: 693; unique sequences: 621.
  - Spearman Δscore versus ΔpKd: -0.1071.
  - Direction accuracy: 42.49%.
  - Correlation 95% CI: [-0.20219641390743148, -0.017549489786270734].
- Exact 21-aa JSON: `prototype3/results/p3_21aa_affinity_audit/native_c_score_21aa/metrics.json`
  - Eligible rows: 129; distinct eligible peptides: 67.
  - Within-context mean Spearman: -0.4872; informative pairs: 9.
  - Informative-pair direction accuracy: 11.11%.
- Position-validation JSON: `prototype3/results/arise_experimental_position_audit/metrics.json`
  - Positions tested: [1, 2, 3]; centered effect: -0.0115; two-sided permutation p: 0.9218.
- Grouped Ridge JSON: `prototype3/results/sequence_only_grouped_target_validation/run_20261027T164512Z/metrics.json`
  - Fold-mean MAE / RMSE / rho: 0.8182 / 1.3656 / 0.0459.
  - Ridge MAE / RMSE / rho: 0.8643 / 1.3957 / 0.1087.

## 5. Rerun status

- broad_native_c_score: **PASS** (`prototype3/results/p3_six_area_final_20261027T164246Z/broad_native_c_score.log`)
- exact_21aa: **PASS** (`prototype3/results/p3_six_area_final_20261027T164246Z/exact_21aa.log`)
- mutation_position: **PASS** (`prototype3/results/p3_six_area_final_20261027T164246Z/mutation_position.log`)
- ridge_vs_fold_mean: **PASS** (`prototype3/results/p3_six_area_final_20261027T164246Z/ridge_vs_fold_mean.log`)
- pepbi_strict_matched: **PASS** (`prototype3/results/p3_six_area_final_20261027T164246Z/pepbi_strict_matched.log`)
- receptor_context_control: **PASS** (`prototype3/results/p3_six_area_final_20261027T164246Z/receptor_context_control.log`)

## 6. Important scientific boundary

The external PPIKB, exact-21-aa, PEPBI, receptor-control, and structural benchmarks
are shared tests of scoring methods or position hypotheses. They cannot produce a
valid adaptive-versus-baseline experimental-accuracy comparison when the generated
candidate sequences have no matching measured affinity labels. The candidate-level
exact-match counts and missing structural poses are therefore reported as missing
evidence, not as zero predictive performance.

Do not claim this native C-score optimizes or predicts binding to HERV-K Env 9SYA
or ARG203/ARG205. It is a sequence-level score without receptor input.

