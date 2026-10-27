# Silicon Virovore P3 — evidence status and claim boundaries

**Prepared:** 26 October 2026
**Purpose:** freeze a defensible interpretation of P3 and define the final reproducible grouped validation. This report is not a claim that every listed dataset was freshly re-run by the author of this document.

## Executive conclusion

P3 has a tested software pipeline and several distinct validation studies, but it does **not** yet establish reliable biological binding prediction. Keep four questions separate:

1. **Native `c_score`:** a sequence-level heuristic with no receptor input. The existing PPIKB Kd benchmark was evaluated within target/assay groups and yielded mean Spearman rho around -0.082 and median around -0.225 across 254 target/assay groups with at least three peptide entries; the group-bootstrap 95% interval for the mean was approximately [-0.153, -0.009]. This does not support affinity-prediction utility.
2. **Adaptive optimizer:** the existing 100-matched-seed comparison reported mean best-across-generations score 0.8338 for adaptive versus 0.7104 for baseline, paired difference +0.1234 (95% paired-bootstrap interval [+0.1189, +0.1279], 100/100 wins). This is an end-to-end sequence-objective result, not proof of better target binding or a clean mutation-guidance causal effect: adaptive preserves/carries forward the best parent, while the baseline restarts mutations from the original seed. Other matched runs (including seeds 616 and 618 in the later setup) favored the baseline and showed greater diversity. Report the setup-specific findings together; do not pool incompatible experiments or say adaptive always wins.
3. **Target-aware structural score:** this is a separate pose-geometry heuristic. The PPIKB experimental mutation analysis scored 18 structure-mapped variant pairs; rho of score change versus pKd change was about 0.059, with 50% directional agreement in the 18-pair analysis. The stricter same-observed-receptor subset had only 10 pairs across five target groups and a very wide bootstrap interval. In an independent-pose sanity check on the 1YCR MDM2/p53 complex, the incorrect docked pose (peptide C-alpha RMSD about 9.97 Å from the crystal peptide after receptor alignment) scored higher than the crystal pose. These results do **not** validate affinity or native-pose ranking.
4. **Separate sequence-only pair model:** previous out-of-fold analysis reported a modest ranking signal (Ridge Spearman rho about 0.243, bootstrap interval about [0.147, 0.337]) but worse point-estimate MAE than the fold-mean baseline; directional improvement was uncertain. This model is not the production `c_score`, is not used by the adaptive optimizer, and is not the structural score. Previous folds must not be described as held-out-target validation unless their split audit proves it.

## Added final audit: target/pair-grouped sequence-only validation

Run `prototype3/scripts/run_grouped_sequence_pair_validation.py` after `prototype3/results/ppikb_kd_benchmark/filtered_kd_measurements.csv` exists. The script:

- collapses measurements by exact target sequence, assay, and peptide before pairing;
- builds pairs only within the same target sequence and assay, same peptide length, and exactly one canonical amino-acid substitution;
- uses a separate fixed-alpha Ridge comparator with only substitution identity, peptide length, mutation position, and flanking-residue features;
- uses GroupKFold over connected leakage components, where rows are joined if they share an exact target sequence, an individual peptide sequence, or an exact unordered peptide-pair signature;
- asserts zero exact target-sequence overlap, individual peptide-sequence overlap, and pair-signature overlap between training and test folds;
- outputs one held-out prediction per pair, a per-fold split audit, split manifest, zero-change baseline, training-fold-mean baseline, Ridge metrics, and cluster-bootstrap paired differences;
- fails rather than silently falling back to a row-random split if there are too few independent components.

This is a stronger *retrospective grouped generalization check* than an ordinary row split. It is **not** a completely untouched prospective test; prior exploratory analyses have already examined this data collection. Exact sequence grouping also does not guarantee that all homologous target families are separated, so do not call it homology-cluster holdout.

### Run command

```bash
python prototype3/scripts/run_grouped_sequence_pair_validation.py \
  --input prototype3/results/ppikb_kd_benchmark/filtered_kd_measurements.csv \
  --out-root prototype3/results/sequence_only_grouped_target_validation \
  --folds 5 --bootstrap-replicates 2000 --seed 20261026
```

Then run:

```bash
python -m pytest -q tests
git diff --check
```

Review `REPORT.md`, `metrics.json`, `fold_split_audit.json`, and `out_of_fold_predictions.csv` in the newly printed run directory. The new run's metrics are not yet embedded here because this patch was tested only with synthetic data in the packaging environment; it must be executed against the actual local PPIKB output before its numerical result is reported.

## Evidence that must remain explicit

- **143 tests passed in 33.56 s** was the user's earlier test result before this patch; this document does not claim a fresh run of the complete repository suite. The three new tests passed on a synthetic dataset in the packaging environment.
- Structural sample counts remain 18 broad mapped pairs and 10 in the strict subset. Do not inflate these counts by accepting fuzzy chain matches without separately validating construct identity and residue mapping.
- The target-aware structural score is not integrated into the candidate optimizer and does not change native `c_score`.
- Docking, contact counts, interface size, and far-separated pose controls are geometric/computational evidence, not experimental binding measurements.
- No claims of biological efficacy, therapeutic utility, or HERV-K-specific predictive accuracy are justified by current data.

## Recommended research claim

> Silicon Virovore is a modular computational peptide-search pipeline whose search behavior and component scores were evaluated against matched sequence-search baselines and retrospective experimental affinity data. A separate sequence-only mutation-pair comparator showed a preliminary ranking association under out-of-fold evaluation; the legacy sequence heuristic did not demonstrate useful within-target affinity association, and the current target-aware structural geometry score did not reliably predict measured affinity changes or consistently rank a native pose above an incorrect pose. The work therefore identifies both software capabilities and specific validation limits that must be addressed before biological predictive utility can be claimed.

This statement is intentionally narrower than “the system reliably predicts binding.” The latter claim is unsupported.
