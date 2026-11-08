# P4 target-context external validation

## What this experiment tests

This experiment asks whether adding a complex-specific structural interaction feature improves prediction of a directly reported experimental peptide–protein binding free energy, compared with the same sequence-only score on the exact same records.

- **P3 comparator:** native sequence-only C fitness score.
- **P4 target-context model:** the same C fitness score plus the matched complex's Rosetta Interface Analyzer `dG_separated` feature.
- **Sensitivity model:** P4 primary features plus `dSASA_int` and `nres_int`.
- **Experimental endpoint:** directly reported PEPBI ITC `ΔG (kcal/mol)`.
- **Split:** leave-one-`Binding Group`-out; all records from a held-out group are excluded from training. Ridge scaling and fitting are done within each training fold.
- **Inferential unit:** binding group, not individual measurement rows.

This is an external transfer experiment using PEPBI. It is not a direct experiment on HERV-K Env 9SYA, and it does not measure HERV-K binding or inhibitor efficacy.

## Fresh result

The filtered complete-case 9–21 residue dataset contained 46 records, 12 unique peptide sequences, and 10 binding groups.

| Model | Macro binding-group MAE (kcal/mol) | Pooled MAE (kcal/mol) | Pooled Spearman ρ |
|---|---:|---:|---:|
| P3 sequence-only | 0.834208 | 1.094243 | -0.3745 |
| P4 sequence + target-context RIA | 0.887456 | 1.257150 | -0.4160 |
| P4 geometry sensitivity | 0.937270 | 1.265934 | -0.4395 |
| Training-group mean | 0.814605 | 0.944033 | -0.5059 |

Paired primary endpoint is `P3 group MAE - P4 group MAE`; positive values favor P4. Observed macro benefit was **−0.053247 kcal/mol**. The group-bootstrap 95% CI was **[−0.122562, −0.009319]**. P4 improved 3/10 groups and worsened 7/10. Exact paired sign-flip p was 0.9893 for the prespecified one-sided alternative that P4 improves, 0.0117 for the opposite direction, and 0.0234 two-sided.

The pooled raw feature correlations were descriptive only: native sequence score vs experimental ΔG, Spearman ρ = −0.2427; RIA `dG_separated` vs experimental ΔG, ρ = −0.0246. Inferential claims use held-out predictions and group-level paired errors, not these pooled correlations.

## Interpretation

On this PEPBI subset, target-context RIA augmentation did **not** improve out-of-group agreement with experimental measurements. The geometry sensitivity model did worse still. This is a negative transfer result, not evidence that all structural evidence is useless. Ten groups and 46 rows remain a small external panel, and these results do not directly adjudicate a 9SYA-specific feature set.

## Target-selection safety gate

P4 target-aware Pareto ordering is now fail-closed behind this external-validation gate. It is authorized only when all predeclared conditions are met:

1. At least 10 independent held-out binding groups.
2. Positive P3-minus-P4 macro-MAE benefit.
3. Group-bootstrap 95% CI lower bound above zero.
4. One-sided exact paired sign-flip p < 0.05 for P4 improvement.
5. More binding groups improve than worsen.

The current result fails four of these five conditions (group count is exactly 10 and passes). Therefore the P4 target-aware parent-ordering mode is blocked before candidate generation. The structural shortlist and diagnostic evidence can still be generated in diagnostic mode, but the target geometry must not alter evolutionary ranking until an external transfer experiment passes.

When authorized by a qualifying future report, P4 target optimization is budget-limited to a deterministic, diverse shortlist; it reorders only original shortlist slots. Every candidate in that shortlist must have exact-sequence 9SYA structural evidence and at least two distinct-seed docking runs under a shared protocol hash. Reliability remains an eligibility check, not a ranking weight. Neither geometry nor docking is labeled affinity.

## Reproduction

From the repository root, with the project's Python environment and the C native library built for the host platform:

```bash
python prototype4/target_context_external_validation.py \
  --input prototype4/data/pepbi_target_context_external_validation.csv \
  --out prototype4/results/target_context_external_validation_fresh \
  --alpha 1.0 --bootstrap-reps 20000 --seed 20261103

pytest -q tests/test_p4_pipeline_integration.py \
  tests/test_p4_target_shortlist.py \
  tests/test_p4_external_validation_gate.py \
  tests/test_p4_target_optimization_config.py \
  tests/test_p4_selection_reporting.py
```

To attempt P4-TA after a future external-validation pass, generate a fresh report first, then run with the explicit flag. The pipeline will stop before candidate generation if the external gate fails:

```bash
python run_pipeline.py --p4-target-optimization \
  --population-size 100 --generations 5 --seed 631 \
  --output-dir results/p4_target_optimization
```

The current checked-in PEPBI result intentionally blocks this command. Do not replace the summary with a fabricated or unrelated positive result.

## Output files

- `prototype4/data/pepbi_target_context_external_validation.csv`: exact analysis cohort.
- `prototype4/results/target_context_external_validation_fresh/summary.json`: metrics, split protocol, and gate decision.
- `out_of_group_predictions.csv`: every held-out prediction and source row.
- `binding_group_metrics.csv`: group-level paired errors.
- `paired_group_effects.csv`: primary group-level contrasts.
- `model_summary.csv`: pooled and macro metrics.
- `external_validation_gate.json`: selection authorization decision.
- `target_context_p3_vs_p4_group_mae.png`: poster figure.

Source workbook: PEPBI, Torres et al. (2025), *Scientific Data*, 12, 1455, doi:10.1038/s41597-025-05754-7. Keep the data/license attribution required by the source dataset.
