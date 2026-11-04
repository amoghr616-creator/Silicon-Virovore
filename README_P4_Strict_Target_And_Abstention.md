# AmphiMize strict P4 target mode and family-matched abstention analysis

## What this change does

1. Adds an explicit, default-off `--p4-target-optimization` mode.
2. In that mode, bypasses the old structure shortlist and requests peptide structures and exact-sequence HADDOCK complexes for the **entire generation**.
3. Requires at least two distinct HADDOCK seeds, at least two sequence-verified poses per seed, a common protocol SHA-256, and complete candidate-level target geometry before target-aware parent reordering is allowed.
4. Fails closed: if any candidate is missing required evidence, target-aware parent ordering is withheld and the generation stops with the evidence-report path. A shortlist-only run is still a diagnostic, not P4-TA.
5. Preserves live candidate-specific poses ahead of historical manifest fallback.
6. Adds a same-generation/same-candidate-object-cohort selection-effect report. It describes whether target ordering changed the winner; it is not affinity validation.
7. Adds an external-benchmark abstention analysis that compares selective retention with random retention **within each held-out family**, preserving the per-family retained counts and bootstrapping complete families.

## Scientific boundary

The cross-run contact-fingerprint metric is an **uncalibrated reproducibility diagnostic**. It is an eligibility gate only; it is not a probability of pose accuracy, binding, or inhibition, and it is not a ranking weight. A candidate passing the gate is only supported by repeatable computational geometry under the recorded protocol.

The selective-risk script evaluates family-held-out SKEMPI-derived protein-protein mutation-cycle predictions. It does not validate HERV-K Env peptide binding, validate the 9SYA pocket hypothesis, or establish therapeutic activity. Report the result only after it has been run on genuinely out-of-fold `split_unit=family` predictions.

## Apply and verify

From the Silicon-Virovore repository root, first inspect the working tree and back it up/commit any unrelated changes. Apply this rebased patch only after its check succeeds:

```bash
git status --short
git apply --check /path/to/AmphiMize_P4_Strict_Target_And_Abstention_Rebase.patch
git apply /path/to/AmphiMize_P4_Strict_Target_And_Abstention_Rebase.patch
python -m py_compile run_pipeline.py src/config.py src/p4_pipeline_integration.py prototype4/selective_abstention.py
python -m pytest -q \
  tests/test_p4_pipeline_integration.py \
  tests/test_p4_pipeline.py \
  tests/test_p4_live_gate.py \
  tests/test_p4_integrated_gate.py \
  tests/test_p4_generation_audit.py \
  tests/test_p4_pose_source_groups.py \
  tests/test_selective_abstention.py \
  -W error::RuntimeWarning
```

Then run the full root suite separately:

```bash
python -m pytest -q tests -W error::RuntimeWarning
```

Do not recursively collect the separate `hybridock-pep` subproject as part of the root suite.

## P4 smoke-test protocol

Use a small generation first; do not start the 100-candidate × 5-generation comparison until the small run has shown all-candidate coverage and independent-run pose counts.

```bash
export SILICON_VIROVORE_P4_DOCK_SEEDS=1917,2917
export SILICON_VIROVORE_P4_TARGET_PARETO=1
export SILICON_VIROVORE_P4_DOCK_NCORES=2

python run_pipeline.py --help
python run_pipeline.py \
  --p4-target-optimization \
  --population-size 4 \
  --generations 1 \
  --seed 631 \
  --run-id p4ta_full_coverage_smoke_631
```

This is a four-candidate software/protocol smoke test, not the preregistered P1–P4 comparison. Inspect the generation report before scaling. The gate should only reorder if all four candidates have exact sequence-matched target complexes, at least two scored poses in each of two distinct-seed runs, and a common protocol hash. Otherwise, expected behavior is an explicit abstention / stop.

Only after the smoke test succeeds should the full treatment be considered. For the planned study, every one of 100 candidates in each generation needs replicated target evidence. Five generations therefore entail as many as 500 candidate structures and 1,000 independent HADDOCK executions per P4-TA run (two seeds per candidate), before retries. Do not describe a run with a withheld Pareto gate as P4-TA.

## Family-matched abstention experiment

After completing family-held-out validation with `prototype4/p4.py validate`, run the new analysis on its out-of-fold output:

```bash
python prototype4/selective_abstention.py \
  --predictions prototype4/results/YOUR_LOFO_RUN/loco_oof_predictions_with_errors.csv \
  --out prototype4/results/YOUR_LOFO_RUN/selective_abstention_family_matched \
  --random-repetitions 1000 \
  --family-bootstrap-repetitions 10000 \
  --seed 20261103
```

The selective arm retains the lowest calibrated expected-error predictions within each family. The random control retains the same number of records within each family. The reported risk is macro-family MAE; the 95% interval resamples whole families, not individual cycles. Positive `random_minus_selective_mae` means selective retention has lower error than the matched random control. It is evidence of abstention utility on this benchmark only if its family-bootstrap interval and randomization comparison support that conclusion.

## Status/reporting rules

- Passing regression tests is software verification, not evidence of biological validity.
- `p4_selection_comparison` is descriptive selection-effect telemetry, not a validation endpoint.
- A strict P4-TA execution that aborts on incomplete evidence is a successful abstention behavior, **not** a completed P4-TA optimization experiment.
- The matched 40-execution P1/P2/P3/P4-TA comparison remains pending until all 40 slots have complete, valid execution records and P4-TA evidence gates pass.
