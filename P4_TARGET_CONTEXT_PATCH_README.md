# AmphiMize P4 target-context validation patch

Apply from the Silicon-Virovore project root:

```bash
patch -p1 < AmphiMize_P4_TargetContext_ExternalValidation.patch
```

The patch adds a group-held-out P3-vs-P4 external-transfer experiment using PEPBI direct ITC ΔG, adds a fail-closed authorization gate to target-aware Pareto ranking, changes P4 live target scoring to a bounded shortlist, and adds tests for the gate and shortlist slot integrity.

The included fresh benchmark result is negative for target-context augmentation on this PEPBI subset. Thus this patch **does not enable an unvalidated target score**. It makes the implementation testable and prevents unsupported target features from changing candidate selection. HERV-K Env 9SYA target-specific efficacy remains unproven.

## Validation performed while preparing this patch

- External PEPBI run: 46 records, 10 binding groups, leave-one-Binding-Group-out.
- Focused P4 test command: 30 passed.
- `python -m py_compile` passed for modified production modules.
- `python run_pipeline.py --pilot --p4-target-optimization ...` failed at the intended preflight gate before candidate generation because the external validation did not pass.
- Full `pytest -q tests` was attempted but did not finish within the execution window; do not claim the complete suite passed.
