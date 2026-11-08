# AmphiMize P4 three-milestone implementation patch

## What this patch does

This patch makes the remaining P4 work executable and auditable. It does **not** manufacture successful structures, prove binding, or certify target optimization.

### 1. Structural coverage

- Generates candidate-specific 9SYA chain-A / peptide chain-B HADDOCK evaluations for the configured, bounded P4 shortlist.
- Evaluates each selected candidate with two separate HADDOCK processes using distinct docking seeds.
- Records failures per candidate so one failed structure does not silently terminate the remaining shortlist.
- Writes execution manifests with exact input, protocol/configuration, output archive, and verified-pose hashes.
- Reuses an existing execution only when its manifest, return code, pose hashes, inputs, protocol, and archive hash verify.
- Adds `scripts/p4_target_coverage_audit.py` to report per-generation candidate coverage and repeatability.

**Important:** coverage is for the explicitly configured structural shortlist (default maximum 4 in the pilot plan), not automatically all 100 population members. The matched study compares equal shortlist sizes and equal structural budgets in both arms. P4-TA reorders only within the shortlist; candidates outside it retain their positions. If target evidence is missing or incomplete, P4 abstains rather than claiming a successful ranking.

### 2. Repeatability

The integration gate now requires at least two independent run IDs, distinct seeds, a shared protocol hash, at least two poses per run, and distinct execution IDs backed by successful execution manifests. A directory name or repeated copy of a pose is not considered an independent execution.

### 3. Matched P3 versus P4-TA study

`scripts/p4_matched_experiment.py plan` creates a deterministic initial FASTA pool for each GA seed and uses that exact pool for both P3 and P4-TA. The generated launcher holds constant the GA seeds, population size, generation count, mutation-guidance mode, candidate shortlist size, HADDOCK seed pair, and structural-evaluation budget.

- **P3:** receives structural scoring as a blinded evaluation/control; target-aware Pareto parent ordering stays off.
- **P4-TA:** can use target-aware Pareto ordering only with complete exact-sequence structural evidence and strict repeatability evidence.
- Both arms receive the same candidate-specific structural evaluation budget.
- The analyzer excludes incomplete or unmatched pairs; it does not fill in missing scores.

The main endpoint is a computational target-geometry measure in the evaluated shortlist. It is **not** an independent experimental binding measurement and is not a direct affinity/efficacy endpoint.

## Run it on the project Mac

Use the project's own `.venv` and native macOS `c/libsafety.so`. Do not replace that library with a Linux build.

### Small pilot (defaults)

From the repository root:

```bash
source .venv/bin/activate
python scripts/p4_matched_experiment.py plan --experiment-dir results/p4_matched_pilot
bash results/p4_matched_pilot/launch_matched_p3_p4ta.sh
python scripts/p4_matched_experiment.py analyze --experiment-dir results/p4_matched_pilot
```

Defaults: GA seeds 631 and 632, population 100, 3 generations, shortlist 4, and independent HADDOCK seeds 1917 and 2917. That is a pilot, not a powered final comparison. It schedules up to 96 separate HADDOCK executions across both arms if all candidates reach structural evaluation (2 GA seeds × 2 arms × 3 generations × 4 candidates × 2 HADDOCK seeds). Actual execution count can be lower when a candidate fails earlier; failures are recorded.

The launcher refuses to overwrite an existing experiment directory. To rerun, choose a new `--experiment-dir` rather than deleting evidence.

### Larger preregistered design

After the pilot pipeline is confirmed operational, create a larger design explicitly:

```bash
python scripts/p4_matched_experiment.py plan \
  --experiment-dir results/p4_matched_full \
  --seeds 631-640 \
  --population-size 100 \
  --generations 5 \
  --shortlist-size 12 \
  --dock-seeds 1917,2917
```

This can schedule up to 2,400 HADDOCK executions. Do not launch it until the small pilot confirms that HADDOCK, peptide structure generation, disk capacity, and evidence reporting all work on the Mac.

### Coverage audit for an individual run

```bash
python scripts/p4_target_coverage_audit.py --run-dir results/p4_matched_pilot/runs/P4_TA_seed631
```

The matched analyzer writes `matched_run_metrics.csv`, `matched_generation_metrics.csv`, `paired_seed_effects.csv`, and `matched_experiment_summary.json`.

## Validation and interpretation

The PEPBI external-transfer comparison in the prior patch did not show improvement from target-context augmentation. This patch therefore adds an explicit `--p4-exploratory-target-optimization` option for a **research pilot only**. It bypasses the failed PEPBI authorization gate while retaining exact-sequence structural and repeatability requirements. Generated reports must remain marked `EXPLORATORY_NOT_VALIDATED`.

A favorable matched geometry result would establish, at most, that this computational objective changed the geometry metric under this experiment. It would not demonstrate HERV-K Env binding affinity, inhibition, safety, or therapeutic efficacy. Direct independent HERV-K experimental measurements are not present in this computational project.

## Test status at patch construction

The focused P4 integration, configuration, gate, shortlist, reporting, generation-audit, pose-source-group, matched-experiment, and execution-manifest tests passed (45 tests). The full repository suite did not finish within the available execution window, so this patch does not claim a full-suite pass. No real HADDOCK matched experiment was run while constructing this patch.
