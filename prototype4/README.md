# Silicon Virovore P4 — empirical validation backbone

P4 is an isolated validation layer. It does not alter the ARISE/ALE optimizer,
relabel its sequence score as affinity, or launch a docking program.

## What this implementation does

1. Reads the original SKEMPI 2.0 table and detects comma, semicolon, or tab delimiter
   from the actual file header.
2. Preserves every source row and every original column in a provenance database
   (`skempi_raw.raw_json`) with a row UID and source-file SHA-256.
3. Extracts only exact, uncensored four-state cycles from a single exact context:
   complex/PDB chain field, Protein 1, Protein 2, reference, temperature, and
   experimental method. Each cycle must have P0, Pi, Pj, and Pij measurements.
4. Rejects censored affinities, conflicting alternative wild-type baselines, and
   conflicting repeated mutant-state values rather than silently averaging them.
5. Stores all four node UIDs on every feature row, source row links on every node,
   and exact feature UIDs plus training-family lists/hashes on every prediction.
6. Runs a nested grouped validation: an outer leave-one-family-out split, inner
   group-held-out predictions for uncertainty calibration, and a fixed-alpha Ridge
   mutation-descriptor model compared with the training-fold-mean baseline.
7. Evaluates whether the calibrated expected absolute error predicts error and
   whether selective retention beats repeated random retention at matched coverage.
8. Imports already-completed structure-run manifests and reports elapsed wall time,
   CPU core-hours, and GPU-hours separately. It never starts a docking engine.

## What it does not do

- It does not determine peptide binding to HERV-K Env or establish that the
  Arg203/Arg205 site is a peptide-binding site.
- SKEMPI 2.0 is a protein–protein interaction mutation benchmark, not a
  peptide-specific test set. The default validation is therefore a general PPI
  mutation-cycle benchmark. A peptide-like subset is available only when exact
  partner-chain sequence lengths and partner identity have been independently
  verified and entered in the family map.
- It does not yet use CABS-dock, HADDOCK, RAPiDock, Rosetta, or another structural
  engine as a validated predictor. `structure_runs` is a provenance/cost ledger,
  not evidence that an engine is accurate.
- It does not calculate EVSI or automatically allocate simulations. EVSI must wait
  until a structural endpoint has passed independent pose-recovery and control
  benchmarks and actual marginal compute costs are measured. A heuristic score
  divided by runtime would not be EVSI.
- It does not use SKEMPI `Hold_out_type` as a substitute for sequence-derived
  receptor-family grouping. The user-supplied family map must be defined before
  examining held-out outcomes.

## Dependencies

From the repository root:

```bash
python -m pip install -r prototype4/requirements.txt
```

Python 3.10+ is required. The existing P3 scripts already use the same scientific
Python stack.

## Step 1 — audit the original SKEMPI file

Place the original, unmodified SKEMPI 2.0 table somewhere outside the generated
P4 results directory. Do not open/save it in a spreadsheet application before
hashing it.

```bash
python prototype4/p4.py audit \
  --input /absolute/path/to/skempi_v2.csv \
  --out prototype4/results/skempi_audit_20261029 \
  --overwrite
```

Inspect these outputs before doing anything else:

- `p4_audit_summary.json`: raw row count, columns, detected delimiter, parse statuses,
  cycle rejection counts, source SHA-256, and actual eligible cycle count.
- `eligible_cycles.csv`: one row per complete four-state cycle.
- `family_map_template.csv`: one row per complex requiring family annotation.
- `p4_evidence.sqlite`: normalized six-table provenance database.

The parser treats `Temperature` as kelvin, consistent with the SKEMPI 2.0
publication. Only positive exact numeric Kd values are eligible. A raw value such
as `>1e-8`, `<1e-8`, `no binding`, or a range is excluded even if a companion
parsed field happens to contain a number.

## Step 2 — assign receptor-family groups

Fill `family_id` in `family_map_template.csv` using the actual sequences of the
receptor-side chains in the PDB complexes. Complexes with homologous receptor
partners belong to the same family group. Preserve chain identity and construct
context. Record the sequence source, clustering method, identity/coverage
threshold, and manual decisions in `family_notes`.

Do not infer family membership from PDB code, protein name alone, SKEMPI
`Hold_out_type`, or `Hold_out_proteins`. If a complex contains a peptide partner,
`family_id` should normally group by the non-peptide receptor chain family so
related receptors cannot appear in both training and test folds.

For a general PPI test, only `family_id` is required. For the optional peptide-like
subset, also fill `protein1_length`, `protein2_length`, and `peptide_partner`
(`protein1`, `protein2`, or `none`) after checking the actual chain sequences in
the cited PDB. `--peptide-only` selects complexes where the indicated peptide
partner is 5–30 residues by default. This length window is an operational filter,
not a universal biological definition of a peptide.

## Step 3 — run nested family-held-out validation

```bash
python prototype4/p4.py validate \
  --db prototype4/results/skempi_audit_20261029/p4_evidence.sqlite \
  --family-map prototype4/results/skempi_audit_20261029/family_map_template.csv \
  --split-unit family \
  --out prototype4/results/p4_validation_20261029
```

For a verified peptide-partner subset only:

```bash
python prototype4/p4.py validate \
  --db prototype4/results/skempi_audit_20261029/p4_evidence.sqlite \
  --family-map prototype4/results/skempi_audit_20261029/family_map_template.csv \
  --split-unit family --peptide-only --max-peptide-length 30 \
  --out prototype4/results/p4_peptide_subset_20261029
```

The model's preprocessing and Ridge fit occur inside each training fold. The
uncertainty-to-error mapping is trained only on inner group-held-out predictions
from that outer training partition. The held-out family is not used to tune the
model, calibrator, or abstention threshold. Ridge alpha is fixed at 10.0; no
hyperparameter search is performed on outer test outcomes.

The code requires at least 20 eligible cycles and five independent groups. Fewer
than ten groups are explicitly labelled exploratory. Five groups are only a
minimum software/statistical gate, not a claim of adequate power. The family-level
bootstrap intervals are conditional on the already-generated out-of-fold
predictions; they do not include retraining variability.

Key outputs:

- `loco_oof_predictions.csv`: one truly outer-held-out prediction per eligible
  cycle, exact Feature_UID, held-out family, full training-family list, and hash.
- `outer_fold_manifest.csv`: outer split membership and training-group hash.
- `outer_group_metrics.csv`: results by held-out family.
- `uncertainty_calibration_bins.csv`: predicted versus observed absolute error.
- `risk_coverage.csv`: calibrated-error selective retention versus repeated random
  retention at matched coverage.
- `p4_validation_metrics.json`: model/baseline, uncertainty, calibration,
  risk–coverage, and conditional group-bootstrap metrics.

Interpretation rules:

- `epsilon_exp` is nonadditivity on the binding-free-energy scale. Nonzero
  epistasis does not prove a direct physical contact between the mutated sites.
- Lower MAE is better. `model_minus_baseline_mae < 0` favours the model.
- Lower risk–coverage AURC is better. A nonpositive improvement versus random
  means uncertainty-guided abstention did not demonstrate an advantage.
- A successful SKEMPI benchmark does not establish transfer to HERV-K or a
  designed 21-aa peptide. That requires a separately defined target-relevant test.

## Step 4 — import structural runs only after independent validation

This code does not call a docking engine. To import completed runs, provide a CSV
with at least these fields:

```text
cycle_uid,node_uid,engine,engine_version,command,config_hash,random_seed,input_structure_hash,output_artifact_hash,status,elapsed_seconds,cpu_core_seconds,gpu_seconds,l_rmsd_angstrom
```

`config_hash`, `input_structure_hash`, and (for successful runs)
`output_artifact_hash` must be 64-character SHA-256 hex digests. A failed run may
have a blank output hash, but its failure status remains in the ledger.

```bash
python prototype4/p4.py import-structures \
  --db prototype4/results/skempi_audit_20261029/p4_evidence.sqlite \
  --manifest /absolute/path/to/completed_structure_runs.csv

python prototype4/p4.py cost-report \
  --db prototype4/results/skempi_audit_20261029/p4_evidence.sqlite \
  --out prototype4/results/p4_structure_costs.json
```

The accounting keeps wall time, CPU core-hours, and GPU-hours separate; it never
adds unlike units. `L_RMSD` should be populated only for genuine native-pose
recovery benchmarks with a known experimental peptide pose. For novel HERV-K
candidates without a native peptide complex, leave it blank.

## Step 5 — run the tests

```bash
python -m unittest tests.test_p4_pipeline -v
pytest -q tests/test_grouped_sequence_pair_validation.py
```

These tests use synthetic data. They verify the software mechanics, not biological
accuracy, data feasibility, novelty, or HERV-K binding.
