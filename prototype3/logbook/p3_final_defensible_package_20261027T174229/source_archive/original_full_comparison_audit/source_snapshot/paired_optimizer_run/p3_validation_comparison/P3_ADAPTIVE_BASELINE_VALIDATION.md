# P3 adaptive-versus-baseline validation

Optimizer runs: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/p3_paired_arm_20261027T022810Z`
Report directory: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/p3_paired_arm_20261027T022810Z/p3_validation_comparison`

## 1. Optimizer outputs

| Arm | Seed | Candidates read | Unique sequences | Mean native fitness | Best | Final-generation mean |
|---|---:|---:|---:|---:|---:|---:|
| adaptive | 2601 | 160 | 151 | 0.7593303900985404 | 0.8432968734829264 | 0.8051635273010702 |
| adaptive | 2602 | 160 | 151 | 0.7531530875283788 | 0.8455310975436344 | 0.8044782434377395 |
| adaptive | 2603 | 160 | 151 | 0.7354922050586472 | 0.8332672856265247 | 0.7919967888833692 |
| baseline | 2601 | 160 | 138 | 0.6510945074517631 | 0.7055481253563461 | 0.6410184150422924 |
| baseline | 2602 | 160 | 135 | 0.6526790917715971 | 0.7145354444591977 | 0.6484208881823028 |
| baseline | 2603 | 160 | 125 | 0.6529478635783754 | 0.7071985064014479 | 0.6569524982230924 |

## 2. Adaptive-selected mutation positions

Source: `adaptive_hotspot_provenance.csv`. Position-window tests use the same adaptive-selected positions for the matched baseline. If positions were not logged, the test is marked blocked; observed mutation frequency is not relabeled as ARISE output.

- Seed 2601: no explicit selected hotspot recovered; positions=; source=.
- Seed 2602: no explicit selected hotspot recovered; positions=; source=.
- Seed 2603: no explicit selected hotspot recovered; positions=; source=.

Mutation statistics: `mutation_position_validation.csv`.

## 3. Native score versus experimental affinity

| Arm | Exact candidate/pKd matches | Exact-match rho | External pairs at observed mutation positions | External rho | External direction accuracy |
|---|---:|---:|---:|---:|---:|
| adaptive | 0 | None | 0 | None | None |
| baseline | 0 | None | 0 | None | None |

Exact candidate matches are direct but may be too sparse to estimate accuracy. The external pair analysis filters a retrospective benchmark by positions observed in each arm; it is not experimental measurement of newly generated candidates.

## 4. PEPBI

- adaptive: 0 exact candidate-sequence matches out of 450 unique sequences.
- baseline: 0 exact candidate-sequence matches out of 275 unique sequences.

Detailed matches: `pepbi_exact_candidate_overlaps.csv`. An overlap count alone is not evidence of predictive accuracy; measured outcomes must be available for the same sequence.

## 5. New target-aware structural scoring

- Candidate-run PDBs inventoried: 0
- Structural scoring attempts: 0
- Eligible attempt metric files: 0

Files: `structural_complex_inventory.csv`, `new_structural_scoring_attempts.csv`. A generated candidate without a valid receptor–peptide complex pose cannot be scored structurally.

## 6. Grouped Ridge versus training-fold mean

- Fresh metrics: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/p3_paired_arm_20261027T022810Z/grouped_sequence_validation/run_20261027T023029Z/metrics.json`
- This is a separate sequence-only prediction benchmark, not a direct optimizer-arm comparison.

## Interpretation

These analyses distinguish internal optimization, retrospective position associations, exact-sequence affinity evidence, PEPBI overlap, and structural pose evidence. No experimental affinity is imputed for an unmeasured candidate. A higher native fitness score must not be described as stronger real-world binding without experimental validation.
