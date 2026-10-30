# Existing-candidate structural translation test

**Generated:** 2026-10-27T18:00:15-04:00
**Target reference:** `/Users/centurion616/Desktop/Silicon-Virovore/prototype2/target_aware_audit/9SYA_chainA_repaired.pdb`
**Target reference sequence loaded:** yes

## Question

Does the adaptive arm's advantage remain visible when existing candidate structures are evaluated with the same target-aware geometric criteria?

This is a structural-proxy test. It is not an experimental measurement of binding affinity and cannot by itself prove stronger real-world binding.

## Predefined criteria

1. Minimum receptor–peptide heavy-atom distance >= 2.0 Å.
2. Minimum peptide-to-ARG203 heavy-atom distance <= 5.0 Å.
3. Minimum peptide-to-ARG205 heavy-atom distance <= 5.0 Å.

A model passes only when it meets all three criteria. Model-level pass rates are reported descriptively; where cluster labels can be identified, the exploratory comparison uses clusters as the resampling/permutation units.

## Scan summary

- Arm-labeled PDB files/archive members discovered: 306.
- PDB models parsed: 306.
- Identical model duplicates skipped within the same arm/study/cluster: 0.
- Models eligible after receptor-sequence verification: 304.
- Model rows recorded, including ineligible structures: 306.

## Comparable study groups

| Existing study group | Baseline pass rate | Adaptive pass rate | Adaptive minus baseline cluster rate | 95% cluster-bootstrap CI | Cluster permutation p |
|---|---:|---:|---:|---:|---:|
| `prototype2/structural_validation/haddock_audit_20261024_184403/{ARM}_models` | 92.5% (37/40) | 86.1% (31/36) | -0.064 | [-0.228, 0.100] | 0.5726 |
| `prototype2/structural_validation/haddock_audit_20261024_184440/{ARM}_models` | 92.5% (37/40) | 86.1% (31/36) | -0.064 | [-0.228, 0.100] | 0.5726 |
| `prototype2/structural_validation/haddock_audit_fixed_20261024_191309/{ARM}_models` | 92.5% (37/40) | 86.1% (31/36) | -0.064 | [-0.228, 0.100] | 0.5726 |
| `prototype2/structural_validation/haddock_audit_retry_20261024_185024/{ARM}_models` | 92.5% (37/40) | 86.1% (31/36) | -0.064 | [-0.228, 0.100] | 0.5726 |

## Interpretation limits

- A higher geometric pass rate would support greater structural plausibility under these predefined criteria, not prove greater experimental affinity.
- A similar or lower pass rate would indicate that the internal-score advantage was not mirrored by this structural proxy.
- Docking models from the same cluster are related observations; model-level percentages must not be treated as independent biological replicates.
- The source path and peptide sequence should be checked against the optimizer's actual candidate inventory before claiming that the tested poses represent the relevant candidate set.
- The existing P3 paired optimizer run was previously reported to contain no candidate PDBs. If these structures come from a different P2 run, label the result as validation of that P2 structural set, not as a direct scoring of all 480 current P3 records per arm.

## Files

- `models_scored.csv`: per-model geometry measurements and eligibility.
- `study_group_summary.csv`: same-study adaptive-versus-baseline comparisons.
- `cluster_summary.csv`: per-cluster pass rates.
- `source_inventory.csv`: exact source paths and scan coverage.
- `structural_translation_pass_rates.png`: comparison graph, when comparable groups exist.
