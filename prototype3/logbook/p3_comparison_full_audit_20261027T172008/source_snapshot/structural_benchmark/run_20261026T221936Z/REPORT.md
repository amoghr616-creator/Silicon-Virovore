# Matched single-substitution experimental benchmark

- Input: `prototype3/results/ppikb_kd_benchmark/filtered_kd_measurements.csv`
- Rows in input: 3439
- Eligible measurement rows: 3439
- Unique peptide/target/assay records: 2496
- Matched single-substitution pairs: 693
- Target-assay groups represented: 147
- Independent target sequences represented: 123
- Pairs where both entries have PDB IDs: 195
- Sequence-score column: not found

## Directional score test

null

## Limitations

- Pairs can share sequences; the pooled pair count is not the independent sample size.
- PDB identifiers need separate structure and chain verification before structural scoring.
- The 0.3 pKd effect threshold is an analysis rule, not a significance test.
- This retrospective test uses experimental affinity records; it does not demonstrate prospective prediction or biological efficacy.

## Reproducibility

`run_config.json` records the eligibility rules, threshold, score column, and random seed.
`benchmark_metrics.json` records counts and metrics; CSV outputs preserve the individual pairs.
