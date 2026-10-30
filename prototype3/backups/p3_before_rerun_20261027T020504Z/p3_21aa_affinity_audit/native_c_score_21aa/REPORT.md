# Native sequence-score audit on PPIKB 21-aa peptides

- Candidate rows examined: **247**
- Rows passing strict filters: **129**
- Distinct scored sequences: **67**
- Independent target/assay/publication contexts: **73**
- Contexts with valid within-context Spearman correlations: **7**
- Mean within-context Spearman rho: **-0.48719876297629444**
- Median within-context Spearman rho: **-0.7999999999999999**
- Group-bootstrap 95% interval for mean rho: **[-0.7714285714285715, -0.030055905833437393]**
- Same-context single-substitution pairs: **20**
- Informative single-substitution pairs (|delta pKd| >= 0.3): **9**
- Mutation direction accuracy, ties counted as half: **0.1111111111111111**
- Exploratory Spearman of score change vs experimental pKd change: **0.6779661016949153**

## Filter counts
- not_exact_canonical_21aa: 0
- not_numeric_uncensored_Kd: 81
- explicit_modification_or_not_linear: 37
- missing_target_sequence: 0

## Interpretation limits

This evaluates a sequence-only score against measured Kd records; it does not evaluate the adaptive optimizer or the target-aware structural score. Scores are compared within target/assay/publication contexts to reduce cross-target confounding. Small context counts or wide uncertainty mean the result remains preliminary. A positive association does not establish HERV-K-specific binding.

## Output files

- `eligible_scored_measurements.csv`
- `context_peptide_observations.csv`
- `within_context_correlations.csv`
- `single_substitution_pairs.csv`
- `mutation_position_summary.csv`
- `metrics.json`
