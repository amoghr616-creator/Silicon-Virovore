# Experimental mutation benchmark: structural-score comparison

- Pair input rows: 693
- Pairs with verified exact-sequence mapped structures for both variants: 18
- Pairs scored successfully: 18
- Unique structure/chain combinations scored: 28
- Results directory: `prototype3/results/experimental_mutation_benchmark/structural_scoring_20261026T223525Z`

## Interpretation rules

- `delta_pKd_B_minus_A > 0` means variant B has the higher measured pKd.
- `delta_score_B_minus_A > 0` means variant B has the higher geometry score.
- Direction agreement is counted only for absolute experimental differences of at least 0.3 pKd units.
- The primary subset requires both variants to have the same observed receptor-chain amino-acid sequence.
- These are experimentally annotated structures scored retrospectively. They are not an independent prospective test, and crystallographic pose quality or construct differences may bias comparisons.
- A high score–affinity correlation or direction agreement in this small selected set would not establish that the scorer predicts affinity generally. Low or negative agreement is meaningful evidence against the current scoring formula.
- The scorer is not wired into the candidate-generation optimizer.

## Metrics

See `benchmark_metrics.json` for counts, associations, direction agreement, and the cluster-bootstrap uncertainty estimate where enough target groups are available.

