# Native sequence score versus experimental mutation effects

- Matched single-substitution pairs: **693**
- Distinct peptide sequences freshly scored: **621**
- Independent connected components: **111**
- Peptide lengths: **9–20**
- Exact 21-aa pairs in this broader set: **0**

## Main results

```json
{
  "n_pairs": 693,
  "spearman_delta_score_vs_delta_pKd": -0.10713922442596775,
  "n_informative_pairs": 466,
  "direction_accuracy_ties_half": 0.4248927038626609,
  "spearman_abs_delta_score_vs_abs_delta_pKd": -0.041827830390897214
}
```

## Cluster-bootstrap 95% intervals

```json
{
  "spearman_delta_score_vs_delta_pKd": [
    -0.20219641390743148,
    -0.017549489786270734
  ],
  "direction_accuracy_ties_half": [
    0.3654991147711051,
    0.4777496065745759
  ],
  "spearman_abs_delta_score_vs_abs_delta_pKd": [
    -0.12619585069276026,
    0.05186891862131836
  ]
}
```

## Position-level shuffle control

```json
{
  "eligible_length_position_cells": 71,
  "minimum_pairs_per_cell": 3,
  "position_profile_spearman": -0.13543037752301498,
  "position_shuffle_pvalue_two_sided": 0.22315536892621476,
  "permutations": 5000
}
```

## Interpretation

The native fitness score is a sequence-level heuristic with no receptor input. This analysis tests whether score differences align with experimental affinity differences, not whether the adaptive optimizer improves or whether HERV-K Env binding is predicted. The positional analysis is an exploratory sensitivity proxy, not direct validation of ARISE/ALE hotspot rankings.
