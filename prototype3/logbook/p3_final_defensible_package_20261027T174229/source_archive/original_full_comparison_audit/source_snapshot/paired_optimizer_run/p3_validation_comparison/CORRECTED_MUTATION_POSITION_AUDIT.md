# Corrected P3 mutation-position audit

Experimental pair dataset: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/native_score_all_mutation_pairs/scored_mutation_pairs.csv`
Usable experimental single-substitution pairs: **693**
Adaptive mutation events successfully extracted: **510**
Baseline mutation events successfully extracted: **480**
Adaptive-observed proxy positions used: **[14, 3, 8]**

These are the most frequently mutated adaptive positions that are represented in the external benchmark. They are not described as ARISE-selected positions.

- Context-adjusted external effect (inside minus outside): **-0.08282590255609804**
- Random equal-size position-set permutation p-value: **0.5494835054981673**
- Pairs inside/outside the window: **179 / 514**

## Mutation-event hit rates

| Arm | Mutation events | Fraction inside adaptive proxy positions |
|---|---:|---:|
| adaptive | 510 | 0.1961 |
| baseline | 480 | 0.1208 |

The external positional statistic is shared by design: the same adaptive-derived positions are tested for both arms. The only arm-specific positional quantity is how often each policy mutates within those positions. This is retrospective and exploratory, not direct experimental validation of generated candidates.