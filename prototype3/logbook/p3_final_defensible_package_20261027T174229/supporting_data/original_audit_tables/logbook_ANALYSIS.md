# Logbook Entry — P3 Comparison: Analysis and Interpretation

**Date:** 27 October 2026  
**Experiment:** Interpretation of the consolidated P3 optimizer, mutation, affinity,
and structural evidence  
**Corresponding data entry:** `logbook_DATA_RESULTS.md`  
**Audit package:** `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/logbook/p3_comparison_full_audit_20261027T172008`

## 1. Central interpretation

This project has two distinct validation questions. The first asks whether adaptive guidance
optimizes the pipeline's internal objective better than the matched uniform baseline. The
second asks whether that objective and the chosen candidates predict real binding affinity.
The current evidence is more favorable for the first question than for the second.

## 2. Hypothesis A — Does adaptive guidance improve optimization?

The three matched seeds showed higher mean native score, best native score, and final-generation
mean score for the adaptive arm. The recorded mean native score across candidate records was
also higher for adaptive (about 0.7493) than baseline (about 0.6522). The candidate inventory
reported 450 unique adaptive sequences versus 275 baseline sequences.

This supports the narrower statement that adaptive guidance improved optimization against
the internal score in this three-seed experiment and yielded a more diverse observed sequence
set. The score difference was directionally consistent across the three seeds.

It does **not** establish that the adaptive candidates bind the target more strongly. If the
native `c_score` does not incorporate the receptor structure, optimizing that score can reward
properties of the sequence-level objective without reliably rewarding receptor-specific
binding. A higher optimizer score is evidence of optimization against that objective, not
evidence of a successful binder by itself.

The mutation audit supports that adaptive and uniform arms do not behave identically: the
adaptive arm has more unique parent-child transitions and different per-position mutation
frequencies. However, position-frequency differences vary by seed and remain descriptive.
The hotspot-guidance flag should not be interpreted as proof that every mutation was selected
because of a hotspot or lies inside the hotspot window.

## 3. Hypothesis B — Do scores or mutation positions track measured affinity?

The broad native-score benchmark provides unfavorable evidence for treating the native score
as a reliable affinity predictor: the score-change/affinity-change correlation was weakly
negative and direction accuracy was below 50%. The broad position 1–3 audit did not show a
general all-pair effect. These observations do not support a confident claim that the native
score or the tested mutation-position rule predicts experimentally favorable substitutions.

The grouped Ridge benchmark gives a mixed model comparison rather than a production-ready
predictor. Ridge ranked affinity changes somewhat better by Spearman correlation, but its
MAE, RMSE, and direction accuracy were worse than the training-fold mean. This illustrates
why model quality should not be declared from one favorable metric while ignoring the others.

The structural benchmark is underpowered: only 18 of 693 candidate mutation pairs were
successfully mapped and scored. Its overall correlation was near zero, its directional
accuracy was at chance, and the same-receptor subset of ten pairs is too small for a robust
claim. The pose sanity check also raises a scoring-validity concern because the recorded score
preferred the unfavorable pose. These results do not validate the structural score as a
reliable experimental-affinity proxy.

The PEPBI strict-matched benchmark has low-to-moderate concordance estimates with broad
confidence intervals and high reported sign-test p-values. Its sample structure is useful
for auditing methods, but the evidence does not establish dependable affinity ranking from
the native score. The comparison with Rosetta does not itself validate either score.

## 4. Critical issue — the new PEPBI exact-match count

The newly supplied worksheet reports 125 exact sequence matches for adaptive and zero for
baseline. This differs from the earlier audit that reported zero matches for both arms. It
is important not to dismiss the new result, but it is equally important not to overstate it.

Before using the count as project evidence, reconcile the audits at the sequence level:

- Export the exact sequence lists from both optimizer arms and the exact peptide sequences
  in the PEPBI source tables.
- Normalize only defensible formatting differences, such as whitespace and letter case;
  do not silently remove amino acids or alter sequences to force matches.
- Produce one row per matched candidate sequence showing the candidate sequence, arm,
  PEPBI source file, matching row or identifier, target/receptor context, experimental
  affinity field, units, and duplicate count.
- Report both the number of distinct matched sequences and the number of independent
  affinity contexts. Multiple rows for the same sequence or receptor context are not
  equivalent to independent experimental validation.
- Confirm whether matches exist because the adaptive candidates overlap the benchmark more
  often, because the previous audit scanned a different set of files, or because the two
  audits used different matching rules.

Until that reconciliation is complete, the accurate interpretation is: *a newer audit reports
potential PEPBI sequence overlap for 125 adaptive candidates, but the count is not yet
reconciled or verified at the sequence-and-affinity-record level.*

If these matches can be reproduced and linked to usable affinity data, they may enable a more
direct retrospective comparison of the two candidate sets. Even then, sequence overlap alone
does not prove better binding; the matched affinity distributions must be compared under
consistent target/receptor contexts, with duplicate handling and uncertainty reported.

## 5. Does the present evidence show translation to real-world binding?

No. It supports internal optimizer improvement and a change in candidate diversity, but
candidate-level experimental superiority has not yet been established. The external benchmarks
mainly test whether the scoring or mutation rules track measured affinity in retrospective
datasets. They are not equivalent to experimentally testing the candidates produced by the
adaptive and baseline arms.

The current data support the following conclusion:

> Across three matched seeds, adaptive guidance improved the optimizer's internal-score
> outcomes and increased the number of unique candidate sequences relative to the uniform
> baseline. However, retrospective experimental-affinity benchmarks did not establish that
> the native sequence score reliably predicts binding affinity, and the structural benchmark
> was limited by a small mapped sample and a pose-scoring concern. A newly supplied PEPBI
> audit reports 125 adaptive exact-sequence matches but conflicts with an earlier audit and
> lacks a visible matched-sequence list. Therefore, whether the optimizer's observed advantage
> translates into improved experimental binding remains unresolved.

## 6. What should happen next

Do not rerun the three-seed optimizer comparison merely to generate another version of the
same result. First reconcile the PEPBI match count and establish which exact candidate sequences
map to independent affinity records. This is the most immediate evidence gap raised by the
new spreadsheet.

After that, report the arm-specific matched affinity results only if the source records support
a valid comparison. Keep internal optimizer scores, retrospective mutation-prediction metrics,
structural proxy scores, and actual experimental binding measurements in separate result
categories. If no suitable affinity records can be linked to generated candidates, state that
candidate-level experimental validation is unavailable rather than treating a sequence match as
a successful binding result.
