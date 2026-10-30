# P3 ANALYSIS LOGBOOK — Scientific Interpretation and Judge-Facing Conclusions

**Document date:** 26 October 2026, local time.  
**Purpose:** interpret completed P3 evidence without conflating optimization performance with biological binding accuracy.

## 1. Overall assessment

The most defensible positive result is that the adaptive search configuration improved the pipeline's **internal optimization score** across all three tested matched seeds. Mean final-generation fitness was 0.8005 for adaptive and 0.6488 for baseline, a descriptive difference of +0.1517. Mean best score per run was 0.8407 and 0.7091, respectively. Adaptive also averaged 151.00 distinct sequences per seed versus 132.67 for baseline.

Those observations do not establish greater experimental binding affinity. The native C score is a sequence-level heuristic that does not evaluate the 9SYA receptor. External PPIKB and PEPBI validations did not establish reliable affinity prediction, and the new optimizer candidates did not exactly match the measured-affinity data or the scanned PEPBI workbook.

## 2. Hypothesis-level conclusions

### Hypothesis A: adaptive search improves the internal objective

**Evidence:** each of the three matched seeds favored adaptive on final-generation mean fitness; the mean paired difference was +0.1517. Adaptive had more recorded unique sequences per run.

**Interpretation:** encouraging, consistent preliminary evidence for the two tested overall policies.

**Limitation:** `uniform_random_reset` changes the parent/reset behavior relative to the adaptive arm, so this comparison does not isolate mutation guidance as the sole cause. The final-generation endpoint is designated as the primary reporting metric during consolidation; it was not preregistered.

### Hypothesis B: adaptively selected positions show larger measured affinity effects

The original all-pairs ARISE test was inconclusive (positions 1–3: effect −0.0115, 95% CI [−0.1690, +0.1543], p = 0.9218). The informative-only subset's positive effect (+0.1716, p = 0.0528) is exploratory and cannot replace the primary result.

The new candidate artifacts did not preserve recoverable explicit ARISE-selected positions. The new proxy used observed mutation frequencies at positions 14, 3 and 8 and produced an external effect of −0.0828, p = 0.5495. Adaptive's higher hit rate at these positions (19.61% vs 12.08%) reflects different mutation behavior, not demonstrated biological importance.

**Conclusion:** the current evidence does not validate the adaptive position-selection mechanism against experimental affinity effects.

### Hypothesis C: the native C score predicts measured affinity

The broad benchmark reported Δscore-versus-ΔpKd Spearman ρ = −0.1071 (95% CI [−0.2022, −0.0176]) and informative direction accuracy of 42.49%. In the small exact-21-aa subset, the mean within-context association was −0.4872, with only nine informative pairs and one correctly directed pair in the recorded analysis.

No generated sequence exactly matched the measured-pKd dataset. Direct experimental affinity accuracy for the newly generated candidates is consequently **not estimable**.

**Conclusion:** the present data do not support a claim of reliable affinity prediction by native C score.

### Hypothesis D: receptor-specific PEPBI and structural validation support binding prediction

The strict PEPBI concordance point estimates were 0.3762 for native score and 0.3927 for Rosetta over 11 groups; both confidence intervals included 0.50. The receptor-context invariance control found zero native-score changes across 24 receptor-variant pairs, demonstrating that the score lacks receptor-specific response when peptide sequence is held fixed.

The existing structural benchmark had 18 mapped pairs, Spearman ρ = +0.0588 and 50% direction accuracy. A pose sanity check showed that the score preferred a poor docked pose over the crystal reference. No fresh candidate structures were available in the six optimizer runs, so there is no new candidate-arm structural comparison.

**Conclusion:** receptor-specific affinity performance and HERV-K Env binding have not been established.

### Hypothesis E: Ridge beats a training-fold mean baseline

The strict grouped comparison prevented target, peptide and exact-pair overlap across folds. Ridge improved the Spearman point estimate (0.1087 vs 0.0459) but was worse on MAE (0.8643 vs 0.8182), RMSE (1.3957 vs 1.3656), and direction accuracy (58.96% vs 62.66%).

**Conclusion:** Ridge did not consistently outperform the simpler mean baseline. This is a separate sequence-only predictive benchmark, not evidence of optimizer superiority.

## 3. Are there too many dependent variables?

The problem is not the raw count of metrics. It is the risk of treating many outcomes as equally important evidence for one broad claim.

For final reporting, organize the work into three experiments:

1. **Optimization experiment:** primary outcome is final-generation mean internal fitness; secondary outcomes are best fitness, unique-sequence fraction, mutation distribution and runtime. Seed is a blocking/matching factor. Preserve the qualification that guidance and parent-reset behavior were not isolated.
2. **Affinity prediction experiment:** experimental ΔpKd is the outcome; native score, Ridge and other scoring methods are predictors. Use identical eligible cases, grouped splits, and a small predefined set of metrics such as MAE, Spearman correlation and direction accuracy.
3. **Mutation-position experiment:** position membership under a frozen position-selection rule is the explanatory factor; absolute experimental ΔpKd change is the outcome, adjusted for assay/target context. The random equal-width window is the control. Keep informative-only filtering secondary.

The algorithm's entropy strength, hotspot weighting, mutation pressure and other configuration settings are independent variables only when they are deliberately manipulated. Changing several at once prevents attributing an observed effect to one component.

## 4. What is the best ISEF narrative?

The strongest defensible narrative is:

> Across three matched-seed computational runs, the adaptive search configuration achieved higher internal native-fitness scores and generated more distinct sequences than a uniform-reset baseline. Retrospective validation did not establish that the native score reliably tracked measured affinity changes or that the adaptive-observed mutation positions had larger experimental effects. Direct affinity accuracy for the generated candidates was unavailable because there were no exact matches to the measured-affinity datasets, and the run folders contained no candidate complex structures for fresh structural scoring. The work therefore provides preliminary evidence about optimization behavior, not demonstrated improvement in HERV-K Env binding.

This narrative communicates an engineering contribution while respecting the experimental evidence. Do not claim that HERV-K binding, ARG203/ARG205 engagement, or real-world affinity improvement is proven.

## 5. Is this all the data needed?

**It is sufficient to document the P3 work completed so far. It is not sufficient to establish biological binding accuracy.**

More analysis of the same reused retrospective data will not fill the central gap. The highest-value follow-up, if time and appropriate resources permit, is direct, independently measured affinity evidence for comparable candidate variants or a properly validated structural benchmark in which valid complex poses exist for both arms. A controlled optimizer ablation with identical parent-retention behavior would separately resolve whether mutation guidance itself caused the internal-score improvement.

Until those data exist, preserve the current result and report missing measurements honestly rather than generating extra metrics that do not test the claim.
