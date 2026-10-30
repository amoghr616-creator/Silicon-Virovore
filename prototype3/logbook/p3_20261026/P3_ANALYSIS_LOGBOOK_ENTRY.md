# P3 Analysis Logbook Entry — Interpretation, Limitations and Next Decisions

**Work date:** 26 October 2026 (America/New_York).  
**Project:** Silicon Virovore — computational peptide optimization and retrospective validation.  
**Purpose:** interpret the combined P3 optimizer and external-validation evidence. This entry is separate from the factual data record.

## 1. Overall conclusion

The best-supported positive result is about **internal optimization**: under the tested configurations, the adaptive policy produced higher values of the pipeline's native fitness objective than the uniform-reset baseline in all three matched seeds. Mean final-generation fitness was 0.8005 for adaptive and 0.6488 for baseline, a descriptive difference of +0.1517. Mean best score per run was 0.8407 versus 0.7091.

That result is encouraging for the search policy but does **not** establish improved experimental binding affinity. The native C score calls the fitness function with `target=NULL`; it is a sequence-property heuristic, not a validated receptor-specific affinity function. The external P3 tests, overall, do not establish reliable prediction of affinity changes or HERV-K Env binding at ARG203/ARG205.

The defensible research narrative is therefore two-part:

1. Adaptive search behavior was promising on the internal optimization objective in the three matched runs.
2. External validation remains a separate and unresolved problem; most direct retrospective score/affinity checks were weak, negative, inconclusive or too small, and the newly generated candidates had no exact measured-affinity matches.

## 2. Hypothesis-by-hypothesis assessment

### Hypothesis A: Adaptive search improves the internal objective

**Observation:** all three matched seeds had higher adaptive final-generation mean scores. The mean paired increase was +0.1517; the difference in mean best scores was +0.1316. Adaptive averaged 151.00 unique sequences per run, while baseline averaged 132.67.

**Assessment:** promising, consistent descriptive evidence for these three runs. It is not a high-powered estimate of general performance.

**Causal limitation:** the arms differ in the overall search policy, particularly parent evolution versus returning to the initial seed. The results cannot isolate the effect of mutation-guidance logic alone. A fair ablation would hold parent-retention/reset behavior fixed and vary only mutation selection.

![Figure 1: Matched optimizer fitness](graphs/p3_fig1_matched_optimizer_fitness.png)

### Hypothesis B: Adaptive-selected mutation positions are experimentally important

**Primary older test:** positions 1–3 did not beat random equal-width windows across all 693 pairs: context-adjusted effect −0.0115, 95% CI [−0.1690, +0.1543], p = 0.9218. The informative-only subset had an effect +0.1716 and p = 0.0528, but it is exploratory and was not the primary all-pairs result.

**New run-specific audit:** the selected hotspots were not explicitly recoverable in any of the three new adaptive run artifacts. The follow-up therefore used positions 14, 3 and 8 because they were frequently mutated by adaptive and represented in the external benchmark. This is a proxy, not an ARISE selection record. The context-adjusted effect was −0.0828 (179 external pairs inside, 514 outside), p = 0.5495. Adaptive mutation events hit those proxy positions more often than baseline events (19.61% vs 12.08%), but that shows a policy-frequency difference, not a biological advantage.

**Assessment:** the current data do not show that adaptively emphasized positions have larger experimental affinity effects. Do not present the proxy positions as verified ARISE hotspots.

![Figure 3: Mutation-position effects](graphs/p3_fig3_mutation_position_effects.png)

### Hypothesis C: Native C score predicts experimental affinity

**All mutation pairs:** Δscore versus ΔpKd had Spearman ρ = −0.1071, with 95% CI [−0.2022, −0.0176]. Informative-pair direction accuracy was 42.49%, 95% CI [36.55%, 47.78%], below the 50% reference.

**Exact 21-aa subset:** the mean within-context Spearman estimate was −0.4872, 95% CI [−0.7714, −0.0301], but only 9 informative pairs were available, with 1 correctly directed pair in the recorded analysis. This sample is too small to generalize confidently.

**Candidate-specific test:** none of the distinct generated sequences from the new optimizer runs matched the pKd benchmark exactly. Consequently, the direct adaptive-versus-baseline experimental-affinity accuracy of generated sequences is **not estimable from the current data**.

**Assessment:** current evidence does not support claiming that higher native C score predicts stronger real-world binding.

![Figure 2: Association estimates](graphs/p3_fig2_affinity_association_estimates.png)

### Hypothesis D: PEPBI and structural scoring support receptor-specific prediction

**PEPBI:** mean group concordance for calculated ΔG from K_D was 0.3762 for native score and 0.3927 for Rosetta, both below the 0.50 reference. Their reported 95% intervals included 0.50, and the independent group count was 11. The reported experimental-ΔG comparison had only one group, too few for a useful general inference.

**Receptor-context test:** the native score stayed unchanged in 24/24 matched receptor-variant pairs, even though 7 pairs had experimental ΔG differences greater than 0.5 kcal/mol. This demonstrates a limitation of using peptide-only score values for receptor-specific prediction.

**Structural scoring:** the existing 18-pair benchmark had ρ = +0.0588, p = 0.8167 and 50% direction accuracy. The pose sanity check also exposed a ranking problem: the score preferred the poor docked pose over the crystal reference. No new complex PDB files were present in the adaptive/baseline optimizer output folders, so fresh structural scores for their candidates were not generated.

**Assessment:** these tests do not establish receptor-specific affinity accuracy or improved HERV-K Env binding.

![Figure 5: PEPBI concordance](graphs/p3_fig5_pepbi_concordance.png)

### Hypothesis E: A learned sequence-only model adds value over a simple baseline

The grouped Ridge model used connected components to ensure that target sequences, peptide sequences and exact unordered pairs did not overlap between folds. It achieved a higher Spearman estimate than the training-fold mean, 0.1087 versus 0.0459, but worse MAE (0.8643 vs 0.8182), worse RMSE (1.3957 vs 1.3656) and lower direction accuracy (58.96% vs 62.66%). The reported confidence interval for the paired rho difference included zero.

**Assessment:** the Ridge model did not consistently outperform the simple training-fold mean. This is a valid negative comparator result.

![Figure 4: Grouped Ridge versus training-fold mean](graphs/p3_fig4_grouped_ridge_vs_mean.png)

## 3. What is the strongest ISEF claim supported?

A defensible statement is:

> In three matched-seed computational runs, the adaptive search configuration achieved higher internal native-fitness scores and yielded more distinct candidate sequences than a uniform-reset baseline. Retrospective affinity benchmarks, however, did not establish that the native score reliably tracks measured affinity changes, and the mutation-position tests did not establish that the new adaptive-observed positions were biologically more consequential. Because generated candidates had no exact matches to the external measured-affinity data and no candidate-run complex PDB files for fresh structural scoring, these results support a preliminary optimization-behavior claim—not a claim of improved experimental binding to HERV-K Env or the ARG203/ARG205 site.

This is stronger scientifically than claiming the pipeline “works” in the biological sense, because it separates a result the data actually support from the unvalidated downstream claim.

## 4. Is this all the data needed?

It is a sufficiently complete package to document **the P3 work performed to date**. It is **not enough to conclude that the full pipeline predicts real-world binding accurately**.

The biggest unresolved gaps are:

1. **Direct candidate-specific affinity evidence:** no exact generated sequence had a measured pKd in the scanned benchmark, and none exactly matched PEPBI. Without such outcomes, direct accuracy for the generated candidates cannot be calculated.
2. **Valid candidate–receptor structures:** there were no candidate-run complex PDB files to score. The old structural benchmark is independent and was carried forward, not a fresh adaptive-versus-baseline comparison.
3. **Traceable hotspot provenance:** the new run outputs did not preserve explicitly recoverable ARISE-selected positions. The position-frequency proxy cannot substitute for that record.
4. **Isolated adaptive effect:** the present optimizer comparison combines guidance-policy differences with different parent/reset behavior. A future controlled ablation should vary mutation guidance while holding parent-selection policy fixed.
5. **Independent confirmation:** the optimizer comparison used three seeds. More seeds would improve precision, but repeating seeds alone would not solve the missing experimental-affinity and structure evidence.

Given the current deadline, the recommendation is to preserve this bundle and document these limits rather than run another poorly matched analysis solely to produce a more favorable number. If later work is possible, the highest-value addition is a valid, same-case evaluation using exact experimental measurements or verified complex structures, coupled with a controlled ablation.

## 5. Data integrity and reproducibility

- Original run results and earlier P3 reports were preserved; this bundle writes to a new logbook folder.
- The three matched optimizer seeds are explicitly listed in the data entry and CSV files.
- No affinity values were imputed for candidates without measurements.
- Proxy positions are not labeled as ARISE-selected hotspots.
- Carried-forward structural results are clearly marked as not rerun.
- Test status `PASS` is interpreted as successful execution, not automatic confirmation of a hypothesis.
- The evidence register distinguishes completed, exploratory, carried-forward and not-estimable analyses.

## 6. Recommended final reporting language

Use “internal native-fitness improvement,” “retrospective affinity association,” “exploratory mutation-position analysis,” “grouped sequence-only comparator,” and “not yet established experimentally” where appropriate.

Avoid “optimized binding affinity,” “validated HERV-K binder,” “confirmed ARG203/ARG205 interaction,” or “proven real-world accuracy” because the current tests do not support those claims.
