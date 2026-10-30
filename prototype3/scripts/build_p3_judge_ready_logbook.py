from pathlib import Path
import csv
import statistics
import zipfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path.cwd()
OUT = ROOT / "prototype3/logbook/p3_20261026_judge_ready"
G = OUT / "graphs"
OUT.mkdir(parents=True, exist_ok=True)
G.mkdir(parents=True, exist_ok=True)

# Values are transcribed from the P3 run outputs and reports presented in
# this conversation. This script organizes existing results; it does not
# rerun analyses or infer missing observations.

runs = [
    dict(arm="Baseline", seed=2601, n=160, unique=138,
         mean=0.6510945075, best=0.7055481254, final=0.6410184150),
    dict(arm="Baseline", seed=2602, n=160, unique=135,
         mean=0.6526790918, best=0.7145354445, final=0.6484208882),
    dict(arm="Baseline", seed=2603, n=160, unique=125,
         mean=0.6529478636, best=0.7071985064, final=0.6569524982),
    dict(arm="Adaptive", seed=2601, n=160, unique=151,
         mean=0.7593303901, best=0.8432968735, final=0.8051635273),
    dict(arm="Adaptive", seed=2602, n=160, unique=151,
         mean=0.7531530875, best=0.8455310975, final=0.8044782434),
    dict(arm="Adaptive", seed=2603, n=160, unique=151,
         mean=0.7354922051, best=0.8332672856, final=0.7919967889),
]

def avg(arm, field):
    return statistics.mean(r[field] for r in runs if r["arm"] == arm)

paired = []
for seed in (2601, 2602, 2603):
    b = next(r for r in runs if r["arm"] == "Baseline" and r["seed"] == seed)
    a = next(r for r in runs if r["arm"] == "Adaptive" and r["seed"] == seed)
    paired.append({
        "seed": seed,
        "adaptive_final_mean": a["final"],
        "baseline_final_mean": b["final"],
        "final_mean_difference": a["final"] - b["final"],
        "adaptive_best": a["best"],
        "baseline_best": b["best"],
        "best_difference": a["best"] - b["best"],
        "adaptive_unique_sequences": a["unique"],
        "baseline_unique_sequences": b["unique"],
    })

final_b, final_a = avg("Baseline", "final"), avg("Adaptive", "final")
best_b, best_a = avg("Baseline", "best"), avg("Adaptive", "best")
mean_b, mean_a = avg("Baseline", "mean"), avg("Adaptive", "mean")
uniq_b, uniq_a = avg("Baseline", "unique"), avg("Adaptive", "unique")

def csv_write(path, rows):
    if not rows:
        path.write_text("status\nno_records\n", encoding="utf-8")
        return
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

# ---------------------- Figure 1: matched optimizer runs -------------------

fig, ax = plt.subplots(figsize=(8.5, 5.5))
x = [0, 1]
for seed in (2601, 2602, 2603):
    b = next(r for r in runs if r["arm"] == "Baseline" and r["seed"] == seed)
    a = next(r for r in runs if r["arm"] == "Adaptive" and r["seed"] == seed)
    ax.plot(x, [b["final"], a["final"]], marker="o",
            linewidth=1.7, label=f"Seed {seed}")
ax.plot(x, [final_b, final_a], marker="D", markersize=8,
        linewidth=3, label="Mean")
ax.set_xticks(x, ["Baseline: uniform reset", "Adaptive: combined guidance"])
ax.set_ylabel("Final-generation mean native fitness")
ax.set_title("P3 optimizer comparison: three matched seeds")
ax.grid(axis="y", alpha=.25)
ax.legend(frameon=False)
fig.text(.5, .012,
         "160 candidate records per run. Native fitness is an internal sequence heuristic, "
         "not measured binding affinity.", ha="center", fontsize=8, wrap=True)
fig.tight_layout(rect=[0, .06, 1, 1])
fig.savefig(G / "01_matched_optimizer_results.png", dpi=200)
plt.close(fig)

# ---------------------- Figure 2: affinity associations --------------------

association = [
    ("Native C: Δscore vs ΔpKd\n693 pairs", -.107139, -.20220, -.01755),
    ("Native C: absolute changes\ninformative pairs", -.04183, -.1262, .0519),
    ("Exact 21-aa within-context\nmean Spearman", -.4872, -.7714, -.0301),
    ("Structural score\n18 mapped pairs", .0588235, None, None),
    ("Training-fold mean\ngrouped predictions", .0458706, None, None),
    ("Ridge\ngrouped predictions", .1087403, None, None),
]
fig, ax = plt.subplots(figsize=(9.8, 6))
ys = np.arange(len(association))[::-1]
for y, (label, rho, lo, hi) in zip(ys, association):
    if lo is None:
        ax.plot(rho, y, "o", markersize=6)
    else:
        ax.errorbar(rho, y, xerr=[[rho-lo], [hi-rho]],
                    fmt="o", capsize=3, markersize=6)
ax.axvline(0, color="black", linestyle="--", linewidth=1)
ax.set_yticks(ys, [r[0] for r in association])
ax.set_xlabel("Spearman correlation (ρ)")
ax.set_title("P3 score/affinity associations across distinct tests")
ax.grid(axis="x", alpha=.25)
fig.text(.5, .012,
         "Intervals shown where reported. Cohorts and endpoints differ; these estimates "
         "are not a single head-to-head ranking.", ha="center", fontsize=8, wrap=True)
fig.tight_layout(rect=[0, .06, 1, 1])
fig.savefig(G / "02_affinity_associations.png", dpi=200)
plt.close(fig)

# ---------------------- Figure 3: position-window results ------------------

positions = [
    ("Positions 1–3, all pairs", -.0115015, -.1690, .1543, "p=.9218"),
    ("Positions 1–3, informative-only", .1715996, .0148125, .3418984, "p=.0528"),
    ("Previous positions 2–3", -.018426, -.2026, .1684, "p=.8454"),
    ("New proxy: positions 14, 3, 8", -.0828259, None, None, "p=.5495"),
]
fig, ax = plt.subplots(figsize=(10, 5.7))
ys = np.arange(len(positions))[::-1]
for y, (label, effect, lo, hi, p) in zip(ys, positions):
    if lo is None:
        ax.plot(effect, y, "o", markersize=7)
    else:
        ax.errorbar(effect, y, xerr=[[effect-lo], [hi-effect]],
                    fmt="o", capsize=3, markersize=6)
    ax.annotate(p, (effect, y), xytext=(7, 5), textcoords="offset points", fontsize=8)
ax.axvline(0, color="black", linestyle="--", linewidth=1)
ax.set_yticks(ys, [r[0] for r in positions])
ax.set_xlabel("Context-adjusted difference in absolute ΔpKd")
ax.set_title("P3 mutation-position validation")
ax.grid(axis="x", alpha=.25)
fig.text(.5, .012,
         "Positive values mean larger experimental effects inside the window. The new proxy "
         "positions were not recovered as logged ARISE selections; no CI is available for that proxy.",
         ha="center", fontsize=8, wrap=True)
fig.tight_layout(rect=[0, .07, 1, 1])
fig.savefig(G / "03_mutation_position_validation.png", dpi=200)
plt.close(fig)

# ---------------------- Figure 4: grouped Ridge comparison -----------------

models = ["Training-fold mean", "Ridge"]
mae = [.8182125546, .8643157903]
rmse = [1.3655814586, 1.3956728700]
direction = [.6266094421, .5896226415]

fig, ax = plt.subplots(figsize=(7.8, 5))
x = np.arange(2)
w = .34
ax.bar(x-w/2, mae, w, label="MAE")
ax.bar(x+w/2, rmse, w, label="RMSE")
ax.set_xticks(x, models)
ax.set_ylabel("Prediction error (pKd-difference units)")
ax.set_title("P3 strict grouped validation")
ax.legend(frameon=False)
ax.grid(axis="y", alpha=.25)
fig.text(.5, .012,
         f"Direction accuracy: mean={direction[0]:.1%}; Ridge={direction[1]:.1%}. "
         "Lower MAE/RMSE is better.", ha="center", fontsize=8, wrap=True)
fig.tight_layout(rect=[0, .06, 1, 1])
fig.savefig(G / "04_grouped_ridge_vs_mean.png", dpi=200)
plt.close(fig)

# ---------------------- Figure 5: strict PEPBI concordance -----------------

pepbi = [
    ("Native sequence score", .3762, .2110, .5315),
    ("Rosetta dG", .3927, .2345, .5475),
]
fig, ax = plt.subplots(figsize=(7.8, 4.4))
for y, (label, value, lo, hi) in zip([1, 0], pepbi):
    ax.errorbar(value, y, xerr=[[value-lo], [hi-value]],
                fmt="o", capsize=4, markersize=7)
ax.axvline(.5, color="black", linestyle="--", linewidth=1.3,
           label="Chance reference = 0.50")
ax.set_yticks([1, 0], [p[0] for p in pepbi])
ax.set_xlim(0, .75)
ax.set_xlabel("Mean group concordance")
ax.set_title("P3 strict PEPBI comparison")
ax.grid(axis="x", alpha=.25)
ax.legend(frameon=False)
fig.text(.5, .012,
         "Calculated ΔG from KD; 11 independent binding groups. Both CIs include 0.50. "
         "The separately reported experimental-ΔG analysis had one group only.",
         ha="center", fontsize=8, wrap=True)
fig.tight_layout(rect=[0, .07, 1, 1])
fig.savefig(G / "05_pepbi_concordance.png", dpi=200)
plt.close(fig)

# ---------------------- Figure 6: observed mutation proxy ------------------

fig, ax = plt.subplots(figsize=(6.8, 4.5))
labels = ["Adaptive", "Baseline"]
values = [.1961, .1208]
bars = ax.bar(labels, values)
ax.set_ylim(0, .25)
ax.set_ylabel("Fraction of extracted mutation events")
ax.set_title("Mutation events at adaptive-observed proxy positions")
ax.grid(axis="y", alpha=.25)
for bar, value in zip(bars, values):
    ax.text(bar.get_x()+bar.get_width()/2, value+.006,
            f"{value:.2%}", ha="center")
fig.text(.5, .012,
         "Proxy positions 14, 3, 8 were inferred from mutation frequency, not recovered "
         "from explicit ARISE hotspot logs.", ha="center", fontsize=8, wrap=True)
fig.tight_layout(rect=[0, .07, 1, 1])
fig.savefig(G / "06_mutation_proxy_frequency.png", dpi=200)
plt.close(fig)

# ---------------------- Shared evidence table ------------------------------

evidence = [
    dict(area="Optimizer comparison", status="completed exploratory",
         result="6 runs; 3 matched seeds; all seeds favored adaptive on final-generation mean fitness",
         source="prototype3/results/p3_paired_arm_20261027T022810Z/"),
    dict(area="PPIKB native score", status="completed retrospective",
         result="693 pairs; rho=-0.1071; direction accuracy=42.49%",
         source="prototype3/results/native_score_all_mutation_pairs/"),
    dict(area="Exact 21-aa score benchmark", status="completed exploratory",
         result="129 rows; 67 peptides; 9 informative pairs; 1/9 correct-direction pair",
         source="prototype3/results/p3_21aa_affinity_audit/native_c_score_21aa/"),
    dict(area="ARISE positions 1-3", status="completed retrospective",
         result="all-pairs effect=-0.0115; CI[-0.1690,0.1543]; p=0.9218",
         source="prototype3/results/arise_experimental_position_audit/"),
    dict(area="New mutation position proxy", status="completed exploratory",
         result="proxy positions 14,3,8; effect=-0.0828; p=0.5495",
         source="prototype3/results/p3_paired_arm_20261027T022810Z/p3_validation_comparison/"),
    dict(area="PEPBI strict matched", status="completed retrospective",
         result="native concordance=0.3762; Rosetta=0.3927; 11 binding groups",
         source="prototype3/scripts/pepbi_strict_matched_benchmark.py"),
    dict(area="PEPBI receptor-context control", status="completed retrospective",
         result="0 native-score changes in 24 matched receptor-variant pairs",
         source="prototype3/scripts/pepbi_target_mutation_invariance.py"),
    dict(area="Structural mutation benchmark", status="carried forward, not rerun",
         result="18 pairs; rho=0.0588; direction accuracy=50%",
         source="prototype3/results/experimental_mutation_benchmark/"),
    dict(area="Grouped Ridge vs training-fold mean", status="rerun completed",
         result="Ridge worse on MAE and direction accuracy; grouped overlap constraints applied",
         source="prototype3/results/p3_paired_arm_20261027T022810Z/grouped_sequence_validation/"),
    dict(area="Generated-candidate pKd accuracy", status="not estimable",
         result="0 exact generated candidate sequence matches in the measured-pKd audit",
         source="P3 candidate sequence audit"),
    dict(area="Generated-candidate PEPBI accuracy", status="not estimable",
         result="0 exact candidate sequence matches in scanned PEPBI workbook",
         source="P3 candidate sequence audit"),
    dict(area="New candidate structural scoring", status="not estimable",
         result="0 complex PDB files found in the six optimizer run folders",
         source="P3 candidate structure audit"),
    dict(area="Explicit hotspot provenance", status="not recovered",
         result="0 new seeds with recoverable explicit selected hotspot positions",
         source="P3 candidate mutation audit"),
]
csv_write(OUT / "P3_EVIDENCE_REGISTER.csv", evidence)
csv_write(OUT / "OPTIMIZER_SEED_RESULTS.csv", runs)
csv_write(OUT / "MATCHED_SEED_DIFFERENCES.csv", paired)
csv_write(OUT / "ASSOCIATION_PLOT_DATA.csv", [
    dict(label=label, rho=rho, ci_low=lo, ci_high=hi)
    for label, rho, lo, hi in association
])
csv_write(OUT / "POSITION_PLOT_DATA.csv", [
    dict(label=label, effect=eff, ci_low=lo, ci_high=hi, p_value=p)
    for label, eff, lo, hi, p in positions
])
csv_write(OUT / "GROUPED_MODEL_METRICS.csv", [
    dict(model="Training-fold mean", mae=mae[0], rmse=rmse[0],
         spearman=0.0458706, direction_accuracy=direction[0], coverage=1.0),
    dict(model="Ridge", mae=mae[1], rmse=rmse[1],
         spearman=0.1087403, direction_accuracy=direction[1],
         coverage=0.9098712),
])

# ---------------------- Data logbook ---------------------------------------

data = f"""# P3 DATA LOGBOOK — Silicon Virovore

**Research period documented:** 26 October 2026 (local time).  
**Scope:** completed P3 retrospective validations and the new 3-seed matched optimizer comparison.  
**Purpose:** preserve numerical observations, configurations, data provenance, and missing evidence.

## 1. Research questions and reporting hierarchy

P3 contains two distinct scientific questions:

1. Does an adaptive search policy improve the pipeline's internal optimization objective compared with a baseline?
2. Do the pipeline's scores and selected positions correspond to experimentally measured affinity changes or validated receptor–peptide structures?

These questions are not interchangeable. For this consolidated record, final-generation mean native fitness is the **primary reported endpoint for the optimizer comparison**. This reporting hierarchy is being formalized during documentation after the runs; it was not preregistered before them. Other optimizer metrics are secondary/exploratory.

For affinity validation, the primary reported retrospective measures are the broad native-score association with measured mutation effects, the all-pairs ARISE window test, and the strict matched PEPBI group-concordance result. The exact-21-aa and new position-proxy results are exploratory or limited by small samples.

## 2. New optimizer comparison: methods and measurements

**Runs:** `prototype3/results/p3_paired_arm_20261027T022810Z/`

- Three matched seeds: 2601, 2602, 2603.
- 16 candidates × 10 generations = 160 recorded candidate records per run.
- Shared seed peptide: `MKLAVFALLVFFAGSSDLIRR` (21 aa).
- Mutation rate: 0.05.
- Elite count: 2.
- Tournament size: 3.
- Adaptive mode: `hotspot_entropy_substitution_epistasis_exploratory_v2`.
- Baseline mode: `uniform_random_reset`.

| Arm | Seed | Records | Unique sequences | Mean fitness over records | Best fitness | Final-generation mean |
|---|---:|---:|---:|---:|---:|---:|
"""
for r in runs:
    data += (
        f"| {r['arm']} | {r['seed']} | {r['n']} | {r['unique']} | "
        f"{r['mean']:.5f} | {r['best']:.5f} | {r['final']:.5f} |\\n"
    )

data += f"""
### Across-seed summary

| Metric | Baseline mean | Adaptive mean | Adaptive − baseline |
|---|---:|---:|---:|
| Mean candidate-record fitness | {mean_b:.5f} | {mean_a:.5f} | {mean_a-mean_b:+.5f} |
| Best score per run | {best_b:.5f} | {best_a:.5f} | {best_a-best_b:+.5f} |
| Final-generation mean fitness | {final_b:.5f} | {final_a:.5f} | {final_a-final_b:+.5f} |
| Unique sequences per run | {uniq_b:.2f} | {uniq_a:.2f} | {uniq_a-uniq_b:+.2f} |

All three matched seeds favored adaptive for final-generation mean native fitness. This is an encouraging descriptive result for the tested search configurations.

**Design qualification:** `uniform_random_reset` returns to the initial seed each generation while the adaptive arm follows an evolving parent. Therefore, this test compares the overall policies and does not isolate mutation guidance from parent/reset behavior. Three seeds provide preliminary evidence, not a high-powered general claim.

![Matched optimizer results](graphs/01_matched_optimizer_results.png)

## 3. Broad native C score versus experimental mutation effects

Source: `prototype3/results/native_score_all_mutation_pairs/scored_mutation_pairs.csv`.

- 693 matched single-substitution pairs.
- 621 unique peptides and 123 target sequences.
- Peptide lengths 9–20 aa.
- 466 informative pairs with |ΔpKd| ≥ 0.3.

| Measurement | Result |
|---|---:|
| Spearman ρ, Δnative score vs ΔpKd | −0.107139 |
| 95% cluster-bootstrap CI | [−0.20220, −0.01755] |
| Informative direction accuracy, ties half | 42.49% |
| Direction-accuracy 95% CI | [36.55%, 47.78%] |
| Spearman ρ, absolute score change vs absolute ΔpKd | −0.04183 |
| 95% CI for absolute-change correlation | [−0.1262, +0.0519] |
| Position-profile association | ρ = −0.1354; shuffle p = 0.223 |

The C score uses a target-null sequence-level fitness calculation. These results do not demonstrate target-specific binding affinity and do not directly measure the newly generated candidates.

![Affinity associations](graphs/02_affinity_associations.png)

## 4. Exact 21-aa affinity subset

Source: `prototype3/results/p3_21aa_affinity_audit/native_c_score_21aa/`.

- 129 eligible measurement rows.
- 67 unique peptides and 73 target/assay/publication contexts.
- Seven contexts had at least three peptides with variable score and affinity.
- 20 single-substitution context rows; nine informative pairs; seven unique peptide pairs and five mutation positions.
- Mean within-context Spearman ρ = −0.4872; 95% group-bootstrap CI [−0.7714, −0.0301].
- Informative direction accuracy = 1/9 = 11.11%.

This is relevant to the 21-aa optimizer seed length, but it is a small exploratory subset and cannot support a precise general estimate.

## 5. ARISE mutation-position validation

Source: `prototype3/results/arise_experimental_position_audit/`.

| Analysis | Context-adjusted difference | 95% CI | Two-sided permutation p |
|---|---:|---:|---:|
| Positions 1–3; all pairs (693) | −0.0115 | [−0.1690, +0.1543] | 0.9218 |
| Positions 2–3; sensitivity check | −0.0184 | [−0.2026, +0.1684] | 0.8454 |
| Positions 1–3; informative-only (466) | +0.1716 | [+0.0148, +0.3419] | 0.0528 |

The all-pairs test is the primary result and found no reliable advantage. The filtered informative-only result is secondary/exploratory; it is not independent confirmation. The external peptides differ from the optimizer seed and are not HERV-K Env-specific measurements.

![Mutation-position results](graphs/03_mutation_position_validation.png)

## 6. New mutation-pattern proxy audit

The new run artifacts did not reveal explicitly logged selected ARISE hotspot positions. The follow-up instead used positions 14, 3 and 8, identified from the most frequent adaptive mutation events that were represented in the experimental dataset. These are **observed-mutation proxy positions**, not confirmed ARISE selections.

- Adaptive parent–child mutation events extracted: 510.
- Baseline events extracted: 480.
- Fraction of events at positions 14, 3 or 8: adaptive 19.61%; baseline 12.08%.
- External pairs inside/outside those positions: 179 / 514.
- Context-adjusted effect: −0.0828.
- Random equal-size position-set permutation p = 0.5495.

The proxy position set did not show larger external affinity effects than the other positions. More frequent adaptive mutation at a position is not equivalent to demonstrating that the position is biologically important.

![Mutation proxy frequency](graphs/06_mutation_proxy_frequency.png)

## 7. Strict PEPBI comparison and receptor-context control

For calculated ΔG from K_D:

| Predictor | Mean group concordance | Independent binding groups | 95% group-bootstrap CI |
|---|---:|---:|---:|
| Native sequence score | 0.3762 | 11 | [0.2110, 0.5315] |
| Rosetta dG | 0.3927 | 11 | [0.2345, 0.5475] |

Both point estimates were below the 0.50 reference, but both CIs included 0.50. The separate reported-ΔG test had one independent group and is insufficient for general conclusions.

Receptor-context invariance analysis: four matched peptide contexts, 24 receptor-variant pairs; the native score changed in 0 pairs, while 7 pairs had measured ΔG differences >0.5 kcal/mol. This illustrates the limitation of a peptide-only score for receptor-specific prediction.

The new candidate sequence scan found zero exact PEPBI sequence matches; therefore, the PEPBI result is a retrospective method benchmark and not separate accuracy evidence for each new optimizer arm.

![PEPBI concordance](graphs/05_pepbi_concordance.png)

## 8. Structural-scoring evidence

The existing structural mutation benchmark was **carried forward, not rerun** in the latest P3 batch.

- 18 mapped pairs; 7 target groups; 28 unique complex-chain combinations.
- Overall Spearman ρ = +0.0588; exploratory p = 0.8167.
- Direction accuracy = 50%.
- Informative subset: 14 pairs, ρ = +0.0593, p = 0.8403, direction accuracy 50%.
- Same observed receptor-sequence subset: 10 pairs, 5 targets, ρ = +0.2849, p = 0.425, direction accuracy 60%; reported intervals were broad.

The 1YCR MDM2/p53 pose sanity check reported 9.965 Å aligned peptide RMSD for the docked pose; the scorer gave that docked pose 0.97751 versus 0.96537 for the crystal pose, preferring the poor pose. This is a limitation of the scorer's current validation.

The six new optimizer output folders had zero candidate complex PDB files available for fresh target-aware scoring. Thus there is no new structural accuracy result for adaptive-generated versus baseline-generated candidates.

## 9. Strict grouped Ridge versus training-fold mean

Source: `{ 'prototype3/results/p3_paired_arm_20261027T022810Z/grouped_sequence_validation/run_20261027T023029Z/metrics.json' }`.

The validation used five folds over 111 connected components, covering 693 pairs and 123 target sequences. Connected-component grouping kept each target sequence, peptide sequence and exact unordered peptide pair within one fold; the audit reported no overlap across folds.

| Predictor | MAE | RMSE | Spearman ρ | Direction accuracy | Coverage |
|---|---:|---:|---:|---:|---:|
| Training-fold mean | 0.8182 | 1.3656 | 0.0459 | 62.66% | 100.00% |
| Ridge | 0.8643 | 1.3957 | 0.1087 | 58.96% | 90.99% |

Ridge had a higher Spearman estimate but worse MAE/RMSE and lower direction accuracy. The reported paired-bootstrap interval for the rho difference included zero. This is a separate sequence-only predictor evaluation—not an optimizer-arm comparison.

![Grouped Ridge versus mean](graphs/04_grouped_ridge_vs_mean.png)

## 10. Candidate-specific evidence and execution status

- Exact generated-sequence matches to the scanned measured-pKd dataset: **0**.
- Exact generated-sequence matches to the PEPBI workbook: **0**.
- New-run complex PDB files available for fresh structural scoring: **0**.
- Recoverable explicit ARISE-selected hotspot sets in the six new runs: **0 seeds**.

The previous retrospective P3 batch reported that all ten listed scripts completed execution. The six new optimizer runs and fresh grouped Ridge evaluation also completed. **Successful execution does not mean the hypothesis passed.** The old structural result remains explicitly carried forward.

## 11. Provenance

See `P3_EVIDENCE_REGISTER.csv` and the CSV inputs to the graphs. Original result files and the earlier logbook bundle remain untouched.

This entry records the P3 work performed to date. It does not claim prospective affinity accuracy or demonstrated HERV-K Env binding.
"""

(OUT / "P3_DATA_LOGBOOK_ENTRY.md").write_text(data, encoding="utf-8")

# ---------------------- Analysis logbook ----------------------------------

analysis = f"""# P3 ANALYSIS LOGBOOK — Scientific Interpretation and Judge-Facing Conclusions

**Document date:** 26 October 2026, local time.  
**Purpose:** interpret completed P3 evidence without conflating optimization performance with biological binding accuracy.

## 1. Overall assessment

The most defensible positive result is that the adaptive search configuration improved the pipeline's **internal optimization score** across all three tested matched seeds. Mean final-generation fitness was {final_a:.4f} for adaptive and {final_b:.4f} for baseline, a descriptive difference of {final_a-final_b:+.4f}. Mean best score per run was {best_a:.4f} and {best_b:.4f}, respectively. Adaptive also averaged {uniq_a:.2f} distinct sequences per seed versus {uniq_b:.2f} for baseline.

Those observations do not establish greater experimental binding affinity. The native C score is a sequence-level heuristic that does not evaluate the 9SYA receptor. External PPIKB and PEPBI validations did not establish reliable affinity prediction, and the new optimizer candidates did not exactly match the measured-affinity data or the scanned PEPBI workbook.

## 2. Hypothesis-level conclusions

### Hypothesis A: adaptive search improves the internal objective

**Evidence:** each of the three matched seeds favored adaptive on final-generation mean fitness; the mean paired difference was +{final_a-final_b:.4f}. Adaptive had more recorded unique sequences per run.

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
"""
(OUT / "P3_ANALYSIS_LOGBOOK_ENTRY.md").write_text(analysis, encoding="utf-8")

# ---------------------- Readme and archive -------------------------------

(OUT / "README.md").write_text("""# P3 judge-ready logbook bundle

- `P3_DATA_LOGBOOK_ENTRY.md`: complete factual data record, methods, metrics and provenance.
- `P3_ANALYSIS_LOGBOOK_ENTRY.md`: interpretation, hypothesis status, dependent-variable structure and next evidence.
- `P3_EVIDENCE_REGISTER.csv`: completed, carried-forward and not-estimable analyses.
- `OPTIMIZER_SEED_RESULTS.csv`: six matched run outcomes.
- `MATCHED_SEED_DIFFERENCES.csv`: paired seed differences.
- `ASSOCIATION_PLOT_DATA.csv`, `POSITION_PLOT_DATA.csv`, `GROUPED_MODEL_METRICS.csv`: graph data.
- `graphs/`: six PNG figures used by the entries.

The figures are generated from the documented metrics. They distinguish internal optimization outcomes from retrospective affinity validation. Existing results are preserved; this package does not claim that all scientific hypotheses passed.
""", encoding="utf-8")

zip_path = OUT / "P3_JUDGE_READY_LOGBOOK_BUNDLE.zip"
with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p != zip_path:
            z.write(p, p.relative_to(OUT))

print("Created judge-ready P3 logbook bundle:")
print(OUT / "P3_DATA_LOGBOOK_ENTRY.md")
print(OUT / "P3_ANALYSIS_LOGBOOK_ENTRY.md")
print(OUT / "P3_EVIDENCE_REGISTER.csv")
print(G)
print(zip_path)
print()
print("Optimizer summary:")
print(f"  Baseline final-generation mean: {final_b:.5f}")
print(f"  Adaptive final-generation mean: {final_a:.5f}")
print(f"  Paired mean difference: {final_a-final_b:+.5f}")
print(f"  Baseline unique sequences/run: {uniq_b:.2f}")
print(f"  Adaptive unique sequences/run: {uniq_a:.2f}")
