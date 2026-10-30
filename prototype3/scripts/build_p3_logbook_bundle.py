from __future__ import annotations

import csv
import json
import math
import statistics
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path.cwd()
OUT = ROOT / "prototype3/logbook/p3_20261026"
GRAPHS = OUT / "graphs"
OUT.mkdir(parents=True, exist_ok=True)
GRAPHS.mkdir(parents=True, exist_ok=True)

# All values below are transcribed from the completed P3 reports/runs
# supplied during this work. Existing source result files are not modified.

RUNS = [
    {"arm": "baseline", "seed": 2601, "records": 160, "unique": 138,
     "mean": 0.6510945074517631, "best": 0.7055481253563461,
     "final": 0.6410184150422924},
    {"arm": "baseline", "seed": 2602, "records": 160, "unique": 135,
     "mean": 0.6526790917715971, "best": 0.7145354444591977,
     "final": 0.6484208881823028},
    {"arm": "baseline", "seed": 2603, "records": 160, "unique": 125,
     "mean": 0.6529478635783754, "best": 0.7071985064014479,
     "final": 0.6569524982230924},
    {"arm": "adaptive", "seed": 2601, "records": 160, "unique": 151,
     "mean": 0.7593303900985404, "best": 0.8432968734829264,
     "final": 0.8051635273010702},
    {"arm": "adaptive", "seed": 2602, "records": 160, "unique": 151,
     "mean": 0.7531530875283788, "best": 0.8455310975436344,
     "final": 0.8044782434377395},
    {"arm": "adaptive", "seed": 2603, "records": 160, "unique": 151,
     "mean": 0.7354922050586472, "best": 0.8332672856265247,
     "final": 0.7919967888833692},
]

P2_REPORT = (
    "prototype3/reports/"
    "p3_baseline_rerun_20261027T020504Z/P3_BASELINE_COMPARISON.md"
)
P2_STATUS = (
    "prototype3/reports/"
    "p3_baseline_rerun_20261027T020504Z/run_status.txt"
)
OPT_ROOT = "prototype3/results/p3_paired_arm_20261027T022810Z"
PAIR_DATA = (
    "prototype3/results/native_score_all_mutation_pairs/"
    "scored_mutation_pairs.csv"
)
EXACT21 = (
    "prototype3/results/p3_21aa_affinity_audit/native_c_score_21aa/"
)
STRUCTURAL = (
    "prototype3/results/experimental_mutation_benchmark/"
    "structural_scoring_20261026T223525Z/benchmark_metrics.json"
)
GROUPED = (
    "prototype3/results/p3_paired_arm_20261027T022810Z/"
    "grouped_sequence_validation/run_20261027T023029Z/metrics.json"
)
MUTATION_AUDIT = (
    "prototype3/results/p3_paired_arm_20261027T022810Z/"
    "p3_validation_comparison/CORRECTED_MUTATION_POSITION_AUDIT.md"
)
FINAL_AUDIT = (
    "prototype3/results/p3_paired_arm_20261027T022810Z/"
    "p3_validation_comparison/P3_ADAPTIVE_BASELINE_VALIDATION.md"
)

# ----------------------------- Summary metrics -----------------------------

arms = {}
for arm in ("baseline", "adaptive"):
    rows = [r for r in RUNS if r["arm"] == arm]
    arms[arm] = {
        "mean_all_candidate_fitness": statistics.mean(r["mean"] for r in rows),
        "mean_best_fitness": statistics.mean(r["best"] for r in rows),
        "mean_final_generation_fitness": statistics.mean(r["final"] for r in rows),
        "mean_unique_sequences_per_seed": statistics.mean(r["unique"] for r in rows),
    }

paired = []
for seed in (2601, 2602, 2603):
    b = next(r for r in RUNS if r["arm"] == "baseline" and r["seed"] == seed)
    a = next(r for r in RUNS if r["arm"] == "adaptive" and r["seed"] == seed)
    paired.append({
        "seed": seed,
        "adaptive_final": a["final"],
        "baseline_final": b["final"],
        "final_difference": a["final"] - b["final"],
        "adaptive_best": a["best"],
        "baseline_best": b["best"],
        "best_difference": a["best"] - b["best"],
        "adaptive_unique": a["unique"],
        "baseline_unique": b["unique"],
    })

mean_final_diff = (
    arms["adaptive"]["mean_final_generation_fitness"]
    - arms["baseline"]["mean_final_generation_fitness"]
)
mean_best_diff = arms["adaptive"]["mean_best_fitness"] - arms["baseline"]["mean_best_fitness"]
mean_all_diff = arms["adaptive"]["mean_all_candidate_fitness"] - arms["baseline"]["mean_all_candidate_fitness"]
unique_diff = (
    arms["adaptive"]["mean_unique_sequences_per_seed"]
    - arms["baseline"]["mean_unique_sequences_per_seed"]
)

# ----------------------------- Graph 1 -------------------------------------

fig, ax = plt.subplots(figsize=(8.5, 5.5))
x = np.array([0, 1])

for seed in (2601, 2602, 2603):
    b = next(r for r in RUNS if r["arm"] == "baseline" and r["seed"] == seed)
    a = next(r for r in RUNS if r["arm"] == "adaptive" and r["seed"] == seed)
    ax.plot(x, [b["final"], a["final"]], marker="o", linewidth=1.6,
            alpha=0.8, label=f"Seed {seed}")

ax.plot(
    x,
    [arms["baseline"]["mean_final_generation_fitness"],
     arms["adaptive"]["mean_final_generation_fitness"]],
    marker="D", linewidth=3, markersize=8, label="Mean across 3 seeds"
)
ax.set_xticks(x, ["Baseline\nuniform reset", "Adaptive\ncombined guidance"])
ax.set_ylabel("Mean native fitness in final generation")
ax.set_title("P3: Matched-seed internal optimization results")
ax.grid(axis="y", alpha=0.25)
ax.legend(frameon=False)
fig.text(
    0.5, 0.015,
    "16 candidates × 10 generations per run; matched seeds 2601–2603. "
    "Native fitness is a sequence heuristic, not measured binding affinity.",
    ha="center", va="bottom", fontsize=8, wrap=True
)
fig.tight_layout(rect=[0, 0.06, 1, 1])
fig.savefig(GRAPHS / "p3_fig1_matched_optimizer_fitness.png", dpi=200)
plt.close(fig)

# ----------------------------- Graph 2 -------------------------------------
# Association estimates from distinct validation tasks. These are not
# head-to-head estimates on a shared cohort, so they must not be ranked as
# though they arose from one common experiment.

association_rows = [
    {
        "label": "Native C: Δscore vs ΔpKd\n693 pairs",
        "rho": -0.107139, "low": -0.20220, "high": -0.01755,
    },
    {
        "label": "Native C: |Δscore| vs |ΔpKd|\ninformative pairs",
        "rho": -0.04183, "low": -0.1262, "high": 0.0519,
    },
    {
        "label": "Exact 21-aa within-context\nmean Spearman",
        "rho": -0.4872, "low": -0.7714, "high": -0.0301,
    },
    {
        "label": "Structural score: Δscore vs ΔpKd\n18 mapped pairs",
        "rho": 0.0588235, "low": None, "high": None,
    },
    {
        "label": "Grouped training-fold mean\nheld-out predictions",
        "rho": 0.0458706, "low": None, "high": None,
    },
    {
        "label": "Grouped Ridge\nheld-out predictions",
        "rho": 0.1087403, "low": None, "high": None,
    },
]

fig, ax = plt.subplots(figsize=(10, 6.2))
y = np.arange(len(association_rows))[::-1]

for yi, row in zip(y, association_rows):
    if row["low"] is None:
        ax.plot(row["rho"], yi, marker="o", linestyle="none", markersize=7)
    else:
        ax.errorbar(
            row["rho"], yi,
            xerr=[[row["rho"] - row["low"]], [row["high"] - row["rho"]]],
            fmt="o", capsize=3, markersize=6
        )

ax.axvline(0, linewidth=1, linestyle="--")
ax.set_yticks(y, [r["label"] for r in association_rows])
ax.set_xlabel("Spearman correlation (ρ)")
ax.set_title("P3: Association with experimental affinity outcomes")
ax.grid(axis="x", alpha=0.25)
fig.text(
    0.5, 0.012,
    "Intervals shown only where the reported 95% CI was available. "
    "Cohorts, endpoints, and validation designs differ; this is not a common-score leaderboard.",
    ha="center", fontsize=8, wrap=True
)
fig.tight_layout(rect=[0, 0.07, 1, 1])
fig.savefig(GRAPHS / "p3_fig2_affinity_association_estimates.png", dpi=200)
plt.close(fig)

# ----------------------------- Graph 3 -------------------------------------

position_rows = [
    {
        "label": "Original positions 1–3\nall pairs; p=0.9218",
        "effect": -0.0115015, "low": -0.1690, "high": 0.1543,
    },
    {
        "label": "Positions 1–3\ninformative-only subset; p=0.0528",
        "effect": 0.1715996, "low": 0.0148125, "high": 0.3418984,
    },
    {
        "label": "Previous positions 2–3\nsensitivity; p=0.8454",
        "effect": -0.018426, "low": -0.2026, "high": 0.1684,
    },
    {
        "label": "New adaptive-observed proxy\npositions 14, 3, 8; p=0.5495",
        "effect": -0.0828259, "low": None, "high": None,
    },
]

fig, ax = plt.subplots(figsize=(10, 5.8))
y = np.arange(len(position_rows))[::-1]
for yi, row in zip(y, position_rows):
    if row["low"] is None:
        ax.plot(row["effect"], yi, marker="o", linestyle="none", markersize=7)
    else:
        ax.errorbar(
            row["effect"], yi,
            xerr=[[row["effect"] - row["low"]], [row["high"] - row["effect"]]],
            fmt="o", capsize=3, markersize=6
        )
ax.axvline(0, linewidth=1, linestyle="--")
ax.set_yticks(y, [r["label"] for r in position_rows])
ax.set_xlabel("Context-adjusted difference in |ΔpKd|: selected window minus other positions")
ax.set_title("P3: Retrospective mutation-position validation")
ax.grid(axis="x", alpha=0.25)
fig.text(
    0.5, 0.012,
    "Positive values indicate larger absolute affinity effects inside the window. "
    "The proxy positions were inferred from mutation frequency, not logged as ARISE-selected hotspots; "
    "no CI was calculated for that proxy result.",
    ha="center", fontsize=8, wrap=True
)
fig.tight_layout(rect=[0, 0.08, 1, 1])
fig.savefig(GRAPHS / "p3_fig3_mutation_position_effects.png", dpi=200)
plt.close(fig)

# ----------------------------- Graph 4 -------------------------------------

models = ["Training-fold mean", "Ridge"]
mae = [0.818212554563056, 0.8643157903240979]
rmse = [1.3655814586084296, 1.395672869959655]
directions = [0.6266094420600858, 0.589622641509434]

fig, ax = plt.subplots(figsize=(8, 5.3))
xx = np.arange(len(models))
width = 0.34
ax.bar(xx - width/2, mae, width, label="MAE")
ax.bar(xx + width/2, rmse, width, label="RMSE")
ax.set_xticks(xx, models)
ax.set_ylabel("Error (pKd-difference units; lower is better)")
ax.set_title("P3: Strict grouped model validation")
ax.legend(frameon=False)
ax.grid(axis="y", alpha=0.25)
fig.text(
    0.5, 0.012,
    f"Direction accuracy: fold mean {directions[0]:.1%}; Ridge {directions[1]:.1%}. "
    "Ridge had higher Spearman ρ (0.109 vs 0.046), but worse MAE/RMSE and lower direction accuracy.",
    ha="center", fontsize=8, wrap=True
)
fig.tight_layout(rect=[0, 0.08, 1, 1])
fig.savefig(GRAPHS / "p3_fig4_grouped_ridge_vs_mean.png", dpi=200)
plt.close(fig)

# ----------------------------- Graph 5 -------------------------------------

pepbi_rows = [
    {"label": "Native sequence score", "mean": 0.3762, "low": 0.2110, "high": 0.5315},
    {"label": "Rosetta dG", "mean": 0.3927, "low": 0.2345, "high": 0.5475},
]

fig, ax = plt.subplots(figsize=(7.8, 4.6))
yy = np.array([1, 0])
for yi, row in zip(yy, pepbi_rows):
    ax.errorbar(
        row["mean"], yi,
        xerr=[[row["mean"] - row["low"]], [row["high"] - row["mean"]]],
        fmt="o", capsize=4, markersize=7
    )
ax.axvline(0.5, linewidth=1.5, linestyle="--", label="Chance reference = 0.50")
ax.set_yticks(yy, [r["label"] for r in pepbi_rows])
ax.set_xlim(0, 0.75)
ax.set_xlabel("Mean group concordance")
ax.set_title("P3: Strict PEPBI matched comparison")
ax.legend(frameon=False)
ax.grid(axis="x", alpha=0.25)
fig.text(
    0.5, 0.012,
    "Calculated ΔG from KD; 11 independent binding groups. Both mean estimates are below 0.50, "
    "but both 95% CIs include 0.50. Reported experimental ΔG had only one group and is not informative.",
    ha="center", fontsize=8, wrap=True
)
fig.tight_layout(rect=[0, 0.08, 1, 1])
fig.savefig(GRAPHS / "p3_fig5_pepbi_concordance.png", dpi=200)
plt.close(fig)

# ----------------------------- Tables --------------------------------------

def write_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

write_csv(OUT / "optimizer_seed_metrics.csv", RUNS)
write_csv(OUT / "optimizer_paired_differences.csv", paired)
write_csv(OUT / "p3_affinity_association_metrics.csv", association_rows)
write_csv(OUT / "p3_position_validation_metrics.csv", position_rows)
write_csv(OUT / "p3_grouped_model_metrics.csv", [
    {"model": "Training-fold mean", "mae": mae[0], "rmse": rmse[0],
     "spearman_rho": 0.045870601293509276,
     "direction_accuracy": directions[0], "coverage": 1.0},
    {"model": "Ridge", "mae": mae[1], "rmse": rmse[1],
     "spearman_rho": 0.10874027339042665,
     "direction_accuracy": directions[1], "coverage": 0.9098712446351931},
])
write_csv(OUT / "p3_pepbi_metrics.csv", [
    {**r, "independent_binding_groups": 11,
     "target": "calculated ΔG from KD"} for r in pepbi_rows
])

# Evidence availability / provenance register. PASS means the script ran,
# never that the scientific hypothesis succeeded.
evidence = [
    {"area": "PPIKB preprocessing", "status": "completed",
     "result": "693 single-substitution pairs; 123 target sequences",
     "source": "prototype3/scripts/ppikb_kd_benchmark.py"},
    {"area": "Native C score vs mutation affinity", "status": "completed",
     "result": "rho=-0.1071; 95% CI [-0.2022,-0.0176]; direction accuracy=42.49%",
     "source": PAIR_DATA},
    {"area": "Exact 21-aa affinity subset", "status": "completed exploratory",
     "result": "129 rows; 67 peptides; 7 eligible contexts; 9 informative pairs",
     "source": EXACT21},
    {"area": "Original ARISE position validation", "status": "completed retrospective",
     "result": "positions 1-3 all-pairs p=0.9218; no reliable primary advantage",
     "source": "prototype3/results/arise_experimental_position_audit/"},
    {"area": "PEPBI strict matched", "status": "completed",
     "result": "11 groups; native concordance 0.3762; Rosetta 0.3927",
     "source": "prototype3/scripts/pepbi_strict_matched_benchmark.py"},
    {"area": "PEPBI receptor-context invariance", "status": "completed",
     "result": "4 contexts; 24 pairs; 0 score changes; 7 experimental differences >0.5 kcal/mol",
     "source": "prototype3/scripts/pepbi_target_mutation_invariance.py"},
    {"area": "Structural mutation benchmark", "status": "carried forward",
     "result": "18 mapped pairs; rho=0.0588; direction accuracy=50%; not rerun in latest batch",
     "source": STRUCTURAL},
    {"area": "Grouped Ridge versus training-fold mean", "status": "rerun completed",
     "result": "grouped target/peptide/exact-pair overlap=0 across folds; Ridge worse on MAE",
     "source": GROUPED},
    {"area": "New adaptive-vs-baseline optimizer", "status": "completed exploratory",
     "result": "six runs; three matched seeds; 160 records per run",
     "source": OPT_ROOT},
    {"area": "New mutation-position proxy audit", "status": "completed exploratory",
     "result": "proxy positions 14,3,8; effect=-0.0828; permutation p=0.5495",
     "source": MUTATION_AUDIT},
    {"area": "Exact generated-candidate affinity validation", "status": "not estimable",
     "result": "zero exact generated-sequence matches to the measured pKd dataset",
     "source": FINAL_AUDIT},
    {"area": "PEPBI overlap for new optimizer candidates", "status": "not estimable",
     "result": "zero exact sequence matches in the workbook audit",
     "source": FINAL_AUDIT},
    {"area": "New candidate structural scoring", "status": "not estimable",
     "result": "zero candidate-run complex PDB files discovered; 0 new scoring attempts",
     "source": FINAL_AUDIT},
    {"area": "New run's explicitly logged ARISE hotspots", "status": "not recovered",
     "result": "positions 14,3,8 are a mutation-frequency proxy, not recorded ARISE selections",
     "source": MUTATION_AUDIT},
]
write_csv(OUT / "P3_EVIDENCE_REGISTER.csv", evidence)

# ----------------------------- Data logbook -------------------------------

data_path = OUT / "P3_DATA_LOGBOOK_ENTRY.md"
data_text = f"""# P3 Data Logbook Entry — Silicon Virovore

**Work date:** 26 October 2026 (America/New_York).  
**Scope:** P3 external retrospective validation plus a new matched-seed optimizer comparison.  
**Purpose of this entry:** preserve the experimental setup, numerical outputs, provenance, and missing-data status without interpreting or omitting negative results.

## 1. Research questions

1. Under matched run budgets, does the adaptive search policy attain higher values of the pipeline's internal fitness objective than a uniform-reset baseline?
2. Do selected mutation positions correspond to positions with larger experimentally measured affinity effects?
3. Does the native sequence score track experimentally measured affinity changes?
4. Do independent PEPBI and structural-scoring benchmarks support receptor-specific binding-affinity prediction?
5. Does a separate sequence-only Ridge predictor improve on a training-fold mean when validation prevents target, peptide, and exact-pair overlap?

These are separate questions. A better internal fitness value is not itself evidence of higher experimental affinity.

## 2. Optimizer comparison: experimental configuration

Six pipeline runs were completed in `prototype3/results/p3_paired_arm_20261027T022810Z/`.

- Matched seeds: **2601, 2602, 2603**.
- Candidate budget: **16 candidates × 10 generations = 160 records per run**.
- Shared explicit settings: starting sequence `MKLAVFALLVFFAGSSDLIRR` (21 residues), mutation rate 0.05, elite count 2, tournament size 3.
- Adaptive arm: `hotspot_entropy_substitution_epistasis_exploratory_v2`, under adaptive configuration.
- Baseline arm: `uniform_random_reset`, under baseline configuration.
- Both arms used the same seed values and recorded run output files under their respective run directories.

**Design caveat:** this compares two *overall policy configurations*. The reset baseline returns to the initial seed, while the adaptive policy follows an evolving parent. The comparison therefore does not isolate mutation guidance as the only causal difference.

### Optimizer results

| Arm | Seed | Candidate records | Unique sequences | Mean fitness over records | Best fitness | Final-generation mean |
|---|---:|---:|---:|---:|---:|---:|
"""
for r in RUNS:
    data_text += (
        f"| {r['arm']} | {r['seed']} | {r['records']} | {r['unique']} | "
        f"{r['mean']:.5f} | {r['best']:.5f} | {r['final']:.5f} |\n"
    )

data_text += f"""
Across the three seeds:

| Summary | Baseline | Adaptive | Adaptive minus baseline |
|---|---:|---:|---:|
| Mean fitness over all records, averaged by seed | {arms['baseline']['mean_all_candidate_fitness']:.5f} | {arms['adaptive']['mean_all_candidate_fitness']:.5f} | {mean_all_diff:+.5f} |
| Best fitness per run, averaged over seeds | {arms['baseline']['mean_best_fitness']:.5f} | {arms['adaptive']['mean_best_fitness']:.5f} | {mean_best_diff:+.5f} |
| Final-generation mean, averaged over seeds | {arms['baseline']['mean_final_generation_fitness']:.5f} | {arms['adaptive']['mean_final_generation_fitness']:.5f} | {mean_final_diff:+.5f} |
| Unique sequences per run, averaged over seeds | {arms['baseline']['mean_unique_sequences_per_seed']:.2f} | {arms['adaptive']['mean_unique_sequences_per_seed']:.2f} | {unique_diff:+.2f} |

All three matched seeds favored the adaptive policy on final-generation mean fitness. Three seeds are an exploratory sample; these descriptive results alone do not establish generalizable statistical superiority.

![Figure 1: Matched-seed optimizer fitness](graphs/p3_fig1_matched_optimizer_fitness.png)

Per-seed paired differences are saved in `optimizer_paired_differences.csv`; the six source runs and summary values are in `optimizer_seed_metrics.csv`.

## 3. PPIKB/native-score retrospective mutation benchmark

Input evidence: `prototype3/results/native_score_all_mutation_pairs/scored_mutation_pairs.csv`.

- 693 matched single-substitution pairs.
- 621 unique peptides; 123 target sequences.
- Peptide lengths represented: 9–20 amino acids.
- 466 informative pairs had absolute experimental ΔpKd at least 0.3.

| Measure | Result |
|---|---:|
| Spearman correlation, Δnative score vs ΔpKd | −0.1071 |
| 95% cluster-bootstrap CI for that correlation | [−0.2022, −0.0176] |
| Informative-pair direction accuracy, ties counted half | 42.49% |
| 95% CI for direction accuracy | [36.55%, 47.78%] |
| Spearman correlation, absolute score change vs absolute pKd change | −0.0418 |
| 95% CI for absolute-change correlation | [−0.1262, +0.0519] |
| Position-profile association | ρ = −0.1354; shuffle p = 0.223 |

The native C score is a sequence-level heuristic with `target=NULL` in the scoring call. This test is not a target-aware docking or binding-affinity calculation. It measures retrospective score/affinity association on external peptide–protein data, not direct affinity outcomes for newly generated candidates.

![Figure 2: Association estimates](graphs/p3_fig2_affinity_association_estimates.png)

## 4. Exact 21-amino-acid benchmark

Source directory: `{EXACT21}`.

- 129 eligible measurement rows; 67 distinct peptides.
- 73 target/assay/publication contexts.
- Only 7 contexts had at least 3 peptides with variable score and affinity.
- 20 single-substitution context rows; 9 informative pairs; 7 unique peptide pairs and 5 mutation positions.

Results:

- Mean within-context Spearman ρ: **−0.4872**.
- Group-bootstrap 95% CI: **[−0.7714, −0.0301]**.
- Informative-pair direction accuracy: **1/9 = 11.11%**.
- In the recorded analysis, 8 of 9 informative rows showed experimental affinity increasing while the native score decreased.

This is an exact-length subset relevant to the optimizer's 21-residue seed, but its informative-pair count is very small and the group estimate is exploratory.

## 5. ARISE mutation-position validation against experimental data

Source: `prototype3/results/arise_experimental_position_audit/` and the prior P3 report.

The primary retrospective test compared positions 1–3 with random contiguous windows of equal width. It did **not** find a reliable primary advantage.

| Test | Context-adjusted effect | 95% CI | Two-sided permutation p |
|---|---:|---:|---:|
| Positions 1–3, all pairs (693) | −0.0115 | [−0.1690, +0.1543] | 0.9218 |
| Previous positions 2–3 sensitivity test | −0.0184 | [−0.2026, +0.1684] | 0.8454 |
| Positions 1–3, filtered informative subset (466 pairs) | +0.1716 | [+0.0148, +0.3419] | 0.0528 |

The informative-only result is secondary and exploratory because it filters to larger observed affinity changes; it is not the primary all-pairs outcome. These peptides are not the same as the HERV-K Env search candidates.

![Figure 3: Mutation-position effects](graphs/p3_fig3_mutation_position_effects.png)

## 6. New adaptive-observed mutation-position proxy

The newer optimizer artifacts did not expose explicitly logged selected ARISE hotspot positions. A separate audit extracted mutation events from parent/child sequence records and chose the three most frequently mutated adaptive positions that were represented in the external pair benchmark: **14, 3 and 8**.

These are labeled as an *observed-mutation proxy*, not actual logged ARISE-selected hotspots.

- Adaptive parent/child mutation events extracted: **510**.
- Baseline parent/child mutation events extracted: **480**.
- Adaptive events at positions 14, 3 or 8: **19.61%**.
- Baseline events at those same positions: **12.08%**.
- External pairs inside/outside the proxy positions: **179 / 514**.
- Context-adjusted external effect, inside minus outside: **−0.0828**.
- Equal-size random-position permutation p: **0.5495**.

The test did not show that the proxy positions have larger external affinity effects than other positions. The same external positional statistic was used for both arms; the arm-specific measurement is the fraction of each arm's mutation events falling inside the proxy positions.

Full source: `{MUTATION_AUDIT}`.

## 7. PEPBI matched validation

Input workbook: `prototype3/data/raw/PEPBI.xlsx`.

For the strict matched analysis of calculated ΔG from K_D:

| Predictor | Mean group concordance | Independent binding groups | 95% group-bootstrap CI |
|---|---:|---:|---:|
| Native sequence score | 0.3762 | 11 | [0.2110, 0.5315] |
| Rosetta dG | 0.3927 | 11 | [0.2345, 0.5475] |

The point estimates are below the 0.50 chance reference, but both intervals include 0.50. The separately reported experimental-ΔG comparison had only one independent binding group; its zero concordance estimate is not a reliable predictive-quality estimate.

The receptor-context invariance control examined 4 matched peptide contexts and 24 receptor-variant pairs. The native score changed in **0/24** pairs, while **7** pairs had experimental ΔG differences greater than 0.5 kcal/mol. That is expected for a score that depends on the peptide sequence alone; it limits use as a receptor-specific predictor.

The new optimizer-candidate audit found **0 exact sequence matches** between its generated sequences and the PEPBI workbook. Consequently, the current PEPBI metrics are retrospective validation of the fixed scoring method, not separate adaptive-versus-baseline prediction accuracies.

![Figure 5: PEPBI concordance](graphs/p3_fig5_pepbi_concordance.png)

## 8. Structural benchmark

Existing benchmark source: `{STRUCTURAL}`. This result was **carried forward**, not rerun as a fresh optimizer-arm comparison.

- 18 mapped mutation pairs, 7 targets and 28 unique complex-chain combinations.
- Overall Spearman ρ between structural-score change and experimental ΔpKd: **+0.0588**, exploratory p = **0.8167**.
- Direction accuracy: **50%**.
- Informative subset: 14 pairs, ρ = **+0.0593**, p = **0.8403**, 50% direction accuracy.
- Same-observed-receptor-sequence subset: 10 pairs across 5 targets, ρ = **+0.2849**, p = **0.4250**, 60% direction accuracy; bootstrap intervals were wide and included no association.

A structural pose sanity check on the MDM2/p53 1YCR complex found 85 receptor Cα correspondences, 13 peptide Cα matches after alignment, and an aligned docked-pose peptide RMSD of 9.965 Å. The score preferred the poor docked pose over the crystal reference (reported scores 0.97751 vs 0.96537). This raises a validation concern about pose ranking.

The new optimizer-run folders contained **0 candidate-run PDB files** available to score, so there were no new structural scores for candidates from either arm.

## 9. Grouped sequence-only Ridge versus training-fold mean

Fresh run: `{GROUPED}`.

Validation used 5 GroupKFold folds built from 111 connected components across 693 pairs and 123 target sequences. The grouping links rows sharing a target sequence, peptide sequence or exact unordered peptide pair. The split audit reports no target, peptide or exact-pair overlap across the folds.

| Model | MAE | RMSE | Spearman ρ | Direction accuracy | Direction coverage |
|---|---:|---:|---:|---:|---:|
| Training-fold mean | 0.8182 | 1.3656 | 0.0459 | 62.66% | 100.00% |
| Ridge | 0.8643 | 1.3957 | 0.1087 | 58.96% | 90.99% |

Ridge had a larger Spearman estimate, but worse MAE, worse RMSE and lower informative-direction accuracy. The reported paired bootstrap CI for the Ridge-minus-mean rho difference was [−0.0637, +0.1726], including zero. This is a separate sequence-only model comparison, not evidence that the optimizer improves.

![Figure 4: Grouped Ridge vs training-fold mean](graphs/p3_fig4_grouped_ridge_vs_mean.png)

## 10. P3 batch execution status

The prior retrospective P3 batch reported successful execution for all ten listed steps:

- PPIKB preprocessing
- Grouped sequence validation
- Native-score mutation-pair benchmark
- Exact 21-aa benchmark
- ARISE position validation
- Corrected PEPBI benchmark
- Strict PEPBI matched benchmark
- PEPBI receptor-context invariance
- Verified single-substitution PEPBI audit
- Exploratory native-score PEPBI audit

**PASS refers to script execution, not to confirmation of the scientific hypothesis.** The structural mutation benchmark was separately carried forward from its existing results.

## 11. Data availability and exclusions

The candidate-specific tests could not yield direct affinity or new structural accuracy metrics for the six optimizer runs:

- Exact generated-sequence matches to the measured pKd dataset: **0**.
- Exact generated-sequence matches to the PEPBI workbook: **0**.
- Candidate-run complex PDB files available for fresh structural scoring: **0**.
- Explicit adaptive-selected hotspot positions recoverable from the new run artifacts: **0 seeds**.

No measured affinity or structural result has been imputed for a candidate without evidence. The proxy-position analysis is explicitly marked as such. The original data and earlier output directories were preserved.

## 12. Provenance and saved outputs

See `P3_EVIDENCE_REGISTER.csv` for the status and source of each analysis.

Important paths:

- Prior consolidated retrospective report: `{P2_REPORT}`
- Matched optimizer runs: `{OPT_ROOT}`
- Corrected mutation-position audit: `{MUTATION_AUDIT}`
- Consolidated new validation audit: `{FINAL_AUDIT}`
- Grouped model metrics: `{GROUPED}`
- Existing structural benchmark: `{STRUCTURAL}`

The five graphs are stored under `graphs/`. Numeric plot inputs and seed values are preserved as CSV files in this logbook folder.
"""
data_path.write_text(data_text, encoding="utf-8")

# ----------------------------- Analysis logbook ----------------------------

analysis_path = OUT / "P3_ANALYSIS_LOGBOOK_ENTRY.md"

analysis_text = f"""# P3 Analysis Logbook Entry — Interpretation, Limitations and Next Decisions

**Work date:** 26 October 2026 (America/New_York).  
**Project:** Silicon Virovore — computational peptide optimization and retrospective validation.  
**Purpose:** interpret the combined P3 optimizer and external-validation evidence. This entry is separate from the factual data record.

## 1. Overall conclusion

The best-supported positive result is about **internal optimization**: under the tested configurations, the adaptive policy produced higher values of the pipeline's native fitness objective than the uniform-reset baseline in all three matched seeds. Mean final-generation fitness was {arms['adaptive']['mean_final_generation_fitness']:.4f} for adaptive and {arms['baseline']['mean_final_generation_fitness']:.4f} for baseline, a descriptive difference of {mean_final_diff:+.4f}. Mean best score per run was {arms['adaptive']['mean_best_fitness']:.4f} versus {arms['baseline']['mean_best_fitness']:.4f}.

That result is encouraging for the search policy but does **not** establish improved experimental binding affinity. The native C score calls the fitness function with `target=NULL`; it is a sequence-property heuristic, not a validated receptor-specific affinity function. The external P3 tests, overall, do not establish reliable prediction of affinity changes or HERV-K Env binding at ARG203/ARG205.

The defensible research narrative is therefore two-part:

1. Adaptive search behavior was promising on the internal optimization objective in the three matched runs.
2. External validation remains a separate and unresolved problem; most direct retrospective score/affinity checks were weak, negative, inconclusive or too small, and the newly generated candidates had no exact measured-affinity matches.

## 2. Hypothesis-by-hypothesis assessment

### Hypothesis A: Adaptive search improves the internal objective

**Observation:** all three matched seeds had higher adaptive final-generation mean scores. The mean paired increase was {mean_final_diff:+.4f}; the difference in mean best scores was {mean_best_diff:+.4f}. Adaptive averaged {arms['adaptive']['mean_unique_sequences_per_seed']:.2f} unique sequences per run, while baseline averaged {arms['baseline']['mean_unique_sequences_per_seed']:.2f}.

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
"""
analysis_path.write_text(analysis_text, encoding="utf-8")

# ----------------------------- Evidence index ------------------------------

index_path = OUT / "README_P3_LOGBOOK_BUNDLE.md"
index_path.write_text(f"""# P3 Logbook Bundle

Contents:

- `P3_DATA_LOGBOOK_ENTRY.md` — factual record of protocols, metrics, datasets, run outputs and data gaps.
- `P3_ANALYSIS_LOGBOOK_ENTRY.md` — hypothesis-by-hypothesis interpretation, strongest supported ISEF claim, and remaining gaps.
- `P3_EVIDENCE_REGISTER.csv` — source/provenance and completion status for each evidence item.
- `optimizer_seed_metrics.csv` — six optimizer run results.
- `optimizer_paired_differences.csv` — matched-seed differences.
- `p3_affinity_association_metrics.csv` — values plotted in Figure 2.
- `p3_position_validation_metrics.csv` — values plotted in Figure 3.
- `p3_grouped_model_metrics.csv` — grouped prediction metrics.
- `p3_pepbi_metrics.csv` — strict matched concordance metrics.
- `graphs/` — five PNG figures used in both entries.

These documents report the supplied P3 metrics transparently, including negative and unavailable results. They do not replace raw data or claim that all tests succeeded scientifically.
""", encoding="utf-8")

# ----------------------------- Bundle ZIP ----------------------------------

zip_path = OUT / "P3_LOGBOOK_BUNDLE.zip"
with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p != zip_path:
            z.write(p, p.relative_to(OUT))

print("P3 logbook bundle created successfully.")
print(f"Data entry:     {data_path}")
print(f"Analysis entry: {analysis_path}")
print(f"Evidence index: {index_path}")
print(f"Evidence CSV:   {OUT / 'P3_EVIDENCE_REGISTER.csv'}")
print(f"Graphs:         {GRAPHS}")
print(f"ZIP bundle:     {zip_path}")
print()
print("Final-generation mean difference (adaptive - baseline):", round(mean_final_diff, 5))
print("Mean best-score difference (adaptive - baseline):", round(mean_best_diff, 5))
print("Mean unique sequences per run difference:", round(unique_diff, 2))
