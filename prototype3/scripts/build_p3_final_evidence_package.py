from pathlib import Path
from datetime import datetime
import csv
import hashlib
import json
import shutil
import sys

ROOT = Path.cwd().resolve()
if not (ROOT / "prototype3").is_dir():
    raise SystemExit("ERROR: Run this script from the Silicon-Virovore project root.")

STAMP = datetime.now().strftime("%Y%m%dT%H%M%S")
OUT = ROOT / "prototype3/logbook" / f"p3_final_defensible_package_{STAMP}"

FULL = ROOT / "prototype3/logbook/p3_comparison_full_audit_20261027T172008"
DIRECT = ROOT / "prototype3/logbook/p3_direct_armwise_20261027T204305Z"
MUTATION = ROOT / "prototype3/logbook/p3_mutation_position_audit_20261027T205650Z"

BROAD = ROOT / "prototype3/results/native_score_all_mutation_pairs/metrics.json"
AA21 = ROOT / "prototype3/results/p3_21aa_affinity_audit/native_c_score_21aa/metrics.json"
POSITIONS = ROOT / "prototype3/results/arise_experimental_position_audit/metrics.json"
GROUPED = ROOT / "prototype3/results/sequence_only_grouped_target_validation"
STRUCTURAL = ROOT / (
    "prototype3/results/experimental_mutation_benchmark/"
    "structural_scoring_20261026T223525Z/benchmark_metrics.json"
)
PEPBI = ROOT / "prototype3/results/pepbi_strict_matched"
PAIRED = ROOT / "prototype3/results/p3_paired_arm_20261027T022810Z"

if OUT.exists():
    raise SystemExit(f"Output already exists: {OUT}")

for sub in [
    "graphs",
    "poster_graphs",
    "supporting_data",
    "source_reports",
    "source_archive",
    "unresolved_evidence",
]:
    (OUT / sub).mkdir(parents=True, exist_ok=True)

copy_log = []

def copy_file(src, rel):
    src = Path(src)
    if not src.is_file():
        copy_log.append({
            "source": str(src),
            "destination": "",
            "status": "NOT_FOUND",
        })
        return None

    dest = OUT / rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    copy_log.append({
        "source": str(src),
        "destination": str(dest.relative_to(OUT)),
        "status": "COPIED",
    })
    return dest

def copy_folder_contents(src, rel):
    src = Path(src)
    if not src.is_dir():
        copy_log.append({
            "source": str(src),
            "destination": "",
            "status": "DIRECTORY_NOT_FOUND",
        })
        return

    dest = OUT / rel
    dest.mkdir(parents=True, exist_ok=True)
    for p in sorted(src.rglob("*")):
        if not p.is_file():
            continue
        try:
            subpath = p.relative_to(src)
        except ValueError:
            continue
        target = dest / subpath
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, target)
        copy_log.append({
            "source": str(p),
            "destination": str(target.relative_to(OUT)),
            "status": "COPIED",
        })

# Preserve the original audit without modifying its reports, source snapshots,
# or original figures. The curated final figures are stored separately.
if FULL.is_dir():
    archive_dest = OUT / "source_archive/original_full_comparison_audit"
    shutil.copytree(FULL, archive_dest)
    copy_log.append({
        "source": str(FULL),
        "destination": str(archive_dest.relative_to(OUT)),
        "status": "ARCHIVED_UNCHANGED",
    })

copy_file(
    DIRECT / "P3_DIRECT_ARMWISE_REPORT.md",
    "source_reports/P3_DIRECT_ARMWISE_REPORT.md",
)
copy_folder_contents(DIRECT, "source_reports/direct_armwise_files")
copy_folder_contents(MUTATION, "supporting_data/candidate_mutation_position_audit")
copy_folder_contents(PAIRED / "paired_analysis", "supporting_data/existing_paired_optimizer_summary")
copy_folder_contents(PEPBI, "supporting_data/pepbi_strict_matched_source")

# Copy the main extraction tables from the previous audit.
if FULL.is_dir():
    for name in [
        "source_file_inventory.csv",
        "candidate_records_extracted.csv",
        "candidate_extraction_summary.csv",
        "paired_optimizer_seed_metrics.csv",
        "external_validation_summary.csv",
        "ridge_vs_mean_summary.csv",
        "pepbi_strict_matched_summary.csv",
        "pepbi_exact_match_audit_provisional.csv",
        "README.md",
        "logbook_DATA_RESULTS.md",
        "logbook_ANALYSIS.md",
    ]:
        copy_file(FULL / name, Path("supporting_data/original_audit_tables") / name)

# Copy the exact metrics and important underlying benchmark data.
for src, rel in [
    (BROAD, "supporting_data/metrics/native_c_score_broad_metrics.json"),
    (AA21, "supporting_data/metrics/native_c_score_exact_21aa_metrics.json"),
    (POSITIONS, "supporting_data/metrics/experimental_position_audit_metrics.json"),
    (STRUCTURAL, "supporting_data/metrics/structural_scoring_metrics.json"),
]:
    copy_file(src, rel)

# Save both grouped-validation metric files if present. The earlier audit
# contains a path discrepancy, so do not silently overwrite one with the other.
copy_file(
    GROUPED / "run_20261027T164512Z/metrics.json",
    "supporting_data/metrics/grouped_validation_run_164512Z_metrics.json",
)
copy_file(
    GROUPED / "run_20261027T205129Z/metrics.json",
    "supporting_data/metrics/grouped_validation_run_205129Z_metrics.json",
)

copy_file(
    PAIRED / "paired_analysis/PAIRED_OPTIMIZER_SUMMARY.md",
    "source_reports/PAIRED_OPTIMIZER_SUMMARY.md",
)
copy_file(
    PAIRED / "paired_analysis/optimizer_run_summary.csv",
    "supporting_data/existing_paired_optimizer_summary/optimizer_run_summary.csv",
)
copy_file(
    PAIRED / "paired_analysis/paired_seed_differences.csv",
    "supporting_data/existing_paired_optimizer_summary/paired_seed_differences.csv",
)

# Include the existing external mutation-pair input if available.
copy_file(
    ROOT / (
        "prototype3/results/sequence_only_mutation_baseline/"
        "run_20261026T234823Z/single_substitution_pairs.csv"
    ),
    "supporting_data/underlying_data/single_substitution_pairs.csv",
)

# Clearly quarantine the unresolved PEPBI claim; do not use it for final figures.
copy_file(
    FULL / "pepbi_exact_match_audit_provisional.csv",
    "unresolved_evidence/pepbi_exact_match_audit_provisional.csv",
)

(OUT / "unresolved_evidence/PEPBI_MATCH_CONFLICT.md").write_text(
    """# Unresolved PEPBI exact-sequence match discrepancy

The prior direct armwise scan reported zero exact candidate-sequence matches
for both arms: baseline 0/275 and adaptive 0/450.

A separately supplied Numbers worksheet was summarized as baseline 0/275 and
adaptive 125/450. The 125 adaptive matches were not accompanied by a sequence-
level match list in the audit output.

The reported adaptive match rate of 27.78% is therefore provisional. Sequence
overlap alone would not prove stronger binding even if the matches were
reproduced. Each match must be linked to the exact source record, experimental
affinity measurement, and receptor context.

The disputed count is intentionally excluded from the curated final graphs
and is not treated as evidence of adaptive biological superiority.
""",
    encoding="utf-8",
)

(OUT / "unresolved_evidence/EXACT_21AA_SUBSET_NOTE.md").write_text(
    """# Exact-21-amino-acid native-score audit: subset caution

The reported mean within-context Spearman correlation is -0.487 across only
seven eligible contexts, with a group-bootstrap 95% interval of approximately
[-0.771, -0.030]. Only nine same-context single-substitution pairs met the
informative-change threshold, with direction accuracy of 1/9 (11.1%).

The same source also reports an exploratory correlation of +0.678 on a
different subset/summary. Those two summaries should not be treated as
equivalent. The positive subset result is not used in the curated poster
figures until the pair-selection and sign conventions are reconciled.
""",
    encoding="utf-8",
)

(OUT / "unresolved_evidence/GROUPED_VALIDATION_PROVENANCE_NOTE.md").write_text(
    """# Grouped-validation provenance note

The original data-results entry names
`run_20261027T205129Z`, while the displayed metrics in the earlier audit were
read from `run_20261027T164512Z`.

This package copies both metric files if available. The grouped-validation
numbers in the curated logbook match the metrics explicitly displayed from
`run_20261027T164512Z`. Verify the source path before making a final
submission. The packaging script does not rerun the model or cross-validation.
""",
    encoding="utf-8",
)

# ---------- Curated logbook: data and results ----------
DATA = f"""# P3 Logbook Entry — Data and Results

**Date:** 27 October 2026  
**Experiment:** Consolidated P3 optimization and retrospective validation audit  
**Purpose:** Preserve the defensible existing results in one reproducible evidence package.  
**Package:** `{OUT}`

## Design and objective

P3 addresses two distinct questions:

1. Does adaptive guidance improve the optimizer's internal objective relative to the matched uniform baseline?
2. Do the internal score, mutation-position signal, or structural-scoring method correspond to independently measured binding-affinity outcomes?

These questions are not interchangeable. An optimizer can improve its assigned score without improving real binding.

This package reuses the existing paired run, candidate records, metric files, logs, and completed validation results. **It does not launch an optimizer, add seeds, perform another seed-by-seed comparison, or rerun the validation scripts.** The figures are summaries of existing results, not new biological experiments.

## A. Existing optimizer and candidate inventory

The existing paired run contains 480 candidate records per arm.

| Measure | Uniform baseline | Adaptive guidance |
|---|---:|---:|
| Candidate records | 480 | 480 |
| Unique candidate sequences | 275 | 450 |
| Mean native score | 0.65224 | 0.74933 |
| Best native score | 0.71454 | 0.84553 |

The armwise score difference is approximately +0.0971 for adaptive minus baseline. The existing paired-run summary reports the same direction for mean score, best score, and final-generation mean across its three matched runs. The underlying run and original summary files are preserved in `supporting_data/` and `source_reports/`.

These are internal objective scores. They are not measured binding affinities.

## B. Mutation behavior audit

The existing parent-child mutation audit reported:

| Measure | Baseline | Adaptive |
|---|---:|---:|
| Parent-child records checked | 480 | 480 |
| Records with an observed sequence change | 480 | 450 |
| Substitution events counted | 480 | 510 |
| Unique parent-child transitions | 275 | 473 |
| Mutation-count consistency | 480/480 | 480/480 |

The metadata marked 432/480 adaptive records (90%) as hotspot-guided; the remaining 48 were not marked guided. All 480 baseline records used the `uniform_matched` mode.

The selected-position frequency differences are descriptive. The presence of a guidance flag does not demonstrate that a particular mutation was caused by a hotspot or that it improves binding.

## C. Broad native C-score versus experimental affinity

The retrospective benchmark included 693 matched single-substitution pairs, 621 unique peptides, peptide lengths of 9–20 amino acids, and 466 informative pairs.

- Spearman correlation between native score change and measured pKd change: **-0.1071**.
- Cluster-bootstrap 95% confidence interval: **[-0.2022, -0.0175]**.
- Direction accuracy on the reported informative comparison: **42.49%**, with a reported 95% interval of approximately **[36.55%, 47.77%]**.

The result does not support treating the native sequence-only score as a reliable predictor of experimental affinity changes. The native score has no receptor input; this benchmark tests the scoring function, not whether the adaptive optimizer produces experimentally better candidates.

## D. Exact-21-amino-acid native-score audit

The filtered audit retained 129 rows, 67 distinct peptide sequences, and 73 measurement contexts. Only seven contexts had enough distinct peptides and variability for the reported within-context correlation analysis.

- Mean within-context Spearman correlation: **-0.4872**.
- Group-bootstrap 95% interval: **[-0.7714, -0.0301]**.
- Informative single-substitution pairs: **9**.
- Direction accuracy on those pairs: **1/9 (11.11%)**.

The small, retrospective sample limits generalization. The source also reports a positive exploratory correlation on a different subset. That result is documented in `unresolved_evidence/EXACT_21AA_SUBSET_NOTE.md` and is not used as a headline finding until the differing subsets and sign conventions are reconciled.

## E. Mutation-position external validation

Across all 693 external mutation pairs, the positions 1–3 comparison reported a centered effect of **-0.0115**, with a target-cluster bootstrap interval of **[-0.1690, +0.1543]** and two-sided permutation p = **0.9218**. This does not show a general position 1–3 effect across all pairs.

Within the 466 informative pairs, the centered difference was **+0.1716**, interval **[+0.0148, +0.3419]**. The one-sided permutation p-value was **0.0200**, and the two-sided value was **0.0528**. This is a narrower exploratory signal and is not independent confirmation of a universal hotspot effect. The data are from a retrospective peptide benchmark, not direct HERV-K Env-specific experimental testing.

## F. Grouped sequence-only Ridge comparison

The recorded grouped validation used 693 pairs, 123 targets, 111 leakage components, and five folds. The audit reported zero target overlap, zero exact-pair overlap, zero peptide-sequence overlap, and one out-of-fold prediction per row.

| Metric | Training-fold mean | Ridge |
|---|---:|---:|
| MAE (lower is better) | 0.8182 | 0.8643 |
| RMSE (lower is better) | 1.3656 | 1.3957 |
| Spearman correlation | 0.0459 | 0.1087 |
| Direction accuracy | 62.66% | 58.96% |

Ridge ranked changes somewhat better by correlation, but the mean baseline had lower error and higher direction accuracy. This sequence-only benchmark is separate from the native C-score and the optimizer. The provenance path discrepancy is documented in `unresolved_evidence/GROUPED_VALIDATION_PROVENANCE_NOTE.md`.

## G. PEPBI strict-matched benchmark

The strict-matched audit reported 290 eligible rows, 69 distinct sequences, 18 matched contexts, and 382 pairwise comparisons.

Across 11 independent groups, mean concordance was:

- Native score: **0.3762**, 95% interval **[0.2110, 0.5315]**, reported one-sided sign-test p = **0.8867**.
- Rosetta delta G: **0.3927**, 95% interval **[0.2345, 0.5476]**, reported one-sided sign-test p = **0.7256**.

The target-reported experimental-delta-G subset had only one independent group, which is insufficient for a dependable general conclusion. These retrospective concordance values do not establish that either score reliably ranks the generated candidates.

The exact candidate-match count remains unresolved: one audit reported zero matches in both arms, while a separate Numbers worksheet was summarized as 125 adaptive matches and zero baseline matches. The latter count is not treated as confirmed evidence because a sequence-level match list and linked affinity records were not supplied in the audit. See `unresolved_evidence/PEPBI_MATCH_CONFLICT.md`.

## H. Structural-scoring benchmark

Of 693 considered experimental mutation pairs, only **18** were successfully mapped and scored.

- Overall Spearman correlation: **0.0588**, p = **0.8167**.
- Overall direction accuracy: **50%**.
- Informative subset: 14 pairs, correlation **0.0593**, p = **0.8403**.
- Same-receptor subset: 10 pairs across five targets, correlation **0.2849**, p = **0.425**, direction accuracy **60%**.

The pose sanity check reported a crystal-pose score of 0.96537 and a docked-pose score of 0.97751; under the recorded score orientation, the unfavorable docked pose received the preferred score. This indicates that the scoring setup needs more validation. The small mapped sample is not sufficient to claim that the structural score reliably predicts affinity.

## I. Receptor-context control

The recorded sequence-score receptor-context audit compared 24 matched receptor-variant pairs. The sequence-only native score was unchanged for zero pairs, while seven pairs had reported experimental energy differences greater than 0.5 kcal/mol.

This illustrates that a score without receptor input cannot directly encode receptor-specific differences. It is not a direct test of adaptive versus baseline binding.

## Files and figures

- `graphs/`: curated, fully labeled evidence plots.
- `poster_graphs/`: three suggested core poster figures.
- `supporting_data/`: source metrics, existing candidate tables, and mutation-audit outputs.
- `source_reports/`: existing report and execution artifacts.
- `source_archive/`: unchanged copy of the prior full audit.
- `unresolved_evidence/`: conflicts that must not be presented as confirmed findings.
- `evidence_manifest.csv`: checksums and paths for files in this package.
"""

(OUT / "logbook_DATA_RESULTS.md").write_text(DATA, encoding="utf-8")

# ---------- Curated logbook: analysis ----------
ANALYSIS = f"""# P3 Logbook Entry — Analysis and Interpretation

**Date:** 27 October 2026  
**Corresponding data entry:** `logbook_DATA_RESULTS.md`  
**Evidence package:** `{OUT}`

## 1. Main conclusion

The current evidence does **not** establish that the adaptive optimizer's advantage translates into stronger real-world binding.

The evidence supports a narrower result: adaptive guidance improved the optimizer's internal-score outcomes and increased candidate-sequence diversity in the existing paired run. The independent retrospective tests do not validate the native sequence score as a dependable experimental-affinity predictor, and there is no confirmed direct experimental comparison of the generated adaptive and baseline candidates.

## 2. Hypothesis A — Does adaptive guidance improve optimization?

The existing armwise inventory reported a mean native score of 0.74933 for adaptive versus 0.65224 for baseline, and 450 unique adaptive sequences versus 275 unique baseline sequences. The paired-run summary also reported the same direction of internal-score advantage across its three existing matched runs.

This supports the claim that adaptive guidance improved optimization against the current internal objective in the tested runs. The mutation audit also shows that adaptive guidance alters mutation behavior and yields more unique parent-child transitions.

The limitation is that these outcomes are based on the pipeline's objective and sequence inventory. They do not show that the generated candidates bind better. The native score is sequence-only rather than receptor-aware, so improving it does not necessarily improve receptor-specific interactions.

## 3. Hypothesis B — Does the scoring or mutation signal predict measured affinity?

The broad native-score benchmark gives unfavorable evidence: the change-score correlation with measured affinity change was weakly negative (-0.107), and direction accuracy was 42.49%. The exact-21-amino-acid audit was also unfavorable in its small within-context and informative-pair summaries, although the source includes a different exploratory subset correlation that needs reconciliation.

The mutation-position test was essentially null across all 693 pairs. The informative subset showed a possible positive association, but the two-sided permutation p-value was about 0.053 and the subset is exploratory. It cannot establish that the same positions are biologically important for the target used by this optimizer.

The Ridge comparison was mixed. Ridge had better rank correlation than the training-fold mean, but worse MAE, worse RMSE, and lower direction accuracy. The PEPBI strict-matched concordance estimates were limited, with broad intervals and high reported sign-test p-values. The structural benchmark had only 18 scored pairs, near-zero overall correlation, and a pose sanity check in which the score favored the unfavorable docked pose.

Taken together, these tests do not establish a dependable link between the existing sequence/structural score outputs and real binding affinity.

## 4. Candidate-level validation is still missing

A retrospective benchmark of external mutations asks whether a scoring rule tracks measurements in historical data. It is not equivalent to testing the actual sequences generated by both optimizer arms.

The PEPBI exact-match evidence is unresolved. The prior direct audit reported zero matches in both arms, but another worksheet was summarized as 125 adaptive matches. Until those matches are reproduced at the sequence level and linked to valid affinity measurements and receptor contexts, neither number should support a candidate-level binding claim.

The package therefore does not claim that adaptive candidates have better measured affinities, does not label an exact sequence match as a successful binder, and does not turn structural-score differences into experimental results.

## 5. Defensible conclusion for the ISEF report

> Across the existing matched optimizer runs, adaptive guidance improved the internal native-score objective and increased the number of unique candidate sequences relative to the uniform baseline. However, retrospective experimental-affinity benchmarks did not establish that the native sequence score reliably predicts binding affinity, and the structural benchmark was limited by a small mapped sample and a pose-scoring concern. Because candidate-level experimental superiority has not been established, whether the observed optimization advantage translates into stronger real-world binding remains unresolved.

## 6. Recommended poster figures

1. **Internal score by arm:** demonstrates the optimization result only. Label it as an internal score, not binding affinity.
2. **Native score versus experimental affinity changes:** shows the critical validation limitation.
3. **Mutation-position validation:** shows the all-pairs result alongside the narrower informative-subset result, making the uncertainty visible.

The unique-candidate-sequence chart is available as a supplementary figure. The disputed PEPBI match count is deliberately absent from final figures.

## 7. Remaining evidence-integrity checks

Before final submission, reconcile the PEPBI match discrepancy at the sequence-and-measurement-record level, confirm whether the grouped sequence-validation metrics came from the 164512Z or 205129Z run, and audit the differing exact-21-amino-acid subset summaries. Do not rerun the optimizer merely to address these source-provenance issues.

A stronger claim about actual binding would require a valid comparison of the generated candidates using an independently validated receptor-aware method and, ideally, matched experimental binding measurements.
"""

(OUT / "logbook_ANALYSIS.md").write_text(ANALYSIS, encoding="utf-8")

# ---------- Figures based only on existing reported results ----------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:
    raise SystemExit(
        "The existing results were not modified, but matplotlib could not be "
        f"loaded for figure creation: {exc}. Activate the project venv and install "
        "the project's existing plotting dependencies before rerunning this packaging script."
    )

G = OUT / "graphs"

def finish(fig, name, note):
    fig.text(0.02, 0.015, note, ha="left", va="bottom", fontsize=8, wrap=True)
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    path = G / name
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path

# 1. Aggregate scores from existing armwise candidate inventory only.
fig, ax = plt.subplots(figsize=(7.6, 5.0))
labels = ["Uniform baseline", "Adaptive guidance"]
values = [0.6522404876, 0.7493252276]
bars = ax.bar(labels, values)
ax.set_ylabel("Mean native score")
ax.set_title("Existing P3 candidate inventory: internal objective")
ax.set_ylim(0, max(values) * 1.18)
for bar, value in zip(bars, values):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 0.012,
        f"{value:.4f}",
        ha="center",
        va="bottom",
    )
finish(
    fig,
    "01_optimizer_internal_score.png",
    "Existing 480-record-per-arm inventory. Internal sequence score only; not measured binding affinity. No optimizer rerun.",
)

# 2. Sequence diversity; no provisional PEPBI values.
fig, ax = plt.subplots(figsize=(7.6, 5.0))
labels = ["Uniform baseline", "Adaptive guidance"]
values = [275, 450]
bars = ax.bar(labels, values)
ax.set_ylabel("Unique candidate sequences")
ax.set_title("Unique sequences in existing candidate inventories")
ax.set_ylim(0, 500)
for bar, value in zip(bars, values):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        value + 8,
        str(value),
        ha="center",
        va="bottom",
    )
finish(
    fig,
    "02_unique_candidate_sequences.png",
    "Unique sequences from the existing inventory. Sequence diversity does not establish stronger binding.",
)

# 3. Broad native score validation, using its reported CI.
fig, ax = plt.subplots(figsize=(8.2, 4.8))
rho = -0.1071392244
lo, hi = -0.2021964139, -0.0175494898
ax.errorbar(
    [rho], [0],
    xerr=[[rho - lo], [hi - rho]],
    fmt="o",
    capsize=5,
)
ax.axvline(0, linestyle="--")
ax.set_xlim(-0.30, 0.15)
ax.set_ylim(-0.8, 0.8)
ax.set_yticks([0])
ax.set_yticklabels(["Native C-score"])
ax.set_xlabel("Spearman correlation: score change vs. measured pKd change")
ax.set_title("Broad retrospective affinity validation")
ax.text(0.02, 0.78, "n = 693 matched mutation pairs", transform=ax.transAxes)
ax.text(
    0.02, 0.67,
    "Direction accuracy on informative pairs: 42.5%",
    transform=ax.transAxes,
)
finish(
    fig,
    "03_broad_native_score_affinity_validation.png",
    "95% cluster-bootstrap CI shown. This benchmark evaluates the native sequence-only score, not the adaptive optimizer's generated candidates.",
)

# 4. Exact 21-aa within-context summary; omit the conflicting +0.678 subset summary.
fig, ax = plt.subplots(figsize=(8.2, 4.8))
rho = -0.4871987630
lo, hi = -0.7714285714, -0.0300559058
ax.errorbar(
    [rho], [0],
    xerr=[[rho - lo], [hi - rho]],
    fmt="o",
    capsize=5,
)
ax.axvline(0, linestyle="--")
ax.set_xlim(-1, 1)
ax.set_ylim(-0.8, 0.8)
ax.set_yticks([0])
ax.set_yticklabels(["Mean within-context rho"])
ax.set_xlabel("Spearman correlation")
ax.set_title("Exact-21-amino-acid audit")
ax.text(0.02, 0.78, "Only 7 contexts qualified for this analysis", transform=ax.transAxes)
ax.text(
    0.02, 0.67,
    "Informative mutation-direction accuracy: 1/9 (11.1%)",
    transform=ax.transAxes,
)
finish(
    fig,
    "04_exact_21aa_score_validation.png",
    "Group-bootstrap interval shown. Small retrospective sample; a different exploratory subset result remains unreconciled and is excluded from this figure.",
)

# 5. Mutation-position window external validation.
fig, ax = plt.subplots(figsize=(8.6, 5.4))
labels = ["All 693 pairs", "466 informative pairs"]
effects = [-0.0115015304, 0.1715996057]
los = [-0.1690128691, 0.0148125001]
his = [0.1543257170, 0.3418983806]
y = [1, 0]
ax.errorbar(
    effects, y,
    xerr=[
        [effects[i] - los[i] for i in range(2)],
        [his[i] - effects[i] for i in range(2)],
    ],
    fmt="o",
    capsize=5,
)
ax.axvline(0, linestyle="--")
ax.set_yticks(y)
ax.set_yticklabels(labels)
ax.set_xlabel("Centered absolute-affinity-effect difference (positions 1–3 minus comparison)")
ax.set_title("Retrospective mutation-position validation")
ax.text(
    0.02, 0.86,
    "All pairs: two-sided permutation p = 0.922",
    transform=ax.transAxes,
)
ax.text(
    0.02, 0.77,
    "Informative subset: two-sided p = 0.0528; one-sided p = 0.020",
    transform=ax.transAxes,
)
finish(
    fig,
    "05_mutation_position_validation.png",
    "Target-cluster bootstrap intervals shown. The informative subset is exploratory and does not establish a target-specific biological hotspot.",
)

# 6. Grouped model validation: error measures.
fig, ax = plt.subplots(figsize=(7.8, 5.0))
labels = ["Training-fold mean", "Sequence-only Ridge"]
mae = [0.8182125546, 0.8643157903]
bars = ax.bar(labels, mae)
ax.set_ylabel("Mean absolute error (lower is better)")
ax.set_title("Grouped sequence-only prediction benchmark")
ax.set_ylim(0, 1.0)
for bar, value in zip(bars, mae):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        value + 0.015,
        f"{value:.3f}",
        ha="center",
        va="bottom",
    )
finish(
    fig,
    "06_grouped_ridge_vs_mean_mae.png",
    "Ridge rho = 0.109 vs mean rho = 0.046; Ridge direction accuracy = 59.0% vs mean = 62.7%. Different metrics favor different methods.",
)

# 7. PEPBI strict-matched concordance with intervals.
fig, ax = plt.subplots(figsize=(8.2, 5.0))
vals = [0.376197, 0.392711]
los = [0.21097, 0.23445]
his = [0.53149, 0.54755]
y = [1, 0]
ax.errorbar(
    vals, y,
    xerr=[
        [vals[i] - los[i] for i in range(2)],
        [his[i] - vals[i] for i in range(2)],
    ],
    fmt="o",
    capsize=5,
)
ax.set_xlim(0, 1)
ax.set_yticks(y)
ax.set_yticklabels(["Native score", "Rosetta delta G"])
ax.set_xlabel("Mean concordance")
ax.set_title("PEPBI strict-matched retrospective benchmark")
ax.text(0.02, 0.88, "11 independent groups per summary", transform=ax.transAxes)
ax.text(
    0.02, 0.78,
    "Reported one-sided sign tests: p = 0.887 and p = 0.726",
    transform=ax.transAxes,
)
finish(
    fig,
    "07_pepbi_strict_matched_concordance.png",
    "Intervals from the existing strict-matched report. These results do not establish dependable ranking of generated candidates.",
)

# 8. Structural score validation.
fig, ax = plt.subplots(figsize=(8.2, 5.0))
vals = [0.0588235, 0.28485]
y = [1, 0]
ax.plot(vals, y, "o")
ax.axvline(0, linestyle="--")
ax.set_xlim(-0.2, 0.6)
ax.set_yticks(y)
ax.set_yticklabels(["Overall: 18 pairs", "Same receptor: 10 pairs"])
ax.set_xlabel("Spearman correlation: structural score change vs. measured affinity change")
ax.set_title("Preliminary structural-scoring benchmark")
ax.text(0.02, 0.88, "Overall p = 0.817; direction accuracy = 50%", transform=ax.transAxes)
ax.text(0.02, 0.78, "Same-receptor p = 0.425; direction accuracy = 60%", transform=ax.transAxes)
finish(
    fig,
    "08_structural_score_validation.png",
    "Only 18 of 693 considered pairs were successfully mapped and scored. Pose sanity check also raised a scoring-validity concern.",
)

# 9. Contextual display of direction accuracies. Different datasets and
# subsets are identified explicitly; this is not a head-to-head model ranking.
fig, ax = plt.subplots(figsize=(9.0, 5.8))
labels = [
    "Native C-score\n(broad, n=466 informative)",
    "Native C-score\n(exact 21 aa, n=9)",
    "Training-fold mean\n(n=466 informative)",
    "Ridge\n(n=466 informative)",
    "Structural score\n(overall, n=18)",
]
values = [42.4893, 11.1111, 62.6609, 58.9623, 50.0]
bars = ax.barh(labels, values)
ax.axvline(50, linestyle="--")
ax.set_xlim(0, 100)
ax.set_xlabel("Reported direction accuracy (%)")
ax.set_title("Direction-accuracy results across distinct retrospective tests")
for bar, value in zip(bars, values):
    ax.text(
        min(value + 1.5, 94),
        bar.get_y() + bar.get_height() / 2,
        f"{value:.1f}%",
        va="center",
    )
finish(
    fig,
    "09_direction_accuracy_across_tests.png",
    "Tests use different datasets, subsets, and methods. Do not interpret this as a single controlled head-to-head comparison.",
)

# Three recommended poster figures.
for src, dest in [
    ("01_optimizer_internal_score.png", "01_optimizer_internal_score.png"),
    ("03_broad_native_score_affinity_validation.png", "02_native_score_affinity_validation.png"),
    ("05_mutation_position_validation.png", "03_mutation_position_validation.png"),
]:
    shutil.copy2(G / src, OUT / "poster_graphs" / dest)

# Package instructions and an explicit conclusion.
README = f"""# P3 Final Defensible Evidence Package

Created: {datetime.now().astimezone().isoformat(timespec="seconds")}

## Main conclusion

The existing evidence supports improved optimization against the internal score.
It does not establish improved real-world or experimentally measured binding.

## Start here

1. `logbook_DATA_RESULTS.md` — curated data/results logbook entry.
2. `logbook_ANALYSIS.md` — corresponding interpretation and limitations.
3. `poster_graphs/` — three suggested core poster figures.
4. `graphs/` — nine curated plots for poster selection or supporting material.
5. `supporting_data/` — source metrics, existing candidate records, and audit tables.
6. `source_reports/` — reports and logs from the completed run and validation scripts.
7. `source_archive/original_full_comparison_audit/` — unmodified copy of the earlier full audit.
8. `unresolved_evidence/` — source conflicts that must not be presented as confirmed results.
9. `evidence_manifest.csv` — checksums and package locations.

## Scope

This was a packaging and figure-generation step. It did not launch the optimizer,
add random seeds, rerun a seed-by-seed comparison, or rerun validation tests.

The aggregate optimizer figure uses the existing inventory. The curated figures
exclude the provisional PEPBI 125-match count and the unreconciled positive
exact-21-aa subset correlation. Original outputs are preserved in the archive.

## Important scientific distinction

Internal score improvement is not equivalent to experimentally stronger binding.
The current retrospective benchmarks do not demonstrate that the advantage
translates into better experimental affinity for the generated candidates.
"""
(OUT / "README_FINAL.md").write_text(README, encoding="utf-8")

# Machine-readable list of copied source artifacts and final package checksums.
with (OUT / "copy_log.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["source", "destination", "status"])
    w.writeheader()
    w.writerows(copy_log)

manifest = []
for p in sorted(OUT.rglob("*")):
    if not p.is_file() or p.name == "evidence_manifest.csv":
        continue
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    manifest.append({
        "relative_path": str(p.relative_to(OUT)),
        "size_bytes": p.stat().st_size,
        "sha256": h.hexdigest(),
    })

with (OUT / "evidence_manifest.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(
        f, fieldnames=["relative_path", "size_bytes", "sha256"]
    )
    w.writeheader()
    w.writerows(manifest)

print("\nP3 FINAL DEFENSIBLE EVIDENCE PACKAGE CREATED")
print("Output folder:", OUT)
print("Data logbook:", OUT / "logbook_DATA_RESULTS.md")
print("Analysis logbook:", OUT / "logbook_ANALYSIS.md")
print("Poster figures:", OUT / "poster_graphs")
print("All curated graphs:", OUT / "graphs")
print("Supporting evidence:", OUT / "supporting_data")
print("Unresolved items:", OUT / "unresolved_evidence")
print("Manifest:", OUT / "evidence_manifest.csv")
print("Files hashed:", len(manifest))
print("\nNo optimizer, seed comparison, or validation experiment was launched.")
