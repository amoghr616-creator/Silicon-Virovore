from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import re
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


ROOT = Path(".")
POSTER = ROOT / "results/ale_isolation_multiseed_20261101/poster_analysis"
FIGURES = POSTER / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

METRICS = POSTER / "comparative_experiment_metrics.csv"
CHECKS = POSTER / "comparative_experiment_run_checks.csv"
SEED_SUMMARY = (
    ROOT / "prototype4/results/teacher_appendix/"
    "ale_10_seed_run_summary.csv"
)
RECON = (
    ROOT / "prototype4/results/p4_final_reconciliation_20261101/"
    "P4_FINAL_RECONCILIATION.json"
)
CLOSEOUT_AUDIT = (
    ROOT / "prototype4/results/p4_closeout_20261101/"
    "p4_closeout_audit.json"
)
SMOKE = (
    ROOT / "results/p4_live_gate_smoke_20261102/"
    "p4_live_gate_smoke_20261102/p4_integrated_assessment/"
    "generation_0001.json"
)


def read_csv(path: Path) -> list[dict]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def nums(value: str) -> list[float]:
    return [
        float(x) for x in re.findall(
            r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",
            str(value),
        )
    ]


def fraction(value: str) -> tuple[int, int]:
    match = re.search(r"(\d+)\s*/\s*(\d+)", str(value))
    if not match:
        raise ValueError(f"Could not parse fraction: {value}")
    return int(match.group(1)), int(match.group(2))


def find_metric(rows: list[dict], text: str) -> dict:
    matches = [r for r in rows if text.lower() in r["outcome"].lower()]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one metric matching {text!r}; found {len(matches)}"
        )
    return matches[0]


def save_figure(fig, stem: str, pdf: PdfPages) -> None:
    fig.tight_layout()
    fig.savefig(FIGURES / f"{stem}.png", dpi=320, bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.svg", bbox_inches="tight")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


# ------------------------------------------------------------
# 1. Validate paired comparison source records.
# ------------------------------------------------------------
metrics = read_csv(METRICS)
checks = read_csv(CHECKS)

if len(checks) != 20:
    raise RuntimeError(f"Expected 20 run checks, found {len(checks)}")

seeds = sorted({r["seed"] for r in checks}, key=int)
arms = sorted({r["arm"] for r in checks})
if len(seeds) != 10 or len(arms) != 2:
    raise RuntimeError(f"Unexpected paired design: seeds={seeds}, arms={arms}")

by_seed = {}
for row in checks:
    by_seed.setdefault(row["seed"], set()).add(row["arm"])
if any(len(x) != 2 for x in by_seed.values()):
    raise RuntimeError("At least one seed is missing a comparison arm")

artifact_columns = [
    "candidate_history.json",
    "experiment_audit.json",
    "pipeline_report.json",
    "top_candidates.csv",
]
artifact_passes = {
    col: sum(row.get(col) == "True" for row in checks)
    for col in artifact_columns
}
p4_counts = sorted({row["p4_generation_reports"] for row in checks})

# The appendix summary includes large serialized JSON cells.
# Increase the standard-library CSV field limit before reading it.
limit = sys.maxsize
while True:
    try:
        csv.field_size_limit(limit)
        break
    except OverflowError:
        limit //= 10

seed_summary_rows = read_csv(SEED_SUMMARY)
summary_seeds = sorted(
    {r.get("seed", "") for r in seed_summary_rows if r.get("seed")},
    key=lambda x: int(x) if x.isdigit() else x,
)
summary_arms = sorted(
    {r.get("arm", "") for r in seed_summary_rows if r.get("arm")}
)

print("========== COMPARATIVE DATA INTEGRITY ==========")
print("Run records:", len(checks))
print("Distinct paired seeds:", len(seeds), seeds)
print("Arms:", arms)
print("Artifact checks:", artifact_passes)
print("P4 reports per run:", p4_counts)
print("Large seed-summary rows:", len(seed_summary_rows))
print("Seed-summary seeds:", len(summary_seeds))
print("Seed-summary arms:", summary_arms)

if any(value != 20 for value in artifact_passes.values()):
    raise RuntimeError("One or more required comparison artifacts are missing")
if p4_counts != ["4"]:
    raise RuntimeError(f"Expected four P4 reports per run, found {p4_counts}")

# ------------------------------------------------------------
# 2. Read aggregate metrics and extract effect estimates.
# ------------------------------------------------------------
mutation = find_metric(metrics, "mutation localization")
best_score = find_metric(metrics, "best overall score")
mean_score = find_metric(metrics, "mean overall score")
identity = find_metric(metrics, "mean pairwise sequence identity")
unique = find_metric(metrics, "unique sequences per generation")
jaccard = find_metric(metrics, "cross-run jaccard")

control_mut, control_total = fraction(mutation["control"])
ale_mut, ale_total = fraction(mutation["ale_guided"])
control_pct = 100 * control_mut / control_total
ale_pct = 100 * ale_mut / ale_total

def effect_and_ci(row: dict) -> tuple[float, float, float]:
    effect_values = nums(row["paired_effect"])
    ci_values = nums(row["ci_95"])
    if len(effect_values) < 1 or len(ci_values) != 2:
        raise ValueError(f"Cannot parse paired effect/CI: {row}")
    return effect_values[0], ci_values[0], ci_values[1]


best_effect, best_lo, best_hi = effect_and_ci(best_score)
mean_effect, mean_lo, mean_hi = effect_and_ci(mean_score)
identity_effect, identity_lo, identity_hi = effect_and_ci(identity)

# Identity values are stored as fractional sequence identity.
# Convert them to percentage points for an intuitive figure axis.
identity_effect *= 100
identity_lo *= 100
identity_hi *= 100

control_unique = nums(unique["control"])[0]
ale_unique = nums(unique["ale_guided"])[0]
jaccard_values = nums(jaccard["ale_guided"])
jaccard_value = jaccard_values[-1]

# ------------------------------------------------------------
# 3. Generate four independent, poster-ready figures.
# ------------------------------------------------------------
pdf_path = FIGURES / "comparative_experiment_figures.pdf"
with PdfPages(pdf_path) as pdf:

    # Figure 1: mutation localization
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    bars = ax.bar(
        ["Control", "ALE-guided"],
        [control_pct, ale_pct],
        width=0.58,
    )
    ax.set_ylabel("Mutations inside selected window (%)")
    ax.set_title("ALE changes where mutations are proposed")
    ax.set_ylim(0, max(70, ale_pct + 12))
    ax.grid(axis="y", alpha=0.25)
    for bar, value in zip(bars, [control_pct, ale_pct]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 1.2,
            f"{value:.2f}%",
            ha="center",
            va="bottom",
        )
    ax.text(
        0.02, 0.97,
        f"Paired difference: {ale_pct-control_pct:+.2f} pp\n"
        f"95% CI: {mutation['ci_95']}\n"
        f"Exact p = {mutation['exact_p']}; "
        f"{mutation['positive_seed_pairs']} seed pairs positive",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
    )
    save_figure(fig, "figure_1_mutation_localization", pdf)

    # Figure 2: paired sequence-fitness effects and confidence intervals
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    entries = [
        ("Mean overall score", mean_effect, mean_lo, mean_hi,
         mean_score["exact_p"]),
        ("Best overall score", best_effect, best_lo, best_hi,
         best_score["exact_p"]),
    ]
    for y, (label, effect, lo, hi, pvalue) in enumerate(entries):
        ax.errorbar(
            effect,
            y,
            xerr=[[effect - lo], [hi - effect]],
            fmt="o",
            capsize=5,
            markersize=7,
        )
        ax.text(
            hi + 0.0007,
            y,
            f"p={pvalue}",
            va="center",
            fontsize=9,
        )
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_yticks([0, 1], [x[0] for x in entries])
    ax.set_xlabel("ALE-guided minus control (sequence-fitness units)")
    ax.set_title("No detected sequence-fitness improvement")
    ax.grid(axis="x", alpha=0.25)
    ax.invert_yaxis()
    save_figure(fig, "figure_2_sequence_fitness_effects", pdf)

    # Figure 3: pairwise sequence-identity effect
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.errorbar(
        identity_effect,
        0,
        xerr=[[
            identity_effect - identity_lo
        ], [
            identity_hi - identity_effect
        ]],
        fmt="o",
        capsize=6,
        markersize=8,
    )
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_yticks([0], ["ALE-guided minus control"])
    ax.set_xlabel("Change in mean pairwise identity (percentage points)")
    ax.set_title("ALE-guided runs show modest sequence convergence")
    ax.grid(axis="x", alpha=0.25)
    ax.text(
        0.02, 0.97,
        f"Effect: +{identity_effect:.3f} pp\n"
        f"95% CI: [{identity_lo:.3f}, {identity_hi:.3f}] pp\n"
        f"Exact p = {identity['exact_p']}; "
        f"{identity['positive_seed_pairs']} pairs positive",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
    )
    save_figure(fig, "figure_3_pairwise_identity_effect", pdf)

    # Figure 4: unique sequences per generation
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    bars = ax.bar(
        ["Control", "ALE-guided"],
        [control_unique, ale_unique],
        width=0.58,
    )
    ax.set_ylabel("Unique sequences per generation")
    ax.set_title("Unique-sequence count was unchanged")
    ax.set_ylim(0, max(control_unique, ale_unique) + 4)
    ax.grid(axis="y", alpha=0.25)
    for bar, value in zip(bars, [control_unique, ale_unique]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.3,
            f"{value:g}",
            ha="center",
        )
    ax.text(
        0.02, 0.97,
        f"ALE final-window mean cross-run Jaccard: {jaccard_value:.4f}\n"
        "Descriptive ALE-only statistic; no control contrast reported.",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
    )
    save_figure(fig, "figure_4_unique_sequence_count", pdf)

# ------------------------------------------------------------
# 4. Write figure captions and results summary.
# ------------------------------------------------------------
captions = f"""# Comparative experiment: figure captions

**Design:** 10 matched seeds (621–630), two arms, 20 total runs.
All 20 runs contain candidate history, experiment audit, pipeline report,
and top-candidate CSV; each has four P4 generation reports.

## Figure 1 — Mutation localization

Control: {control_mut}/{control_total} ({control_pct:.2f}%).
ALE-guided: {ale_mut}/{ale_total} ({ale_pct:.2f}%).
Paired difference: {mutation['paired_effect']}; 95% CI {mutation['ci_95']};
exact paired p = {mutation['exact_p']}; positive seed pairs:
{mutation['positive_seed_pairs']}.

Interpretation: ALE changes the mutation proposal distribution toward the
selected sequence window. This does not establish improved binding or
experimental affinity.

## Figure 2 — Sequence-fitness effects

Best overall score, generations 2–4: ALE minus control =
{best_score['paired_effect']}; 95% CI {best_score['ci_95']};
exact p = {best_score['exact_p']}.

Mean overall score, generations 2–4: ALE minus control =
{mean_score['paired_effect']}; 95% CI {mean_score['ci_95']};
exact p = {mean_score['exact_p']}.

Interpretation: both intervals include zero. The experiment did not detect
a statistically significant sequence-fitness improvement.

## Figure 3 — Pairwise sequence identity

ALE minus control: {identity['paired_effect']}; 95% CI
{identity['ci_95']}; exact p = {identity['exact_p']};
positive seed pairs: {identity['positive_seed_pairs']}.

Interpretation: ALE-guided populations became modestly more similar.
This is consistent with a possible convergence/diversity trade-off.

## Figure 4 — Unique sequences per generation

Control: {unique['control']}; ALE-guided: {unique['ale_guided']}.
The number of unique sequences was unchanged.

ALE final-window mean cross-run Jaccard similarity was {jaccard_value:.4f};
this is descriptive and does not have a separately reported control contrast.

## Overall conclusion

The tested adaptive policy altered mutation localization and modestly
increased pairwise sequence identity, but did not demonstrate improved
sequence-level fitness. No structural binding or experimental affinity
conclusion follows from this comparison.
"""
(POSTER / "COMPARATIVE_EXPERIMENT_CAPTIONS.md").write_text(
    captions, encoding="utf-8"
)

# ------------------------------------------------------------
# 5. Verify the live smoke-test report and P4 closeout evidence.
# ------------------------------------------------------------
smoke = read_json(SMOKE)
recon = read_json(RECON)
closeout = read_json(CLOSEOUT_AUDIT)

gate = smoke.get("integrated_gate", {})
reliability = smoke.get("reliability", {})
epistasis = smoke.get("epistasis_live_loop", {})
risk = recon.get("risk_coverage", {}).get("validation_metrics", {})
md_inputs = recon.get("md_inputs", {})
closeout_benchmark = closeout.get("benchmark", {})

model_rmse = (
    risk.get("sequence_mutation_descriptor_model", {})
    .get("rmse_kcal_mol")
)
baseline_rmse = (
    risk.get("baseline_training_fold_mean", {})
    .get("rmse_kcal_mol")
)
uncertainty = risk.get("uncertainty_error_association", {})
trajectory_entropy = md_inputs.get("trajectory_entropy_calculated")
entropy_rmse = recon.get("entropy_rmse_results", [])

family_count_closeout = closeout_benchmark.get("unique_family_ids")

closeout_text = f"""# P4 evidence closeout — 2026-11-02

## Verified implementation

- Candidate-level gate is called in the live generation loop.
- One-generation smoke test exited successfully.
- Candidates adjudicated: {smoke.get('candidate_count')}.
- Target evidence coverage: {smoke.get('target_evidence_coverage')}.
- Integrated-gate eligible: {gate.get('eligible_count')}.
- Integrated-gate abstentions: {gate.get('abstain_count')}.
- Candidate target evidence valid: {gate.get('target_valid_count')}.
- Candidate reliability evidence valid: {gate.get('reliability_valid_count')}.
- Candidate epistasis evidence valid: {gate.get('epistasis_valid_count')}.
- Pairwise epistasis operator active: {epistasis.get('active_candidate_count')}.
- Candidate-level validated delta-epsilon predictions:
  {epistasis.get('candidate_level_delta_epsilon_validated_count')}.
- Reliability calibrated in smoke test: {reliability.get('calibrated')}.
- Scalar fitness generated by the gate: {gate.get('scalar_fitness_generated')}.

The smoke test used a historical 9SYA manifest with a limited sequence
coverage. Pareto ordering was correctly withheld because only a subset
of candidates had valid target evidence. Abstention did not stop the
legacy sequence-level search.

## Reliability benchmark result

- Family-held-out sequence/mutation model RMSE:
  {model_rmse} kcal/mol.
- Training-fold-mean baseline RMSE:
  {baseline_rmse} kcal/mol.
- Uncertainty/error Spearman rho:
  {uncertainty.get('calibrated_error_abs_error_spearman_rho')}.
- Uncertainty/error association p-value:
  {uncertainty.get('calibrated_error_abs_error_p')}.

The sequence/mutation model RMSE is worse than the baseline. The current
evidence does not justify using candidate-level reliability as an
optimizer objective or feedback weight.

## Remaining P4 scientific requirements

1. **Trajectory contact entropy:** reconciliation records
   `trajectory_entropy_calculated = {trajectory_entropy}`.
   There are {len(md_inputs.get('trajectory_files', []))} trajectory input
   records and {len(md_inputs.get('topology_candidates', []))}
   topology candidates, but the required trajectory entropy result has
   not been produced in the reconciliation artifact.

2. **Entropy versus predictive error:** recorded result rows:
   {len(entropy_rmse)}. The incremental entropy/RMSE comparison remains
   unreported.

3. **Candidate-level reliability calibration:** not demonstrated by the
   live smoke test. The repeatability diagnostic is not a calibration
   of structural prediction error.

4. **Candidate-level epistasis:** the pairwise epistasis/entropy operator
   is active in proposal generation, but the smoke test has
   {epistasis.get('candidate_level_delta_epsilon_validated_count')}
   validated candidate-specific delta-epsilon predictions.

5. **EVSI and pruning:** computational scores alone do not satisfy the
   independent experimental-evidence requirement for candidate pruning.
   The current gate correctly abstains.

6. **Family identity audit:** the closeout artifact reports
   `unique_family_ids = {family_count_closeout}`. Reconcile this field
   against the family identifiers used for the 242-cycle LOFO benchmark
   before making a family-count claim in the final presentation.

7. **Biological validation:** no experimental HERV-K Env binding
   validation is established by these computational analyses.

## Comparative experiment

The 10-seed, 20-run comparison is complete. Its verified outputs and
poster figures are stored in:
`results/ale_isolation_multiseed_20261101/poster_analysis/`.

## Final interpretation

P4 has a functioning live candidate-evidence adjudication path and
reproducible generation-level reporting. It is not a validated
HERV-K binder-discovery pipeline. The outstanding requirements concern
trajectory-derived uncertainty, reliability validation, candidate-level
epistasis evidence, and experimental decision support.
"""
(ROOT / "P4_EVIDENCE_CLOSEOUT_20261102.md").write_text(
    closeout_text, encoding="utf-8"
)

# ------------------------------------------------------------
# 6. Save a compact provenance manifest for generated outputs.
# ------------------------------------------------------------
inputs = [METRICS, CHECKS, SEED_SUMMARY, RECON, CLOSEOUT_AUDIT, SMOKE]
outputs = [
    FIGURES / "figure_1_mutation_localization.png",
    FIGURES / "figure_2_sequence_fitness_effects.png",
    FIGURES / "figure_3_pairwise_identity_effect.png",
    FIGURES / "figure_4_unique_sequence_count.png",
    pdf_path,
    POSTER / "COMPARATIVE_EXPERIMENT_CAPTIONS.md",
    ROOT / "P4_EVIDENCE_CLOSEOUT_20261102.md",
]

manifest = {
    "generated_utc": __import__("datetime").datetime.now(
        __import__("datetime").timezone.utc
    ).isoformat(),
    "paired_seeds": seeds,
    "arms": arms,
    "run_records": len(checks),
    "artifact_pass_counts": artifact_passes,
    "p4_reports_per_run": p4_counts,
    "input_sha256": {
        str(p): sha256(p) for p in inputs
    },
    "outputs": [str(p) for p in outputs],
}
(POSTER / "FINAL_COMPARISON_OUTPUT_MANIFEST.json").write_text(
    json.dumps(manifest, indent=2) + "\n",
    encoding="utf-8",
)

print("\n========== FIGURES CREATED ==========")
for path in outputs:
    print(path, f"({path.stat().st_size:,} bytes)" if path.exists() else "MISSING")

print("\n========== P4 OPEN REQUIREMENTS ==========")
print("trajectory_entropy_calculated:", trajectory_entropy)
print("entropy_vs_RMSE_result_rows:", len(entropy_rmse))
print("reliability model RMSE:", model_rmse)
print("baseline RMSE:", baseline_rmse)
print("uncertainty/error rho:", uncertainty.get(
    "calibrated_error_abs_error_spearman_rho"
))
print("candidate-level validated epistasis:",
      epistasis.get("candidate_level_delta_epsilon_validated_count"))
print("closeout unique_family_ids:", family_count_closeout)
print("gmx executable:", shutil.which("gmx"))
print("MDAnalysis installed:", importlib.util.find_spec("MDAnalysis") is not None)
print("mdtraj installed:", importlib.util.find_spec("mdtraj") is not None)

print("\nBUILD COMPLETE")
