from pathlib import Path
import json
from datetime import datetime, timezone

ROOT = Path.cwd()
OUT = ROOT / "prototype3/results/P3_FINAL_20261026"
OUT.mkdir(parents=True, exist_ok=True)

def load_json(path):
    try:
        return json.loads(path.read_text())
    except Exception as e:
        return {"_read_error": str(e), "_path": str(path)}

def newest(pattern):
    matches = list(ROOT.glob(pattern))
    matches = [p for p in matches if p.is_file()]
    return max(matches, key=lambda p: p.stat().st_mtime) if matches else None

def fmt(x):
    if isinstance(x, float):
        return f"{x:.4f}"
    return str(x)

position_path = (
    ROOT / "prototype3/results/arise_experimental_position_audit/metrics.json"
)
sequence_path = newest(
    "prototype3/results/sequence_only_mutation_baseline/run_*/metrics.json"
)
structural_path = newest(
    "prototype3/results/experimental_mutation_benchmark/"
    "structural_scoring_*/benchmark_metrics.json"
)

position = load_json(position_path) if position_path.exists() else {}
sequence = load_json(sequence_path) if sequence_path else {}
structural = load_json(structural_path) if structural_path else {}

files_used = {
    "position_audit": str(position_path) if position_path.exists() else None,
    "sequence_affinity_metrics": str(sequence_path) if sequence_path else None,
    "structural_metrics": str(structural_path) if structural_path else None,
}

lines = [
    "# Silicon Virovore — P3 Evidence Report",
    "",
    f"Generated: {datetime.now(timezone.utc).isoformat()}",
    "",
    "## 1. Research question",
    "",
    "Does the computational pipeline provide evidence of useful mutation-position "
    "signals and agreement between computational scores and experimentally "
    "measured peptide affinity?",
    "",
    "## 2. Experimental mutation-position test",
    "",
]

primary = position.get("primary_window_test", {})
filtered = position.get("informative_subset_sensitivity", {})
previous = position.get("previous_generation_window_sensitivity", {})

if primary:
    lines += [
        "The latest logged hotspot maps to peptide positions 1–3. The preceding "
        "generation's hotspot maps to positions 2–3.",
        "",
        "| Analysis | Centered effect difference | Cluster-bootstrap 95% CI | "
        "One-sided window permutation p | Two-sided p |",
        "|---|---:|---:|---:|---:|",
    ]
    for label, result in [
        ("All eligible mutation pairs; positions 1–3", primary),
        ("Larger-effect subset; positions 1–3", filtered),
        ("All eligible pairs; positions 2–3 sensitivity check", previous),
    ]:
        if result.get("centered_effect_difference_window_minus_elsewhere") is None:
            continue
        ci = result.get("target_cluster_bootstrap_95ci")
        ci_text = f"[{fmt(ci[0])}, {fmt(ci[1])}]" if ci else "Unavailable"
        lines.append(
            f"| {label} | "
            f"{fmt(result['centered_effect_difference_window_minus_elsewhere'])} | "
            f"{ci_text} | "
            f"{fmt(result.get('random_contiguous_window_permutation_p_greater'))} | "
            f"{fmt(result.get('random_contiguous_window_permutation_p_two_sided'))} |"
        )
    lines += [
        "",
        f"- Total single-substitution pairs: {position.get('n_single_substitution_pairs', 'unknown')}.",
        f"- Pairs with |ΔpKd| ≥ 0.3: {position.get('n_informative_pairs_abs_delta_pKd_ge_0_3', 'unknown')}.",
        "",
        "**Interpretation:** The full-pair analysis does not show that the selected "
        "positions outperform the comparison positions. The larger-effect subset "
        "has a positive exploratory association, but the two-sided window "
        "permutation result is borderline. Both results must be reported. The "
        "subset result is not independent confirmation, and this cross-target "
        "analysis does not establish HERV-K-specific biological importance.",
        "",
    ]
else:
    lines += [
        "Position-window metrics were not found at the expected path. Do not "
        "claim that this test is complete until its metrics are located.",
        "",
    ]

lines += [
    "## 3. Sequence-based affinity benchmark",
    "",
]

if sequence:
    lines += [
        f"Metrics source: `{sequence_path.relative_to(ROOT)}`",
        "",
        "The sequence-based model is a separate benchmark model. Its performance "
        "must not be attributed to the native C `c_score`, ARISE/ALE, or the "
        "target-aware structural score.",
        "",
        "```json",
        json.dumps(sequence, indent=2),
        "```",
        "",
        "**Interpretation requirement:** Verify that validation folds separate "
        "related targets and that no duplicate or closely related peptide "
        "measurements leak across train and test partitions before describing "
        "the result as generalization to unseen targets.",
        "",
    ]
else:
    lines += [
        "Sequence-affinity metrics were not found. No claim is made about this "
        "benchmark in this report.",
        "",
    ]

lines += [
    "## 4. Experimental structural-scoring benchmark",
    "",
]

if structural:
    lines += [
        f"Metrics source: `{structural_path.relative_to(ROOT)}`",
        "",
        "```json",
        json.dumps(structural, indent=2),
        "```",
        "",
        "**Interpretation:** Treat this as a small structural-score benchmark. "
        "A geometric score is not an experimentally measured affinity, and "
        "correlation with ΔpKd must be demonstrated on held-out complexes "
        "before claiming predictive accuracy.",
        "",
    ]
else:
    lines += [
        "Structural benchmark metrics were not found. Do not imply that the "
        "structural score has been validated against experimental affinity.",
        "",
    ]

lines += [
    "## 5. What the evidence supports",
    "",
    "- The project has implemented and tested computational search and scoring components.",
    "- A retrospective experimental mutation dataset allows an exploratory test of position-level signals.",
    "- Sequence-based and structural scores can be evaluated against measured affinity, but their validity must be established separately.",
    "",
    "## 6. What the evidence does not yet support",
    "",
    "- That the native `c_score` predicts experimental binding affinity.",
    "- That the current geometric structural score reliably predicts affinity or binding specificity.",
    "- That ARISE/ALE improves affinity optimization over a properly matched baseline.",
    "- That a hotspot inferred from sequence-score patterns is biologically important for HERV-K Env residues ARG203/ARG205.",
    "- That the pipeline is accurate or reliable enough for real-world biological use.",
    "",
    "## 7. Required reporting safeguards",
    "",
    "1. Preserve negative and inconclusive results alongside positive exploratory findings.",
    "2. Do not pad or truncate peptides to satisfy the native 21-aa `c_score` interface.",
    "3. Do not select score weights or thresholds on the same data later described as validation.",
    "4. Distinguish internal optimization performance from external predictive validity.",
    "5. State that retrospective cross-target analyses are not prospective HERV-K validation.",
    "",
    "## 8. Source files",
    "",
]
for label, path in files_used.items():
    lines.append(f"- **{label}:** `{path or 'NOT FOUND'}`")

report = "\n".join(lines) + "\n"
(OUT / "P3_FINAL_EVIDENCE_REPORT.md").write_text(report)
(OUT / "source_manifest.json").write_text(json.dumps(files_used, indent=2) + "\n")

print(report)
print(f"\nSaved final report to: {OUT / 'P3_FINAL_EVIDENCE_REPORT.md'}")
