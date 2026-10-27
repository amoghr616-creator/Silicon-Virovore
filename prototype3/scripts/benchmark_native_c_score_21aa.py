from pathlib import Path
import csv
import hashlib
import importlib.util
import inspect
import json
import math
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path.cwd()
INPUT = ROOT / "prototype3/results/p3_21aa_affinity_audit/ppi_kb_exact21_candidates.csv"
OUT = ROOT / "prototype3/results/p3_21aa_affinity_audit/native_c_score_21aa"
OUT.mkdir(parents=True, exist_ok=True)

AA = set("ACDEFGHIKLMNPQRSTVWY")
AFFINITY_RE = re.compile(
    r"^\s*Kd\s*=\s*([0-9]*\.?[0-9]+(?:[eE][+-]?\d+)?)"
    r"\s*([fpnumµμ]?M)\s*$",
    re.IGNORECASE,
)
UNIT_TO_MOLAR = {
    "m": 1.0,
    "mm": 1e-3,
    "um": 1e-6,
    "nm": 1e-9,
    "pm": 1e-12,
    "fm": 1e-15,
}
EMPTY = {"", "nan", "none", "null", "na", "n/a"}
UNMODIFIED = {"", "nan", "none", "null", "na", "n/a",
              "no", "false", "0", "unmodified"}
FREE_TERMINUS = UNMODIFIED | {"free"}

def text(row, column):
    value = row.get(column, "")
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value.lower() in EMPTY else value

def canonical_sequence(value):
    seq = re.sub(r"\s+", "", str(value or "")).upper()
    return seq if len(seq) == 21 and set(seq) <= AA else None

def parse_exact_kd(affinity):
    """Convert an uncensored, explicitly reported Kd to pKd."""
    match = AFFINITY_RE.fullmatch(str(affinity or "").strip())
    if not match:
        return None

    value = float(match.group(1))
    unit = (
        match.group(2).strip()
        .replace("μ", "u")
        .replace("µ", "u")
        .lower()
    )
    if value <= 0 or unit not in UNIT_TO_MOLAR:
        return None

    molar = value * UNIT_TO_MOLAR[unit]
    return -math.log10(molar)

def no_explicit_modification(row):
    # The sequence-only scorer cannot represent these chemical differences.
    for column in ("Chemical_Modification", "Residue_Modification"):
        value = text(row, column).lower()
        if value not in UNMODIFIED:
            return False

    for column in ("N-terminal_Modification", "C-terminal_Modification"):
        value = text(row, column).lower()
        if value not in FREE_TERMINUS:
            return False

    status = text(row, "Linear/Cyclic").lower()
    if status != "linear":
        return False

    for column in ("Cyclization_Method", "Cyclization_Way"):
        value = text(row, column).lower()
        if value not in EMPTY and value not in {"no", "none", "linear"}:
            return False

    return True

def stable_id(value):
    return hashlib.sha256(str(value).encode()).hexdigest()[:16]

if not INPUT.exists():
    raise SystemExit(f"Candidate CSV not found: {INPUT}")

# Import the established helper without running the benchmark script's main.
script_path = ROOT / "prototype3/scripts/ppikb_kd_benchmark.py"
spec = importlib.util.spec_from_file_location("project_ppikb_benchmark", script_path)
if spec is None or spec.loader is None:
    raise SystemExit(f"Cannot load scorer helper: {script_path}")

module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)

score_function = getattr(module, "native_score", None)
if not callable(score_function):
    raise SystemExit(
        "No callable native_score(sequence) found in "
        "prototype3/scripts/ppikb_kd_benchmark.py"
    )

print("Scoring helper:", script_path)
print("Function signature:", inspect.signature(score_function))
print("Function implementation:")
print(inspect.getsource(score_function))

raw = pd.read_csv(INPUT, dtype=str, keep_default_na=False)
required = {"Peptide_Sequence", "Protein_Sequence", "Affinity"}
missing = required - set(raw.columns)
if missing:
    raise SystemExit(f"Missing required CSV columns: {sorted(missing)}")

eligible = []
rejection = {
    "not_exact_canonical_21aa": 0,
    "not_numeric_uncensored_Kd": 0,
    "explicit_modification_or_not_linear": 0,
    "missing_target_sequence": 0,
}

for row in raw.to_dict("records"):
    seq = canonical_sequence(row.get("Peptide_Sequence", ""))
    if seq is None:
        rejection["not_exact_canonical_21aa"] += 1
        continue

    pkd = parse_exact_kd(row.get("Affinity", ""))
    if pkd is None:
        rejection["not_numeric_uncensored_Kd"] += 1
        continue

    if not no_explicit_modification(row):
        rejection["explicit_modification_or_not_linear"] += 1
        continue

    target = re.sub(
        r"\s+", "", text(row, "Protein_Sequence")
    ).upper()
    if not target:
        rejection["missing_target_sequence"] += 1
        continue

    assay = text(row, "Affinity_Assay") or "assay_not_reported"
    doi = text(row, "DOI")
    pmid = text(row, "PMID")
    title = text(row, "Title")
    source = text(row, "audit_source_file") or "source_not_recorded"
    paper = doi or (f"PMID:{pmid}" if pmid else "") or title or source

    context_payload = "\x1f".join([target, assay, paper])
    context_id = stable_id(context_payload)

    eligible.append({
        "source_id": text(row, "ID"),
        "source_file": source,
        "peptide_sequence": seq,
        "target_sequence": target,
        "target_group_id": stable_id(target),
        "assay": assay,
        "study_id": paper,
        "context_id": context_id,
        "affinity_raw": text(row, "Affinity"),
        "pKd": pkd,
        "DOI": doi,
        "PMID": pmid,
        "Title": title,
    })

if not eligible:
    raise SystemExit(
        "No records passed the strict exact-21-aa, numeric-Kd, "
        "linear/no-explicit-modification filters. Review the source annotations."
    )

eligible_df = pd.DataFrame(eligible)

# Cache one deterministic score per distinct peptide sequence.
score_cache = {}
score_errors = {}
for seq in sorted(eligible_df["peptide_sequence"].unique()):
    try:
        value = score_function(seq)
        if value is None:
            raise ValueError("native_score returned None")
        value = float(value)
        if not math.isfinite(value):
            raise ValueError(f"non-finite score: {value}")
        score_cache[seq] = value
    except Exception as exc:
        score_errors[seq] = f"{type(exc).__name__}: {exc}"

eligible_df["native_score"] = eligible_df["peptide_sequence"].map(score_cache)
eligible_df.to_csv(OUT / "eligible_scored_measurements.csv", index=False)

scored = eligible_df.dropna(subset=["native_score"]).copy()
if scored.empty:
    raise SystemExit(
        "Scoring failed for every eligible sequence. "
        "Check the printed scorer implementation and native-library setup."
    )

context_cols = ["context_id", "target_group_id", "assay", "study_id"]

# Collapse repeated measurements for the same peptide within one context.
observations = (
    scored.groupby(context_cols + ["peptide_sequence"], as_index=False)
    .agg(
        pKd=("pKd", "median"),
        native_score=("native_score", "first"),
        source_records=("source_id", "size"),
    )
)
observations.to_csv(OUT / "context_peptide_observations.csv", index=False)

# Primary analysis: rank agreement within independent target/assay/study contexts.
group_results = []
for keys, group in observations.groupby(context_cols, dropna=False):
    if len(group) < 3:
        continue
    if group["pKd"].nunique() < 2 or group["native_score"].nunique() < 2:
        continue

    rho, pvalue = spearmanr(group["native_score"], group["pKd"])
    if not math.isfinite(float(rho)):
        continue

    context_id, target_id, assay, study = keys
    group_results.append({
        "context_id": context_id,
        "target_group_id": target_id,
        "assay": assay,
        "study_id": study,
        "n_distinct_peptides": len(group),
        "spearman_rho": float(rho),
        "exploratory_pvalue": float(pvalue),
    })

groups_df = pd.DataFrame(group_results)
groups_df.to_csv(OUT / "within_context_correlations.csv", index=False)

rng = np.random.default_rng(20261027)
if len(groups_df):
    rhos = groups_df["spearman_rho"].to_numpy()
    boot = np.array([
        rng.choice(rhos, size=len(rhos), replace=True).mean()
        for _ in range(2000)
    ])
    mean_rho = float(rhos.mean())
    median_rho = float(np.median(rhos))
    rho_ci = [float(np.quantile(boot, 0.025)),
              float(np.quantile(boot, 0.975))]
else:
    mean_rho = median_rho = None
    rho_ci = None

# Secondary analysis: single substitutions within the same target/assay/study.
pair_rows = []
for keys, group in observations.groupby(context_cols, dropna=False):
    records = group.to_dict("records")
    context_id, target_id, assay, study = keys

    for i in range(len(records)):
        a = records[i]
        seq_a = a["peptide_sequence"]
        for j in range(i + 1, len(records)):
            b = records[j]
            seq_b = b["peptide_sequence"]

            differences = [
                pos for pos in range(21) if seq_a[pos] != seq_b[pos]
            ]
            if len(differences) != 1:
                continue

            pos = differences[0]
            delta_pkd = float(b["pKd"] - a["pKd"])
            delta_score = float(b["native_score"] - a["native_score"])

            if abs(delta_pkd) < 0.3:
                direction = None
            elif abs(delta_score) < 1e-12:
                direction = 0.5
            else:
                direction = float(
                    (delta_pkd > 0) == (delta_score > 0)
                )

            pair_rows.append({
                "context_id": context_id,
                "target_group_id": target_id,
                "assay": assay,
                "study_id": study,
                "sequence_a": seq_a,
                "sequence_b": seq_b,
                "pKd_a": a["pKd"],
                "pKd_b": b["pKd"],
                "delta_pKd_b_minus_a": delta_pkd,
                "score_a": a["native_score"],
                "score_b": b["native_score"],
                "delta_score_b_minus_a": delta_score,
                "mutation_position_1based": pos + 1,
                "mutation_A_to_B": f"{seq_a[pos]}>{seq_b[pos]}",
                "informative_abs_delta_pKd_ge_0_3": abs(delta_pkd) >= 0.3,
                "direction_accuracy_ties_half": direction,
            })

pairs = pd.DataFrame(pair_rows)
pairs.to_csv(OUT / "single_substitution_pairs.csv", index=False)

informative = (
    pairs[pairs["informative_abs_delta_pKd_ge_0_3"]].copy()
    if len(pairs) else pd.DataFrame()
)

if len(informative):
    direction_values = informative["direction_accuracy_ties_half"].dropna()
    direction_accuracy = (
        float(direction_values.mean()) if len(direction_values) else None
    )
    if len(informative) >= 3 and informative["delta_score_b_minus_a"].nunique() > 1:
        mutation_rho = float(spearmanr(
            informative["delta_score_b_minus_a"],
            informative["delta_pKd_b_minus_a"],
        ).statistic)
    else:
        mutation_rho = None
else:
    direction_accuracy = mutation_rho = None

if len(pairs):
    position_summary = (
        pairs.groupby("mutation_position_1based")
        .agg(
            n_single_substitution_pairs=("delta_pKd_b_minus_a", "size"),
            median_abs_experimental_delta_pKd=(
                "delta_pKd_b_minus_a", lambda x: float(np.median(np.abs(x)))
            ),
            median_abs_score_change=(
                "delta_score_b_minus_a", lambda x: float(np.median(np.abs(x)))
            ),
        )
        .reset_index()
    )
    position_summary.to_csv(OUT / "mutation_position_summary.csv", index=False)
else:
    position_summary = pd.DataFrame(columns=[
        "mutation_position_1based", "n_single_substitution_pairs",
        "median_abs_experimental_delta_pKd", "median_abs_score_change"
    ])
    position_summary.to_csv(OUT / "mutation_position_summary.csv", index=False)

metrics = {
    "benchmark": "PPIKB exact canonical 21-aa native sequence-score audit",
    "input_candidate_rows": int(len(raw)),
    "rows_passing_exact21_kd_linear_no_explicit_modification_filters": int(len(eligible_df)),
    "unique_peptide_sequences_eligible": int(eligible_df["peptide_sequence"].nunique()),
    "unique_peptides_scored_successfully": int(len(score_cache)),
    "score_failures": score_errors,
    "filter_rejection_counts": rejection,
    "measurement_contexts": int(eligible_df["context_id"].nunique()),
    "contexts_with_at_least_3_distinct_scored_peptides_and_variable_values": int(len(groups_df)),
    "within_context_spearman_mean": mean_rho,
    "within_context_spearman_median": median_rho,
    "within_context_spearman_group_bootstrap_95ci_for_mean": rho_ci,
    "single_substitution_pairs_same_context": int(len(pairs)),
    "informative_single_substitution_pairs_abs_delta_pKd_ge_0_3": int(len(informative)),
    "informative_mutation_direction_accuracy_ties_half": direction_accuracy,
    "informative_delta_score_vs_delta_pKd_spearman_exploratory": mutation_rho,
    "interpretation_warning": (
        "These are retrospective associations, not prospective prediction. "
        "The native score is sequence-only and has no receptor input. "
        "Pairs overlap in peptides and are not independent observations. "
        "The position summary describes this filtered dataset, not a validated "
        "biological hotspot selector. The grouped correlations are exploratory."
    ),
}

(OUT / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

report = [
    "# Native sequence-score audit on PPIKB 21-aa peptides",
    "",
    f"- Candidate rows examined: **{len(raw)}**",
    f"- Rows passing strict filters: **{len(eligible_df)}**",
    f"- Distinct scored sequences: **{len(score_cache)}**",
    f"- Independent target/assay/publication contexts: **{eligible_df['context_id'].nunique()}**",
    f"- Contexts with valid within-context Spearman correlations: **{len(groups_df)}**",
    f"- Mean within-context Spearman rho: **{mean_rho}**",
    f"- Median within-context Spearman rho: **{median_rho}**",
    f"- Group-bootstrap 95% interval for mean rho: **{rho_ci}**",
    f"- Same-context single-substitution pairs: **{len(pairs)}**",
    f"- Informative single-substitution pairs (|delta pKd| >= 0.3): **{len(informative)}**",
    f"- Mutation direction accuracy, ties counted as half: **{direction_accuracy}**",
    f"- Exploratory Spearman of score change vs experimental pKd change: **{mutation_rho}**",
    "",
    "## Filter counts",
]
report.extend([f"- {key}: {value}" for key, value in rejection.items()])
report.extend([
    "",
    "## Interpretation limits",
    "",
    "This evaluates a sequence-only score against measured Kd records; it does not "
    "evaluate the adaptive optimizer or the target-aware structural score. Scores "
    "are compared within target/assay/publication contexts to reduce cross-target "
    "confounding. Small context counts or wide uncertainty mean the result remains "
    "preliminary. A positive association does not establish HERV-K-specific binding.",
    "",
    "## Output files",
    "",
    "- `eligible_scored_measurements.csv`",
    "- `context_peptide_observations.csv`",
    "- `within_context_correlations.csv`",
    "- `single_substitution_pairs.csv`",
    "- `mutation_position_summary.csv`",
    "- `metrics.json`",
])
(OUT / "REPORT.md").write_text("\n".join(report) + "\n")

print("\n========== NATIVE 21-AA BENCHMARK SUMMARY ==========")
print(json.dumps(metrics, indent=2))
print("\nSaved outputs in:", OUT)
