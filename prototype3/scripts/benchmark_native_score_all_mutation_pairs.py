from pathlib import Path
import hashlib
import importlib.util
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path.cwd()
PAIR_FILE = (
    ROOT / "prototype3/results/sequence_only_mutation_baseline/"
    "run_20261026T234823Z/single_substitution_pairs.csv"
)
SCORER_FILE = ROOT / "prototype3/scripts/ppikb_kd_benchmark.py"
OUT = ROOT / "prototype3/results/native_score_all_mutation_pairs"
OUT.mkdir(parents=True, exist_ok=True)

BOOTSTRAPS = 2000
PERMUTATIONS = 5000
SEED = 20261027
AA = set("ACDEFGHIKLMNPQRSTVWY")
rng = np.random.default_rng(SEED)


def rho(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3 or len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return np.nan
    result = float(spearmanr(x, y).statistic)
    return result if np.isfinite(result) else np.nan


def direction_accuracy(frame):
    # Only assess direction when experimental effect is meaningful.
    f = frame.loc[frame["abs_delta_pKd"] >= 0.3]
    if f.empty:
        return np.nan
    score_delta = f["delta_score"].to_numpy()
    exp_delta = f["delta_pKd"].to_numpy()
    result = np.where(
        np.abs(score_delta) <= 1e-12,
        0.5,
        ((score_delta > 0) == (exp_delta > 0)).astype(float),
    )
    return float(result.mean())


def metric_bundle(frame):
    informative = frame.loc[frame["abs_delta_pKd"] >= 0.3]
    return {
        "n_pairs": int(len(frame)),
        "spearman_delta_score_vs_delta_pKd": rho(
            frame["delta_score"], frame["delta_pKd"]
        ),
        "n_informative_pairs": int(len(informative)),
        "direction_accuracy_ties_half": direction_accuracy(frame),
        "spearman_abs_delta_score_vs_abs_delta_pKd": rho(
            informative["delta_score"].abs(),
            informative["delta_pKd"].abs(),
        ),
    }


if not PAIR_FILE.exists():
    raise SystemExit(f"Missing pair file: {PAIR_FILE}")
if not SCORER_FILE.exists():
    raise SystemExit(f"Missing native scorer helper: {SCORER_FILE}")

# Import the exact native_score(sequence) helper used in the project.
spec = importlib.util.spec_from_file_location(
    "project_ppikb_benchmark_for_pair_validation", SCORER_FILE
)
if spec is None or spec.loader is None:
    raise SystemExit("Could not import the project scorer helper.")
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
native_score = module.native_score

pairs = pd.read_csv(PAIR_FILE)
required = {
    "target_group_id", "seq_len", "sequence_a", "sequence_b",
    "pKd_a", "pKd_b",
}
missing = required - set(pairs.columns)
if missing:
    raise SystemExit(f"Pair file is missing columns: {sorted(missing)}")

pairs["sequence_a"] = pairs["sequence_a"].astype(str).str.strip().str.upper()
pairs["sequence_b"] = pairs["sequence_b"].astype(str).str.strip().str.upper()
pairs["pKd_a"] = pd.to_numeric(pairs["pKd_a"], errors="coerce")
pairs["pKd_b"] = pd.to_numeric(pairs["pKd_b"], errors="coerce")
pairs["delta_pKd"] = pairs["pKd_b"] - pairs["pKd_a"]
pairs["abs_delta_pKd"] = pairs["delta_pKd"].abs()

# Verify that these really are canonical, equal-length, single-substitution pairs.
positions = []
for row in pairs.itertuples(index=False):
    a, b = row.sequence_a, row.sequence_b
    if (
        len(a) != len(b)
        or len(a) < 9
        or not set(a) <= AA
        or not set(b) <= AA
    ):
        raise SystemExit(f"Invalid peptide pair: {a!r}, {b!r}")
    diff = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    if len(diff) != 1:
        raise SystemExit(
            f"Expected exactly one substitution; found {len(diff)}: {a}, {b}"
        )
    positions.append(diff[0] + 1)

pairs["mutation_position_1based"] = positions
pairs["seq_len"] = pairs["sequence_a"].str.len()
pairs["mutation_fractional_position"] = (
    (pairs["mutation_position_1based"] - 1) /
    (pairs["seq_len"] - 1)
)

# Score every distinct peptide fresh through the native C API.
unique_sequences = sorted(
    set(pairs["sequence_a"]) | set(pairs["sequence_b"])
)
score_map = {}
for i, sequence in enumerate(unique_sequences, 1):
    score = float(native_score(sequence))
    if not np.isfinite(score) or not 0.0 <= score <= 1.0:
        raise SystemExit(f"Invalid native score {score} for {sequence}")
    score_map[sequence] = score
    if i % 200 == 0:
        print(f"Scored {i}/{len(unique_sequences)} unique sequences")

pairs["score_a"] = pairs["sequence_a"].map(score_map)
pairs["score_b"] = pairs["sequence_b"].map(score_map)
pairs["delta_score"] = pairs["score_b"] - pairs["score_a"]

if pairs[["score_a", "score_b"]].isna().any().any():
    raise SystemExit("At least one paired peptide is missing a native score.")

# Build independent leakage components: connect rows sharing a target sequence
# group or either individual peptide. Resample whole components in the bootstrap.
n = len(pairs)
parent = list(range(n))
size = [1] * n


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(a, b):
    ra, rb = find(a), find(b)
    if ra == rb:
        return
    if size[ra] < size[rb]:
        ra, rb = rb, ra
    parent[rb] = ra
    size[ra] += size[rb]


owners = {}
for i, row in enumerate(pairs.itertuples(index=False)):
    target = str(row.target_group_id)
    tokens = [
        ("target", target),
        ("peptide", row.sequence_a),
        ("peptide", row.sequence_b),
    ]
    for token in tokens:
        if token in owners:
            union(i, owners[token])
        else:
            owners[token] = i

components = {}
for i in range(n):
    components.setdefault(find(i), []).append(i)
component_rows = list(components.values())

# Primary observed metrics.
observed = metric_bundle(pairs)

# Cluster bootstrap CIs, keeping related peptide/target rows together.
boot = {key: [] for key in (
    "spearman_delta_score_vs_delta_pKd",
    "direction_accuracy_ties_half",
    "spearman_abs_delta_score_vs_abs_delta_pKd",
)}
for _ in range(BOOTSTRAPS):
    chosen = rng.integers(0, len(component_rows), size=len(component_rows))
    idx = [i for c in chosen for i in component_rows[c]]
    m = metric_bundle(pairs.iloc[idx])
    for key in boot:
        value = m[key]
        if value is not None and np.isfinite(value):
            boot[key].append(float(value))


def interval(values):
    if len(values) < 100:
        return None
    return [
        float(np.quantile(values, 0.025)),
        float(np.quantile(values, 0.975)),
    ]


bootstrap_ci = {key: interval(values) for key, values in boot.items()}

# Position-level proxy:
# Compare median absolute native-score changes and median absolute experimental
# affinity changes for observed mutation positions, stratified by peptide length.
# This is NOT the adaptive hotspot selector itself.
position_rows = []
for (length, pos), g in pairs.groupby(
    ["seq_len", "mutation_position_1based"]
):
    informative = g[g["abs_delta_pKd"] >= 0.3]
    position_rows.append({
        "seq_len": int(length),
        "mutation_position_1based": int(pos),
        "n_pairs": int(len(g)),
        "n_informative_pairs": int(len(informative)),
        "median_abs_delta_score": float(g["delta_score"].abs().median()),
        "median_abs_delta_pKd": float(g["delta_pKd"].abs().median()),
        "direction_accuracy_ties_half": direction_accuracy(g),
    })

position_summary = pd.DataFrame(position_rows)
position_summary.to_csv(OUT / "position_sensitivity.csv", index=False)

# Position-label shuffle control:
# Within each peptide length, shuffle the score-sensitivity profile across
# eligible positions. Keep only length-position cells with at least 3 pairs.
eligible_pos = position_summary[position_summary["n_pairs"] >= 3].copy()
eligible_pos = eligible_pos.groupby("seq_len").filter(
    lambda g: len(g) >= 2
).reset_index(drop=True)

position_test = {
    "eligible_length_position_cells": int(len(eligible_pos)),
    "minimum_pairs_per_cell": 3,
    "position_profile_spearman": None,
    "position_shuffle_pvalue_two_sided": None,
    "permutations": PERMUTATIONS,
}

if (
    len(eligible_pos) >= 8
    and eligible_pos["median_abs_delta_score"].nunique() > 1
    and eligible_pos["median_abs_delta_pKd"].nunique() > 1
):
    observed_pos_rho = rho(
        eligible_pos["median_abs_delta_score"],
        eligible_pos["median_abs_delta_pKd"],
    )
    position_test["position_profile_spearman"] = observed_pos_rho

    score_profiles = eligible_pos["median_abs_delta_score"].to_numpy()
    exp_profiles = eligible_pos["median_abs_delta_pKd"].to_numpy()
    length_values = eligible_pos["seq_len"].to_numpy()
    grouped_indices = [
        np.flatnonzero(length_values == length)
        for length in np.unique(length_values)
    ]
    null = []
    for _ in range(PERMUTATIONS):
        shuffled = score_profiles.copy()
        for idx in grouped_indices:
            shuffled[idx] = rng.permutation(shuffled[idx])
        val = rho(shuffled, exp_profiles)
        if np.isfinite(val):
            null.append(val)

    if null and np.isfinite(observed_pos_rho):
        pvalue = (
            1 + sum(abs(v) >= abs(observed_pos_rho) for v in null)
        ) / (1 + len(null))
        position_test["position_shuffle_pvalue_two_sided"] = float(pvalue)

metrics = {
    "benchmark": "Native sequence fitness versus experimental mutation effects",
    "pair_file": str(PAIR_FILE),
    "score_function": "native_score(sequence) -> c_evaluate_sequence -> fitness_score",
    "n_matched_single_substitution_pairs": int(len(pairs)),
    "n_unique_peptide_sequences_freshly_scored": int(len(score_map)),
    "peptide_length_min": int(pairs["seq_len"].min()),
    "peptide_length_max": int(pairs["seq_len"].max()),
    "n_exact_21aa_pairs": int((pairs["seq_len"] == 21).sum()),
    "bootstrap_unit": (
        "Connected components formed by shared target groups or individual "
        "peptide sequences"
    ),
    "n_independent_components": int(len(component_rows)),
    "bootstrap_replicates_requested": BOOTSTRAPS,
    "observed_metrics": observed,
    "cluster_bootstrap_95ci": bootstrap_ci,
    "position_profile_test": position_test,
    "interpretation_limits": [
        "The C score is a sequence-property heuristic with no receptor input.",
        "This benchmark tests the native score, not the adaptive optimizer.",
        "The broad dataset includes peptide lengths 9-20; it is not an exact-21-aa-only result.",
        "Position profiling is exploratory and is not validation of ARISE/ALE hotspot selection.",
        "Pairs overlap in peptides and target groups; connected-component bootstrap partly "
        "addresses this dependence.",
    ],
}

pairs.to_csv(OUT / "scored_mutation_pairs.csv", index=False)
(OUT / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

report = [
    "# Native sequence score versus experimental mutation effects",
    "",
    f"- Matched single-substitution pairs: **{len(pairs)}**",
    f"- Distinct peptide sequences freshly scored: **{len(score_map)}**",
    f"- Independent connected components: **{len(component_rows)}**",
    f"- Peptide lengths: **{pairs['seq_len'].min()}–{pairs['seq_len'].max()}**",
    f"- Exact 21-aa pairs in this broader set: **{int((pairs['seq_len'] == 21).sum())}**",
    "",
    "## Main results",
    "",
    "```json",
    json.dumps(observed, indent=2),
    "```",
    "",
    "## Cluster-bootstrap 95% intervals",
    "",
    "```json",
    json.dumps(bootstrap_ci, indent=2),
    "```",
    "",
    "## Position-level shuffle control",
    "",
    "```json",
    json.dumps(position_test, indent=2),
    "```",
    "",
    "## Interpretation",
    "",
    "The native fitness score is a sequence-level heuristic with no receptor input. "
    "This analysis tests whether score differences align with experimental affinity "
    "differences, not whether the adaptive optimizer improves or whether HERV-K Env "
    "binding is predicted. The positional analysis is an exploratory sensitivity proxy, "
    "not direct validation of ARISE/ALE hotspot rankings.",
    "",
]
(OUT / "REPORT.md").write_text("\n".join(report))

print("\n========== NATIVE SCORE MUTATION BENCHMARK ==========")
print(json.dumps(metrics, indent=2))
print("\nFiles written to:", OUT)
