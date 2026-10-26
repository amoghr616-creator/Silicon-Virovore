import ctypes
import json
import math
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path.cwd()
WORKBOOK = ROOT / "prototype3/data/raw/PEPBI.xlsx"
OUT = ROOT / "prototype3/results/pepbi_native_score_exploratory"
OUT.mkdir(parents=True, exist_ok=True)

if not WORKBOOK.is_file():
    raise SystemExit(f"Workbook missing: {WORKBOOK}")

# Import the existing native bridge, but call its C function directly.
# Do NOT change or monkey-patch the production length validator.
sys.path.insert(0, str(ROOT))
import c.bridge as bridge

df = pd.read_excel(WORKBOOK, header=1)
df.columns = [str(c).strip() for c in df.columns]

required = [
    "Binding Group",
    "PEPBI Complex Name",
    "Peptide Sequence",
    "ΔG (kcal/mol)",
    "T (K)",
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise SystemExit(f"Missing columns: {missing}")

df["_seq"] = (
    df["Peptide Sequence"].fillna("").astype(str)
    .str.upper().str.replace(r"\s+", "", regex=True)
)
df["_length"] = df["_seq"].str.len()
df["_dg"] = pd.to_numeric(df["ΔG (kcal/mol)"], errors="coerce")
df["_temp"] = pd.to_numeric(df["T (K)"], errors="coerce")
df["_group"] = df["Binding Group"].fillna("").astype(str).str.strip()
df["_valid"] = df["_seq"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+")

def native_score(sequence):
    """Call c_evaluate_sequence directly; do not bypass the C API itself."""
    metrics = bridge.NativeFitnessMetrics()
    status = bridge.lib.c_evaluate_sequence(
        sequence.encode("ascii"),
        ctypes.byref(metrics),
    )
    if status != 0:
        raise RuntimeError(
            f"C scorer returned status={status} for length={len(sequence)}"
        )

    score = float(metrics.fitness_score)
    if not math.isfinite(score) or not 0.0 <= score <= 1.0:
        raise RuntimeError(
            f"Unexpected native score {score!r} for length={len(sequence)}"
        )
    return score

# Controlled smoke tests: one real workbook sequence at every length >= 9.
# Stop at the first failure; do not assume arbitrary-length safety.
sequences = sorted(
    set(df.loc[df["_valid"], "_seq"]),
    key=lambda s: (len(s), s),
)
representatives = {}
for seq in sequences:
    if len(seq) >= 9:
        representatives.setdefault(len(seq), seq)

print("===== DIRECT C API LENGTH SMOKE TEST =====")
for length, seq in sorted(representatives.items()):
    score = native_score(seq)
    print(f"length={length:2d} status=0 score={score:.8f}")

# If the smoke tests all pass, score each distinct canonical sequence
# accepted by the C function. Shorter sequences are excluded because
# the native function rejects inputs shorter than nine residues.
score_cache = {}
for seq in sequences:
    if len(seq) >= 9:
        score_cache[seq] = native_score(seq)

df["_c_score"] = df["_seq"].map(score_cache)

scored = df[
    df["_valid"]
    & df["_dg"].notna()
    & df["_c_score"].notna()
].copy()

# Preserve the individual measurement rows and their contexts.
rows_out = scored[
    [
        "Binding Group",
        "PEPBI Complex Name",
        "Peptide Sequence",
        "_length",
        "_temp",
        "_dg",
        "_c_score",
    ]
].copy()
rows_out.insert(0, "dataframe_row_1based", scored.index + 1)
rows_out = rows_out.rename(columns={
    "_length": "peptide_length",
    "_temp": "temperature_K",
    "_dg": "experimental_delta_g_kcal_mol",
    "_c_score": "native_sequence_prior_score",
})
rows_out.to_csv(OUT / "pepbi_native_score_rows.csv", index=False)

# Descriptive sequence-level summaries prevent duplicated measurements
# from being mistaken for independent peptide sequences.
seq_summary = (
    scored.groupby("_seq")
    .agg(
        native_sequence_prior_score=("_c_score", "first"),
        mean_experimental_delta_g=("_dg", "mean"),
        sd_experimental_delta_g=("_dg", "std"),
        measurement_rows=("_dg", "count"),
        binding_groups=("_group", "nunique"),
        temperatures=("_temp", "nunique"),
    )
    .reset_index()
    .rename(columns={"_seq": "peptide_sequence"})
)
seq_summary["peptide_length"] = seq_summary["peptide_sequence"].str.len()
seq_summary.to_csv(
    OUT / "pepbi_native_score_sequence_summary.csv", index=False
)

def descriptive_corr(data, x, y, method):
    pair = data[[x, y]].dropna()
    if len(pair) < 3 or pair[x].nunique() < 2 or pair[y].nunique() < 2:
        return None
    return float(pair[x].corr(pair[y], method=method))

summary = {
    "workbook_rows": int(len(df)),
    "unique_canonical_sequences": int(df.loc[df["_valid"], "_seq"].nunique()),
    "unique_sequences_scored_by_direct_C_api": int(len(score_cache)),
    "native_C_smoke_test_lengths": sorted(representatives),
    "minimum_native_sequence_length": 9,
    "production_bridge_modified": False,
    "experimental_delta_g_rows_total": int(df["_dg"].notna().sum()),
    "experimental_rows_scored": int(len(scored)),
    "experimental_rows_excluded_due_to_short_sequence": int(
        (df["_valid"] & df["_dg"].notna() & (df["_length"] < 9)).sum()
    ),
    "distinct_measured_sequences_scored": int(
        scored["_seq"].nunique()
    ),
    "row_level_spearman_rho_descriptive": descriptive_corr(
        scored, "_c_score", "_dg", "spearman"
    ),
    "sequence_mean_spearman_rho_descriptive": descriptive_corr(
        seq_summary,
        "native_sequence_prior_score",
        "mean_experimental_delta_g",
        "spearman",
    ),
    "sequence_mean_pearson_r_descriptive": descriptive_corr(
        seq_summary,
        "native_sequence_prior_score",
        "mean_experimental_delta_g",
        "pearson",
    ),
    "interpretation": (
        "Exploratory association only. The native score passes target=NULL; "
        "measurements span multiple protein targets, groups, mutations, and "
        "conditions. Row-level correlation is not independent validation. "
        "Sequence-level mean Delta G can also confound distinct targets."
    ),
}

report_path = OUT / "pepbi_native_score_audit_summary.json"
report_path.write_text(json.dumps(summary, indent=2))

print("\n===== EXPLORATORY RESULTS =====")
print(json.dumps(summary, indent=2))
print("\n===== SEQUENCE-LEVEL SUMMARY =====")
print(seq_summary.to_string(index=False))
print("\n===== OUTPUT FILES =====")
for path in [
    OUT / "pepbi_native_score_rows.csv",
    OUT / "pepbi_native_score_sequence_summary.csv",
    report_path,
]:
    print(f"{path} ({path.stat().st_size:,} bytes)")
