import ctypes
import json
import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path.cwd()
WORKBOOK = ROOT / "prototype3/data/raw/PEPBI.xlsx"
OUT = ROOT / "prototype3/results/pepbi_corrected_benchmark"
OUT.mkdir(parents=True, exist_ok=True)

if not WORKBOOK.is_file():
    raise SystemExit(f"Workbook not found: {WORKBOOK}")

sys.path.insert(0, str(ROOT))
import c.bridge as bridge

df = pd.read_excel(WORKBOOK, header=1)
df.columns = [str(c).strip() for c in df.columns]

required = [
    "Binding Group",
    "PEPBI Complex Name",
    "Peptide Sequence",
    "Peptide Length",
    "ΔG (kcal/mol)",
    "Calculated ΔG (kcal/mol)",
    "dG_separated",
]
missing = [c for c in required if c not in df.columns]
if missing:
    print("Available potentially relevant columns:")
    print([c for c in df.columns
           if "g" in c.lower() or "delta" in c.lower()
           or "dg" in c.lower()])
    raise SystemExit(f"Missing expected columns: {missing}")

df["_seq"] = (
    df["Peptide Sequence"].fillna("").astype(str)
    .str.upper().str.replace(r"\s+", "", regex=True)
)
df["_length"] = df["_seq"].str.len()
df["_reported_length"] = pd.to_numeric(
    df["Peptide Length"], errors="coerce"
)
df["_direct_dg"] = pd.to_numeric(df["ΔG (kcal/mol)"], errors="coerce")
df["_calc_dg"] = pd.to_numeric(
    df["Calculated ΔG (kcal/mol)"], errors="coerce"
)
df["_ria_dg"] = pd.to_numeric(df["dG_separated"], errors="coerce")
df["_group"] = df["Binding Group"].fillna("").astype(str).str.strip()

canonical = df["_seq"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+")
df["_canonical"] = canonical.fillna(False)

def native_score(sequence):
    metrics = bridge.NativeFitnessMetrics()
    status = bridge.lib.c_evaluate_sequence(
        sequence.encode("ascii"),
        ctypes.byref(metrics),
    )
    if status != 0:
        raise RuntimeError(
            f"C API status={status}; peptide length={len(sequence)}"
        )
    score = float(metrics.fitness_score)
    if not math.isfinite(score) or not 0 <= score <= 1:
        raise RuntimeError(f"Invalid score {score!r} for {sequence}")
    return score

# The C API has a minimum length of 9. Do not pad or truncate.
eligible_sequences = sorted({
    seq for seq in df.loc[df["_canonical"], "_seq"].unique()
    if len(seq) >= 9
})
scores = {}
failures = []

for seq in eligible_sequences:
    try:
        scores[seq] = native_score(seq)
    except Exception as exc:
        failures.append({"sequence": seq, "error": str(exc)})

df["_c_score"] = df["_seq"].map(scores)
df["_primary_eligible"] = (
    df["_canonical"]
    & df["_length"].between(9, 20)
    & df["_c_score"].notna()
)

# Repeat evaluations of selected sequences to test observed determinism.
stability = []
for seq in eligible_sequences[:min(10, len(eligible_sequences))]:
    vals = [native_score(seq) for _ in range(3)]
    stability.append({
        "sequence": seq,
        "length": len(seq),
        "score_1": vals[0],
        "score_2": vals[1],
        "score_3": vals[2],
        "max_minus_min": max(vals) - min(vals),
    })

stability_df = pd.DataFrame(stability)
stable = (
    stability_df.empty
    or bool((stability_df["max_minus_min"] <= 1e-12).all())
)
stability_df.to_csv(OUT / "score_repeatability.csv", index=False)

df.to_csv(OUT / "pepbi_all_rows_with_native_scores.csv", index=False)

def summarize_corr(label, data, target, predictor, transform=1.0):
    z = data[[target, predictor]].copy().dropna()
    result = {
        "analysis": label,
        "target": target,
        "predictor": predictor,
        "n_rows": int(len(z)),
        "n_unique_sequences": int(
            data.loc[z.index, "_seq"].nunique()
        ) if len(z) else 0,
        "n_binding_groups": int(
            data.loc[z.index, "_group"].nunique()
        ) if len(z) else 0,
        "spearman_rho": None,
        "pearson_r": None,
    }
    if len(z) >= 3 and z[target].nunique() > 1 and z[predictor].nunique() > 1:
        # Both predictors below are oriented so higher means
        # stronger predicted binding; experimental ΔG is lower
        # for more favorable binding. A negative association is
        # therefore the nominally concordant direction.
        predicted_strength = transform * z[predictor]
        result["spearman_rho"] = float(
            predicted_strength.corr(z[target], method="spearman")
        )
        result["pearson_r"] = float(predicted_strength.corr(z[target]))
    return result

primary = df[df["_primary_eligible"]].copy()
corrs = []

for target, label in [
    ("_calc_dg", "calculated_DG_from_KD"),
    ("_direct_dg", "reported_experimental_DG"),
]:
    corrs.append(summarize_corr(
        f"{label}: native sequence prior",
        primary, target, "_c_score",
    ))
    corrs.append(summarize_corr(
        f"{label}: Rosetta RIA",
        primary, target, "_ria_dg", transform=-1.0,
    ))

# Independent of differences in mean affinity between targets, examine
# relative differences between peptide sequences within each binding group.
for target, label in [
    ("_calc_dg", "calculated_DG_from_KD"),
    ("_direct_dg", "reported_experimental_DG"),
]:
    available = primary[
        primary[target].notna() & (primary["_group"] != "")
    ].copy()

    grouped = (
        available.groupby(["_group", "_seq"], as_index=False)
        .agg(
            target_dg=(target, "mean"),
            native_score=("_c_score", "first"),
            ria_dg=("_ria_dg", "mean"),
            n_complex_rows=("_seq", "size"),
        )
    )

    if grouped.empty:
        continue

    counts = grouped.groupby("_group")["_seq"].transform("nunique")
    within = grouped[counts >= 2].copy()

    if not within.empty:
        within["target_centered"] = (
            within["target_dg"]
            - within.groupby("_group")["target_dg"].transform("mean")
        )
        within["native_centered"] = (
            within["native_score"]
            - within.groupby("_group")["native_score"].transform("mean")
        )
        within["ria_strength"] = -within["ria_dg"]
        within["ria_centered"] = (
            within["ria_strength"]
            - within.groupby("_group")["ria_strength"].transform("mean")
        )

        for predictor, pred_label in [
            ("native_centered", "native sequence prior"),
            ("ria_centered", "Rosetta RIA"),
        ]:
            corrs.append(summarize_corr(
                f"{label}: within-group {pred_label}",
                within.rename(columns={
                    "target_centered": "_target_centered",
                }),
                "_target_centered", predictor,
            ))

        within.to_csv(
            OUT / f"within_group_{label}.csv", index=False
        )

corr_df = pd.DataFrame(corrs)
corr_df.to_csv(OUT / "correlation_results.csv", index=False)

summary = {
    "workbook_rows": int(len(df)),
    "canonical_sequence_rows": int(df["_canonical"].sum()),
    "unique_canonical_sequences": int(df.loc[df["_canonical"], "_seq"].nunique()),
    "unique_sequences_scored": int(len(scores)),
    "score_failures": failures,
    "primary_length_range": [9, 20],
    "primary_eligible_rows": int(df["_primary_eligible"].sum()),
    "primary_unique_sequences": int(
        df.loc[df["_primary_eligible"], "_seq"].nunique()
    ),
    "primary_groups": int(
        df.loc[df["_primary_eligible"], "_group"].nunique()
    ),
    "reported_delta_g_numeric_rows": int(df["_direct_dg"].notna().sum()),
    "calculated_delta_g_numeric_rows": int(df["_calc_dg"].notna().sum()),
    "ria_dg_separated_numeric_rows": int(df["_ria_dg"].notna().sum()),
    "score_repeatability_passed_for_tested_sequences": stable,
    "n_repeatability_sequences_tested": int(len(stability_df)),
    "production_bridge_modified": False,
    "target_specific_binding_objective": False,
    "correlations": corrs,
    "interpretation": (
        "Exploratory benchmark only. Native C fitness passes target=NULL. "
        "Rows within binding groups are related, and sequence-only scoring "
        "cannot respond to protein mutations for a fixed peptide. Correlations "
        "do not validate HERV-K Env binding or establish affinity prediction."
    ),
}

summary_path = OUT / "benchmark_summary.json"
summary_path.write_text(json.dumps(summary, indent=2))

print("===== PEPBI CORRECTED BENCHMARK =====")
print(json.dumps(summary, indent=2))
print("\n===== CORRELATIONS =====")
print(corr_df.to_string(index=False))
print("\n===== OUTPUT FILES =====")
for path in [
    OUT / "pepbi_all_rows_with_native_scores.csv",
    OUT / "score_repeatability.csv",
    OUT / "correlation_results.csv",
    summary_path,
]:
    print(f"{path} ({path.stat().st_size:,} bytes)")
