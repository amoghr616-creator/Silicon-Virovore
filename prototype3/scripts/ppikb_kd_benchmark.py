from pathlib import Path
import ctypes
import json
import math
import re
import sys

import numpy as np
import pandas as pd

ROOT = Path.cwd()
SRC = ROOT / "prototype3/data/raw/PPIKB/PPIKB_Research.xlsx"
OUT = ROOT / "prototype3/results/ppikb_kd_benchmark"
OUT.mkdir(parents=True, exist_ok=True)

if not SRC.is_file():
    raise SystemExit(f"Workbook missing: {SRC}")

sys.path.insert(0, str(ROOT))
import c.bridge as bridge

df = pd.read_excel(SRC, header=0)
df.columns = [str(c).strip() for c in df.columns]

required = [
    "Protein_Name", "Protein_Sequence", "Peptide_Sequence",
    "Peptide_Length", "Chemical_Modification", "Cyclization_Method",
    "Linear/Cyclic", "N-terminal_Modification",
    "C-terminal_Modification", "Affinity", "Affinity_Assay",
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise SystemExit(
        f"Missing fields: {missing}\nAvailable fields: {df.columns.tolist()}"
    )

def norm_text(series):
    return series.fillna("").astype(str).str.strip()

def is_unmodified(series, allowed):
    return norm_text(series).str.casefold().isin(allowed)

# Preserve sequence case: lowercase letters may encode stereochemical
# information. Do not silently uppercase D-amino-acid annotations.
df["_peptide"] = norm_text(df["Peptide_Sequence"]).str.replace(
    r"\s+", "", regex=True
)
df["_protein"] = norm_text(df["Protein_Sequence"]).str.upper().str.replace(
    r"\s+", "", regex=True
)
df["_assay"] = norm_text(df["Affinity_Assay"])
df["_name"] = norm_text(df["Protein_Name"])
df["_length"] = df["_peptide"].str.len()
df["_reported_length"] = pd.to_numeric(df["Peptide_Length"], errors="coerce")

canonical = df["_peptide"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+")
df["_canonical"] = canonical.fillna(False)
df["_length_ok"] = df["_length"].between(9, 20)
df["_linear"] = norm_text(df["Linear/Cyclic"]).str.casefold().eq("linear")

none_values = {"", "none", "no", "n/a", "na", "false", "0", "unmodified"}
free_values = none_values | {"free"}

df["_unmodified"] = (
    is_unmodified(df["Chemical_Modification"], none_values)
    & is_unmodified(df["Cyclization_Method"], none_values)
    & is_unmodified(df["N-terminal_Modification"], free_values)
    & is_unmodified(df["C-terminal_Modification"], free_values)
)

# Include residue modifications if that field exists in this workbook.
if "Residue_Modification" in df.columns:
    df["_unmodified"] &= is_unmodified(
        df["Residue_Modification"], none_values
    )

# Parse explicit affinity types and units. Censored values such as Kd < 5 nM
# are retained in the audit counts but excluded from exact-value analysis.
affinity_re = re.compile(
    r"\b(Kd|Ki|IC50)\s*=\s*([<>≤≥]?)\s*"
    r"(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*"
    r"(pM|nM|uM|µM|μM|mM|M)\b",
    re.IGNORECASE,
)
unit_to_molar = {
    "pm": 1e-12, "nm": 1e-9, "um": 1e-6,
    "µm": 1e-6, "μm": 1e-6, "mm": 1e-3, "m": 1.0,
}

def parse_affinity(value):
    s = str(value).strip()
    m = affinity_re.search(s)
    if not m:
        return pd.Series({
            "_affinity_type": None,
            "_censored": None,
            "_value": np.nan,
            "_unit": None,
            "_molar": np.nan,
            "_p_affinity": np.nan,
        })

    kind, comparator, number, unit = m.groups()
    val = float(number)
    molar = val * unit_to_molar[unit.casefold()]

    return pd.Series({
        "_affinity_type": kind.upper(),
        "_censored": bool(comparator),
        "_value": val,
        "_unit": unit,
        "_molar": molar,
        "_p_affinity": -math.log10(molar) if molar > 0 else np.nan,
    })

parsed = df["Affinity"].apply(parse_affinity)
df = pd.concat([df, parsed], axis=1)

# Primary benchmark: exact, uncensored Kd values only.
df["_target_valid"] = df["_protein"].str.fullmatch(
    r"[ACDEFGHIKLMNPQRSTVWY]+"
).fillna(False)

df["_eligible"] = (
    df["_canonical"]
    & df["_length_ok"]
    & df["_linear"]
    & df["_unmodified"]
    & df["_target_valid"]
    & df["_affinity_type"].eq("KD")
    & df["_censored"].eq(False)
    & df["_p_affinity"].notna()
    & df["_assay"].ne("")
)

eligible = df.loc[df["_eligible"]].copy()

# Call the native C API directly so that the production bridge's fixed
# 21-residue Python validation remains unchanged. Only score peptides within
# the preselected 9–20 residue range.
def native_score(sequence):
    metrics = bridge.NativeFitnessMetrics()
    status = bridge.lib.c_evaluate_sequence(
        sequence.encode("ascii"), ctypes.byref(metrics)
    )
    if status != 0:
        raise RuntimeError(
            f"C API returned status {status} for sequence length {len(sequence)}"
        )
    score = float(metrics.fitness_score)
    if not math.isfinite(score) or not 0.0 <= score <= 1.0:
        raise RuntimeError(f"Invalid score {score!r}")
    return score

seqs = sorted(eligible["_peptide"].unique())
scores = {}
failures = []

for seq in seqs:
    try:
        scores[seq] = native_score(seq)
    except Exception as exc:
        failures.append({"sequence": seq, "error": str(exc)})

if failures:
    print("Scoring failures:", json.dumps(failures[:10], indent=2))

eligible["_native_score"] = eligible["_peptide"].map(scores)
eligible = eligible[eligible["_native_score"].notna()].copy()

# Use target sequence + assay as the comparison context. Replicate
# measurements are combined for each target/assay/peptide combination.
group_cols = ["_protein", "_assay", "_peptide"]
collapsed = (
    eligible.groupby(group_cols, as_index=False)
    .agg(
        pKd=("_p_affinity", "median"),
        native_score=("_native_score", "first"),
        measurement_rows=("_p_affinity", "size"),
        protein_name=("_name", "first"),
    )
)

group_results = []
for (protein, assay), block in collapsed.groupby(["_protein", "_assay"]):
    block = block.dropna(subset=["pKd", "native_score"])
    n_seq = int(block["_peptide"].nunique())

    if n_seq < 3 or block["pKd"].nunique() < 2 or block["native_score"].nunique() < 2:
        continue

    group_results.append({
        "protein_sequence": protein,
        "protein_name": block["protein_name"].iloc[0],
        "assay": assay,
        "distinct_peptides": n_seq,
        "measurement_rows": int(block["measurement_rows"].sum()),
        "spearman_rho": float(
            block["native_score"].corr(block["pKd"], method="spearman")
        ),
        "pearson_r": float(block["native_score"].corr(block["pKd"])),
    })

groups = pd.DataFrame(group_results)
groups.to_csv(OUT / "target_assay_correlations.csv", index=False)
collapsed.to_csv(OUT / "target_peptide_kd_summary.csv", index=False)
eligible.to_csv(OUT / "filtered_kd_measurements.csv", index=False)

def group_summary(min_peptides, seed):
    g = groups[groups["distinct_peptides"] >= min_peptides]
    values = g["spearman_rho"].dropna().to_numpy()
    if not len(values):
        return {
            "minimum_distinct_peptides_per_group": min_peptides,
            "independent_target_assay_groups": 0,
        }

    rng = np.random.default_rng(seed)
    boot = np.array([
        rng.choice(values, size=len(values), replace=True).mean()
        for _ in range(20000)
    ])
    return {
        "minimum_distinct_peptides_per_group": min_peptides,
        "independent_target_assay_groups": int(len(values)),
        "mean_group_spearman": float(values.mean()),
        "median_group_spearman": float(np.median(values)),
        "groups_positive": int((values > 0).sum()),
        "groups_negative": int((values < 0).sum()),
        "bootstrap_95ci_low": float(np.quantile(boot, 0.025)),
        "bootstrap_95ci_high": float(np.quantile(boot, 0.975)),
    }

kind_counts = (
    df["_affinity_type"].fillna("UNPARSED").value_counts().to_dict()
)

summary = {
    "input_file": str(SRC),
    "input_rows": int(len(df)),
    "rows_with_canonical_peptide": int(df["_canonical"].sum()),
    "rows_with_9_to_20_residue_peptides": int(
        (df["_canonical"] & df["_length_ok"]).sum()
    ),
    "rows_unmodified_and_linear": int(
        (df["_canonical"] & df["_length_ok"] & df["_linear"] & df["_unmodified"]).sum()
    ),
    "parsed_affinity_types": {str(k): int(v) for k, v in kind_counts.items()},
    "exact_uncensored_Kd_rows_before_scoring": int(
        (
            df["_canonical"] & df["_length_ok"] & df["_linear"]
            & df["_unmodified"] & df["_target_valid"]
            & df["_affinity_type"].eq("KD")
            & df["_censored"].eq(False)
            & df["_p_affinity"].notna()
        ).sum()
    ),
    "eligible_rows_after_successful_native_scoring": int(len(eligible)),
    "distinct_peptides_scored": int(len(scores) - len(failures)),
    "native_score_failures": failures,
    "distinct_target_assay_groups_with_correlations": int(len(groups)),
    "group_level_results_min_3_peptides": group_summary(3, 20261026),
    "group_level_results_min_5_peptides": group_summary(5, 20261027),
    "score_semantics": (
        "Existing sequence-property heuristic; C API receives peptide sequence "
        "only and passes target=NULL. Higher pKd means stronger measured affinity."
    ),
    "interpretation": (
        "Exploratory external benchmark. Correlations are computed within target "
        "sequence and assay, with replicate measurements collapsed per peptide. "
        "Bootstrap intervals resample target-assay groups, not raw rows. This does "
        "not validate target-specific binding prediction."
    ),
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2))

print("===== PPIKB KD BENCHMARK =====")
print(json.dumps(summary, indent=2))
print("\n===== TARGET / ASSAY CORRELATIONS =====")
print(groups.to_string(index=False) if len(groups) else "No eligible groups.")
print("\n===== OUTPUT FILES =====")
for p in [
    OUT / "filtered_kd_measurements.csv",
    OUT / "target_peptide_kd_summary.csv",
    OUT / "target_assay_correlations.csv",
    OUT / "summary.json",
]:
    print(f"{p} ({p.stat().st_size:,} bytes)")
