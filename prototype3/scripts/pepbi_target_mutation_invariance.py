from pathlib import Path
from itertools import combinations
import json

import pandas as pd

ROOT = Path.cwd()
SRC = ROOT / "prototype3/results/pepbi_corrected_benchmark/pepbi_all_rows_with_native_scores.csv"
OUT = ROOT / "prototype3/results/pepbi_target_mutation_invariance"
OUT.mkdir(parents=True, exist_ok=True)

if not SRC.exists():
    raise SystemExit(f"Missing scored data: {SRC}")

df = pd.read_csv(SRC)

def numeric_column(candidates, output):
    for name in candidates:
        if name in df.columns:
            df[output] = pd.to_numeric(df[name], errors="coerce")
            return
    raise SystemExit(f"Could not find any of {candidates}")

def text_column(candidates, output, default="<missing>"):
    for name in candidates:
        if name in df.columns:
            df[output] = df[name].fillna(default).astype(str).str.strip()
            return
    df[output] = default

text_column(["_seq", "Peptide Sequence"], "seq")
text_column(["_group", "Binding Group"], "group")
text_column(["Protein Mutation(s)"], "protein_mut")
text_column(["Protein Sequence"], "protein_seq")
text_column(["Binding Site"], "binding_site")
numeric_column(["_temp", "T (K)"], "temperature_K")
numeric_column(["_direct_dg", "ΔG (kcal/mol)"], "experimental_dg")
numeric_column(["_calc_dg", "Calculated ΔG (kcal/mol)"], "calculated_dg")
numeric_column(["_c_score"], "native_score")
numeric_column(["_ria_dg", "dG_separated"], "rosetta_dg")

df["seq"] = df["seq"].str.upper().str.replace(r"\s+", "", regex=True)
df["protein_seq"] = df["protein_seq"].str.upper().str.replace(r"\s+", "", regex=True)

valid = (
    df["seq"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+").fillna(False)
    & df["experimental_dg"].notna()
    & df["native_score"].notna()
    & df["temperature_K"].notna()
    & df["group"].ne("<missing>")
)
df = df.loc[valid].copy()

# Match the same peptide, binding system, temperature and site.
# Protein mutation/sequence is deliberately allowed to differ.
keys = ["group", "seq", "temperature_K", "binding_site"]
context = (
    df.groupby(keys + ["protein_mut", "protein_seq"], dropna=False, as_index=False)
    .agg(
        experimental_dg=("experimental_dg", "mean"),
        calculated_dg=("calculated_dg", "mean"),
        native_score=("native_score", "first"),
        rosetta_dg=("rosetta_dg", "mean"),
        source_rows=("seq", "size"),
    )
)

pairs = []
invariance = []

for key, block in context.groupby(keys, dropna=False):
    variants = block.to_dict("records")
    if len(variants) < 2:
        continue

    scores = [r["native_score"] for r in variants]
    dg_values = [r["experimental_dg"] for r in variants]
    spread = max(scores) - min(scores)

    invariance.append({
        "binding_group": key[0],
        "peptide_sequence": key[1],
        "temperature_K": key[2],
        "binding_site": key[3],
        "length": len(key[1]),
        "protein_contexts": len(variants),
        "experimental_dg_min": min(dg_values),
        "experimental_dg_max": max(dg_values),
        "experimental_dg_range": max(dg_values) - min(dg_values),
        "native_score_min": min(scores),
        "native_score_max": max(scores),
        "native_score_range": spread,
        "score_invariant_to_protein_context": spread <= 1e-12,
    })

    for a, b in combinations(variants, 2):
        if a["protein_mut"] == b["protein_mut"] and a["protein_seq"] == b["protein_seq"]:
            continue
        pairs.append({
            "binding_group": key[0],
            "peptide_sequence": key[1],
            "peptide_length": len(key[1]),
            "temperature_K": key[2],
            "binding_site": key[3],
            "protein_mutation_A": a["protein_mut"],
            "protein_mutation_B": b["protein_mut"],
            "protein_sequence_A": a["protein_seq"],
            "protein_sequence_B": b["protein_seq"],
            "experimental_dg_A": a["experimental_dg"],
            "experimental_dg_B": b["experimental_dg"],
            "delta_experimental_dg_B_minus_A": (
                b["experimental_dg"] - a["experimental_dg"]
            ),
            "calculated_dg_A": a["calculated_dg"],
            "calculated_dg_B": b["calculated_dg"],
            "native_score_A": a["native_score"],
            "native_score_B": b["native_score"],
            "native_score_difference_B_minus_A": (
                b["native_score"] - a["native_score"]
            ),
            "rosetta_dg_A": a["rosetta_dg"],
            "rosetta_dg_B": b["rosetta_dg"],
            "source_rows_A": a["source_rows"],
            "source_rows_B": b["source_rows"],
        })

pairs_df = pd.DataFrame(pairs)
invariance_df = pd.DataFrame(invariance)

pairs_df.to_csv(OUT / "same_peptide_target_mutation_pairs.csv", index=False)
invariance_df.to_csv(OUT / "target_context_invariance_summary.csv", index=False)

summary = {
    "input": str(SRC),
    "measured_rows_examined": int(len(df)),
    "matched_peptide_contexts": int(len(invariance_df)),
    "matched_target_variant_pairs": int(len(pairs_df)),
    "pairs_with_native_score_difference_over_1e-12": (
        int((pairs_df["native_score_difference_B_minus_A"].abs() > 1e-12).sum())
        if not pairs_df.empty else 0
    ),
    "pairs_with_experimental_delta_g_difference_over_0_5_kcal_mol": (
        int((pairs_df["delta_experimental_dg_B_minus_A"].abs() > 0.5).sum())
        if not pairs_df.empty else 0
    ),
    "largest_observed_experimental_delta_g_difference": (
        float(pairs_df["delta_experimental_dg_B_minus_A"].abs().max())
        if not pairs_df.empty else None
    ),
    "interpretation": (
        "This measures whether a sequence-only score is invariant to receptor "
        "context. It is not a test of causal mutation effects, binding prediction "
        "accuracy, or HERV-K Env interaction."
    ),
}

(OUT / "summary.json").write_text(json.dumps(summary, indent=2))

print("===== TARGET-MUTATION INVARIANCE AUDIT =====")
print(json.dumps(summary, indent=2))
print("\n===== MATCHED PAIRS =====")
print(pairs_df.to_string(index=False) if not pairs_df.empty else "No matched pairs found.")
print("\n===== OUTPUTS =====")
for p in [
    OUT / "same_peptide_target_mutation_pairs.csv",
    OUT / "target_context_invariance_summary.csv",
    OUT / "summary.json",
]:
    print(f"{p} ({p.stat().st_size:,} bytes)")
