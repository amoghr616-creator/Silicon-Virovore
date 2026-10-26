from pathlib import Path
from itertools import combinations
from math import comb
import json

import numpy as np
import pandas as pd

ROOT = Path.cwd()
SRC = ROOT / "prototype3/results/pepbi_corrected_benchmark/pepbi_all_rows_with_native_scores.csv"
OUT = ROOT / "prototype3/results/pepbi_strict_matched"
OUT.mkdir(parents=True, exist_ok=True)

if not SRC.exists():
    raise SystemExit(f"Missing input: {SRC}")

df = pd.read_csv(SRC)

# The scored CSV retains the original workbook temperature column.
# Restore the normalized field expected by the matching code.
if "_temp" not in df.columns and "T (K)" in df.columns:
    df["_temp"] = pd.to_numeric(df["T (K)"], errors="coerce")

needed = [
    "_seq", "_length", "_group", "_temp", "_c_score",
    "_calc_dg", "_direct_dg", "_ria_dg",
    "Protein Sequence", "Protein Mutation(s)", "Binding Site",
]
missing = [c for c in needed if c not in df.columns]
if missing:
    raise SystemExit(f"Missing columns: {missing}")

for col in ["_length", "_temp", "_c_score", "_calc_dg",
            "_direct_dg", "_ria_dg"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# Normalize context annotations without modifying source values on disk.
for col in ["_group", "Protein Sequence", "Protein Mutation(s)", "Binding Site"]:
    df[col] = df[col].fillna("<missing>").astype(str).str.strip()

df["Protein Sequence"] = df["Protein Sequence"].str.replace(
    r"\s+", "", regex=True
).str.upper()

valid = (
    df["_seq"].fillna("").str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+")
    & df["_length"].between(9, 20)
    & df["_c_score"].notna()
    & df["_temp"].notna()
    & df["_group"].ne("<missing>")
    & df["Protein Sequence"].ne("<missing>")
    & df["Protein Sequence"].ne("X")
)
df = df.loc[valid].copy()

# Exact matching variables. Peptide sequence is intentionally NOT a key,
# because different peptide sequences are the variants being compared.
keys = [
    "_group", "Protein Sequence", "Protein Mutation(s)",
    "Binding Site", "_temp", "_length",
]

seq = (
    df.groupby(keys + ["_seq"], dropna=False, as_index=False)
      .agg(
          native_score=("_c_score", "first"),
          calculated_dg=("_calc_dg", "mean"),
          experimental_dg=("_direct_dg", "mean"),
          rosetta_dg=("_ria_dg", "mean"),
          source_rows=("_seq", "size"),
      )
)
seq["_context_id"] = seq.groupby(keys, dropna=False, sort=False).ngroup()

pair_rows = []
target_cols = {
    "calculated_dg": "Calculated ΔG from K_D",
    "experimental_dg": "Reported experimental ΔG",
}
predictors = ["native_score", "rosetta_dg"]
target_ties = {name: 0 for name in target_cols}
contexts_with_pairs = {name: set() for name in target_cols}

for context_id, block in seq.groupby("_context_id", sort=False):
    records = block.to_dict("records")
    if len(records) < 2:
        continue

    for a, b in combinations(records, 2):
        for target, target_label in target_cols.items():
            ta, tb = a[target], b[target]
            if pd.isna(ta) or pd.isna(tb):
                continue

            delta_target = float(ta - tb)
            if abs(delta_target) < 1e-9:
                target_ties[target] += 1
                continue

            for predictor in predictors:
                pa, pb = a[predictor], b[predictor]
                if pd.isna(pa) or pd.isna(pb):
                    continue

                delta_predictor = float(pa - pb)
                if abs(delta_predictor) < 1e-12:
                    concordance = 0.5
                elif predictor == "native_score":
                    # Higher native score should correspond to lower ΔG.
                    concordance = float(delta_predictor * delta_target < 0)
                else:
                    # Lower Rosetta dG_separated should correspond to lower ΔG.
                    concordance = float(delta_predictor * delta_target > 0)

                row = {
                    "context_id": int(context_id),
                    "binding_group": a["_group"],
                    "protein_sequence": a["Protein Sequence"],
                    "protein_mutations": a["Protein Mutation(s)"],
                    "binding_site": a["Binding Site"],
                    "temperature_K": a["_temp"],
                    "peptide_length": int(a["_length"]),
                    "sequence_A": a["_seq"],
                    "sequence_B": b["_seq"],
                    "target": target_label,
                    "predictor": predictor,
                    "delta_target_A_minus_B": delta_target,
                    "delta_predictor_A_minus_B": delta_predictor,
                    "pairwise_concordance": concordance,
                }
                pair_rows.append(row)
                contexts_with_pairs[target].add(int(context_id))

pair_cols = [
    "context_id", "binding_group", "protein_sequence",
    "protein_mutations", "binding_site", "temperature_K",
    "peptide_length", "sequence_A", "sequence_B", "target",
    "predictor", "delta_target_A_minus_B",
    "delta_predictor_A_minus_B", "pairwise_concordance",
]
pairs = pd.DataFrame(pair_rows, columns=pair_cols)
pairs.to_csv(OUT / "matched_pair_predictions.csv", index=False)

summary_rows = []
context_rows = []

if not pairs.empty:
    context = (
        pairs.groupby(
            ["target", "predictor", "context_id", "binding_group"],
            as_index=False,
        )
        .agg(
            concordance=("pairwise_concordance", "mean"),
            evaluated_pairs=("pairwise_concordance", "size"),
        )
    )
    context.to_csv(OUT / "matched_context_summary.csv", index=False)

    # First average contexts within a binding group. Then bootstrap whole
    # binding groups, so groups—not individual pairwise comparisons—are
    # the uncertainty units.
    group = (
        context.groupby(["target", "predictor", "binding_group"], as_index=False)
        .agg(
            group_concordance=("concordance", "mean"),
            matched_contexts=("context_id", "nunique"),
            evaluated_pairs=("evaluated_pairs", "sum"),
        )
    )
    group.to_csv(OUT / "matched_group_summary.csv", index=False)

    rng = np.random.default_rng(20261026)
    for (target, predictor), block in group.groupby(["target", "predictor"]):
        values = block["group_concordance"].to_numpy()
        n = len(values)
        boot = np.array([
            rng.choice(values, size=n, replace=True).mean()
            for _ in range(20000)
        ])
        non_ties = values[np.abs(values - 0.5) > 1e-12]
        wins = int((non_ties > 0.5).sum())
        losses = int((non_ties < 0.5).sum())
        m = wins + losses
        p_one_sided = (
            sum(comb(m, k) for k in range(wins, m + 1)) / (2 ** m)
            if m else None
        )

        summary_rows.append({
            "target": target,
            "predictor": predictor,
            "independent_binding_groups": n,
            "mean_group_concordance": float(values.mean()),
            "median_group_concordance": float(np.median(values)),
            "groups_above_chance": wins,
            "groups_below_chance": losses,
            "group_sign_test_one_sided_p": p_one_sided,
            "group_bootstrap_95ci_low": float(np.quantile(boot, 0.025)),
            "group_bootstrap_95ci_high": float(np.quantile(boot, 0.975)),
            "total_matched_pairs": int(block["evaluated_pairs"].sum()),
        })

summary = {
    "source": str(SRC),
    "filter": "canonical peptides length 9-20; exact match on binding group, protein sequence, protein mutation annotation, binding site, temperature, and peptide length",
    "eligible_rows": int(len(df)),
    "distinct_sequences": int(df["_seq"].nunique()),
    "strict_matched_contexts_with_pairwise_data": int(
        pairs["context_id"].nunique() if not pairs.empty else 0
    ),
    "pairwise_comparisons": int(len(pairs)),
    "target_ties_excluded": target_ties,
    "results": summary_rows,
    "interpretation": (
        "Pairwise ranking diagnostic, not target-binding validation. "
        "A lack of matched contexts means this dataset cannot estimate "
        "within-context sequence-ranking performance under these strict rules."
    ),
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2))

print("===== STRICT MATCHED BENCHMARK =====")
print(json.dumps(summary, indent=2))
print("\nFiles:", OUT.resolve())
