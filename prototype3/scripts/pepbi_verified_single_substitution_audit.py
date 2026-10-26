from pathlib import Path
from itertools import combinations
import json
import pandas as pd
import numpy as np

ROOT = Path.cwd()
SRC = ROOT / "prototype3/results/pepbi_corrected_benchmark/pepbi_all_rows_with_native_scores.csv"
OUT = ROOT / "prototype3/results/pepbi_verified_single_substitution"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(SRC)

def text(names, out, default="X"):
    for n in names:
        if n in df.columns:
            df[out] = df[n].fillna(default).astype(str).str.strip()
            return
    df[out] = default

def numeric(names, out):
    for n in names:
        if n in df.columns:
            df[out] = pd.to_numeric(df[n], errors="coerce")
            return
    df[out] = np.nan

text(["_seq", "Peptide Sequence"], "peptide")
text(["_group", "Binding Group"], "group", "<missing>")
text(["Protein Sequence"], "protein")
text(["Protein Mutation(s)"], "mutation")
text(["Binding Site"], "site")
numeric(["_temp", "T (K)"], "temperature_K")
numeric(["_direct_dg", "ΔG (kcal/mol)"], "experimental_dg")
numeric(["_calc_dg", "Calculated ΔG (kcal/mol)"], "calculated_dg")
numeric(["_c_score"], "native_score")
numeric(["_ria_dg", "dG_separated"], "rosetta_dg")

df["peptide"] = df["peptide"].str.upper().str.replace(r"\s+", "", regex=True)
df["protein"] = df["protein"].str.upper().str.replace(r"\s+", "", regex=True)

df = df[
    df["peptide"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+").fillna(False)
    & df["protein"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+").fillna(False)
    & df["experimental_dg"].notna()
    & df["native_score"].notna()
    & df["temperature_K"].notna()
    & df["group"].ne("<missing>")
].copy()

# Collapse duplicate records of the same peptide/protein context.
keys = [
    "group", "peptide", "temperature_K", "site",
    "protein", "mutation",
]
contexts = (
    df.groupby(keys, dropna=False, as_index=False)
    .agg(
        experimental_dg=("experimental_dg", "mean"),
        calculated_dg=("calculated_dg", "mean"),
        native_score=("native_score", "first"),
        rosetta_dg=("rosetta_dg", "mean"),
        source_rows=("peptide", "size"),
    )
)

# Combine annotation duplicates that map to the identical protein sequence.
seq_keys = ["group", "peptide", "temperature_K", "site", "protein"]
contexts = (
    contexts.groupby(seq_keys, dropna=False, as_index=False)
    .agg(
        mutations=("mutation", lambda x: "; ".join(sorted(set(x)))),
        experimental_dg=("experimental_dg", "mean"),
        calculated_dg=("calculated_dg", "mean"),
        native_score=("native_score", "first"),
        rosetta_dg=("rosetta_dg", "mean"),
        source_rows=("source_rows", "sum"),
    )
)

pairs = []
for key, block in contexts.groupby(
    ["group", "peptide", "temperature_K", "site"], dropna=False
):
    records = block.to_dict("records")
    for a, b in combinations(records, 2):
        sa, sb = a["protein"], b["protein"]
        if len(sa) != len(sb):
            continue

        changes = [
            (i + 1, x, y)
            for i, (x, y) in enumerate(zip(sa, sb)) if x != y
        ]
        # Exact single-substitution filter: reject homologs and multi-mutants.
        if len(changes) != 1:
            continue

        pos, aa, bb = changes[0]
        dg_delta = a["experimental_dg"] - b["experimental_dg"]
        ns_delta = a["native_score"] - b["native_score"]
        r_delta = a["rosetta_dg"] - b["rosetta_dg"]

        pairs.append({
            "binding_group": key[0],
            "peptide_sequence": key[1],
            "peptide_length": len(key[1]),
            "within_documented_5_to_20_aa_range": 5 <= len(key[1]) <= 20,
            "temperature_K": key[2],
            "binding_site_annotation": key[3],
            "protein_mutation_annotation_A": a["mutations"],
            "protein_mutation_annotation_B": b["mutations"],
            "protein_sequence_A": sa,
            "protein_sequence_B": sb,
            "protein_sequence_difference_position_1based": pos,
            "protein_residue_A": aa,
            "protein_residue_B": bb,
            "protein_mutation_annotation_available": (
                a["mutations"] != "X" or b["mutations"] != "X"
            ),
            "experimental_dg_A": a["experimental_dg"],
            "experimental_dg_B": b["experimental_dg"],
            "experimental_dg_A_minus_B": dg_delta,
            "native_score_A": a["native_score"],
            "native_score_B": b["native_score"],
            "native_score_A_minus_B": ns_delta,
            "native_score_rank_concordant": (
                np.nan if abs(dg_delta) < 1e-12 or abs(ns_delta) < 1e-12
                else bool(dg_delta * ns_delta < 0)
            ),
            "rosetta_dg_A": a["rosetta_dg"],
            "rosetta_dg_B": b["rosetta_dg"],
            "rosetta_dg_A_minus_B": r_delta,
            "rosetta_rank_concordant": (
                np.nan if pd.isna(r_delta) or abs(dg_delta) < 1e-12
                or abs(r_delta) < 1e-12
                else bool(dg_delta * r_delta > 0)
            ),
        })

pair_df = pd.DataFrame(pairs)
pair_path = OUT / "single_residue_receptor_variant_pairs.csv"
pair_df.to_csv(pair_path, index=False)

def rate(s):
    s = s.dropna()
    if s.empty:
        return None
    return float(s.astype(float).mean())

summary = {
    "source_rows_with_experimental_dg_and_protein_sequence": int(len(df)),
    "distinct_peptide_protein_contexts": int(len(contexts)),
    "exact_single_residue_receptor_variant_pairs": int(len(pair_df)),
    "pairs_with_explicit_mutation_annotation": (
        int(pair_df["protein_mutation_annotation_available"].sum())
        if not pair_df.empty else 0
    ),
    "pairs_within_documented_peptide_length_range": (
        int(pair_df["within_documented_5_to_20_aa_range"].sum())
        if not pair_df.empty else 0
    ),
    "distinct_binding_groups": (
        int(pair_df["binding_group"].nunique()) if not pair_df.empty else 0
    ),
    "experimental_dg_non_tied_pairs": (
        int((pair_df["experimental_dg_A_minus_B"].abs() > 1e-12).sum())
        if not pair_df.empty else 0
    ),
    "native_score_identical_pair_fraction": (
        float((pair_df["native_score_A_minus_B"].abs() < 1e-12).mean())
        if not pair_df.empty else None
    ),
    "native_pairwise_rank_concordance": (
        rate(pair_df["native_score_rank_concordant"])
        if not pair_df.empty else None
    ),
    "rosetta_pairwise_rank_concordance": (
        rate(pair_df["rosetta_rank_concordant"])
        if not pair_df.empty else None
    ),
    "warning": (
        "Small exploratory audit, not a generalizable affinity-validation study. "
        "Protein sequence indexing may refer to a fragment; verify annotations "
        "against source metadata. A peptide-only score is expected to be identical "
        "when the peptide sequence is held fixed."
    ),
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2))

print("===== VERIFIED SINGLE-SUBSTITUTION AUDIT =====")
print(json.dumps(summary, indent=2))
print("\n===== PAIRS =====")
if pair_df.empty:
    print("No exact single-residue receptor variants found.")
else:
    columns = [
        "binding_group", "peptide_sequence", "peptide_length",
        "protein_mutation_annotation_A", "protein_mutation_annotation_B",
        "protein_sequence_difference_position_1based",
        "protein_residue_A", "protein_residue_B",
        "experimental_dg_A", "experimental_dg_B",
        "native_score_A", "native_score_B",
        "rosetta_dg_A", "rosetta_dg_B",
    ]
    print(pair_df[columns].to_string(index=False))
print(f"\nSaved: {pair_path}")
print(f"Saved: {OUT / 'summary.json'}")
