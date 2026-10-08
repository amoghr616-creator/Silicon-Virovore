from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("prototype2/p2_rapidock")

RUNS = {
    "P2-01": ROOT / "P2-01-clean",
    "P2-02": ROOT / "P2-02",
    "P2-03": ROOT / "P2-03",
    "P2-04": ROOT / "P2-04",
    "P2-05": ROOT / "P2-05",
    "P2-06": ROOT / "P2-06",
}

OUT = Path("prototype2")
rows = []

for candidate, run_dir in RUNS.items():

    ranked_path = run_dir / "ranked_poses.csv"
    cluster_path = run_dir / "cluster_summary.csv"

    if not ranked_path.exists():
        print(f"WARNING: missing {ranked_path}")
        continue

    df = pd.read_csv(ranked_path)

    print(f"\n{candidate}")
    print("Columns:", ", ".join(df.columns))

    row = {
        "candidate": candidate,
        "n_rows": len(df),
    }

    # Number of scored poses
    row["n_scored"] = len(df)

    # ----------------------------
    # ΔG
    # ----------------------------
    if "delta_g" in df.columns:
        dg = pd.to_numeric(df["delta_g"], errors="coerce").dropna()
        row["best_delta_g"] = dg.min() if len(dg) else np.nan
        row["median_top5_delta_g"] = (
            dg.nsmallest(min(5, len(dg))).median()
            if len(dg) else np.nan
        )
        row["median_delta_g_all"] = dg.median() if len(dg) else np.nan

    # ----------------------------
    # rank score
    # ----------------------------
    if "rank_score" in df.columns:
        rs = pd.to_numeric(df["rank_score"], errors="coerce").dropna()
        row["best_rank_score"] = rs.min() if len(rs) else np.nan
        row["median_top5_rank_score"] = (
            rs.nsmallest(min(5, len(rs))).median()
            if len(rs) else np.nan
        )
        row["median_rank_score_all"] = rs.median() if len(rs) else np.nan

    # ----------------------------
    # Interface metrics
    # ----------------------------
    for col in ["bsa", "contacts", "charged_frustration", "n_clash"]:
        if col in df.columns:
            x = pd.to_numeric(df[col], errors="coerce").dropna()

            if len(x):
                row[f"median_{col}"] = x.median()
                row[f"max_{col}"] = x.max()

    # ----------------------------
    # Charged confidence
    # ----------------------------
    if "charged_confidence" in df.columns:
        vals = df["charged_confidence"].dropna().astype(str)
        if len(vals):
            row["charged_confidence_mode"] = vals.mode().iloc[0]

    # ----------------------------
    # Initial clashed flag
    # ----------------------------
    if "is_clashed" in df.columns:
        vals = df["is_clashed"].fillna(False).astype(bool)
        row["fraction_initially_clashed"] = vals.mean()

    # ----------------------------
    # Cluster information
    # ----------------------------
    if cluster_path.exists():
        try:
            cdf = pd.read_csv(cluster_path)
            row["n_clusters"] = len(cdf)
        except Exception:
            row["n_clusters"] = np.nan

    rows.append(row)

summary = pd.DataFrame(rows)

# Numerical formatting
numeric_cols = summary.select_dtypes(include=[np.number]).columns
summary[numeric_cols] = summary[numeric_cols].round(4)

# Sort by robust top-5 ΔG, not the single most negative pose
if "median_top5_delta_g" in summary.columns:
    summary = summary.sort_values("median_top5_delta_g", ascending=True)

out_file = OUT / "p2_docking_summary.csv"
summary.to_csv(out_file, index=False)

print("\n" + "=" * 80)
print("P2 DOCKING SUMMARY")
print("=" * 80)
print(summary.to_string(index=False))
print("\nSaved to:")
print(out_file)
