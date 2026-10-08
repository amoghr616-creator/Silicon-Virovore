from pathlib import Path
import pandas as pd


ROOT = Path("prototype2/p2_rapidock")

rows = []

for candidate_dir in sorted(ROOT.iterdir()):

    if not candidate_dir.is_dir():
        continue

    csv_path = candidate_dir / "ranked_poses.csv"

    if not csv_path.exists():
        continue

    df = pd.read_csv(csv_path)

    # Reject obviously unusable poses first.
    usable = df[
        (df["is_clashed"] == False) &
        (df["is_clipped"] == False)
    ].copy()

    if usable.empty:
        print(f"{candidate_dir.name}: no usable poses")
        continue

    # Do NOT use predicted ΔG alone.
    usable = usable.sort_values(
        ["rank_score", "n_contact_residues"],
        ascending=[True, False]
    )

    best = usable.iloc[0]

    rows.append({
        "candidate": candidate_dir.name,
        "pose": best["pose_filename"],
        "delta_g_kcal_mol": best["delta_g"],
        "rank_score": best["rank_score"],
        "contacts": best["n_contact_residues"],
        "clashed": best["is_clashed"],
        "confidence": best["charged_confidence"],
    })

out = pd.DataFrame(rows)

print("\n=== P2 candidate ranking ===")
print(out.to_string(index=False))

out.to_csv(
    ROOT / "p2_candidate_ranking.csv",
    index=False
)

print(
    f"\nSaved: {ROOT / 'p2_candidate_ranking.csv'}"
)
