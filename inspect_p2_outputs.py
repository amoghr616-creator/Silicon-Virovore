from pathlib import Path
import pandas as pd

ROOT = Path("prototype2/p2_rapidock")

runs = {
    "P2-01": "P2-01-clean",
    "P2-02": "P2-02",
    "P2-03": "P2-03",
    "P2-04": "P2-04",
    "P2-05": "P2-05",
    "P2-06": "P2-06",
}

for candidate, folder in runs.items():
    print("\n" + "=" * 80)
    print(candidate)
    print("=" * 80)

    for filename in ["cluster_summary.csv", "ranked_poses.csv"]:
        path = ROOT / folder / filename

        if not path.exists():
            print(f"\nMISSING: {path}")
            continue

        df = pd.read_csv(path)

        print(f"\n--- {filename} ---")
        print("Rows:", len(df))
        print("Columns:")
        for col in df.columns:
            print(" ", col)

        print("\nFirst 3 rows:")
        print(df.head(3).to_string(index=False))

