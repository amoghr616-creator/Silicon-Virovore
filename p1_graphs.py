import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

runs = {
    "Baseline": "results/ablation-baseline-619/candidate_history.json",
    "ARISE-only": "results/ablation-arise-only-619/candidate_history.json",
    "Hotspot-only": "results/ablation-hotspot-only-619/candidate_history.json",
    "Adaptive": "results/ablation-adaptive-619/candidate_history.json",
}

out = Path("results/p1_graphs")
out.mkdir(parents=True, exist_ok=True)

data = {}

for name, path in runs.items():
    with open(path) as f:
        rows = json.load(f)

    df = pd.DataFrame(rows)

    if "posthoc_run_normalized_score" not in df.columns:
        raise RuntimeError(f"{name}: missing posthoc_run_normalized_score")

    data[name] = df.sort_values("generation")

# --------------------------------------------------
# Graph 1: ablation trajectories
# --------------------------------------------------

plt.figure(figsize=(10, 6))

for name, df in data.items():
    plt.plot(
        df["generation"],
        df["posthoc_run_normalized_score"],
        linewidth=2,
        label=name
    )

plt.xlabel("Generation")
plt.ylabel("Run-normalized score")
plt.title("Prototype 1: Ablation Search Trajectories")
plt.legend()
plt.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(out / "P1_ablation_trajectory.png", dpi=300)
plt.close()

# --------------------------------------------------
# Graph 2: final ablation score
# --------------------------------------------------

summary = []

for name, df in data.items():
    final = df.iloc[-1]

    summary.append({
        "method": name,
        "final_score": final["posthoc_run_normalized_score"],
        "max_score": df["posthoc_run_normalized_score"].max(),
        "max_score_generation": int(
            df.loc[
                df["posthoc_run_normalized_score"].idxmax(),
                "generation"
            ]
        ),
        "final_native_fitness": final["native_fitness"],
        "final_sequence": final["sequence"],
    })

summary_df = pd.DataFrame(summary)
summary_df.to_csv(out / "P1_ablation_summary.csv", index=False)

plt.figure(figsize=(9, 5))

plt.bar(
    summary_df["method"],
    summary_df["final_score"]
)

plt.xlabel("Method")
plt.ylabel("Final run-normalized score")
plt.title("Prototype 1: Final Ablation Performance")
plt.grid(axis="y", alpha=0.25)
plt.tight_layout()
plt.savefig(out / "P1_ablation_final_score.png", dpi=300)
plt.close()

print("\nDONE")
print("Created:")
print(out / "P1_ablation_trajectory.png")
print(out / "P1_ablation_final_score.png")
print(out / "P1_ablation_summary.csv")
print("\nSummary:")
print(summary_df.to_string(index=False))
