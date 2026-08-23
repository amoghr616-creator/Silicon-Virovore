# plot_dashboard.py
import pandas as pd
import matplotlib.pyplot as plt

# Load exported metrics from C run
df = pd.read_csv("metrics.csv")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

# 1. Trajectory of Peak Fitness
ax1.plot(df["Generation"], df["FitnessScore"], color="#D32F2F", linewidth=2, label="Peak Fitness")
ax1.set_ylabel("C Fitness Score", fontsize=11, fontweight="bold")
ax1.set_title("10,000-Generation Standalone C Engine Evolution", fontsize=13, fontweight="bold")
ax1.grid(True, linestyle="--", alpha=0.5)
ax1.legend(loc="lower right")

# 2. Structural Metrics (Helix Propensity & Eisenberg Moment)
ax2.plot(df["Generation"], df["HelixPropensity"], color="#00796B", linewidth=1.5, label="α-Helix Propensity")
ax2.plot(df["Generation"], df["Eisenberg"], color="#FFA000", linewidth=1.5, label="Hydrophobic Moment")
ax2.set_xlabel("Generation", fontsize=11, fontweight="bold")
ax2.set_ylabel("Biophysical Score", fontsize=11, fontweight="bold")
ax2.grid(True, linestyle="--", alpha=0.5)
ax2.legend(loc="lower right")

plt.tight_layout()
plt.savefig("c_engine_trajectory.png", dpi=300)
print("[EXPORT] Evolutionary dashboard saved to c_engine_trajectory.png")