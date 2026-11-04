from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path.cwd()
OUT = (
    ROOT / "prototype4" / "results"
    / "p4_pose_reproducibility_20261103"
)

source = OUT / "pose_similarity_matrix.csv"
destination = OUT / "interface_similarity_heatmap.png"

with source.open("r", encoding="utf-8", newline="") as handle:
    rows = list(csv.reader(handle))

labels = rows[0][1:]
matrix = np.array([
    [float(value) if value else np.nan for value in row[1:]]
    for row in rows[1:]
], dtype=float)

if matrix.shape[0] != matrix.shape[1] or matrix.shape[0] != len(labels):
    raise ValueError("Similarity matrix is not square or labels do not match.")

fig, ax = plt.subplots(figsize=(10, 8))
image = ax.imshow(matrix, vmin=0, vmax=1, aspect="equal")

midpoint = sum(label.startswith("1917:") for label in labels)
if 0 < midpoint < len(labels):
    ax.axvline(midpoint - 0.5, linewidth=2)
    ax.axhline(midpoint - 0.5, linewidth=2)

ax.set_xticks(range(len(labels)))
ax.set_yticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
ax.set_yticklabels(labels, fontsize=8)

ax.set_title("HADDOCK Interface Fingerprint Similarity")
ax.set_xlabel("Pose representative")
ax.set_ylabel("Pose representative")

colorbar = fig.colorbar(image, ax=ax)
colorbar.set_label("Jaccard similarity of residue-pair contacts")

fig.tight_layout()
fig.savefig(destination, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {destination}")
