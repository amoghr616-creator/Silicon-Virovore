import matplotlib.pyplot as plt

# Initialize lists to hold our data
time = []
rmsd = []

# Parse the GROMACS xvg file
with open("rmsd.xvg", "r") as f:
    for line in f:
        # Skip comment and metadata lines
        if line.startswith(("@", "#")):
            continue
        parts = line.split()
        if len(parts) == 2:
            time.append(float(parts[0]))
            rmsd.append(float(parts[1]))

# Create a highly professional, publication-ready plot
plt.figure(figsize=(8, 5))
plt.plot(time, rmsd, color='#1f77b4', linewidth=1.5, label='Backbone RMSD')

# Style the chart
plt.title("Virovore Peptide Prototype 1.0 - Structural Stability Test", fontsize=14, fontweight='bold', pad=15)
plt.xlabel("Simulation Time (ns)", fontsize=12, fontweight='bold')
plt.ylabel("RMSD (nm)", fontsize=12, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.6)
plt.xlim(0, 10)
plt.ylim(0, max(rmsd) * 1.2 if rmsd else 1.0)
plt.legend(loc='upper left', frameon=True, facecolor='white', edgecolor='none')
plt.tight_layout()

# Save the plot for your engineering notebook
plt.savefig("rmsd_plot.png", dpi=300)
print("Success! Plot saved as 'rmsd_plot.png'")
plt.show()