from pathlib import Path
import csv

master_path = Path("prototype2/p2_master_results.csv")
metrics_path = Path("prototype2/p2_sequence_metrics.tsv")

# Read existing master table
with master_path.open() as f:
    master_rows = list(csv.DictReader(f))

# Read the sequence metrics you already calculated
with metrics_path.open() as f:
    metric_rows = list(csv.DictReader(f, delimiter="\t"))

metrics_by_candidate = {
    row["candidate"]: row
    for row in metric_rows
}

# Map existing metric names into the master table
mapping = {
    "length_aa": "length",
    "molecular_weight_da": "molecular_weight_da",
    "theoretical_pI": "pi",
    "charge_pH_5": "charge_pH5",
    "charge_pH_6": "charge_pH6",
    "charge_pH_7": "charge_pH7",
    "charge_pH_7_4": "charge_pH7_4",
    "charge_pH_8": "charge_pH8",
    "charge_pH_9": "charge_pH9",
    "gravy_hydropathy": "gravy",
    "eisenberg_hydrophobic_moment": "hydrophobic_moment",
    "hydrophobic_fraction": "hydrophobic_fraction",
    "max_hydrophobic_run": "max_hydrophobic_run",
    "all_L_canonical": "stereochemistry",
    "old_best_score": "ga_fitness",
    "runs_seen": "sequence_diversity",
}

for master in master_rows:
    candidate = master["candidate"]
    source = metrics_by_candidate[candidate]

    for source_field, master_field in mapping.items():
        if source_field in source and source[source_field] != "":
            value = source[source_field]

            if source_field == "all_L_canonical":
                value = "L-canonical amino acids" if value.lower() == "yes" else value

            master[master_field] = value

# Rewrite master table
fields = list(master_rows[0].keys())

with master_path.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(master_rows)

print(f"Updated {master_path}")
print(f"Merged {len(master_rows)} candidates")
