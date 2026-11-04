import csv
import os
from collections import defaultdict
from pathlib import Path

ROOT = Path.cwd()
FAMILY_MAP = Path(
    "prototype4/results/skempi_audit_20261029/family_map.csv"
)
OUT = Path("prototype4/results/data_isolation_20261101")

with FAMILY_MAP.open(newline="", encoding="utf-8-sig") as f:
    complexes = list(csv.DictReader(f))

pdb_ids = {
    row["pdb_id"].strip().upper()
    for row in complexes
    if row.get("pdb_id", "").strip()
}

extensions = (".pdb", ".pdb.gz", ".cif", ".cif.gz",
              ".mmcif", ".mmcif.gz", ".ent", ".ent.gz")

skip = {
    ".git", ".venv", "venv", "__pycache__", "node_modules",
    "hybridock-pep", ".pytest_cache",
}

matches = defaultdict(list)
all_counts = defaultdict(int)

for current, dirs, files in os.walk(ROOT):
    dirs[:] = [
        d for d in dirs
        if d not in skip
        and not (
            d.startswith("haddock3_chain_")
            or d == "haddock_pose_audit_20261031"
        )
    ]

    for filename in files:
        lower = filename.lower()
        if not lower.endswith(extensions):
            continue

        path = Path(current) / filename
        relative = path.relative_to(ROOT)
        all_counts[relative.parts[0]] += 1

        stem = lower
        for suffix in (".gz", ".pdb", ".cif", ".mmcif", ".ent"):
            if stem.endswith(suffix):
                stem = stem[:-len(suffix)]

        # Match PDB identifiers at the beginning of structure filenames.
        for pdb_id in pdb_ids:
            if stem == pdb_id.lower() or stem.startswith(pdb_id.lower() + "_"):
                matches[pdb_id].append(str(relative))

OUT.mkdir(parents=True, exist_ok=True)

report = []
for pdb_id in sorted(pdb_ids):
    paths = sorted(set(matches[pdb_id]))
    report.append({
        "pdb_id": pdb_id,
        "available": bool(paths),
        "matching_file_count": len(paths),
        "paths": ";".join(paths),
    })

with (OUT / "skempi_structure_inventory.tsv").open(
    "w", newline="", encoding="utf-8"
) as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["pdb_id", "available", "matching_file_count", "paths"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(report)

available = sum(bool(x["available"]) for x in report)

print("SKEMPI STRUCTURE INVENTORY")
print("Complexes in family map:", len(complexes))
print("Unique PDB identifiers:", len(pdb_ids))
print("PDBs with local structure files:", available)
print("PDBs missing local structure files:", len(pdb_ids) - available)

print("\nSTRUCTURE FILE COUNTS BY TOP-LEVEL DIRECTORY")
for name, count in sorted(all_counts.items()):
    print(f"{name}: {count}")

print("\nAVAILABLE SKEMPI STRUCTURES")
for row in report:
    if row["available"]:
        print(row["pdb_id"], row["paths"])

print("\nMISSING PDB IDENTIFIERS")
missing = [x["pdb_id"] for x in report if not x["available"]]
print(", ".join(missing) if missing else "None")

print("\nInventory saved:", OUT / "skempi_structure_inventory.tsv")
