import csv
from pathlib import Path

BASE = Path(
    "prototype4/results/p4_contact_occupancy_20261103"
)

def read_csv(name):
    with (BASE / name).open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))

pair_rows = read_csv("residue_pair_occupancy.csv")
residue_rows = read_csv("residue_occupancy.csv")
target_rows = read_csv("target_site_distances.csv")

SEEDS = ("1917", "2917")
CUTOFFS = (5.0, 8.0)

def jaccard(a, b):
    union = a | b
    return len(a & b) / len(union) if union else float("nan")

print("COMBINED CONTACT-SET OVERLAP")
print("These are unions across five representatives per seed.")
for cutoff in CUTOFFS:
    pair_sets = {}
    for seed in SEEDS:
        pair_sets[seed] = {
            (r["receptor_residue"], r["peptide_residue"])
            for r in pair_rows
            if r["seed"] == seed and float(r["cutoff_A"]) == cutoff
        }

    print(f"\nCutoff {cutoff:g} Å")
    print("  Unique receptor-peptide pairs:",
          {s: len(pair_sets[s]) for s in SEEDS})
    print("  Shared exact pairs:", len(pair_sets["1917"] & pair_sets["2917"]))
    print("  Pair-union Jaccard:",
          f"{jaccard(pair_sets['1917'], pair_sets['2917']):.4f}")

    for chain, label in (("A", "receptor"), ("B", "peptide")):
        residue_sets = {}
        for seed in SEEDS:
            residue_sets[seed] = {
                r["residue"]
                for r in residue_rows
                if r["seed"] == seed
                and float(r["cutoff_A"]) == cutoff
                and r["chain"] == chain
            }

        print(
            f"  {label} residue-union Jaccard: "
            f"{jaccard(residue_sets['1917'], residue_sets['2917']):.4f} "
            f"(shared={len(residue_sets['1917'] & residue_sets['2917'])})"
        )

    # Receptor residues recurring in at least 2/5 representatives
    by_seed = {}
    for seed in SEEDS:
        by_seed[seed] = {
            r["residue"]: int(r["poses_with_contact"])
            for r in residue_rows
            if r["seed"] == seed
            and float(r["cutoff_A"]) == cutoff
            and r["chain"] == "A"
            and int(r["poses_with_contact"]) >= 2
        }

    shared = sorted(set(by_seed["1917"]) & set(by_seed["2917"]))
    print("  Shared receptor residues recurring in >=2/5 poses:")
    if not shared:
        print("    None")
    for residue in shared:
        print(
            f"    {residue}: "
            f"1917={by_seed['1917'][residue]}/5, "
            f"2917={by_seed['2917'][residue]}/5"
        )

print("\nTARGET-SITE ENGAGEMENT")
for position in ("203", "205"):
    print(f"\nArg{position}")
    for cutoff in (5.0, 8.0, 10.0):
        counts = {}
        for seed in SEEDS:
            distances = [
                float(r["minimum_distance_to_peptide_A"])
                for r in target_rows
                if r["seed"] == seed
                and r["reference_position_to_check"] == position
                and r["matched_output_residue"].startswith("A:ARG:")
                and r["minimum_distance_to_peptide_A"] != ""
            ]
            counts[seed] = sum(d <= cutoff for d in distances)

        print(
            f"  <= {cutoff:g} Å: "
            f"1917={counts['1917']}/5, "
            f"2917={counts['2917']}/5, "
            f"pooled={sum(counts.values())}/10"
        )

print("\nSIMULTANEOUS ARG203 + ARG205 PROXIMITY")
for cutoff in (5.0, 8.0, 10.0):
    counts = {}
    for seed in SEEDS:
        members = sorted({
            r["member"] for r in target_rows
            if r["seed"] == seed
        })
        count = 0
        for member in members:
            distances = {
                r["reference_position_to_check"]:
                    r["minimum_distance_to_peptide_A"]
                for r in target_rows
                if r["seed"] == seed and r["member"] == member
            }
            a = distances.get("203", "")
            b = distances.get("205", "")
            if a and b and float(a) <= cutoff and float(b) <= cutoff:
                count += 1
        counts[seed] = count

    print(
        f"  <= {cutoff:g} Å: "
        f"1917={counts['1917']}/5, "
        f"2917={counts['2917']}/5, "
        f"pooled={sum(counts.values())}/10"
    )
