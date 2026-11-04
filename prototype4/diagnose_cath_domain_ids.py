import csv
from collections import Counter, defaultdict
from pathlib import Path

DATA = Path("prototype4/data/cath_latest_release")
DOMAIN_LIST = DATA / "cath-domain-list.txt"
BOUNDARIES = DATA / "cath-domain-boundaries.txt"

FAMILY_MAP = Path(
    "prototype4/results/data_isolation_20261101/"
    "family_map_recovered.csv"
)

with FAMILY_MAP.open(newline="", encoding="utf-8-sig") as f:
    pdbs = {
        r["pdb_id"].strip().lower()
        for r in csv.DictReader(f)
    }

classification_ids = set()

with DOMAIN_LIST.open(encoding="utf-8", errors="replace") as f:
    for line in f:
        fields = line.split()
        if fields and not line.lstrip().startswith("#"):
            classification_ids.add(fields[0])

classification_by_pdb = defaultdict(list)

for domain_id in classification_ids:
    classification_by_pdb[domain_id[:4].lower()].append(domain_id)

rows = []
format_counts = Counter()
match_counts = Counter()

with BOUNDARIES.open(encoding="utf-8", errors="replace") as f:
    for line_number, line in enumerate(f, 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue

        fields = line.split()
        if len(fields) < 3:
            continue

        key = fields[0]
        pdb = key[:4].lower()

        if pdb not in pdbs:
            continue

        try:
            n_domains = int(fields[1][1:])
            n_fragments = int(fields[2][1:])
        except (ValueError, IndexError):
            continue

        # Test competing domain-ID conventions against the
        # actual downloaded classification IDs.
        candidates = {}

        for i in range(n_domains):
            candidates[f"zero_based_{i}"] = key + f"{i:02d}"
            candidates[f"one_based_{i+1}"] = key + f"{i+1:02d}"

            candidates[f"zero_based_1digit_{i}"] = key + str(i)
            candidates[f"one_based_1digit_{i+1}"] = key + str(i+1)

        candidates["exact_key"] = key

        matches = {
            label: candidate
            for label, candidate in candidates.items()
            if candidate in classification_ids
        }

        for label in matches:
            match_counts[label] += 1

        format_counts[(len(key), fields[1], fields[2])] += 1

        rows.append({
            "line": line_number,
            "pdb": pdb,
            "key": key,
            "key_length": len(key),
            "domain_count": n_domains,
            "fragment_count": n_fragments,
            "matches": matches,
            "classification_examples": (
                classification_by_pdb[pdb][:12]
            ),
            "raw": line.strip(),
        })

print("\n=== BOUNDARY FORMAT ===")
print("Relevant boundary rows:", len(rows))

print("\nBoundary key lengths:")
print(dict(Counter(r["key_length"] for r in rows)))

print("\n=== DOMAIN ID MATCH COUNTS ===")
print(dict(match_counts))

print("\n=== REPRESENTATIVE RAW ROWS ===")
shown = set()

for row in rows:
    if row["pdb"] in shown:
        continue

    shown.add(row["pdb"])
    print(f"\nPDB: {row['pdb']}")
    print("Raw:", row["raw"])
    print("Boundary key:", row["key"])
    print("Key length:", row["key_length"])
    print("Domain count:", row["domain_count"])
    print("Actual CATH IDs:", row["classification_examples"])
    print("Candidate matches:", row["matches"])

    if len(shown) >= 12:
        break

print("\n=== AMBIGUOUS OR UNMATCHED EXAMPLES ===")

for row in rows:
    if len(row["matches"]) != row["domain_count"]:
        print({
            "pdb": row["pdb"],
            "boundary_key": row["key"],
            "domain_count": row["domain_count"],
            "candidate_matches": row["matches"],
            "classification_ids": row["classification_examples"],
        })

print("\nDIAGNOSTIC COMPLETE")
