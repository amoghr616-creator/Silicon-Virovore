from pathlib import Path
from collections import defaultdict

BASE = Path("prototype4/data/cath_latest_release")
DOMAIN_LIST = BASE / "cath-domain-list.txt"
BOUNDARIES = BASE / "cath-domain-boundaries.txt"

TARGETS = {
    "1ao7A", "1ao7E",
    "1bj1H", "1danH",
    "1yy9A", "4cvwA",
    "4j2lA", "5m2oB",
}

def normalize_chain(value):
    value = value.strip()
    return value[:4].lower() + value[4:5].upper()

def normalize_domain(value):
    value = value.strip()
    return value[:4].lower() + value[4:]

def read_data(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing CATH file: {path}")

    return [
        line.strip()
        for line in path.read_text(
            encoding="utf-8", errors="replace"
        ).splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]

domains = defaultdict(dict)
boundaries = defaultdict(dict)

# Domain list:
# ID, C, A, T, H, S35, S60, S95, S100, S100-count, length, resolution
for line in read_data(DOMAIN_LIST):
    fields = line.split()
    if len(fields) < 12:
        continue

    domain_id = normalize_domain(fields[0])
    chain_key = normalize_chain(domain_id[:5])

    if chain_key not in TARGETS:
        continue

    domains[chain_key][domain_id] = {
        "cath": ".".join(fields[1:5]),
        "s35": fields[5],
        "length": fields[10],
        "resolution": fields[11],
    }

# Boundary file: domain ID followed by residue-range string.
for line in read_data(BOUNDARIES):
    fields = line.split(None, 1)
    if len(fields) != 2:
        continue

    domain_id = normalize_domain(fields[0])
    chain_key = normalize_chain(domain_id[:5])

    if chain_key in TARGETS:
        boundaries[chain_key][domain_id] = fields[1]

print("\n=== EXACT DOMAIN-ID JOIN ===")

for chain in sorted(TARGETS):
    domain_rows = domains[chain]
    boundary_rows = boundaries[chain]

    print(f"\n{chain}")
    print(
        f"  classified domains: {len(domain_rows)} | "
        f"boundary records: {len(boundary_rows)}"
    )

    all_ids = sorted(set(domain_rows) | set(boundary_rows))

    if not all_ids:
        print("  NO MATCHING CATH DOMAIN OR BOUNDARY RECORD")
        continue

    for domain_id in all_ids:
        d = domain_rows.get(domain_id)
        boundary = boundary_rows.get(domain_id)

        if d:
            print(
                f"  {domain_id} | CATH={d['cath']} "
                f"| S35={d['s35']} | length={d['length']} "
                f"| resolution={d['resolution']} "
                f"| ranges={boundary or 'MISSING'}"
            )
        else:
            print(
                f"  {domain_id} | NOT IN DOMAIN LIST "
                f"| ranges={boundary or 'MISSING'}"
            )

    only_domains = set(domain_rows) - set(boundary_rows)
    only_boundaries = set(boundary_rows) - set(domain_rows)

    if only_domains:
        print("  MISSING BOUNDARIES:", sorted(only_domains))
    if only_boundaries:
        print("  BOUNDARY-ONLY IDS:", sorted(only_boundaries))

print("\n=== COMPLETE ===")
