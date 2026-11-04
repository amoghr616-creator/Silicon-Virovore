import csv
import json
from pathlib import Path
from collections import Counter

ROOTS = [
    Path("prototype4/results"),
    Path("results"),
]

MATRIX = Path(
    "prototype4/results/live_target_coverage_audit_20261102/"
    "candidate_sequence_coverage.csv"
)

STRUCTURE_EXTENSIONS = {
    ".pdb", ".ent", ".cif", ".mmcif",
}

TEXT_EXTENSIONS = {
    ".json", ".csv", ".tsv", ".txt", ".yaml", ".yml",
}

EXCLUDED = {
    ".git", ".venv", "__pycache__", "node_modules",
}


def load_candidates():
    with MATRIX.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    candidates = load_candidates()
    sequences = {
        row["sequence"]: row["candidate_id"]
        for row in candidates
    }

    structures = []
    metadata_hits = []

    for root in ROOTS:
        if not root.exists():
            continue

        for path in root.rglob("*"):
            if not path.is_file():
                continue

            if any(part in EXCLUDED for part in path.parts):
                continue

            suffix = path.suffix.lower()

            if suffix in STRUCTURE_EXTENSIONS:
                structures.append(path)

            if suffix not in TEXT_EXTENSIONS:
                continue

            try:
                if path.stat().st_size > 20_000_000:
                    continue
                content = path.read_text(errors="ignore")
            except OSError:
                continue

            for sequence, candidate_id in sequences.items():
                if sequence in content:
                    metadata_hits.append(
                        (candidate_id, sequence, path)
                    )

    print("========== STRUCTURAL FILE INVENTORY ==========")
    print(f"Structure files found: {len(structures)}")

    by_extension = Counter(
        path.suffix.lower() for path in structures
    )

    for extension, count in sorted(by_extension.items()):
        print(f"{extension}: {count}")

    for path in sorted(structures):
        print(f"{path} ({path.stat().st_size:,} bytes)")

    print("\n========== CANDIDATE SEQUENCE METADATA HITS ==========")

    for candidate_id, sequence, path in metadata_hits:
        print(f"{candidate_id} | {sequence} | {path}")

    print("\n========== MANIFEST FILES ==========")

    manifest_count = 0

    for root in ROOTS:
        if not root.exists():
            continue

        for path in root.rglob("*manifest*.json"):
            if any(part in EXCLUDED for part in path.parts):
                continue
            manifest_count += 1
            print(path)

    print(f"\nTotal manifest files: {manifest_count}")
    print("Inventory completed. No files were modified.")


if __name__ == "__main__":
    main()
