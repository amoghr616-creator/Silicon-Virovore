import csv
import json
import hashlib
from pathlib import Path
from collections import OrderedDict, Counter

MATRIX = Path(
    "prototype4/results/live_target_coverage_audit_20261102/"
    "candidate_sequence_coverage.csv"
)

ROOTS = [
    Path("results"),
    Path("prototype4/results"),
]

OUT = Path(
    "prototype4/results/pdb_chain_audit_20261102"
)

EXTENSIONS = {".pdb", ".ent"}

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N",
    "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}

EXCLUDED = {
    ".git", ".venv", "__pycache__", "node_modules"
}


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_pdb_chains(path):
    chains = OrderedDict()
    current_model = 1

    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                try:
                    current_model = int(line[10:14].strip())
                except ValueError:
                    pass

                if current_model > 1:
                    break

            if not line.startswith(("ATOM  ", "HETATM")):
                continue

            if len(line) < 27:
                continue

            altloc = line[16]
            if altloc not in (" ", "A"):
                continue

            residue_name = line[17:20].strip().upper()
            if residue_name not in AA3:
                continue

            chain = line[21].strip() or "_"
            residue_number = line[22:26].strip()
            insertion_code = line[26].strip()
            residue_key = (residue_number, insertion_code)

            if chain not in chains:
                chains[chain] = OrderedDict()

            if residue_key not in chains[chain]:
                chains[chain][residue_key] = AA3[residue_name]

    return {
        chain: "".join(residues.values())
        for chain, residues in chains.items()
        if residues
    }


def main():
    if not MATRIX.is_file():
        raise FileNotFoundError(MATRIX)

    with MATRIX.open(newline="") as handle:
        candidates = list(csv.DictReader(handle))

    candidate_sequences = {
        row["candidate_id"]: row["sequence"].strip().upper()
        for row in candidates
    }

    files = set()
    counts_by_root = Counter()

    for root in ROOTS:
        if not root.is_dir():
            print(f"WARNING: directory missing: {root}")
            continue

        for path in root.rglob("*"):
            if not path.is_file():
                continue

            if any(part in EXCLUDED for part in path.parts):
                continue

            if path.suffix.lower() in EXTENSIONS:
                resolved = path.resolve()
                files.add(resolved)

    records = []
    hits_by_candidate = {
        candidate_id: []
        for candidate_id in candidate_sequences
    }

    for path in sorted(files):
        try:
            chains = read_pdb_chains(path)
            digest = sha256_file(path)
        except OSError as exc:
            print(f"READ ERROR: {path}: {exc}")
            continue

        if str(path).startswith(str(Path("prototype4/results").resolve())):
            root_label = "prototype4/results"
        else:
            root_label = "results"

        counts_by_root[root_label] += 1

        exact_matches = []
        substring_matches = []

        for candidate_id, sequence in candidate_sequences.items():
            for chain_id, chain_sequence in chains.items():

                if chain_sequence == sequence:
                    exact_matches.append({
                        "candidate_id": candidate_id,
                        "chain": chain_id,
                        "match_type": "exact_chain_sequence",
                    })

                elif sequence in chain_sequence:
                    substring_matches.append({
                        "candidate_id": candidate_id,
                        "chain": chain_id,
                        "match_type": "candidate_is_subsequence",
                    })

        if not exact_matches and not substring_matches:
            continue

        chain_summary = {
            chain_id: {
                "length": len(sequence),
                "sequence": sequence,
            }
            for chain_id, sequence in chains.items()
        }

        record = {
            "path": str(path),
            "root": root_label,
            "sha256": digest,
            "chain_count": len(chains),
            "chains": chain_summary,
            "exact_matches": exact_matches,
            "substring_matches": substring_matches,
        }

        records.append(record)

        for match in exact_matches + substring_matches:
            candidate_id = match["candidate_id"]
            hits_by_candidate[candidate_id].append({
                "path": str(path),
                "sha256": digest,
                "root": root_label,
                "chain_count": len(chains),
                "chains": chain_summary,
                **match,
            })

    OUT.mkdir(parents=True, exist_ok=True)

    unique_hashes = {
        record["sha256"] for record in records
    }

    result = {
        "roots_scanned": [str(root) for root in ROOTS],
        "pdb_files_scanned": len(files),
        "counts_by_root": dict(counts_by_root),
        "candidate_count": len(candidate_sequences),
        "candidate_hits": hits_by_candidate,
        "matching_structure_records": records,
        "unique_matching_file_hashes": len(unique_hashes),
        "interpretation": (
            "Chain sequence matches are preliminary only. "
            "They do not establish HERV-K receptor identity, "
            "complex provenance, structural quality, or binding."
        ),
    }

    json_path = OUT / "full_repository_chain_inventory.json"
    json_path.write_text(json.dumps(result, indent=2) + "\n")

    print("========== COMPLETE PDB INVENTORY ==========")
    print(f"Total PDB/ENT files: {len(files)}")

    for root, count in counts_by_root.items():
        print(f"{root}: {count}")

    print("\n========== CANDIDATE MATCHES ==========")

    for candidate_id, sequence in candidate_sequences.items():
        hits = hits_by_candidate[candidate_id]

        print(
            f"\n{candidate_id} | sequence={sequence} "
            f"| exact={sum(h['match_type'] == 'exact_chain_sequence' for h in hits)} "
            f"| subsequence={sum(h['match_type'] == 'candidate_is_subsequence' for h in hits)}"
        )

        for hit in hits:
            print(
                f"  {hit['match_type']} | {hit['path']} "
                f"| SHA256={hit['sha256'][:16]}"
            )

            for chain_id, info in hit["chains"].items():
                print(
                    f"    chain {chain_id}: length={info['length']} "
                    f"sequence={info['sequence']}"
                )

    print("\n========== SUMMARY ==========")
    print(f"Matching structure records: {len(records)}")
    print(f"Unique matching file hashes: {len(unique_hashes)}")
    print(f"JSON report: {json_path}")
    print("No PDB files or manifests were modified.")


if __name__ == "__main__":
    main()
