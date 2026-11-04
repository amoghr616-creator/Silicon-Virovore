import csv
import json
from pathlib import Path
from collections import OrderedDict

MATRIX = Path(
    "prototype4/results/live_target_coverage_audit_20261102/"
    "candidate_sequence_coverage.csv"
)

ROOT = Path("results")

OUT = Path(
    "prototype4/results/pdb_chain_audit_20261102"
)

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N",
    "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}


def read_pdb_chains(path):
    chains = OrderedDict()

    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue

            if len(line) < 27:
                continue

            residue_name = line[17:20].strip().upper()

            if residue_name not in AA3:
                continue

            chain = line[21].strip() or "_"
            residue_number = line[22:26].strip()
            insertion_code = line[26].strip()

            # Ignore alternate atom records after the first residue record.
            key = (residue_number, insertion_code)

            if chain not in chains:
                chains[chain] = OrderedDict()

            if key not in chains[chain]:
                chains[chain][key] = AA3[residue_name]

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

    if not ROOT.is_dir():
        raise FileNotFoundError(ROOT)

    pdb_files = sorted(ROOT.rglob("*.pdb"))

    records = []
    candidate_hits = {key: [] for key in candidate_sequences}

    print("========== PDB CHAIN COMPOSITION AUDIT ==========")
    print(f"PDB files discovered: {len(pdb_files)}")

    for path in pdb_files:
        try:
            chains = read_pdb_chains(path)
        except OSError as exc:
            print(f"Could not read {path}: {exc}")
            continue

        chain_lengths = {
            chain: len(sequence)
            for chain, sequence in chains.items()
        }

        exact_matches = []
        substring_matches = []

        for candidate_id, candidate_sequence in candidate_sequences.items():
            for chain, chain_sequence in chains.items():
                if chain_sequence == candidate_sequence:
                    exact_matches.append({
                        "candidate_id": candidate_id,
                        "chain": chain,
                        "match_type": "exact_chain_sequence",
                    })

                elif candidate_sequence in chain_sequence:
                    substring_matches.append({
                        "candidate_id": candidate_id,
                        "chain": chain,
                        "match_type": "candidate_is_subsequence",
                    })

        record = {
            "path": str(path),
            "chain_count": len(chains),
            "chain_lengths": chain_lengths,
            "chain_sequences": chains,
            "exact_candidate_matches": exact_matches,
            "substring_candidate_matches": substring_matches,
        }

        records.append(record)

        for match in exact_matches + substring_matches:
            candidate_hits[match["candidate_id"]].append({
                "path": str(path),
                **match,
                "chain_count": len(chains),
                "chain_lengths": chain_lengths,
            })

    OUT.mkdir(parents=True, exist_ok=True)

    report = {
        "pdb_files_scanned": len(pdb_files),
        "candidate_count": len(candidate_sequences),
        "candidate_hits": candidate_hits,
        "structures": records,
        "interpretation": (
            "Chain sequence matching only. A candidate sequence match "
            "does not establish HERV-K receptor identity, structural "
            "quality, docking validity, or binding affinity."
        ),
    }

    json_path = OUT / "pdb_chain_inventory.json"
    json_path.write_text(json.dumps(report, indent=2) + "\n")

    csv_path = OUT / "candidate_pdb_chain_matches.csv"

    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "candidate_id",
            "sequence",
            "pdb_path",
            "match_type",
            "peptide_chain",
            "total_chains",
            "chain_lengths",
        ])

        for candidate_id, sequence in candidate_sequences.items():
            hits = candidate_hits[candidate_id]

            if not hits:
                writer.writerow([
                    candidate_id, sequence, "NONE",
                    "NO_CHAIN_MATCH", "", 0, "",
                ])
                continue

            for hit in hits:
                writer.writerow([
                    candidate_id,
                    sequence,
                    hit["path"],
                    hit["match_type"],
                    hit["chain"],
                    hit["chain_count"],
                    json.dumps(hit["chain_lengths"]),
                ])

    print("\n========== CANDIDATE MATCH SUMMARY ==========")

    for candidate_id, sequence in candidate_sequences.items():
        hits = candidate_hits[candidate_id]
        exact = sum(
            hit["match_type"] == "exact_chain_sequence"
            for hit in hits
        )
        partial = sum(
            hit["match_type"] == "candidate_is_subsequence"
            for hit in hits
        )

        print(
            f"{candidate_id}: exact={exact}, "
            f"partial={partial}, sequence={sequence}"
        )

    print(f"\nJSON report: {json_path}")
    print(f"CSV report: {csv_path}")
    print("\nNo structural files were modified.")


if __name__ == "__main__":
    main()
