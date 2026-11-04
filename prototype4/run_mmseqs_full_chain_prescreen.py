#!/usr/bin/env python3
import csv
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "prototype4/results/data_isolation_20261101"
QUERY = BASE / "all_partner_group_sequences.fasta"
OUT = BASE / "mmseqs_full_chain_prescreen"

# Both prepared receptor files use PDB chain A.
# The suffix identifies the original 9SYA source chain, not the prepared PDB chain.
PDBS = {
    "9SYA_source_chain_A": (
        ROOT / "prototype4/data/haddock3/receptors/9sya_source_chain_A.pdb",
        "A",
        "A",
    ),
    "9SYA_source_chain_B": (
        ROOT / "prototype4/data/haddock3/receptors/9sya_source_chain_B.pdb",
        "A",
        "B",
    ),
}

AA3 = {
    "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C",
    "GLN":"Q", "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I",
    "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P",
    "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",
}

FIELDS = [
    "query", "target", "pident", "alnlen",
    "qstart", "qend", "qlen",
    "tstart", "tend", "tlen",
    "evalue", "bits", "qcov", "tcov",
]

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def extract_chain(path, expected_chain):
    residues = {}
    observed_chains = set()

    for line in path.read_text(errors="strict").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM  ") or len(line) < 54:
            continue
        if line[16:17] not in {" ", "A"}:
            continue

        chain = line[21:22].strip()
        observed_chains.add(chain)

        if chain != expected_chain:
            continue

        resname = line[17:20].strip().upper()
        if resname not in AA3:
            continue

        try:
            number = int(line[22:26])
        except ValueError:
            raise ValueError(f"Invalid residue number in {path}: {line}")

        insertion = line[26:27].strip()
        key = (number, insertion)

        if key in residues and residues[key] != AA3[resname]:
            raise ValueError(f"Conflicting residue identity at {key} in {path}")

        residues[key] = AA3[resname]

    if not residues:
        raise RuntimeError(
            f"No standard ATOM residues found for chain {expected_chain} in {path}. "
            f"Observed chains: {sorted(observed_chains)}"
        )

    ordered = sorted(residues.items(), key=lambda x: (x[0][0], x[0][1]))
    sequence = "".join(aa for _, aa in ordered)
    residue_map = {key: aa for key, aa in ordered}

    for position in (203, 205):
        if residue_map.get((position, "")) != "R":
            raise RuntimeError(
                f"{path}: expected Arg at author position {position}; "
                f"observed {residue_map.get((position, ''))!r}. "
                "Stop rather than silently using incorrect numbering."
            )

    return sequence, len(ordered)

def read_fasta(path):
    records = {}
    current = None

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            if current in records:
                raise ValueError(f"Duplicate FASTA ID: {current}")
            records[current] = ""
        else:
            if current is None:
                raise ValueError("FASTA sequence encountered before header")
            records[current] += line.upper()

    if not records or any(not seq for seq in records.values()):
        raise ValueError("FASTA is empty or contains an empty sequence")
    return records

def write_fasta(path, records):
    with path.open("w") as f:
        for name, seq in records.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + "\n")

def main():
    mmseqs = shutil.which("mmseqs")
    if not mmseqs:
        raise RuntimeError("MMseqs2 executable not found on PATH")
    if not QUERY.is_file():
        raise FileNotFoundError(QUERY)
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError(f"Refusing to overwrite existing results: {OUT}")

    queries = read_fasta(QUERY)
    OUT.mkdir(parents=True, exist_ok=True)

    targets = {}
    target_metadata = []

    for target_id, (pdb_path, actual_chain_id, source_chain_label) in PDBS.items():
        if not pdb_path.is_file():
            raise FileNotFoundError(pdb_path)

        sequence, residue_count = extract_chain(
            pdb_path, actual_chain_id
        )
        targets[target_id] = sequence

        target_metadata.append({
            "target_id": target_id,
            "source_pdb": str(pdb_path.relative_to(ROOT)),
            "source_sha256": sha256(pdb_path),
            "original_9sya_source_chain": source_chain_label,
            "actual_prepared_pdb_chain_id": actual_chain_id,
            "sequence_length": len(sequence),
            "observed_standard_residues": residue_count,
            "arg203_verified": True,
            "arg205_verified": True,
            "sequence_scope": "full_observed_protein_chain",
        })

    target_fasta = OUT / "9sya_full_chains.fasta"
    write_fasta(target_fasta, targets)

    raw_output = OUT / "mmseqs_raw.tsv"
    tmp = OUT / "tmp"
    threads = min(4, os.cpu_count() or 1)

    command = [
        mmseqs, "easy-search",
        str(QUERY), str(target_fasta), str(raw_output), str(tmp),
        "--min-seq-id", "0.30",
        "-e", "1e-3",
        "--alignment-mode", "3",
        "--max-seqs", "10000",
        "--threads", str(threads),
        "--format-output", ",".join(FIELDS),
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )
    if result.returncode:
        raise RuntimeError(
            f"MMseqs2 failed ({result.returncode}):\n{result.stderr}"
        )

    manifest_path = BASE / "sequence_manifest.tsv"
    source_manifest = {}
    if manifest_path.is_file():
        with manifest_path.open(newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                fid = row.get("fasta_id", "")
                if fid:
                    source_manifest.setdefault(fid, set()).add(
                        row.get("mapping_status", "unknown")
                    )

    all_hits = []
    strict_hits = []
    coverage_hits = []

    with raw_output.open(newline="") as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row:
                continue
            if len(row) != len(FIELDS):
                raise ValueError(
                    f"Expected {len(FIELDS)} output columns, received {len(row)}"
                )

            hit = dict(zip(FIELDS, row))
            pident = float(hit["pident"])
            evalue = float(hit["evalue"])
            qcov = float(hit["qcov"])
            tcov = float(hit["tcov"])
            qlen = int(hit["qlen"])
            tlen = int(hit["tlen"])

            shorter_coverage = qcov if qlen <= tlen else tcov
            passes_identity_evalue = pident >= 30.0 and evalue < 1e-3
            passes_coverage = shorter_coverage >= 0.70

            hit.update({
                "shorter_sequence_coverage": shorter_coverage,
                "passes_30pct_identity_and_evalue": passes_identity_evalue,
                "passes_70pct_shorter_coverage": passes_coverage,
                "passes_30pct_70pct_combined": (
                    passes_identity_evalue and passes_coverage
                ),
                "query_mapping_statuses": ";".join(sorted(
                    source_manifest.get(hit["query"], {"not_in_manifest"})
                )),
                "comparison_scope": "full_chain_prescreen_not_interface_specific",
            })

            all_hits.append(hit)
            if passes_identity_evalue:
                strict_hits.append(hit)
            if passes_identity_evalue and passes_coverage:
                coverage_hits.append(hit)

    output_fields = list(all_hits[0]) if all_hits else (
        FIELDS + [
            "shorter_sequence_coverage",
            "passes_30pct_identity_and_evalue",
            "passes_70pct_shorter_coverage",
            "passes_30pct_70pct_combined",
            "query_mapping_statuses",
            "comparison_scope",
        ]
    )

    def write_tsv(path, records):
        with path.open("w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=output_fields, delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(records)

    write_tsv(OUT / "all_full_chain_hits.tsv", all_hits)
    write_tsv(OUT / "identity_evalue_hits.tsv", strict_hits)
    write_tsv(OUT / "identity_evalue_coverage70_hits.tsv", coverage_hits)

    version = subprocess.run(
        [mmseqs, "version"], capture_output=True, text=True, check=True
    ).stdout.strip()

    summary = {
        "mmseqs_version": version,
        "query_fasta": str(QUERY.relative_to(ROOT)),
        "query_fasta_sha256": sha256(QUERY),
        "target_fasta": str(target_fasta.relative_to(ROOT)),
        "target_fasta_sha256": sha256(target_fasta),
        "query_unique_sequences": len(queries),
        "targets": target_metadata,
        "parameters": {
            "minimum_sequence_identity": 0.30,
            "maximum_evalue": 1e-3,
            "postfilter_evalue_rule": "strictly_less_than_1e-3",
            "identity_mode": 3,
            "coverage_filter_during_search": None,
            "reported_coverage": ["qcov", "tcov"],
            "additional_shorter_sequence_coverage_flag": 0.70,
        },
        "results": {
            "total_reported_hits": len(all_hits),
            "identity_evalue_hits": len(strict_hits),
            "identity_evalue_and_70pct_coverage_hits": len(coverage_hits),
        },
        "interpretation": (
            "Preliminary full-chain homology screen only. "
            "It does not satisfy the interface-specific MMseqs2 requirement."
        ),
        "command": command,
    }

    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    print("MMSEQS2 FULL-CHAIN PRESCREEN COMPLETE")
    print("Version:", version)
    print("Unique SKEMPI sequences:", len(queries))
    print("9SYA targets:", len(targets))
    print("Total hits:", len(all_hits))
    print("Identity >=30% and E-value <1e-3:", len(strict_hits))
    print("Also >=70% shorter-sequence coverage:", len(coverage_hits))
    print("Results:", OUT)
    print("\nScope: full-chain prescreen only; not interface-specific.")

if __name__ == "__main__":
    main()
