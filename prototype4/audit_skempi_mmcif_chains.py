#!/usr/bin/env python3

import csv
import json
import hashlib
from collections import OrderedDict
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

ROOT = Path.cwd()
FAMILY_MAP = Path(
    "prototype4/results/skempi_audit_20261029/family_map.csv"
)
STRUCTURES = Path("prototype4/data/structures/mmcif")
OUT = Path("prototype4/results/data_isolation_20261101")

AA3 = {
    "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C",
    "GLN":"Q", "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I",
    "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P",
    "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",
    "MSE":"M", "SEC":"U", "PYL":"O", "HYP":"P",
}

def digest(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()

def json_list(value, field):
    result = json.loads(value) if value else []
    if not isinstance(result, list):
        raise ValueError(f"{field} must contain a JSON list")
    return result

def cif_column(data, key, n, default="?"):
    values = data.get(key)
    if values is None:
        return [default] * n
    if not isinstance(values, list):
        return [values] * n
    if len(values) != n:
        raise ValueError(
            f"mmCIF column length mismatch: {key}: {len(values)} != {n}"
        )
    return values

def extract_chains(cif_path):
    data = MMCIF2Dict(str(cif_path))

    group = data["_atom_site.group_PDB"]
    n = len(group)

    auth_chain = cif_column(data, "_atom_site.auth_asym_id", n)
    label_chain = cif_column(data, "_atom_site.label_asym_id", n)
    auth_seq = cif_column(data, "_atom_site.auth_seq_id", n)
    ins_code = cif_column(data, "_atom_site.pdbx_PDB_ins_code", n)
    auth_comp = cif_column(data, "_atom_site.auth_comp_id", n)
    label_comp = cif_column(data, "_atom_site.label_comp_id", n)
    model = cif_column(data, "_atom_site.pdbx_PDB_model_num", n, "1")
    altloc = cif_column(data, "_atom_site.label_alt_id", n, ".")

    chains = {}

    for i in range(n):
        if group[i] != "ATOM":
            continue

        if str(model[i]) not in {"1", "1.0"}:
            continue

        if str(altloc[i]) not in {".", "?", "A"}:
            continue

        chain = str(auth_chain[i])
        if chain in {".", "?"}:
            chain = str(label_chain[i])

        if chain in {".", "?"}:
            continue

        comp = str(auth_comp[i])
        if comp in {".", "?"}:
            comp = str(label_comp[i])

        comp = comp.upper()
        if comp not in AA3:
            continue

        residue_number = str(auth_seq[i])
        insertion = str(ins_code[i])
        if insertion in {".", "?"}:
            insertion = ""

        key = (residue_number, insertion)

        chains.setdefault(chain, OrderedDict())

        previous = chains[chain].get(key)
        letter = AA3[comp]

        if previous is not None and previous != letter:
            raise ValueError(
                f"Conflicting residue identity: {cif_path}, "
                f"chain {chain}, residue {key}"
            )

        chains[chain][key] = letter

    return {
        chain: "".join(residues.values())
        for chain, residues in chains.items()
    }

def main():
    with FAMILY_MAP.open(
        newline="", encoding="utf-8-sig"
    ) as f:
        complexes = list(csv.DictReader(f))

    records = []
    verified_fasta = OrderedDict()
    group_results = []

    for row in complexes:
        pdb_id = row["pdb_id"].strip().upper()
        cif_path = STRUCTURES / f"{pdb_id.lower()}.cif"

        if not cif_path.is_file():
            raise FileNotFoundError(cif_path)

        observed = extract_chains(cif_path)

        for partner in ("protein1", "protein2"):
            chain_ids = json_list(
                row[f"{partner}_auth_chain_ids_json"],
                f"{partner}_auth_chain_ids_json",
            )
            stored_sequences = [
                "".join(str(s).split()).upper()
                for s in json_list(
                    row[f"{partner}_sequences_json"],
                    f"{partner}_sequences_json",
                )
            ]

            chain_ids = [str(x) for x in chain_ids]

            missing = [
                c for c in chain_ids if c not in observed
            ]

            present = [
                c for c in chain_ids if c in observed
            ]

            sequences = [
                observed[c] for c in present
            ]

            status = "unresolved"

            if missing:
                status = "missing_author_chain"

            elif len(chain_ids) == len(stored_sequences) and chain_ids:
                exact = all(
                    observed[c] == seq
                    for c, seq in zip(chain_ids, stored_sequences)
                )
                if exact:
                    status = "verified_per_chain_exact"

            elif len(chain_ids) == 1 and len(stored_sequences) == 1:
                status = (
                    "verified_per_chain_exact"
                    if observed[chain_ids[0]] == stored_sequences[0]
                    else "sequence_mismatch"
                )

            elif (
                len(chain_ids) > 1
                and len(stored_sequences) == 1
                and sequences
            ):
                if all(s == stored_sequences[0] for s in sequences):
                    status = "verified_shared_identical_chain_sequence"
                elif "".join(sequences) == stored_sequences[0]:
                    status = "verified_group_concatenation"

            group_results.append({
                "complex_id": row["complex_id"],
                "pdb_id": pdb_id,
                "partner": partner,
                "family_id": row["family_id"],
                "expected_chain_ids": json.dumps(chain_ids),
                "observed_chain_ids": json.dumps(sorted(observed)),
                "missing_chain_ids": json.dumps(missing),
                "stored_sequence_count": len(stored_sequences),
                "observed_sequence_count": len(sequences),
                "mapping_status": status,
            })

            # Export only chains whose correspondence has been verified.
            verified = status.startswith("verified_")

            for chain in chain_ids:
                sequence = observed.get(chain, "")
                chain_verified = verified and bool(sequence)

                fasta_id = ""
                if chain_verified:
                    fasta_id = (
                        f"SK_{pdb_id}_{partner}_{chain}_"
                        f"{digest(sequence)[:10]}"
                    )
                    verified_fasta[fasta_id] = sequence

                records.append({
                    "complex_id": row["complex_id"],
                    "pdb_id": pdb_id,
                    "partner": partner,
                    "family_id": row["family_id"],
                    "family_group_type": "IF30_not_CATH_SCOP",
                    "author_chain_id": chain,
                    "sequence_length": len(sequence),
                    "sequence_sha256": digest(sequence) if sequence else "",
                    "fasta_id": fasta_id,
                    "chain_mapping_status": status,
                    "chain_mapping_verified": chain_verified,
                    "cif_path": str(cif_path),
                })

    OUT.mkdir(parents=True, exist_ok=True)

    group_fields = list(group_results[0])
    with (OUT / "mmcif_partner_mapping.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(
            f, fieldnames=group_fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(group_results)

    record_fields = list(records[0])
    with (OUT / "mmcif_chain_manifest.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(
            f, fieldnames=record_fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(records)

    with (OUT / "verified_partner_chains.fasta").open(
        "w", encoding="utf-8"
    ) as f:
        for name, sequence in verified_fasta.items():
            f.write(f">{name}\n")
            for i in range(0, len(sequence), 70):
                f.write(sequence[i:i + 70] + "\n")

    counts = {}
    for row in group_results:
        counts[row["mapping_status"]] = (
            counts.get(row["mapping_status"], 0) + 1
        )

    summary = {
        "complex_count": len(complexes),
        "partner_groups_audited": len(group_results),
        "chain_records": len(records),
        "verified_chain_records": sum(
            r["chain_mapping_verified"] for r in records
        ),
        "unique_verified_sequences": len(verified_fasta),
        "mapping_status_counts": counts,
        "cath_scop_assignments_created": False,
        "interface_residues_extracted": False,
        "source_family_groups_remain_IF30": True,
    }

    (OUT / "mmcif_mapping_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    print("SKEMPI MMCIF MAPPING AUDIT COMPLETE")
    print(json.dumps(summary, indent=2))
    print("\nVerified chain FASTA:", len(verified_fasta))
    print("Mapping report:", OUT / "mmcif_partner_mapping.tsv")
    print("Chain manifest:", OUT / "mmcif_chain_manifest.tsv")

if __name__ == "__main__":
    main()
