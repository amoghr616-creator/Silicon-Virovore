#!/usr/bin/env python3
"""Build a provenance-preserving FASTA catalog from the populated P4 family map."""

import csv
import hashlib
import json
from collections import Counter, OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "prototype4/results/skempi_audit_20261029/family_map.csv"
OUT = ROOT / "prototype4/results/data_isolation_20261101"

ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWYX")

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def parse_json_list(value, field, row_number):
    try:
        result = json.loads(value) if value else []
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Invalid JSON in {field}, CSV row {row_number}: {exc}"
        ) from exc

    if not isinstance(result, list):
        raise ValueError(
            f"{field} must contain a JSON list at CSV row {row_number}"
        )
    return result

def write_fasta(path, records):
    with path.open("w", encoding="utf-8") as f:
        for sequence_id, sequence in records:
            f.write(f">{sequence_id}\n")
            for start in range(0, len(sequence), 70):
                f.write(sequence[start:start + 70] + "\n")

def main():
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)

    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError(
            f"{OUT} already contains files. Inspect it before rerunning."
        )

    OUT.mkdir(parents=True, exist_ok=True)

    unique = OrderedDict()
    rows_out = []
    mapping_counts = Counter()
    invalid_sequences = []

    with SOURCE.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        required = {
            "complex_id", "pdb_id", "family_id",
            "protein1_auth_chain_ids_json",
            "protein2_auth_chain_ids_json",
            "protein1_sequences_json",
            "protein2_sequences_json",
            "protein1_chain_group", "protein2_chain_group",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {sorted(missing)}")

    for row_number, row in enumerate(rows, start=2):
        for partner in ("protein1", "protein2"):
            chain_ids = parse_json_list(
                row[f"{partner}_auth_chain_ids_json"],
                f"{partner}_auth_chain_ids_json",
                row_number,
            )
            sequences = parse_json_list(
                row[f"{partner}_sequences_json"],
                f"{partner}_sequences_json",
                row_number,
            )

            if not all(isinstance(x, str) for x in chain_ids):
                raise ValueError(f"Non-string chain ID at CSV row {row_number}")
            if not all(isinstance(x, str) for x in sequences):
                raise ValueError(f"Non-string sequence at CSV row {row_number}")

            # Only a single chain and single sequence are unambiguously paired
            # without relying on undocumented ordering of parallel arrays.
            if len(chain_ids) == 1 and len(sequences) == 1:
                mapping_status = "single_chain_single_sequence"
            elif len(chain_ids) == len(sequences) and len(chain_ids) > 1:
                mapping_status = "parallel_arrays_order_unverified"
            elif not chain_ids:
                mapping_status = "no_auth_chain_ids"
            else:
                mapping_status = "chain_sequence_count_mismatch"

            mapping_counts[mapping_status] += len(sequences)

            for index, raw_sequence in enumerate(sequences, start=1):
                sequence = "".join(raw_sequence.split()).upper()
                invalid = sorted(set(sequence) - ALLOWED_AA)

                record = {
                    "complex_id": row["complex_id"],
                    "pdb_id": row["pdb_id"],
                    "partner": partner,
                    "sequence_index": index,
                    "family_id": row["family_id"],
                    "family_group_type": "IF30_interaction_family_not_CATH_SCOP",
                    "chain_group": row[f"{partner}_chain_group"],
                    "auth_chain_ids_json": json.dumps(chain_ids),
                    "auth_chain_id": (
                        chain_ids[0]
                        if mapping_status == "single_chain_single_sequence"
                        else ""
                    ),
                    "mapping_status": mapping_status,
                    "sequence_length": len(sequence),
                    "sequence_sha256": "",
                    "fasta_id": "",
                    "sequence_scope": "full_partner_chain_sequence_not_interface_extracted",
                    "status": "valid",
                    "invalid_characters": ",".join(invalid),
                }

                if not sequence:
                    record["status"] = "excluded_empty_sequence"
                elif invalid:
                    record["status"] = "excluded_noncanonical_alphabet"
                    invalid_sequences.append({
                        "complex_id": row["complex_id"],
                        "pdb_id": row["pdb_id"],
                        "partner": partner,
                        "invalid_characters": invalid,
                    })
                else:
                    digest = hashlib.sha256(
                        sequence.encode("ascii")
                    ).hexdigest()

                    record["sequence_sha256"] = digest

                    if digest not in unique:
                        unique[digest] = {
                            "fasta_id": f"SKSEQ{len(unique) + 1:06d}",
                            "sequence": sequence,
                            "has_single_chain_mapping": False,
                        }

                    record["fasta_id"] = unique[digest]["fasta_id"]

                    if mapping_status == "single_chain_single_sequence":
                        unique[digest]["has_single_chain_mapping"] = True

                rows_out.append(record)

    if not rows_out:
        raise RuntimeError("No sequence records were found.")

    # All unique, valid partner sequences.
    all_fasta = [
        (item["fasta_id"], item["sequence"])
        for item in unique.values()
    ]

    # Restrict this FASTA to sequences with at least one unambiguous
    # one-chain/one-sequence source record.
    mapped_fasta = [
        (item["fasta_id"], item["sequence"])
        for item in unique.values()
        if item["has_single_chain_mapping"]
    ]

    unresolved_fasta = [
        (item["fasta_id"], item["sequence"])
        for item in unique.values()
        if not item["has_single_chain_mapping"]
    ]

    write_fasta(OUT / "all_partner_group_sequences.fasta", all_fasta)
    write_fasta(OUT / "single_chain_mapped_sequences.fasta", mapped_fasta)
    write_fasta(OUT / "unresolved_group_sequences.fasta", unresolved_fasta)

    manifest_path = OUT / "sequence_manifest.tsv"
    fields = list(rows_out[0])

    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    summary = {
        "source_csv": str(SOURCE.relative_to(ROOT)),
        "source_csv_sha256": sha256(SOURCE),
        "source_complex_rows": len(rows),
        "source_family_groups": len({
            r["family_id"] for r in rows if r["family_id"]
        }),
        "family_group_interpretation": (
            "Existing IF30 interaction-family labels; not CATH or SCOP assignments"
        ),
        "sequence_entries": len(rows_out),
        "unique_valid_sequences": len(unique),
        "single_chain_mapped_unique_sequences": len(mapped_fasta),
        "unresolved_group_unique_sequences": len(unresolved_fasta),
        "mapping_status_counts": dict(mapping_counts),
        "invalid_sequences": invalid_sequences,
        "sequence_scope": (
            "Full partner-chain/group sequences from family_map.csv. "
            "Interface residues have not been extracted."
        ),
        "mmseqs2_ready_for_full_sequence_prescreen": bool(all_fasta),
        "mmseqs2_interface_screen_ready": False,
        "cath_scop_assignments_present": False,
    }

    (OUT / "sequence_catalog_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("SEQUENCE CATALOG CREATED")
    print("Output:", OUT)
    print("Source complexes:", len(rows))
    print("Existing IF30 groups:", summary["source_family_groups"])
    print("Sequence entries:", len(rows_out))
    print("Unique valid sequences:", len(unique))
    print("Unambiguously single-chain mapped:", len(mapped_fasta))
    print("Unresolved group sequences:", len(unresolved_fasta))
    print("Invalid sequences:", len(invalid_sequences))
    print("\nMapping categories:")
    for key, value in sorted(mapping_counts.items()):
        print(f"  {key}: {value}")

    print("\nIMPORTANT")
    print("These are full partner-chain/group sequences, not interface-only sequences.")
    print("No CATH/SCOP labels have been assigned.")
    print("No existing family-map files were modified.")

if __name__ == "__main__":
    main()
