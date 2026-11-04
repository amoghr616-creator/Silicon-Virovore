import csv
import json
import re
from collections import Counter, OrderedDict
from pathlib import Path

from Bio.Align import PairwiseAligner
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

FAMILY_MAP = Path(
    "prototype4/results/skempi_audit_20261029/family_map.csv"
)
STRUCTURES = Path("prototype4/data/structures/mmcif")
OUT = Path("prototype4/results/data_isolation_20261101")

OUT.mkdir(parents=True, exist_ok=True)

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O", "HYP": "P",
    "ASX": "B", "GLX": "Z", "UNK": "X",
}


def as_list(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def clean_sequence(sequence):
    if sequence is None:
        return ""
    sequence = str(sequence).upper()
    sequence = re.sub(r"\([^)]*\)", "X", sequence)
    return re.sub(r"[^A-Z]", "", sequence)


def parse_json_list(value):
    if not value:
        return []
    try:
        result = json.loads(value)
    except (TypeError, json.JSONDecodeError):
        return []
    if isinstance(result, list):
        return [x for x in result if x not in (None, "")]
    if result in (None, ""):
        return []
    return [result]


def mmcif_column(data, key, n, default="?"):
    values = data.get(key)
    if values is None:
        return [default] * n
    if not isinstance(values, list):
        values = [values]
    if len(values) == 1 and n > 1:
        return values * n
    if len(values) != n:
        raise ValueError(
            f"{key}: expected {n} rows, found {len(values)}"
        )
    return values


def extract_mmcif_chains(path):
    """
    Returns protein chain records with:
      author chain ID
      label chain ID
      entity ID
      full entity-polymer sequence
      coordinate-observed sequence
    """
    data = MMCIF2Dict(str(path))

    # Full polymer sequences, including residues not resolved in the model.
    entity_ids = as_list(data.get("_entity_poly.entity_id"))
    entity_sequences = as_list(
        data.get("_entity_poly.pdbx_seq_one_letter_code_can")
    )

    if not entity_sequences:
        entity_sequences = as_list(
            data.get("_entity_poly.pdbx_seq_one_letter_code")
        )

    entity_sequence = {}
    for entity_id, sequence in zip(entity_ids, entity_sequences):
        seq = clean_sequence(sequence)
        if seq:
            entity_sequence[str(entity_id)] = seq

    # Label asymmetry ID -> entity ID.
    asym_ids = as_list(data.get("_struct_asym.id"))
    asym_entities = as_list(data.get("_struct_asym.entity_id"))

    entity_by_label = {
        str(label): str(entity)
        for label, entity in zip(asym_ids, asym_entities)
    }

    groups = as_list(data.get("_atom_site.group_PDB"))
    n = len(groups)

    if n == 0:
        return []

    auth = mmcif_column(data, "_atom_site.auth_asym_id", n)
    label = mmcif_column(data, "_atom_site.label_asym_id", n)
    entity = mmcif_column(data, "_atom_site.label_entity_id", n)
    comp = mmcif_column(data, "_atom_site.label_comp_id", n)
    label_seq = mmcif_column(data, "_atom_site.label_seq_id", n)
    auth_seq = mmcif_column(data, "_atom_site.auth_seq_id", n)
    insertion = mmcif_column(data, "_atom_site.pdbx_PDB_ins_code", n)
    model = mmcif_column(data, "_atom_site.pdbx_PDB_model_num", n, "1")
    alt = mmcif_column(data, "_atom_site.label_alt_id", n, ".")

    # One coordinate sequence per (author chain, label chain).
    chain_data = {}

    for i in range(n):
        if str(groups[i]).upper() != "ATOM":
            continue
        if str(model[i]) not in {"1", "1.0"}:
            continue
        if str(alt[i]) not in {".", "?", "A"}:
            continue

        aa = AA3_TO_1.get(str(comp[i]).upper())
        if aa is None:
            continue

        auth_id = str(auth[i])
        label_id = str(label[i])

        if auth_id in {".", "?"}:
            auth_id = label_id
        if label_id in {".", "?"}:
            label_id = auth_id

        if auth_id in {".", "?"} or label_id in {".", "?"}:
            continue

        entity_id = str(entity[i])
        if entity_id in {".", "?"}:
            entity_id = entity_by_label.get(label_id, "")

        key = (auth_id, label_id)

        if key not in chain_data:
            chain_data[key] = {
                "entity_id": entity_id,
                "residues": OrderedDict(),
            }

        # label_seq_id is preferred because it represents polymer order.
        ls = str(label_seq[i])
        if ls not in {".", "?"}:
            residue_key = ("label", ls)
        else:
            residue_key = (
                "auth",
                str(auth_seq[i]),
                str(insertion[i]),
            )

        chain_data[key]["residues"].setdefault(residue_key, aa)

    records = []

    for (auth_id, label_id), item in chain_data.items():
        observed = "".join(item["residues"].values())
        entity_id = item["entity_id"]

        records.append({
            "auth_id": auth_id,
            "label_id": label_id,
            "entity_id": entity_id,
            "observed_seq": observed,
            "entity_seq": entity_sequence.get(entity_id, ""),
        })

    return records


ALIGNER = PairwiseAligner()
ALIGNER.mode = "local"
ALIGNER.match_score = 2.0
ALIGNER.mismatch_score = -1.0
ALIGNER.open_gap_score = -6.0
ALIGNER.extend_gap_score = -0.5


def compare_sequences(query, reference):
    query = clean_sequence(query)
    reference = clean_sequence(reference)

    if not query or not reference:
        return {
            "identity": 0.0,
            "query_coverage": 0.0,
            "reference_coverage": 0.0,
            "exact": False,
        }

    if query == reference:
        return {
            "identity": 1.0,
            "query_coverage": 1.0,
            "reference_coverage": 1.0,
            "exact": True,
        }

    alignment = ALIGNER.align(query, reference)[0]

    query_positions = set()
    reference_positions = set()
    matches = 0
    aligned_pairs = 0

    for qi, ri in zip(
        alignment.indices[0],
        alignment.indices[1],
    ):
        qi = int(qi)
        ri = int(ri)

        if qi < 0 or ri < 0:
            continue

        query_positions.add(qi)
        reference_positions.add(ri)
        aligned_pairs += 1

        if query[qi] == reference[ri]:
            matches += 1

    if aligned_pairs == 0:
        identity = 0.0
    else:
        identity = matches / aligned_pairs

    return {
        "identity": identity,
        "query_coverage": len(query_positions) / len(query),
        "reference_coverage": len(reference_positions) / len(reference),
        "exact": False,
    }


def score_chain_match(query, chain):
    options = []

    for source, sequence in (
        ("entity_poly", chain["entity_seq"]),
        ("observed_coordinates", chain["observed_seq"]),
    ):
        if not sequence:
            continue

        metrics = compare_sequences(query, sequence)
        metrics["source"] = source
        options.append(metrics)

    if not options:
        return None

    return max(
        options,
        key=lambda x: (
            x["exact"],
            x["identity"],
            x["query_coverage"],
            x["reference_coverage"],
        ),
    )


def is_acceptable(metrics):
    if metrics["exact"]:
        return True

    # Conservative fallback for missing residues / construct differences.
    return (
        metrics["identity"] >= 0.95
        and metrics["query_coverage"] >= 0.90
        and metrics["reference_coverage"] >= 0.70
    )


def resolve_partner(row, partner, chains):
    expected_ids = [
        str(x) for x in parse_json_list(
            row.get(f"{partner}_auth_chain_ids_json")
        )
    ]

    stored_sequences = [
        clean_sequence(x)
        for x in parse_json_list(
            row.get(f"{partner}_sequences_json")
        )
    ]
    stored_sequences = [x for x in stored_sequences if x]

    per_sequence = []
    mapped_auth = set()
    mapped_label = set()
    methods = set()
    any_ambiguous = False
    all_exact = bool(stored_sequences)
    expanded_search_used = False

    for sequence_index, query in enumerate(stored_sequences, start=1):
        candidates = []

        for chain in chains:
            metrics = score_chain_match(query, chain)
            if metrics is None:
                continue

            expected_match = (
                chain["auth_id"] in expected_ids
                or chain["label_id"] in expected_ids
            )

            candidates.append({
                **chain,
                **metrics,
                "expected_chain_match": expected_match,
            })

        expected_candidates = [
            c for c in candidates if c["expected_chain_match"]
        ]

        preferred = [
            c for c in expected_candidates if is_acceptable(c)
        ]

        if preferred:
            valid = preferred
        else:
            valid = [c for c in candidates if is_acceptable(c)]
            if valid and expected_ids:
                expanded_search_used = True

        if not valid:
            per_sequence.append({
                "sequence_index": sequence_index,
                "sequence_length": len(query),
                "status": "unresolved",
                "matches": [],
            })
            all_exact = False
            continue

        def rank(c):
            return (
                int(c["exact"]),
                round(c["identity"], 8),
                round(c["query_coverage"], 8),
                round(c["reference_coverage"], 8),
            )

        best_rank = max(rank(c) for c in valid)
        best = [c for c in valid if rank(c) == best_rank]

        # Multiple exact matches are acceptable when they are all
        # among the partner's explicitly listed chain IDs.
        all_best_are_expected = all(
            c["expected_chain_match"] for c in best
        )

        if len(best) > 1 and not all_best_are_expected:
            any_ambiguous = True

        for c in best:
            mapped_auth.add(c["auth_id"])
            mapped_label.add(c["label_id"])
            methods.add(
                "exact_" + c["source"]
                if c["exact"]
                else "aligned_" + c["source"]
            )

        if not all(c["exact"] for c in best):
            all_exact = False

        per_sequence.append({
            "sequence_index": sequence_index,
            "sequence_length": len(query),
            "status": "matched",
            "matches": [
                {
                    "auth_chain": c["auth_id"],
                    "label_chain": c["label_id"],
                    "entity_id": c["entity_id"],
                    "source": c["source"],
                    "identity": round(c["identity"], 5),
                    "query_coverage": round(c["query_coverage"], 5),
                    "reference_coverage": round(
                        c["reference_coverage"], 5
                    ),
                    "exact": c["exact"],
                    "was_expected_chain": c["expected_chain_match"],
                }
                for c in best
            ],
        })

    unmatched_sequences = any(
        item["status"] != "matched" for item in per_sequence
    )

    expected_covered = set(expected_ids).issubset(
        mapped_auth | mapped_label
    ) if expected_ids else True

    if not stored_sequences or not mapped_auth:
        status = "unresolved"
    elif unmatched_sequences:
        status = "partial_mapping"
    elif any_ambiguous:
        status = "ambiguous"
    elif not expected_covered:
        status = "recovered_chain_id_mismatch"
    elif all_exact:
        status = "verified_exact"
    else:
        status = "verified_high_identity"

    return {
        "partner": partner,
        "expected_chain_ids": expected_ids,
        "stored_sequence_count": len(stored_sequences),
        "mapped_auth_chain_ids": sorted(mapped_auth),
        "mapped_label_chain_ids": sorted(mapped_label),
        "status": status,
        "mapping_methods": sorted(methods),
        "expanded_search_used": expanded_search_used,
        "per_sequence": per_sequence,
    }


def main():
    if not FAMILY_MAP.exists():
        raise FileNotFoundError(FAMILY_MAP)

    with FAMILY_MAP.open(
        newline="", encoding="utf-8-sig"
    ) as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError("family_map.csv has no data rows")

    recovery_rows = []
    recovered_family_rows = []
    status_counts = Counter()
    structure_cache = {}

    for row in rows:
        pdb = str(row["pdb_id"]).strip().upper()
        structure_path = STRUCTURES / f"{pdb.lower()}.cif"

        if pdb not in structure_cache:
            if not structure_path.exists():
                structure_cache[pdb] = []
            else:
                structure_cache[pdb] = extract_mmcif_chains(
                    structure_path
                )

        chains = structure_cache[pdb]
        recovered_row = dict(row)

        for partner in ("protein1", "protein2"):
            result = resolve_partner(row, partner, chains)
            status_counts[result["status"]] += 1

            recovered_row[
                f"{partner}_resolved_auth_chain_ids_json"
            ] = json.dumps(result["mapped_auth_chain_ids"])

            recovered_row[
                f"{partner}_resolved_label_chain_ids_json"
            ] = json.dumps(result["mapped_label_chain_ids"])

            recovered_row[
                f"{partner}_recovery_status"
            ] = result["status"]

            recovery_rows.append({
                "pdb_id": pdb,
                "family_id": row.get("family_id", ""),
                **{
                    key: value
                    for key, value in result.items()
                    if key != "per_sequence"
                },
                "per_sequence_json": json.dumps(
                    result["per_sequence"]
                ),
            })

        recovered_family_rows.append(recovered_row)

    # Detailed partner-level mapping report.
    report_path = OUT / "partner_mapping_recovery.tsv"
    with report_path.open("w", newline="", encoding="utf-8") as handle:
        fields = list(recovery_rows[0].keys())
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(recovery_rows)

    # New family map; original is never overwritten.
    family_path = OUT / "family_map_recovered.csv"
    with family_path.open("w", newline="", encoding="utf-8") as handle:
        fields = list(recovered_family_rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(recovered_family_rows)

    summary = {
        "complex_count": len(rows),
        "partner_groups": len(recovery_rows),
        "mmcif_chain_records": sum(
            len(chains) for chains in structure_cache.values()
        ),
        "status_counts": dict(status_counts),
        "acceptance_rule": {
            "exact_sequence_match": True,
            "fallback_min_identity": 0.95,
            "fallback_min_query_coverage": 0.90,
            "fallback_min_reference_coverage": 0.70,
            "ambiguous_matches_auto_accepted": False,
        },
        "outputs": [
            str(report_path),
            str(family_path),
        ],
    }

    summary_path = OUT / "mapping_recovery_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print("\nSEQUENCE MAPPING RECOVERY COMPLETE")
    print(json.dumps(summary, indent=2))
    print("\nNo original family-map or source data was overwritten.")


if __name__ == "__main__":
    main()
