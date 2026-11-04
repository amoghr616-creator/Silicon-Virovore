import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

ROOT = Path.cwd()
FAMILY_MAP = Path(
    "prototype4/results/data_isolation_20261101/"
    "family_map_recovered.csv"
)
STRUCTURES = Path("prototype4/data/structures/mmcif")
OUT = Path(
    "prototype4/results/data_isolation_20261101/"
    "interface_extraction"
)

CUTOFF_A = 5.0

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


def parse_json_list(value):
    if not value:
        return []
    result = json.loads(value)
    if isinstance(result, list):
        return [str(x) for x in result if x not in (None, "")]
    return [str(result)]


def column(data, key, n, default="?"):
    values = data.get(key)
    if values is None:
        return [default] * n
    if not isinstance(values, list):
        values = [values]
    if len(values) == 1 and n > 1:
        values = values * n
    if len(values) != n:
        raise ValueError(
            f"mmCIF column length mismatch: {key}: "
            f"{len(values)} != {n}"
        )
    return values


def numeric_sort_value(value):
    try:
        return (0, int(value))
    except (TypeError, ValueError):
        return (1, str(value))


def residue_label(residue):
    return (
        f'{residue["auth_chain"]}:'
        f'{residue["auth_seq"]}'
        f'{residue["insertion"]}:'
        f'{residue["resname"]}'
    )


def extract_atoms(cif_path):
    data = MMCIF2Dict(str(cif_path))

    groups = as_list(data.get("_atom_site.group_PDB"))
    n = len(groups)
    if n == 0:
        raise ValueError(f"No atom records in {cif_path}")

    fields = {
        "auth_chain": column(data, "_atom_site.auth_asym_id", n),
        "label_chain": column(data, "_atom_site.label_asym_id", n),
        "auth_comp": column(data, "_atom_site.auth_comp_id", n),
        "label_comp": column(data, "_atom_site.label_comp_id", n),
        "auth_seq": column(data, "_atom_site.auth_seq_id", n),
        "label_seq": column(data, "_atom_site.label_seq_id", n),
        "ins": column(data, "_atom_site.pdbx_PDB_ins_code", n, "?"),
        "element": column(data, "_atom_site.type_symbol", n),
        "x": column(data, "_atom_site.Cartn_x", n),
        "y": column(data, "_atom_site.Cartn_y", n),
        "z": column(data, "_atom_site.Cartn_z", n),
        "model": column(data, "_atom_site.pdbx_PDB_model_num", n, "1"),
        "alt": column(data, "_atom_site.label_alt_id", n, "."),
    }

    atoms = []
    residues = {}

    for i in range(n):
        if str(groups[i]).upper() not in {"ATOM", "HETATM"}:
            continue
        if str(fields["model"][i]) not in {"1", "1.0"}:
            continue
        if str(fields["alt"][i]) not in {".", "?", "A"}:
            continue

        element = str(fields["element"][i]).strip().upper()
        if element in {"H", "D"}:
            continue

        resname = str(fields["label_comp"][i]).upper()
        aa = AA3_TO_1.get(resname)

        if aa is None:
            resname = str(fields["auth_comp"][i]).upper()
            aa = AA3_TO_1.get(resname)

        if aa is None:
            continue

        auth_chain = str(fields["auth_chain"][i])
        label_chain = str(fields["label_chain"][i])

        if auth_chain in {".", "?"}:
            auth_chain = label_chain
        if label_chain in {".", "?"}:
            label_chain = auth_chain

        auth_seq = str(fields["auth_seq"][i])
        label_seq = str(fields["label_seq"][i])
        insertion = str(fields["ins"][i])

        if insertion in {".", "?"}:
            insertion = ""

        if label_seq in {".", "?"}:
            label_seq = auth_seq

        key = (
            auth_chain, label_chain, auth_seq,
            insertion, label_seq, resname
        )

        residue = {
            "auth_chain": auth_chain,
            "label_chain": label_chain,
            "auth_seq": auth_seq,
            "label_seq": label_seq,
            "insertion": insertion,
            "resname": resname,
            "aa": aa,
        }
        residues[key] = residue

        try:
            xyz = (
                float(fields["x"][i]),
                float(fields["y"][i]),
                float(fields["z"][i]),
            )
        except (TypeError, ValueError):
            continue

        atoms.append({
            "xyz": xyz,
            "residue_key": key,
            "auth_chain": auth_chain,
        })

    return atoms, residues


def atom_subset(atoms, author_chains):
    selected = [
        atom for atom in atoms
        if atom["auth_chain"] in set(author_chains)
    ]
    found = {atom["auth_chain"] for atom in selected}
    missing = sorted(set(author_chains) - found)

    if missing:
        raise ValueError(
            f"Expected author chains absent from coordinates: {missing}"
        )

    if not selected:
        raise ValueError(
            f"No protein heavy atoms for chains {author_chains}"
        )

    return selected


def get_fragments(interface_residues, pdb, family, partner):
    """
    Split interface residues into contiguous polymer sequence fragments.
    Do not join separated sequence regions.
    """
    grouped = defaultdict(list)

    for residue in interface_residues:
        grouped[
            (residue["auth_chain"], residue["label_chain"])
        ].append(residue)

    fragments = []

    for (auth_chain, label_chain), residues in sorted(grouped.items()):
        residues.sort(
            key=lambda r: (
                numeric_sort_value(r["label_seq"]),
                r["insertion"],
                r["auth_seq"],
            )
        )

        segments = []
        current = []

        for residue in residues:
            if not current:
                current = [residue]
                continue

            previous = current[-1]
            try:
                previous_num = int(previous["label_seq"])
                current_num = int(residue["label_seq"])
                contiguous = current_num <= previous_num + 1
            except (TypeError, ValueError):
                contiguous = False

            if contiguous:
                current.append(residue)
            else:
                segments.append(current)
                current = [residue]

        if current:
            segments.append(current)

        for index, segment in enumerate(segments, start=1):
            sequence = "".join(r["aa"] for r in segment)

            fragments.append({
                "pdb_id": pdb,
                "family_id": family,
                "partner": partner,
                "auth_chain": auth_chain,
                "label_chain": label_chain,
                "fragment_index": index,
                "start_label_seq": segment[0]["label_seq"],
                "end_label_seq": segment[-1]["label_seq"],
                "residue_count": len(segment),
                "sequence": sequence,
                "residues": segment,
            })

    return fragments


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    with FAMILY_MAP.open(
        newline="", encoding="utf-8-sig"
    ) as f:
        rows = list(csv.DictReader(f))

    if len(rows) != 35:
        raise RuntimeError(
            f"Expected 35 complexes; found {len(rows)}. "
            "Inspect the recovered family map before continuing."
        )

    residue_output = []
    fragment_output = []
    complex_summary = []

    for row in rows:
        pdb = row["pdb_id"].strip().upper()
        family = row.get("family_id", "")
        cif = STRUCTURES / f"{pdb.lower()}.cif"

        if not cif.exists():
            raise FileNotFoundError(cif)

        for partner in ("protein1", "protein2"):
            status = row.get(f"{partner}_recovery_status", "")
            if status != "verified_exact":
                raise RuntimeError(
                    f"{pdb} {partner} has status {status!r}"
                )

        chains1 = parse_json_list(
            row["protein1_resolved_auth_chain_ids_json"]
        )
        chains2 = parse_json_list(
            row["protein2_resolved_auth_chain_ids_json"]
        )

        overlap = set(chains1) & set(chains2)
        if overlap:
            raise RuntimeError(
                f"{pdb}: partner chain assignments overlap: {overlap}"
            )

        atoms, residue_map = extract_atoms(cif)

        atoms1 = atom_subset(atoms, chains1)
        atoms2 = atom_subset(atoms, chains2)

        xyz1 = np.asarray([a["xyz"] for a in atoms1], dtype=float)
        xyz2 = np.asarray([a["xyz"] for a in atoms2], dtype=float)

        tree1 = cKDTree(xyz1)
        tree2 = cKDTree(xyz2)

        neighbor_lists = tree1.query_ball_tree(
            tree2, r=CUTOFF_A
        )

        contact_data = {
            "protein1": {},
            "protein2": {},
        }

        contact_pairs = 0

        for i, neighbors in enumerate(neighbor_lists):
            if not neighbors:
                continue

            atom1 = atoms1[i]
            key1 = atom1["residue_key"]

            for j in neighbors:
                atom2 = atoms2[j]
                key2 = atom2["residue_key"]

                distance = float(
                    np.linalg.norm(xyz1[i] - xyz2[j])
                )
                contact_pairs += 1

                for side, own_key, other_key in (
                    ("protein1", key1, key2),
                    ("protein2", key2, key1),
                ):
                    entry = contact_data[side].setdefault(
                        own_key,
                        {
                            "min_distance": float("inf"),
                            "atom_pair_count": 0,
                            "counterparts": set(),
                        },
                    )
                    entry["min_distance"] = min(
                        entry["min_distance"], distance
                    )
                    entry["atom_pair_count"] += 1
                    entry["counterparts"].add(other_key)

        partner_counts = {}

        for side in ("protein1", "protein2"):
            interface_residues = []

            for key, info in contact_data[side].items():
                residue = dict(residue_map[key])
                residue["min_distance"] = info["min_distance"]
                residue["atom_pair_count"] = info["atom_pair_count"]
                residue["counterpart_keys"] = info["counterparts"]
                interface_residues.append(residue)

            partner_counts[side] = len(interface_residues)

            for residue in interface_residues:
                counterpart_labels = sorted(
                    residue_label(residue_map[k])
                    for k in residue["counterpart_keys"]
                )

                residue_output.append({
                    "pdb_id": pdb,
                    "family_id": family,
                    "partner": side,
                    "auth_chain": residue["auth_chain"],
                    "label_chain": residue["label_chain"],
                    "auth_seq_id": residue["auth_seq"],
                    "label_seq_id": residue["label_seq"],
                    "insertion_code": residue["insertion"],
                    "resname": residue["resname"],
                    "one_letter": residue["aa"],
                    "minimum_interpartner_distance_A": round(
                        residue["min_distance"], 3
                    ),
                    "contacting_atom_pair_count": (
                        residue["atom_pair_count"]
                    ),
                    "contacting_partner_residues_json": json.dumps(
                        counterpart_labels
                    ),
                })

            fragments = get_fragments(
                interface_residues, pdb, family, side
            )

            for fragment in fragments:
                fragment_output.append({
                    "pdb_id": pdb,
                    "family_id": family,
                    "partner": side,
                    "auth_chain": fragment["auth_chain"],
                    "label_chain": fragment["label_chain"],
                    "fragment_index": fragment["fragment_index"],
                    "start_label_seq": fragment["start_label_seq"],
                    "end_label_seq": fragment["end_label_seq"],
                    "residue_count": fragment["residue_count"],
                    "sequence": fragment["sequence"],
                })

        complex_summary.append({
            "pdb_id": pdb,
            "family_id": family,
            "protein1_chains": json.dumps(chains1),
            "protein2_chains": json.dumps(chains2),
            "protein1_interface_residues": partner_counts["protein1"],
            "protein2_interface_residues": partner_counts["protein2"],
            "interpartner_atom_pairs_within_cutoff": contact_pairs,
            "cutoff_A": CUTOFF_A,
        })

        if contact_pairs == 0:
            print(
                f"WARNING: {pdb} has no interpartner contacts "
                f"within {CUTOFF_A} A."
            )

    def write_tsv(path, records):
        if not records:
            raise RuntimeError(f"No records to write: {path}")
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=list(records[0].keys()),
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(records)

    residue_path = OUT / "interface_residues.tsv"
    fragment_path = OUT / "interface_fragments.tsv"
    summary_path = OUT / "interface_complex_summary.tsv"
    fasta_path = OUT / "interface_fragments.fasta"

    write_tsv(residue_path, residue_output)
    write_tsv(fragment_path, fragment_output)
    write_tsv(summary_path, complex_summary)

    with fasta_path.open("w", encoding="utf-8") as f:
        for i, fragment in enumerate(fragment_output, start=1):
            header = (
                f"IF{i:05d}|{fragment['pdb_id']}|"
                f"{fragment['family_id']}|{fragment['partner']}|"
                f"auth{fragment['auth_chain']}|"
                f"label{fragment['label_chain']}|"
                f"frag{fragment['fragment_index']}|"
                f"n{fragment['residue_count']}"
            )
            f.write(f">{header}\n")
            f.write(f"{fragment['sequence']}\n")

    summary = {
        "complexes_processed": len(complex_summary),
        "partner_groups_processed": 2 * len(complex_summary),
        "total_interface_residues": len(residue_output),
        "interface_fragments": len(fragment_output),
        "contact_cutoff_A": CUTOFF_A,
        "atom_selection": (
            "Protein heavy atoms; first model; altloc . ? A; "
            "recognized amino-acid residues"
        ),
        "fragment_policy": (
            "Only contiguous polymer sequence segments are emitted; "
            "discontinuous segments are kept separate."
        ),
        "outputs": [
            str(residue_path),
            str(fragment_path),
            str(summary_path),
            str(fasta_path),
        ],
    }

    (OUT / "interface_extraction_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print("\nINTERFACE EXTRACTION COMPLETE")
    print(json.dumps(summary, indent=2))
    print("\nPer-complex summary:")
    for item in complex_summary:
        print(
            item["pdb_id"],
            "P1 residues =", item["protein1_interface_residues"],
            "P2 residues =", item["protein2_interface_residues"],
            "atom pairs =", item[
                "interpartner_atom_pairs_within_cutoff"
            ],
        )


if __name__ == "__main__":
    main()
