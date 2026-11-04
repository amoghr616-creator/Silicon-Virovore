import csv
import json
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

ROOT = Path.cwd()

SKEMPI_INTERFACE = Path(
    "prototype4/results/data_isolation_20261101/"
    "interface_extraction/interface_residues.tsv"
)

STRUCTURES = Path("prototype4/data/structures/mmcif")

OUT = Path(
    "prototype4/results/data_isolation_20261101/"
    "interface_window_mmseqs"
)

TARGET_PDB = "9SYA"
TARGET_AUTHOR_CHAINS = ("A", "B")
TARGET_HOTSPOTS = (203, 205)

TARGET_PATCH_CUTOFF_A = 8.0
WINDOW_FLANK = 15

MIN_IDENTITY = 0.30
MAX_EVALUE = 1e-3
MIN_COVERAGE = 0.70

THREADS = 4

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O", "HYP": "P",
    "ASX": "B", "GLX": "Z", "UNK": "X",
}

FORMAT_FIELDS = [
    "query", "target", "pident", "alnlen",
    "qstart", "qend", "tstart", "tend",
    "qcov", "tcov", "evalue", "bits",
]


def as_list(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def clean_sequence(value):
    return "".join(
        char for char in str(value).upper()
        if char in "ABCDEFGHIKLMNPQRSTVWXYZUO"
    )


def cif_column(data, key, n, default="?"):
    values = data.get(key)

    if values is None:
        return [default] * n

    if not isinstance(values, list):
        values = [values]

    if len(values) == 1 and n > 1:
        values = values * n

    if len(values) != n:
        raise ValueError(
            f"mmCIF column mismatch for {key}: "
            f"{len(values)} values, expected {n}"
        )

    return values


def load_structure(path):
    """
    Load canonical entity-polymer sequences and recognized protein
    heavy atoms. Polymer sequence indices are label_seq_id positions.
    """
    data = MMCIF2Dict(str(path))

    entity_ids = as_list(data.get("_entity_poly.entity_id"))
    entity_sequences = as_list(
        data.get("_entity_poly.pdbx_seq_one_letter_code_can")
    )

    if not entity_sequences:
        entity_sequences = as_list(
            data.get("_entity_poly.pdbx_seq_one_letter_code")
        )

    entity_seq = {
        str(entity): clean_sequence(sequence)
        for entity, sequence in zip(entity_ids, entity_sequences)
        if clean_sequence(sequence)
    }

    asym_ids = as_list(data.get("_struct_asym.id"))
    asym_entities = as_list(data.get("_struct_asym.entity_id"))

    entity_by_label = {
        str(label): str(entity)
        for label, entity in zip(asym_ids, asym_entities)
    }

    groups = as_list(data.get("_atom_site.group_PDB"))
    n = len(groups)

    if n == 0:
        raise ValueError(f"No atom records in {path}")

    fields = {
        "auth_chain": cif_column(data, "_atom_site.auth_asym_id", n),
        "label_chain": cif_column(data, "_atom_site.label_asym_id", n),
        "auth_comp": cif_column(data, "_atom_site.auth_comp_id", n),
        "label_comp": cif_column(data, "_atom_site.label_comp_id", n),
        "auth_seq": cif_column(data, "_atom_site.auth_seq_id", n),
        "label_seq": cif_column(data, "_atom_site.label_seq_id", n),
        "ins": cif_column(data, "_atom_site.pdbx_PDB_ins_code", n),
        "element": cif_column(data, "_atom_site.type_symbol", n),
        "atom_name": cif_column(data, "_atom_site.label_atom_id", n),
        "x": cif_column(data, "_atom_site.Cartn_x", n),
        "y": cif_column(data, "_atom_site.Cartn_y", n),
        "z": cif_column(data, "_atom_site.Cartn_z", n),
        "model": cif_column(
            data, "_atom_site.pdbx_PDB_model_num", n, "1"
        ),
        "alt": cif_column(data, "_atom_site.label_alt_id", n, "."),
    }

    atoms = []
    observed_by_label = defaultdict(dict)
    residue_atoms = defaultdict(list)
    residue_info = {}

    for i in range(n):
        if str(groups[i]).upper() not in {"ATOM", "HETATM"}:
            continue

        if str(fields["model"][i]) not in {"1", "1.0"}:
            continue

        if str(fields["alt"][i]) not in {".", "?", "A"}:
            continue

        element = str(fields["element"][i]).upper().strip()

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

        if auth_chain in {".", "?"} or label_chain in {".", "?"}:
            continue

        auth_seq = str(fields["auth_seq"][i])
        label_seq = str(fields["label_seq"][i])

        try:
            label_position = int(label_seq)
        except (TypeError, ValueError):
            continue

        insertion = str(fields["ins"][i])
        if insertion in {".", "?"}:
            insertion = ""

        atom_name = str(fields["atom_name"][i]).strip().upper()

        try:
            xyz = np.array([
                float(fields["x"][i]),
                float(fields["y"][i]),
                float(fields["z"][i]),
            ], dtype=float)
        except (TypeError, ValueError):
            continue

        key = (
            auth_chain, label_chain, auth_seq,
            insertion, label_seq, resname,
        )

        info = {
            "auth_chain": auth_chain,
            "label_chain": label_chain,
            "auth_seq": auth_seq,
            "label_seq": label_seq,
            "insertion": insertion,
            "resname": resname,
            "aa": aa,
        }

        residue_info[key] = info
        residue_atoms[key].append(xyz)

        observed_by_label[label_chain][label_position] = aa

        atoms.append({
            "key": key,
            "info": info,
            "xyz": xyz,
            "atom_name": atom_name,
            "element": element,
        })

    # Prefer complete canonical polymer sequences, including unresolved
    # positions. Fall back to coordinate sequences with X placeholders.
    sequence_by_label = {}

    for label_chain, entity in entity_by_label.items():
        sequence = entity_seq.get(entity, "")

        if sequence:
            sequence_by_label[label_chain] = sequence

    for label_chain, positions in observed_by_label.items():
        if label_chain in sequence_by_label:
            continue

        max_position = max(positions)
        sequence_by_label[label_chain] = "".join(
            positions.get(position, "X")
            for position in range(1, max_position + 1)
        )

    return {
        "sequence_by_label": sequence_by_label,
        "atoms": atoms,
        "residue_atoms": residue_atoms,
        "residue_info": residue_info,
    }


def make_windows(
    sequence,
    interface_positions,
    window_id_prefix,
    metadata,
):
    """
    Create contiguous sequence windows centered on interface positions.
    Overlapping/adjacent windows are merged. Interface positions remain
    explicitly recorded relative to each window.
    """
    sequence = clean_sequence(sequence)

    positions = sorted({
        int(position)
        for position in interface_positions
        if 1 <= int(position) <= len(sequence)
    })

    if not positions:
        return []

    intervals = []

    for position in positions:
        start = max(1, position - WINDOW_FLANK)
        end = min(len(sequence), position + WINDOW_FLANK)

        if intervals and start <= intervals[-1][1] + 1:
            intervals[-1] = (
                intervals[-1][0],
                max(intervals[-1][1], end),
            )
        else:
            intervals.append((start, end))

    records = []

    for start, end in intervals:
        absolute_interface_positions = [
            position for position in positions
            if start <= position <= end
        ]

        relative_positions = [
            position - start + 1
            for position in absolute_interface_positions
        ]

        subsequence = sequence[start - 1:end]

        if not subsequence:
            continue

        record = {
            "window_id": f"{window_id_prefix}{len(records) + 1:05d}",
            **metadata,
            "start_label_seq": start,
            "end_label_seq": end,
            "sequence_length": len(subsequence),
            "interface_residue_count": len(relative_positions),
            "interface_density": round(
                len(relative_positions) / len(subsequence), 5
            ),
            "interface_positions_absolute_json": json.dumps(
                absolute_interface_positions
            ),
            "interface_positions_relative_json": json.dumps(
                relative_positions
            ),
            "sequence": subsequence,
        }

        records.append(record)

    return records



def load_target_pdb(path):
    """
    Load the verified two-chain 9SYA PDB.

    auth_seq retains the PDB author residue number (e.g. Arg203).
    label_seq is the sequential polymer position used by FASTA windows.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("9SYA_target", str(path))
    model = next(structure.get_models())

    sequence_by_label = {}
    atoms = []
    residue_atoms = defaultdict(list)
    residue_info = {}

    for chain in model:
        chain_sequence = []

        polymer_residues = [
            residue for residue in chain
            if is_aa(residue, standard=False)
            and residue.resname.upper() in AA3_TO_1
        ]

        for polymer_position, residue in enumerate(
            polymer_residues, start=1
        ):
            resname = residue.resname.upper()
            aa = AA3_TO_1[resname]

            auth_seq = str(residue.id[1])
            insertion = str(residue.id[2]).strip()
            label_seq = str(polymer_position)

            key = (
                chain.id,
                chain.id,
                auth_seq,
                insertion,
                label_seq,
                resname,
            )

            info = {
                "auth_chain": chain.id,
                "label_chain": chain.id,
                "auth_seq": auth_seq,
                "label_seq": label_seq,
                "insertion": insertion,
                "resname": resname,
                "aa": aa,
            }

            residue_info[key] = info
            chain_sequence.append(aa)

            for atom in residue:
                atom_name = atom.get_name().strip().upper()
                element = str(
                    atom.element or atom_name[:1]
                ).strip().upper()

                if element in {"H", "D"}:
                    continue

                # Keep blank and A alternate conformations only.
                altloc = str(atom.get_altloc()).strip()
                if altloc not in {"", "A"}:
                    continue

                xyz = np.asarray(atom.coord, dtype=float)

                residue_atoms[key].append(xyz)
                atoms.append({
                    "key": key,
                    "info": info,
                    "xyz": xyz,
                    "atom_name": atom_name,
                    "element": element,
                })

        if chain_sequence:
            sequence_by_label[chain.id] = "".join(chain_sequence)

    missing = [
        chain for chain in TARGET_AUTHOR_CHAINS
        if chain not in sequence_by_label
    ]
    if missing:
        raise RuntimeError(
            f"Expected 9SYA chains missing from {path}: {missing}"
        )

    for chain in TARGET_AUTHOR_CHAINS:
        arg_numbers = {
            int(atom["info"]["auth_seq"])
            for atom in atoms
            if atom["info"]["auth_chain"] == chain
            and atom["info"]["resname"] == "ARG"
            and atom["atom_name"] in {"NE", "NH1", "NH2", "CZ"}
        }

        missing_hotspots = set(TARGET_HOTSPOTS) - arg_numbers

        if missing_hotspots:
            raise RuntimeError(
                f"Chain {chain}: missing Arg hotspot residues "
                f"{sorted(missing_hotspots)}"
            )

    return {
        "sequence_by_label": sequence_by_label,
        "atoms": atoms,
        "residue_atoms": residue_atoms,
        "residue_info": residue_info,
    }


def find_target_cif():
    preferred = STRUCTURES / "9sya.cif"

    if preferred.exists():
        return preferred

    candidates = list(
        Path("prototype4/data/structures").rglob("9sya.cif")
    )

    if len(candidates) == 1:
        return candidates[0]

    if not candidates:
        raise FileNotFoundError(
            "Could not locate 9sya.cif under "
            "prototype4/data/structures"
        )

    raise RuntimeError(
        "Multiple 9sya.cif files found. Inspect these paths:\n"
        + "\n".join(str(p) for p in candidates)
    )


def extract_target_patch(target_structure):
    """
    Identify chain residues within 8 A of the guanidinium-group heavy
    atoms of Arg203 and Arg205, separately for author chains A and B.
    """
    atoms = target_structure["atoms"]
    residue_atoms = target_structure["residue_atoms"]
    residue_info = target_structure["residue_info"]

    patch_positions = defaultdict(set)
    patch_rows = []
    chain_status = {}

    guanidinium_names = {"NE", "NH1", "NH2", "CZ"}

    for author_chain in TARGET_AUTHOR_CHAINS:
        hotspot_atoms = []

        for atom in atoms:
            info = atom["info"]

            if info["auth_chain"] != author_chain:
                continue

            if info["resname"] != "ARG":
                continue

            try:
                residue_number = int(info["auth_seq"])
            except ValueError:
                continue

            if residue_number not in TARGET_HOTSPOTS:
                continue

            if atom["atom_name"] in guanidinium_names:
                hotspot_atoms.append(atom["xyz"])

        if not hotspot_atoms:
            chain_status[author_chain] = "missing_hotspot_atoms"
            continue

        tree = cKDTree(np.asarray(hotspot_atoms))
        chain_residue_count = 0

        for key, coords in residue_atoms.items():
            info = residue_info[key]

            if info["auth_chain"] != author_chain:
                continue

            distances, _ = tree.query(
                np.asarray(coords, dtype=float), k=1
            )
            minimum_distance = float(np.min(distances))

            if minimum_distance > TARGET_PATCH_CUTOFF_A:
                continue

            try:
                label_position = int(info["label_seq"])
            except ValueError:
                continue

            patch_positions[info["label_chain"]].add(label_position)
            chain_residue_count += 1

            patch_rows.append({
                "pdb_id": TARGET_PDB,
                "auth_chain": info["auth_chain"],
                "label_chain": info["label_chain"],
                "auth_seq_id": info["auth_seq"],
                "label_seq_id": info["label_seq"],
                "insertion_code": info["insertion"],
                "resname": info["resname"],
                "one_letter": info["aa"],
                "minimum_distance_to_arg203_arg205_guanidinium_A": round(
                    minimum_distance, 3
                ),
            })

        chain_status[author_chain] = {
            "status": "patch_extracted",
            "patch_residue_count": chain_residue_count,
        }

    if not patch_positions.get(
        next(iter(patch_positions), "")
    ) and not patch_rows:
        raise RuntimeError(
            "No target patch residues found. Check 9SYA author-chain "
            "IDs, Arg203/Arg205 numbering, and atom names."
        )

    return patch_positions, patch_rows, chain_status


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


def write_fasta(path, records):
    with path.open("w", encoding="utf-8") as f:
        for record in records:
            f.write(f'>{record["window_id"]}\n')
            f.write(f'{record["sequence"]}\n')


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    if not SKEMPI_INTERFACE.exists():
        raise FileNotFoundError(SKEMPI_INTERFACE)

    target_path = Path("data/receptor/9SYA.pdb")
    if not target_path.exists():
        raise FileNotFoundError(target_path)

    # ------------------------------------------------------------
    # A. Build SKEMPI interface-centered sequence windows.
    # ------------------------------------------------------------
    with SKEMPI_INTERFACE.open(
        newline="", encoding="utf-8"
    ) as f:
        interface_rows = list(csv.DictReader(f, delimiter="\t"))

    groups = defaultdict(set)

    for row in interface_rows:
        try:
            position = int(row["label_seq_id"])
        except (KeyError, TypeError, ValueError):
            continue

        key = (
            row["pdb_id"].upper(),
            row["family_id"],
            row["partner"],
            row["auth_chain"],
            row["label_chain"],
        )
        groups[key].add(position)

    sk_windows = []
    sk_mapping_rows = []
    sk_counter = 0
    structure_cache = {}

    for key, positions in sorted(groups.items()):
        pdb, family, partner, auth_chain, label_chain = key
        cif_path = STRUCTURES / f"{pdb.lower()}.cif"

        if not cif_path.exists():
            raise FileNotFoundError(cif_path)

        if pdb not in structure_cache:
            structure_cache[pdb] = load_structure(cif_path)

        structure = structure_cache[pdb]
        sequence = structure["sequence_by_label"].get(label_chain)

        if not sequence:
            sk_mapping_rows.append({
                "pdb_id": pdb,
                "family_id": family,
                "partner": partner,
                "auth_chain": auth_chain,
                "label_chain": label_chain,
                "status": "missing_polymer_sequence",
                "interface_positions": len(positions),
            })
            continue

        sk_counter += 1

        metadata = {
            "source": "SKEMPI",
            "pdb_id": pdb,
            "family_id": family,
            "partner": partner,
            "auth_chain": auth_chain,
            "label_chain": label_chain,
        }

        windows = make_windows(
            sequence,
            positions,
            f"S{sk_counter:04d}_",
            metadata,
        )

        sk_windows.extend(windows)

        sk_mapping_rows.append({
            "pdb_id": pdb,
            "family_id": family,
            "partner": partner,
            "auth_chain": auth_chain,
            "label_chain": label_chain,
            "status": "windows_created" if windows else "no_valid_positions",
            "interface_positions": len(positions),
            "sequence_length": len(sequence),
            "windows_created": len(windows),
        })

    # ------------------------------------------------------------
    # B. Build 9SYA Arg203/Arg205-centered target windows.
    # ------------------------------------------------------------
    target_structure = load_target_pdb(target_path)

    patch_positions, patch_rows, chain_status = (
        extract_target_patch(target_structure)
    )

    target_windows = []
    target_counter = 0

    for label_chain, positions in sorted(patch_positions.items()):
        sequence = target_structure["sequence_by_label"].get(label_chain)

        if not sequence:
            continue

        # Determine author chain for this label chain from coordinates.
        auth_ids = sorted({
            atom["info"]["auth_chain"]
            for atom in target_structure["atoms"]
            if atom["info"]["label_chain"] == label_chain
        })

        for auth_chain in auth_ids:
            target_counter += 1

            metadata = {
                "source": "HERVK_9SYA",
                "pdb_id": TARGET_PDB,
                "family_id": "TARGET_REFERENCE",
                "partner": "Env_target_patch",
                "auth_chain": auth_chain,
                "label_chain": label_chain,
            }

            target_windows.extend(
                make_windows(
                    sequence,
                    positions,
                    f"T{target_counter:04d}_",
                    metadata,
                )
            )

    if not sk_windows:
        raise RuntimeError("No SKEMPI sequence windows were generated.")

    if not target_windows:
        raise RuntimeError(
            "No HERV-K target windows were generated. "
            "Inspect the Arg-site extraction."
        )

    # ------------------------------------------------------------
    # C. Write auditable inputs and maps.
    # ------------------------------------------------------------
    sk_fasta = OUT / "skempi_interface_windows.fasta"
    target_fasta = OUT / "hervk_9sya_target_windows.fasta"

    write_fasta(sk_fasta, sk_windows)
    write_fasta(target_fasta, target_windows)

    write_tsv(
        OUT / "skempi_interface_window_map.tsv", sk_windows
    )
    write_tsv(
        OUT / "hervk_target_window_map.tsv", target_windows
    )
    write_tsv(
        OUT / "skempi_window_mapping_audit.tsv", sk_mapping_rows
    )
    write_tsv(
        OUT / "hervk_target_patch_residues.tsv", patch_rows
    )

    # ------------------------------------------------------------
    # D. Run MMseqs2 with the prespecified sequence thresholds.
    # ------------------------------------------------------------
    executable = shutil.which("mmseqs")

    if executable is None:
        raise RuntimeError(
            "MMseqs2 executable not found in PATH."
        )

    raw_path = OUT / "mmseqs_raw.tsv"
    tmp_path = OUT / "mmseqs_tmp"

    if raw_path.exists():
        raw_path.unlink()

    if tmp_path.exists():
        shutil.rmtree(tmp_path)

    command = [
        executable,
        "easy-search",
        str(sk_fasta),
        str(target_fasta),
        str(raw_path),
        str(tmp_path),
        "--min-seq-id", str(MIN_IDENTITY),
        "-e", str(MAX_EVALUE),
        "-c", str(MIN_COVERAGE),
        "--cov-mode", "0",
        "--alignment-mode", "3",
        "--max-seqs", "100000",
        "--threads", str(THREADS),
        "--format-output", ",".join(FORMAT_FIELDS),
    ]

    print("\nRUNNING MMSEQS2")
    print(" ".join(command))
    subprocess.run(command, check=True)

    # ------------------------------------------------------------
    # E. Require alignment overlap with interface positions
    #    on BOTH the SKEMPI query and the HERV-K target window.
    # ------------------------------------------------------------
    query_by_id = {
        row["window_id"]: row for row in sk_windows
    }
    target_by_id = {
        row["window_id"]: row for row in target_windows
    }

    all_hits = []
    interface_hits = []
    multi_interface_hits = []

    if raw_path.exists():
        with raw_path.open(encoding="utf-8") as f:
            raw_lines = [line.rstrip("\n") for line in f if line.strip()]
    else:
        raw_lines = []

    for line_number, line in enumerate(raw_lines, start=1):
        parts = line.split("\t")

        if len(parts) != len(FORMAT_FIELDS):
            raise ValueError(
                f"MMseqs row {line_number}: expected "
                f"{len(FORMAT_FIELDS)} fields, got {len(parts)}"
            )

        hit = dict(zip(FORMAT_FIELDS, parts))

        qmeta = query_by_id.get(hit["query"])
        tmeta = target_by_id.get(hit["target"])

        if qmeta is None or tmeta is None:
            raise KeyError(
                f"MMseqs returned unknown window: "
                f"{hit['query']} -> {hit['target']}"
            )

        qstart = int(hit["qstart"])
        qend = int(hit["qend"])
        tstart = int(hit["tstart"])
        tend = int(hit["tend"])

        q_positions = json.loads(
            qmeta["interface_positions_relative_json"]
        )
        t_positions = json.loads(
            tmeta["interface_positions_relative_json"]
        )

        q_overlap = sum(
            qstart <= p <= qend for p in q_positions
        )
        t_overlap = sum(
            tstart <= p <= tend for p in t_positions
        )

        enriched = {
            **hit,
            "query_pdb": qmeta["pdb_id"],
            "query_family": qmeta["family_id"],
            "query_partner": qmeta["partner"],
            "query_auth_chain": qmeta["auth_chain"],
            "query_label_chain": qmeta["label_chain"],
            "target_auth_chain": tmeta["auth_chain"],
            "target_label_chain": tmeta["label_chain"],
            "query_interface_positions_overlapped": q_overlap,
            "target_interface_positions_overlapped": t_overlap,
            "both_sides_interface_overlap": (
                q_overlap >= 1 and t_overlap >= 1
            ),
            "both_sides_multi_residue_overlap": (
                q_overlap >= 2 and t_overlap >= 2
            ),
        }

        all_hits.append(enriched)

        if q_overlap >= 1 and t_overlap >= 1:
            interface_hits.append(enriched)

        if q_overlap >= 2 and t_overlap >= 2:
            multi_interface_hits.append(enriched)

    output_fields = (
        FORMAT_FIELDS
        + [
            "query_pdb", "query_family", "query_partner",
            "query_auth_chain", "query_label_chain",
            "target_auth_chain", "target_label_chain",
            "query_interface_positions_overlapped",
            "target_interface_positions_overlapped",
            "both_sides_interface_overlap",
            "both_sides_multi_residue_overlap",
        ]
    )

    def write_hit_table(path, records):
        with path.open(
            "w", newline="", encoding="utf-8"
        ) as f:
            writer = csv.DictWriter(
                f,
                fieldnames=output_fields,
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(records)

    write_hit_table(OUT / "mmseqs_hits_all.tsv", all_hits)
    write_hit_table(
        OUT / "mmseqs_hits_interface_overlap.tsv",
        interface_hits,
    )
    write_hit_table(
        OUT / "mmseqs_hits_multi_interface_overlap.tsv",
        multi_interface_hits,
    )

    summary = {
        "skempi_partner_chain_groups": len(groups),
        "skempi_windows": len(sk_windows),
        "skempi_windows_shorter_than_15": sum(
            r["sequence_length"] < 15 for r in sk_windows
        ),
        "skempi_windows_with_unresolved_sequence": sum(
            "X" in r["sequence"] for r in sk_windows
        ),
        "target_structure": str(target_path),
        "target_chain_status": chain_status,
        "target_patch_cutoff_A": TARGET_PATCH_CUTOFF_A,
        "target_patch_hotspots": list(TARGET_HOTSPOTS),
        "target_windows": len(target_windows),
        "window_flank_residues": WINDOW_FLANK,
        "mmseqs_thresholds": {
            "minimum_sequence_identity": MIN_IDENTITY,
            "maximum_evalue": MAX_EVALUE,
            "minimum_query_and_target_coverage": MIN_COVERAGE,
            "coverage_mode": 0,
            "alignment_mode": 3,
        },
        "mmseqs_hits_passing_sequence_thresholds": len(all_hits),
        "hits_overlapping_interface_on_both_sides": len(interface_hits),
        "hits_overlapping_at_least_two_interface_residues_per_side": (
            len(multi_interface_hits)
        ),
        "unique_skempi_complexes_with_interface_hits": len({
            r["query_pdb"] for r in interface_hits
        }),
        "unique_if30_groups_with_interface_hits": len({
            r["query_family"] for r in interface_hits
        }),
        "outputs": [
            str(sk_fasta),
            str(target_fasta),
            str(OUT / "skempi_interface_window_map.tsv"),
            str(OUT / "hervk_target_window_map.tsv"),
            str(OUT / "hervk_target_patch_residues.tsv"),
            str(OUT / "mmseqs_hits_all.tsv"),
            str(OUT / "mmseqs_hits_interface_overlap.tsv"),
            str(OUT / "mmseqs_hits_multi_interface_overlap.tsv"),
        ],
    }

    (OUT / "interface_window_mmseqs_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print("\nINTERFACE WINDOW + MMSEQS2 SCREEN COMPLETE")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
