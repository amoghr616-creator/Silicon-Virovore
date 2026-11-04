from __future__ import annotations

import csv
import io
import math
import re
import tarfile
from pathlib import Path
from itertools import product

import numpy as np
from scipy.spatial import cKDTree


ROOT = Path.cwd()
RUN_ROOT = ROOT / "prototype4/results/p4_reliability_repeats_20261103"
OUT = ROOT / "prototype4/results/p4_pose_reproducibility_20261103"

RUNS = {
    "1917": RUN_ROOT / "HADDOCK3_chain_A_seed1917",
    "2917": RUN_ROOT / "HADDOCK3_chain_A_seed2917",
}

CLUSTERS = [f"cluster_{i}_model_1.pdb" for i in range(1, 6)]
CUTOFFS = [4.0, 5.0, 8.0, 10.0]


def locate_pose_bytes(run_dir: Path, filename: str) -> bytes:
    """Find an extracted PDB or read it directly from summary.tgz."""

    direct = list(run_dir.rglob(filename))
    if direct:
        return direct[0].read_bytes()

    archives = sorted(run_dir.rglob("summary.tgz"))

    for archive in archives:
        try:
            with tarfile.open(archive, "r:*") as tar:
                matches = [
                    member for member in tar.getmembers()
                    if member.isfile()
                    and Path(member.name).name == filename
                ]

                if matches:
                    if len(matches) > 1:
                        print(
                            f"WARNING: {archive} contains "
                            f"{len(matches)} copies of {filename}; "
                            "using the first."
                        )

                    handle = tar.extractfile(matches[0])
                    if handle is not None:
                        return handle.read()
        except (tarfile.TarError, OSError) as exc:
            print(f"WARNING: Could not read {archive}: {exc}")

    raise FileNotFoundError(
        f"Could not find {filename} in {run_dir} "
        "or any summary.tgz archive."
    )


def parse_chain_atoms(pdb_bytes: bytes):
    """
    Parse protein heavy atoms from chains A and B.

    Chain A = receptor
    Chain B = peptide

    Alternate conformations are restricted to blank/A.
    """

    receptor = []
    peptide = []

    text = pdb_bytes.decode("utf-8", errors="replace")

    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue

        if len(line) < 54:
            continue

        chain = line[21:22].strip()
        if chain not in {"A", "B"}:
            continue

        altloc = line[16:17].strip()
        if altloc not in {"", "A"}:
            continue

        atom_name = line[12:16].strip()
        element = line[76:78].strip().upper()

        if not element:
            element = re.sub(r"^[0-9]+", "", atom_name)[:1].upper()

        if element in {"H", "D"}:
            continue

        try:
            xyz = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
        except ValueError:
            continue

        resname = line[17:20].strip()
        resseq = line[22:26].strip()
        icode = line[26:27].strip()

        residue = f"{chain}:{resname}:{resseq}{icode}"
        atom_label = f"{residue}:{atom_name}"

        record = {
            "xyz": xyz,
            "residue": residue,
            "atom": atom_label,
        }

        if chain == "A":
            receptor.append(record)
        else:
            peptide.append(record)

    if not receptor or not peptide:
        raise ValueError(
            "Could not parse both receptor chain A and peptide chain B. "
            f"Parsed A={len(receptor)}, B={len(peptide)}."
        )

    return receptor, peptide


def analyze_pose(pdb_bytes: bytes):
    receptor, peptide = parse_chain_atoms(pdb_bytes)

    xyz_a = np.asarray([a["xyz"] for a in receptor], dtype=float)
    xyz_b = np.asarray([a["xyz"] for a in peptide], dtype=float)

    tree_a = cKDTree(xyz_a)
    tree_b = cKDTree(xyz_b)

    nearest_distances, nearest_b_indices = tree_b.query(xyz_a, k=1)
    min_index_a = int(np.argmin(nearest_distances))
    min_index_b = int(nearest_b_indices[min_index_a])
    minimum_distance = float(nearest_distances[min_index_a])

    closest_pair = (
        receptor[min_index_a]["atom"],
        peptide[min_index_b]["atom"],
    )

    contact_sets = {}
    cutoff_metrics = {}

    for cutoff in CUTOFFS:
        neighbors = tree_a.query_ball_tree(tree_b, r=cutoff)

        residue_pairs = set()
        for i, peptide_indices in enumerate(neighbors):
            for j in peptide_indices:
                residue_pairs.add(
                    (
                        receptor[i]["residue"],
                        peptide[j]["residue"],
                    )
                )

        contact_sets[cutoff] = residue_pairs

        cutoff_metrics[cutoff] = {
            "contact_pairs": len(residue_pairs),
            "receptor_residues": len(
                {pair[0] for pair in residue_pairs}
            ),
            "peptide_residues": len(
                {pair[1] for pair in residue_pairs}
            ),
        }

    return {
        "receptor_atoms": len(receptor),
        "peptide_atoms": len(peptide),
        "minimum_distance": minimum_distance,
        "closest_receptor_atom": closest_pair[0],
        "closest_peptide_atom": closest_pair[1],
        "contact_sets": contact_sets,
        "cutoff_metrics": cutoff_metrics,
    }


def jaccard(set_a: set, set_b: set) -> float:
    union = set_a | set_b
    if not union:
        return float("nan")
    return len(set_a & set_b) / len(union)


def write_csv(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    poses = {}

    print("Loading archived HADDOCK representatives...")

    for seed, run_dir in RUNS.items():
        if not run_dir.is_dir():
            raise FileNotFoundError(f"Missing run directory: {run_dir}")

        for filename in CLUSTERS:
            pdb_bytes = locate_pose_bytes(run_dir, filename)
            key = f"{seed}_{filename}"

            print(f"  {key}: {len(pdb_bytes):,} bytes")
            poses[key] = {
                "seed": seed,
                "member": filename,
                "analysis": analyze_pose(pdb_bytes),
            }

    # Verify atom counts against the previous analysis.
    old_metrics = OUT / "pose_metrics.csv"
    if old_metrics.exists():
        with old_metrics.open(encoding="utf-8", newline="") as handle:
            old_rows = list(csv.DictReader(handle))

        old_lookup = {
            (row["seed"], row["member"]): row
            for row in old_rows
        }

        print("\nAtom-count cross-check:")
        for key, pose in poses.items():
            seed = pose["seed"]
            member = pose["member"]
            parsed = pose["analysis"]
            old = old_lookup.get((seed, member))

            if old is None:
                print(f"  {key}: no previous metrics row")
                continue

            a_old = int(old["receptor_atom_count"])
            b_old = int(old["peptide_atom_count"])

            matches = (
                parsed["receptor_atoms"] == a_old
                and parsed["peptide_atoms"] == b_old
            )

            print(
                f"  {key}: "
                f"A {parsed['receptor_atoms']}/{a_old}, "
                f"B {parsed['peptide_atoms']}/{b_old}, "
                f"match={matches}"
            )

    # Per-pose distance and cutoff metrics.
    pose_rows = []

    for key, pose in poses.items():
        result = pose["analysis"]

        row = {
            "pose": key,
            "seed": pose["seed"],
            "member": pose["member"],
            "receptor_atoms": result["receptor_atoms"],
            "peptide_atoms": result["peptide_atoms"],
            "minimum_distance_A": result["minimum_distance"],
            "closest_receptor_atom": result["closest_receptor_atom"],
            "closest_peptide_atom": result["closest_peptide_atom"],
        }

        for cutoff in CUTOFFS:
            m = result["cutoff_metrics"][cutoff]
            suffix = str(int(cutoff))
            row[f"contacts_{suffix}A"] = m["contact_pairs"]
            row[f"receptor_residues_{suffix}A"] = m["receptor_residues"]
            row[f"peptide_residues_{suffix}A"] = m["peptide_residues"]

        pose_rows.append(row)

    pose_fields = list(pose_rows[0])
    write_csv(
        OUT / "pose_distance_metrics.csv",
        pose_rows,
        pose_fields,
    )

    # Full 10 x 10 Jaccard matrix at 5 Å.
    keys = list(poses)
    matrix_rows = []

    for key_a in keys:
        row = {"pose": key_a}
        set_a = poses[key_a]["analysis"]["contact_sets"][5.0]

        for key_b in keys:
            set_b = poses[key_b]["analysis"]["contact_sets"][5.0]
            value = jaccard(set_a, set_b)
            row[key_b] = "" if math.isnan(value) else value

        matrix_rows.append(row)

    write_csv(
        OUT / "similarity_matrix.csv",
        matrix_rows,
        ["pose"] + keys,
    )

    # Cross-seed comparison across all 25 pairs and all cutoffs.
    seeds = list(RUNS)
    seed_a, seed_b = seeds

    run_a = [
        key for key, pose in poses.items()
        if pose["seed"] == seed_a
    ]
    run_b = [
        key for key, pose in poses.items()
        if pose["seed"] == seed_b
    ]

    cross_rows = []

    for key_a, key_b in product(run_a, run_b):
        a = poses[key_a]["analysis"]
        b = poses[key_b]["analysis"]

        row = {
            "pose_seed_1917": key_a,
            "pose_seed_2917": key_b,
            "minimum_distance_1917_A": a["minimum_distance"],
            "minimum_distance_2917_A": b["minimum_distance"],
        }

        for cutoff in CUTOFFS:
            row[f"jaccard_{int(cutoff)}A"] = jaccard(
                a["contact_sets"][cutoff],
                b["contact_sets"][cutoff],
            )

        cross_rows.append(row)

    write_csv(
        OUT / "cross_run_similarity.csv",
        cross_rows,
        list(cross_rows[0]),
    )

    # Summary statistics.
    print("\nPER-POSE DISTANCE / CONTACT RESULTS")
    print("-" * 100)

    for row in pose_rows:
        print(
            f"{row['pose']:36s} "
            f"minimum={row['minimum_distance_A']:.2f} Å  "
            f"4Å={row['contacts_4A']:>3}  "
            f"5Å={row['contacts_5A']:>3}  "
            f"8Å={row['contacts_8A']:>3}  "
            f"10Å={row['contacts_10A']:>3}"
        )

    print("\nCROSS-SEED JACCARD SUMMARY")
    print("-" * 70)

    for cutoff in CUTOFFS:
        column = f"jaccard_{int(cutoff)}A"
        values = [
            row[column] for row in cross_rows
            if not math.isnan(row[column])
        ]

        undefined = len(cross_rows) - len(values)

        if values:
            print(
                f"{int(cutoff):>2} Å: "
                f"mean={np.mean(values):.5f}, "
                f"median={np.median(values):.5f}, "
                f"defined={len(values)}/25, "
                f"undefined={undefined}"
            )
        else:
            print(f"{int(cutoff):>2} Å: no defined comparisons")

    print("\nRANK-1 COMPARISON")
    rank1_a = f"{seed_a}_cluster_1_model_1.pdb"
    rank1_b = f"{seed_b}_cluster_1_model_1.pdb"

    for cutoff in CUTOFFS:
        value = jaccard(
            poses[rank1_a]["analysis"]["contact_sets"][cutoff],
            poses[rank1_b]["analysis"]["contact_sets"][cutoff],
        )
        print(
            f"{int(cutoff)} Å Jaccard: "
            f"{value:.5f}" if not math.isnan(value)
            else f"{int(cutoff)} Å Jaccard: undefined (empty union)"
        )

    # Best cross-run match at 5 Å, excluding undefined comparisons.
    defined_pairs = [
        row for row in cross_rows
        if not math.isnan(row["jaccard_5A"])
    ]

    print("\nBEST CROSS-RUN MATCH AT 5 Å")
    if defined_pairs:
        best = max(defined_pairs, key=lambda row: row["jaccard_5A"])
        print(best)
    else:
        print("No pair has a nonempty contact-set union.")

    print("\nOUTPUTS")
    for filename in (
        "pose_distance_metrics.csv",
        "similarity_matrix.csv",
        "cross_run_similarity.csv",
    ):
        print(" ", OUT / filename)


if __name__ == "__main__":
    main()
