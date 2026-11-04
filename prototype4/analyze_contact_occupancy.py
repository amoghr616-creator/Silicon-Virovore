from __future__ import annotations

import csv
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


ROOT = Path.cwd()
RUN_ROOT = ROOT / "prototype4/results/p4_reliability_repeats_20261103"
OUT = ROOT / "prototype4/results/p4_contact_occupancy_20261103"

RUNS = {
    "1917": RUN_ROOT / "HADDOCK3_chain_A_seed1917",
    "2917": RUN_ROOT / "HADDOCK3_chain_A_seed2917",
}

CLUSTERS = [f"cluster_{i}_model_1.pdb" for i in range(1, 6)]
CUTOFFS = (5.0, 8.0)
TARGET_POSITIONS = (203, 205)


def locate_pose(run_dir: Path, filename: str) -> bytes:
    direct = list(run_dir.rglob(filename))
    if direct:
        return direct[0].read_bytes()

    for archive in sorted(run_dir.rglob("summary.tgz")):
        with tarfile.open(archive, "r:*") as tar:
            matches = [
                member for member in tar.getmembers()
                if member.isfile()
                and Path(member.name).name == filename
            ]
            if matches:
                handle = tar.extractfile(matches[0])
                if handle is not None:
                    return handle.read()

    raise FileNotFoundError(f"{filename} not found under {run_dir}")


def parse_pose(data: bytes):
    atoms = {"A": [], "B": []}
    residues = {"A": defaultdict(list), "B": defaultdict(list)}
    residue_info = {}

    for line in data.decode("utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM") or len(line) < 54:
            continue

        chain = line[21:22].strip()
        if chain not in atoms:
            continue

        altloc = line[16:17].strip()
        if altloc not in ("", "A"):
            continue

        atom_name = line[12:16].strip()
        element = line[76:78].strip().upper()
        if not element:
            element = atom_name.lstrip("0123456789")[:1].upper()
        if element in {"H", "D"}:
            continue

        try:
            xyz = np.array([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ])
        except ValueError:
            continue

        resname = line[17:20].strip()
        resseq = line[22:26].strip()
        icode = line[26:27].strip()
        label = f"{chain}:{resname}:{resseq}{icode}"

        atoms[chain].append((xyz, label))
        residues[chain][label].append(xyz)
        residue_info[label] = {
            "chain": chain,
            "resname": resname,
            "resseq": resseq,
            "icode": icode,
        }

    if not atoms["A"] or not atoms["B"]:
        raise ValueError(
            f"Missing chain A or B atoms: "
            f"A={len(atoms['A'])}, B={len(atoms['B'])}"
        )

    return atoms, residues, residue_info


def analyze_pose(data: bytes):
    atoms, residues, residue_info = parse_pose(data)

    coords_a = np.array([x[0] for x in atoms["A"]])
    coords_b = np.array([x[0] for x in atoms["B"]])
    labels_a = [x[1] for x in atoms["A"]]
    labels_b = [x[1] for x in atoms["B"]]

    tree_a = cKDTree(coords_a)
    tree_b = cKDTree(coords_b)

    contact_sets = {}
    residue_contact_sets = {}

    for cutoff in CUTOFFS:
        neighbors = tree_a.query_ball_tree(tree_b, r=cutoff)
        pairs = set()

        for i, js in enumerate(neighbors):
            for j in js:
                pairs.add((labels_a[i], labels_b[j]))

        contact_sets[cutoff] = pairs
        residue_contact_sets[cutoff] = {
            "A": {a for a, _ in pairs},
            "B": {b for _, b in pairs},
        }

    # Minimum distance from each receptor residue to any peptide atom.
    target_distances = {}
    for label, xyz_list in residues["A"].items():
        xyz = np.asarray(xyz_list)
        distances, _ = tree_b.query(xyz, k=1)
        target_distances[label] = float(np.min(distances))

    return {
        "contact_sets": contact_sets,
        "residue_contact_sets": residue_contact_sets,
        "target_distances": target_distances,
        "residue_info": residue_info,
        "atom_counts": (len(coords_a), len(coords_b)),
    }


def write_csv(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    poses = {}

    print("Loading pose representatives...")
    for seed, run_dir in RUNS.items():
        for member in CLUSTERS:
            key = (seed, member)
            result = analyze_pose(locate_pose(run_dir, member))
            poses[key] = result
            a_count, b_count = result["atom_counts"]
            print(f"  {seed} {member}: chain A={a_count}, B={b_count}")

    pair_rows = []
    residue_rows = []
    target_rows = []

    for seed in RUNS:
        seed_poses = [
            (member, poses[(seed, member)])
            for member in CLUSTERS
        ]

        # Pair and single-residue occupancy, separately at each cutoff.
        for cutoff in CUTOFFS:
            pair_counts = defaultdict(int)
            side_counts = {"A": defaultdict(int), "B": defaultdict(int)}

            for member, result in seed_poses:
                pairs = result["contact_sets"][cutoff]
                for pair in pairs:
                    pair_counts[pair] += 1

                for side in ("A", "B"):
                    for residue in result["residue_contact_sets"][cutoff][side]:
                        side_counts[side][residue] += 1

            for (res_a, res_b), count in sorted(
                pair_counts.items(),
                key=lambda item: (-item[1], item[0]),
            ):
                pair_rows.append({
                    "seed": seed,
                    "cutoff_A": cutoff,
                    "receptor_residue": res_a,
                    "peptide_residue": res_b,
                    "poses_with_contact": count,
                    "representatives": len(CLUSTERS),
                    "occupancy": count / len(CLUSTERS),
                })

            for side in ("A", "B"):
                for residue, count in sorted(
                    side_counts[side].items(),
                    key=lambda item: (-item[1], item[0]),
                ):
                    residue_rows.append({
                        "seed": seed,
                        "cutoff_A": cutoff,
                        "chain": side,
                        "residue": residue,
                        "poses_with_contact": count,
                        "representatives": len(CLUSTERS),
                        "occupancy": count / len(CLUSTERS),
                    })

        # Target-position distance audit. This checks numbering in the
        # output structures; it does not assume numbering is preserved.
        for member, result in seed_poses:
            for position in TARGET_POSITIONS:
                matches = []
                for label, info in result["residue_info"].items():
                    if (
                        info["chain"] == "A"
                        and info["resseq"].strip() == str(position)
                    ):
                        matches.append(label)

                if matches:
                    nearest_label = min(
                        matches,
                        key=lambda label: result["target_distances"][label],
                    )
                    distance = result["target_distances"][nearest_label]
                    identity = nearest_label
                else:
                    distance = ""
                    identity = "NOT_FOUND"

                target_rows.append({
                    "seed": seed,
                    "member": member,
                    "reference_position_to_check": position,
                    "matched_output_residue": identity,
                    "minimum_distance_to_peptide_A": distance,
                })

    write_csv(
        OUT / "residue_pair_occupancy.csv",
        pair_rows,
        [
            "seed", "cutoff_A", "receptor_residue", "peptide_residue",
            "poses_with_contact", "representatives", "occupancy",
        ],
    )
    write_csv(
        OUT / "residue_occupancy.csv",
        residue_rows,
        [
            "seed", "cutoff_A", "chain", "residue",
            "poses_with_contact", "representatives", "occupancy",
        ],
    )
    write_csv(
        OUT / "target_site_distances.csv",
        target_rows,
        [
            "seed", "member", "reference_position_to_check",
            "matched_output_residue", "minimum_distance_to_peptide_A",
        ],
    )

    print("\nTOP RECEPTOR RESIDUES BY OCCUPANCY")
    for seed in RUNS:
        for cutoff in CUTOFFS:
            rows = [
                row for row in residue_rows
                if row["seed"] == seed
                and row["cutoff_A"] == cutoff
                and row["chain"] == "A"
            ]
            rows.sort(key=lambda row: (-row["occupancy"], row["residue"]))
            print(f"\nSeed {seed}, cutoff {cutoff:g} Å")
            for row in rows[:12]:
                print(
                    f"  {row['residue']:20s} "
                    f"{row['poses_with_contact']}/5 "
                    f"({row['occupancy']:.0%})"
                )

    print("\nTOP PEPTIDE RESIDUES BY OCCUPANCY")
    for seed in RUNS:
        for cutoff in CUTOFFS:
            rows = [
                row for row in residue_rows
                if row["seed"] == seed
                and row["cutoff_A"] == cutoff
                and row["chain"] == "B"
            ]
            rows.sort(key=lambda row: (-row["occupancy"], row["residue"]))
            print(f"\nSeed {seed}, cutoff {cutoff:g} Å")
            for row in rows[:12]:
                print(
                    f"  {row['residue']:20s} "
                    f"{row['poses_with_contact']}/5 "
                    f"({row['occupancy']:.0%})"
                )

    print("\nTARGET POSITION DISTANCES")
    for row in target_rows:
        print(
            f"  seed={row['seed']} {row['member']} "
            f"position={row['reference_position_to_check']} "
            f"residue={row['matched_output_residue']} "
            f"distance={row['minimum_distance_to_peptide_A']}"
        )

    print("\nOUTPUT FILES")
    for name in (
        "residue_pair_occupancy.csv",
        "residue_occupancy.csv",
        "target_site_distances.csv",
    ):
        print(" ", OUT / name)


if __name__ == "__main__":
    main()
