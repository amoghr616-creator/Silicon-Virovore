#!/usr/bin/env python3
"""Compare final HADDOCK interface fingerprints across independent seeds.

Input:
  prototype4/results/p4_reliability_repeats_20261103/
    HADDOCK3_chain_A_seed1917/analysis/6_seletopclusts_analysis/summary.tgz
    HADDOCK3_chain_A_seed2917/analysis/6_seletopclusts_analysis/summary.tgz

Output:
  prototype4/results/p4_pose_reproducibility_20261103/

This is a structural reproducibility diagnostic, NOT an affinity predictor.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import tarfile
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path


ROOT = Path.cwd()
REPEAT_ROOT = (
    ROOT / "prototype4" / "results"
    / "p4_reliability_repeats_20261103"
)
OUTPUT = (
    ROOT / "prototype4" / "results"
    / "p4_pose_reproducibility_20261103"
)

RUNS = {
    "1917": REPEAT_ROOT / "HADDOCK3_chain_A_seed1917",
    "2917": REPEAT_ROOT / "HADDOCK3_chain_A_seed2917",
}

RECEPTOR_CHAIN = "A"
PEPTIDE_CHAIN = "B"
CONTACT_CUTOFF_A = 5.0
FINAL_STAGE = "6_seletopclusts_analysis"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def parse_pdb(data: bytes):
    """Parse heavy-atom coordinates and retain chain/residue identity."""
    atoms = defaultdict(list)
    chains = set()

    for raw_line in data.decode("ascii", errors="replace").splitlines():
        record = raw_line[0:6].strip()
        if record not in {"ATOM", "HETATM"}:
            continue
        if len(raw_line) < 54:
            continue

        atom_name = raw_line[12:16].strip().upper()
        altloc = raw_line[16:17].strip()
        resname = raw_line[17:20].strip().upper()
        chain = raw_line[21:22].strip()
        resseq = raw_line[22:26].strip()
        icode = raw_line[26:27].strip()

        # Ignore hydrogens/deuteriums and alternate conformers other than A.
        if atom_name.startswith(("H", "D")):
            continue
        if altloc not in {"", "A"}:
            continue
        if not chain or not resseq:
            continue
        if resname in {"HOH", "WAT", "DOD"}:
            continue

        try:
            x = float(raw_line[30:38])
            y = float(raw_line[38:46])
            z = float(raw_line[46:54])
        except ValueError:
            continue

        residue_id = f"{chain}:{resname}{resseq}{icode}"
        atoms[chain].append((x, y, z, residue_id))
        chains.add(chain)

    return atoms, sorted(chains)


def contact_fingerprint(atoms):
    """Return receptor-peptide residue pairs with any heavy-atom contact."""
    receptor_atoms = atoms.get(RECEPTOR_CHAIN, [])
    peptide_atoms = atoms.get(PEPTIDE_CHAIN, [])

    if not receptor_atoms or not peptide_atoms:
        raise ValueError(
            "Expected both receptor and peptide chains; "
            f"found chains={sorted(atoms)}"
        )

    cutoff2 = CONTACT_CUTOFF_A ** 2
    cell_width = CONTACT_CUTOFF_A
    grid = defaultdict(list)

    def cell_key(x, y, z):
        return (
            math.floor(x / cell_width),
            math.floor(y / cell_width),
            math.floor(z / cell_width),
        )

    for x, y, z, residue_id in receptor_atoms:
        grid[cell_key(x, y, z)].append((x, y, z, residue_id))

    contacts = set()

    for px, py, pz, peptide_residue in peptide_atoms:
        cx, cy, cz = cell_key(px, py, pz)

        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for rx, ry, rz, receptor_residue in grid.get(
                        (cx + dx, cy + dy, cz + dz), ()
                    ):
                        distance2 = (
                            (px - rx) ** 2
                            + (py - ry) ** 2
                            + (pz - rz) ** 2
                        )
                        if distance2 <= cutoff2:
                            contacts.add(
                                (receptor_residue, peptide_residue)
                            )

    return contacts, len(receptor_atoms), len(peptide_atoms)


def jaccard(left: set, right: set) -> float:
    union = left | right
    if not union:
        return float("nan")
    return len(left & right) / len(union)


def load_final_poses(seed: str, run_dir: Path):
    archive = (
        run_dir / "analysis" / FINAL_STAGE / "summary.tgz"
    )
    if not archive.is_file():
        raise FileNotFoundError(f"Missing final-stage archive: {archive}")

    poses = []
    with tarfile.open(archive, "r:gz") as tar:
        members = [
            m for m in tar.getmembers()
            if m.isfile()
            and m.name.lower().endswith(".pdb")
            and re.search(r"cluster_\d+_model_1\.pdb$", m.name)
        ]

        if not members:
            raise RuntimeError(
                f"No selected cluster representative PDBs in {archive}"
            )

        for member in sorted(members, key=lambda m: m.name):
            handle = tar.extractfile(member)
            if handle is None:
                continue

            pdb_bytes = handle.read()
            atoms, chains = parse_pdb(pdb_bytes)
            contacts, receptor_atom_count, peptide_atom_count = (
                contact_fingerprint(atoms)
            )

            poses.append({
                "seed": seed,
                "member": member.name,
                "archive": str(archive.relative_to(ROOT)),
                "archive_sha256": sha256_bytes(archive.read_bytes()),
                "pdb_sha256": sha256_bytes(pdb_bytes),
                "chains": chains,
                "receptor_atom_count": receptor_atom_count,
                "peptide_atom_count": peptide_atom_count,
                "contact_count": len(contacts),
                "receptor_contact_residue_count": len({
                    pair[0] for pair in contacts
                }),
                "peptide_contact_residue_count": len({
                    pair[1] for pair in contacts
                }),
                "contacts": contacts,
            })

    if not poses:
        raise RuntimeError(f"No usable poses found for seed {seed}")

    return poses


def mean_or_none(values):
    valid = [v for v in values if math.isfinite(v)]
    return sum(valid) / len(valid) if valid else None


def binary_entropy_bits(p: float) -> float:
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(
        p * math.log2(p) + (1.0 - p) * math.log2(1.0 - p)
    )


def compare_stage_configs():
    """Compare like-for-like generated stage configs, ignoring only iniseed."""
    stages = [
        "0_topoaa",
        "1_rigidbody",
        "2_seletop",
        "3_flexref",
        "4_emref",
        "5_clustfcc",
        "6_seletopclusts",
    ]

    results = []

    for stage in stages:
        paths = {
            seed: run_dir / stage / "params.cfg"
            for seed, run_dir in RUNS.items()
        }
        if not all(path.is_file() for path in paths.values()):
            results.append({
                "stage": stage,
                "status": "missing_stage_config",
                "equivalent_after_seed_normalization": None,
            })
            continue

        texts = {
            seed: path.read_text(encoding="utf-8", errors="replace")
            for seed, path in paths.items()
        }

        def normalize(text):
            return re.sub(
                r"(?im)^\s*iniseed\s*=\s*[^\r\n]*\r?\n?",
                "",
                text,
            ).strip()

        normalized = {seed: normalize(text) for seed, text in texts.items()}
        same = normalized["1917"] == normalized["2917"]

        results.append({
            "stage": stage,
            "status": "compared",
            "equivalent_after_seed_normalization": same,
            "seed1917_sha256": sha256_bytes(
                normalized["1917"].encode("utf-8")
            ),
            "seed2917_sha256": sha256_bytes(
                normalized["2917"].encode("utf-8")
            ),
        })

    return results


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)

    all_poses = []
    by_seed = {}

    for seed, run_dir in RUNS.items():
        if not run_dir.is_dir():
            raise FileNotFoundError(f"Missing run directory: {run_dir}")

        poses = load_final_poses(seed, run_dir)
        by_seed[seed] = poses
        all_poses.extend(poses)

    # Pose-level table.
    pose_csv = OUTPUT / "pose_metrics.csv"
    with pose_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "seed", "member", "chains", "pdb_sha256",
                "receptor_atom_count", "peptide_atom_count",
                "contact_count", "receptor_contact_residue_count",
                "peptide_contact_residue_count",
            ],
        )
        writer.writeheader()
        for pose in all_poses:
            writer.writerow({
                key: (
                    ";".join(pose[key])
                    if key == "chains"
                    else pose[key]
                )
                for key in writer.fieldnames
            })

    # Full pairwise Jaccard similarity matrix.
    labels = [
        f"{p['seed']}:{Path(p['member']).name}"
        for p in all_poses
    ]
    matrix = []
    for left in all_poses:
        row = [jaccard(left["contacts"], right["contacts"])
               for right in all_poses]
        matrix.append(row)

    similarity_csv = OUTPUT / "pose_similarity_matrix.csv"
    with similarity_csv.open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["pose", *labels])
        for label, row in zip(labels, matrix):
            writer.writerow([
                label,
                *[
                    "" if not math.isfinite(value) else f"{value:.6f}"
                    for value in row
                ],
            ])

    # Compare cluster representatives within and between seeds.
    within = {}
    for seed, poses in by_seed.items():
        values = [
            jaccard(a["contacts"], b["contacts"])
            for a, b in combinations(poses, 2)
        ]
        within[seed] = {
            "pair_count": len(values),
            "mean_jaccard": mean_or_none(values),
            "values": values,
        }

    cross_values = [
        jaccard(a["contacts"], b["contacts"])
        for a in by_seed["1917"]
        for b in by_seed["2917"]
    ]

    # Consensus contacts: present in at least half the cluster representatives.
    consensus = {}
    for seed, poses in by_seed.items():
        counts = Counter()
        for pose in poses:
            counts.update(pose["contacts"])
        minimum_count = math.ceil(len(poses) / 2)
        consensus[seed] = {
            pair for pair, count in counts.items()
            if count >= minimum_count
        }

    consensus_jaccard = jaccard(
        consensus["1917"], consensus["2917"]
    )

    # Descriptive contact-occupancy entropy across the ten representatives.
    contact_counts = Counter()
    for pose in all_poses:
        contact_counts.update(pose["contacts"])

    total_poses = len(all_poses)
    entropies = [
        binary_entropy_bits(count / total_poses)
        for count in contact_counts.values()
    ]

    summary = {
        "analysis_type": "descriptive_interface_reproducibility",
        "contact_cutoff_angstrom": CONTACT_CUTOFF_A,
        "receptor_chain": RECEPTOR_CHAIN,
        "peptide_chain": PEPTIDE_CHAIN,
        "seeds": sorted(RUNS),
        "poses_per_seed": {
            seed: len(poses) for seed, poses in by_seed.items()
        },
        "total_pose_representatives": total_poses,
        "within_seed_cluster_similarity": {
            seed: {
                "pair_count": values["pair_count"],
                "mean_jaccard": values["mean_jaccard"],
            }
            for seed, values in within.items()
        },
        "between_seed_pose_similarity": {
            "pair_count": len(cross_values),
            "mean_jaccard": mean_or_none(cross_values),
            "median_jaccard": (
                sorted(cross_values)[len(cross_values) // 2]
                if cross_values else None
            ),
        },
        "between_seed_consensus_similarity": {
            "seed1917_consensus_contact_count": len(consensus["1917"]),
            "seed2917_consensus_contact_count": len(consensus["2917"]),
            "jaccard": consensus_jaccard,
        },
        "contact_occupancy_entropy": {
            "observed_contact_pairs": len(contact_counts),
            "mean_binary_entropy_bits": mean_or_none(entropies),
            "interpretation": (
                "Descriptive variability across selected cluster "
                "representatives; not a calibrated reliability score."
            ),
        },
        "stage_config_comparison": compare_stage_configs(),
        "limitations": [
            "Cluster representatives are alternative poses, not independent experiments.",
            "Two random seeds are insufficient to calibrate a general reliability model.",
            "Contact overlap does not establish binding affinity or biological activity.",
            "The 5-angstrom contact cutoff is an analysis convention.",
        ],
    }

    summary_path = OUTPUT / "pose_reproducibility_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nPose metrics: {pose_csv}")
    print(f"Similarity matrix: {similarity_csv}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
