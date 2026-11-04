from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
import json
import re
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "prototype4" / "results"

RUNS = {
    "A_pilot": RESULTS / "haddock3_chain_A_pilot_20261031",
    "A_repeat": RESULTS / "haddock3_chain_A_rep_seed20261031",
    "B_pilot": RESULTS / "haddock3_chain_B_pilot_20261031",
    "B_repeat": RESULTS / "haddock3_chain_B_rep_seed20261031",
}

CUTOFF_A = 5.0
OUT = RESULTS / "haddock_run_reproducibility_20261101"
OUT.mkdir(parents=True, exist_ok=True)


def open_text(path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", errors="replace")
    return path.open("r", errors="replace")


def sha256_file(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_protein_chains(path):
    """Read heavy-atom coordinates, grouped by chain and residue."""
    chains = defaultdict(lambda: defaultdict(list))

    with open_text(path) as f:
        for line in f:
            if not line.startswith("ATOM  "):
                continue

            altloc = line[16:17]
            if altloc not in (" ", "A"):
                continue

            atom_name = line[12:16].strip()
            element = line[76:78].strip().upper()

            if element in {"H", "D"} or atom_name.upper().startswith(("H", "D")):
                continue

            try:
                xyz = np.array([
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ], dtype=float)
            except ValueError:
                continue

            chain = line[21:22].strip() or "_"
            residue = (
                line[22:26].strip(),
                line[26:27].strip(),
                line[17:20].strip(),
            )

            chains[chain][residue].append(xyz)

    return chains


def residue_label(residue):
    number, insertion, name = residue
    return f"{number}{insertion}:{name}"


def contact_fingerprint(path):
    chains = read_protein_chains(path)

    if len(chains) != 2:
        raise ValueError(
            f"{path.name}: expected exactly 2 protein chains, "
            f"found {[(c, len(r)) for c, r in chains.items()]}"
        )

    # The shorter protein chain is treated as the peptide.
    # The longer chain is treated as the receptor.
    ranked = sorted(chains.items(), key=lambda x: len(x[1]))

    peptide_chain, peptide = ranked[0]
    receptor_chain, receptor = ranked[1]

    if len(peptide) > 60:
        raise ValueError(
            f"{path.name}: shorter chain has {len(peptide)} residues; "
            "automatic peptide-chain assignment is uncertain."
        )

    fingerprint = set()
    receptor_interface = set()
    peptide_interface = set()
    arg203_peptide_contacts = set()
    arg205_peptide_contacts = set()

    min_interface_distance = float("inf")

    for rec_residue, rec_atoms in receptor.items():
        rec_xyz = np.asarray(rec_atoms)

        for pep_residue, pep_atoms in peptide.items():
            pep_xyz = np.asarray(pep_atoms)

            delta = rec_xyz[:, None, :] - pep_xyz[None, :, :]
            min_distance = float(
                np.sqrt(np.sum(delta * delta, axis=2).min())
            )

            min_interface_distance = min(
                min_interface_distance, min_distance
            )

            if min_distance <= CUTOFF_A:
                rec_label = residue_label(rec_residue)
                pep_label = residue_label(pep_residue)

                fingerprint.add(f"{rec_label}|{pep_label}")
                receptor_interface.add(rec_label)
                peptide_interface.add(pep_label)

                rec_number, _, rec_name = rec_residue

                if rec_number == "203" and rec_name == "ARG":
                    arg203_peptide_contacts.add(pep_label)

                if rec_number == "205" and rec_name == "ARG":
                    arg205_peptide_contacts.add(pep_label)

    cluster_match = re.search(r"cluster_(\d+)_model", path.name)
    cluster = int(cluster_match.group(1)) if cluster_match else None

    return {
        "run_model": path.name,
        "cluster": cluster,
        "target_chain": receptor_chain,
        "peptide_chain": peptide_chain,
        "target_residues": len(receptor),
        "peptide_residues": len(peptide),
        "contact_pair_count": len(fingerprint),
        "target_interface_residues": len(receptor_interface),
        "peptide_interface_residues": len(peptide_interface),
        "arg203_peptide_contact_count": len(arg203_peptide_contacts),
        "arg205_peptide_contact_count": len(arg205_peptide_contacts),
        "minimum_interface_distance_A": (
            min_interface_distance
            if min_interface_distance != float("inf")
            else None
        ),
        "fingerprint": fingerprint,
        "sha256": sha256_file(path),
    }


def jaccard(a, b):
    union = a | b
    if not union:
        return None
    return len(a & b) / len(union)


def natural_key(path):
    match = re.search(r"cluster_(\d+)_model", path.name)
    return int(match.group(1)) if match else 999999


def find_models(run_dir):
    folder = run_dir / "6_seletopclusts"

    files = list(folder.glob("cluster_*_model_1.pdb.gz"))
    files += list(folder.glob("cluster_*_model_1.pdb"))

    return sorted(set(files), key=natural_key)


def metadata_seed_hits(run_dir):
    hits = []
    allowed = {".cfg", ".toml", ".yaml", ".yml", ".json", ".log", ".txt"}

    for path in run_dir.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in allowed:
            continue
        if "analysis" in path.parts:
            continue

        try:
            if path.stat().st_size > 2_000_000:
                continue
            lines = path.read_text(errors="replace").splitlines()
        except OSError:
            continue

        for line_number, line in enumerate(lines, 1):
            if re.search(
                r"\b(seed|random_seed|randomseed|rng_seed)\b",
                line,
                re.IGNORECASE,
            ):
                hits.append({
                    "file": str(path.relative_to(ROOT)),
                    "line": line_number,
                    "text": line.strip()[:240],
                })

                if len(hits) >= 80:
                    return hits

    return hits


def summarize(values):
    values = [x for x in values if x is not None]

    if not values:
        return {"n": 0}

    return {
        "n": len(values),
        "mean": statistics.mean(values),
        "median": statistics.median(values),
        "minimum": min(values),
        "maximum": max(values),
    }


def write_csv(path, rows, fields):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    poses_by_run = {}
    pose_rows = []
    metadata = {}

    for run_name, run_dir in RUNS.items():
        print(f"\nLoading {run_name}: {run_dir}")

        if not run_dir.exists():
            raise FileNotFoundError(run_dir)

        files = find_models(run_dir)
        print("Representative models:", len(files))

        if not files:
            raise FileNotFoundError(
                f"No representative PDB models found in {run_dir}"
            )

        poses = []

        for path in files:
            result = contact_fingerprint(path)
            result["run"] = run_name
            result["path"] = str(path.relative_to(ROOT))
            poses.append(result)

            row = {
                key: value
                for key, value in result.items()
                if key != "fingerprint"
            }
            row["fingerprint_json"] = json.dumps(
                sorted(result["fingerprint"])
            )
            pose_rows.append(row)

            print(
                f"  {path.name}: "
                f"{result['contact_pair_count']} contact pairs; "
                f"{result['target_residues']} target residues; "
                f"{result['peptide_residues']} peptide residues"
            )

        poses_by_run[run_name] = poses

        metadata[run_name] = {
            "directory": str(run_dir.relative_to(ROOT)),
            "representative_models": len(poses),
            "seed_metadata_hits": metadata_seed_hits(run_dir),
            "cluster_output_present": (
                run_dir / "5_clustfcc" / "cluster.out.gz"
            ).exists(),
        }

    write_csv(
        OUT / "pose_contact_fingerprints.csv",
        pose_rows,
        [
            "run", "run_model", "cluster", "path",
            "target_chain", "peptide_chain",
            "target_residues", "peptide_residues",
            "contact_pair_count",
            "target_interface_residues",
            "peptide_interface_residues",
            "arg203_peptide_contact_count",
            "arg205_peptide_contact_count",
            "minimum_interface_distance_A",
            "sha256", "fingerprint_json",
        ],
    )

    pair_rows = []
    comparison_groups = {}

    def add_comparison(name, left, right):
        values = []

        for a, b in itertools.product(left, right):
            score = jaccard(a["fingerprint"], b["fingerprint"])
            if score is None:
                continue

            values.append(score)
            pair_rows.append({
                "comparison": name,
                "model_1": f"{a['run']}/{a['run_model']}",
                "model_2": f"{b['run']}/{b['run_model']}",
                "jaccard_similarity": score,
                "intersection_contacts": len(
                    a["fingerprint"] & b["fingerprint"]
                ),
                "union_contacts": len(
                    a["fingerprint"] | b["fingerprint"]
                ),
            })

        comparison_groups[name] = summarize(values)

        if values:
            print(f"\n{name}:")
            print(json.dumps(comparison_groups[name], indent=2))

    # Within-run pose heterogeneity.
    for run_name, poses in poses_by_run.items():
        add_comparison(
            f"{run_name}_within_run",
            poses,
            poses,
        )

        # Remove self-comparisons and duplicate symmetric pairs.
        pair_rows[:] = [
            row for row in pair_rows
            if not (
                row["comparison"] == f"{run_name}_within_run"
                and row["model_1"] >= row["model_2"]
            )
        ]

        within_values = [
            row["jaccard_similarity"]
            for row in pair_rows
            if row["comparison"] == f"{run_name}_within_run"
        ]
        comparison_groups[f"{run_name}_within_run"] = summarize(
            within_values
        )

    # Between-run reproducibility, separately by protomer.
    best_match_rows = []

    for protomer in ("A", "B"):
        pilot = poses_by_run[f"{protomer}_pilot"]
        repeat = poses_by_run[f"{protomer}_repeat"]

        comparison_name = f"{protomer}_pilot_vs_repeat"
        add_comparison(comparison_name, pilot, repeat)

        for source_name, source, target in (
            ("pilot", pilot, repeat),
            ("repeat", repeat, pilot),
        ):
            for pose in source:
                scores = [
                    jaccard(pose["fingerprint"], other["fingerprint"])
                    for other in target
                ]
                scores = [x for x in scores if x is not None]

                if scores:
                    best_match_rows.append({
                        "protomer": protomer,
                        "source_run": source_name,
                        "model": pose["run_model"],
                        "best_cross_run_jaccard": max(scores),
                    })

    write_csv(
        OUT / "pairwise_interface_similarity.csv",
        pair_rows,
        [
            "comparison", "model_1", "model_2",
            "jaccard_similarity",
            "intersection_contacts", "union_contacts",
        ],
    )

    write_csv(
        OUT / "best_cross_run_matches.csv",
        best_match_rows,
        [
            "protomer", "source_run", "model",
            "best_cross_run_jaccard",
        ],
    )

    summary = {
        "experiment": "HADDOCK representative-pose interface reproducibility",
        "contact_definition": (
            "Protein residue pairs with minimum heavy-atom distance "
            f"<= {CUTOFF_A} Angstrom"
        ),
        "chain_assignment": (
            "Shorter protein chain treated as peptide; longer chain "
            "treated as receptor. Confirm from residue counts and sequences."
        ),
        "run_metadata": metadata,
        "comparisons": comparison_groups,
        "best_cross_run_match_summary": {
            protomer: summarize([
                row["best_cross_run_jaccard"]
                for row in best_match_rows
                if row["protomer"] == protomer
            ])
            for protomer in ("A", "B")
        },
        "limitations": [
            "Representative cluster models are not independent observations.",
            "Pairwise similarities are descriptive, not independent statistical replicates.",
            "Only two runs per protomer are available in this comparison.",
            "Seed metadata must be checked before describing runs as independent-seed replicates.",
            "Residue numbering must be consistent across the compared structures.",
            "A 5 Angstrom contact fingerprint measures interface similarity, not binding affinity.",
            "This structural pilot does not validate HERV-K inhibition.",
        ],
    }

    (OUT / "reproducibility_summary.json").write_text(
        json.dumps(summary, indent=2)
    )

    print("\n=== REPRODUCIBILITY SUMMARY ===")
    print(json.dumps(summary, indent=2))
    print("\nOUTPUT:", OUT)


if __name__ == "__main__":
    main()
