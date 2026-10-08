from pathlib import Path
from itertools import combinations
import math

import numpy as np
import pandas as pd

ROOT = Path("prototype2/p2_rapidock")

RUNS = {
    "P2-01": ROOT / "P2-01-clean",
    "P2-02": ROOT / "P2-02",
    "P2-03": ROOT / "P2-03",
    "P2-04": ROOT / "P2-04",
    "P2-05": ROOT / "P2-05",
    "P2-06": ROOT / "P2-06",
}

# Your current structural hypothesis.
SITE = np.array([-8.405, -5.823, 1.523], dtype=float)
ANCHOR_RESIDUES = {"203", "205"}

CONTACT_CUTOFF = 4.5


def parse_pdb(path):
    """
    Return:
      chains[chain_id][(resseq, icode, resname)] = list of xyz coordinates
    Only ATOM/HETATM records are read.
    Hydrogen atoms are ignored.
    """
    chains = {}

    with open(path, errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            atom_name = line[12:16].strip()

            # Ignore hydrogens.
            if atom_name.upper().startswith("H"):
                continue

            chain = line[21].strip() or "_"
            resname = line[17:20].strip()
            resseq = line[22:26].strip()
            icode = line[26].strip()

            try:
                xyz = np.array([
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ])
            except ValueError:
                continue

            key = (resseq, icode, resname)

            chains.setdefault(chain, {})
            chains[chain].setdefault(key, [])
            chains[chain][key].append(xyz)

    return chains


def chain_residue_count(chain):
    return len(chain)


def identify_receptor_peptide(chains):
    """
    Assumes receptor is the largest protein-like chain
    and docked peptide is the smallest chain with <=50 residues.
    """
    counts = {
        chain: chain_residue_count(residues)
        for chain, residues in chains.items()
    }

    peptide_candidates = [
        chain for chain, n in counts.items()
        if n <= 50
    ]

    if not peptide_candidates:
        raise RuntimeError(
            f"Could not identify peptide chain from chain sizes: {counts}"
        )

    peptide_chain = min(
        peptide_candidates,
        key=lambda c: counts[c]
    )

    receptor_candidates = [
        chain for chain in chains
        if chain != peptide_chain
    ]

    receptor_chain = max(
        receptor_candidates,
        key=lambda c: counts[c]
    )

    return receptor_chain, peptide_chain, counts


def residue_label(chain, key):
    resseq, icode, resname = key
    suffix = icode if icode else ""
    return f"{chain}:{resseq}{suffix}:{resname}"


def analyze_pose(path):
    chains = parse_pdb(path)

    receptor_chain, peptide_chain, chain_counts = (
        identify_receptor_peptide(chains)
    )

    receptor_residues = chains[receptor_chain]
    peptide_residues = chains[peptide_chain]

    receptor_atoms = []
    receptor_atom_labels = []

    for key, atoms in receptor_residues.items():
        for xyz in atoms:
            receptor_atoms.append(xyz)
            receptor_atom_labels.append(
                residue_label(receptor_chain, key)
            )

    peptide_atoms = []

    for atoms in peptide_residues.values():
        peptide_atoms.extend(atoms)

    receptor_atoms = np.array(receptor_atoms)
    peptide_atoms = np.array(peptide_atoms)

    # Peptide centroid and distance to structural hypothesis.
    centroid = peptide_atoms.mean(axis=0)
    site_distance = float(np.linalg.norm(centroid - SITE))

    # Contact receptor residues.
    contact_residues = set()

    for atom in peptide_atoms:
        distances = np.linalg.norm(
            receptor_atoms - atom,
            axis=1
        )

        hit_indices = np.where(distances <= CONTACT_CUTOFF)[0]

        for idx in hit_indices:
            contact_residues.add(
                receptor_atom_labels[idx]
            )

    return {
        "receptor_chain": receptor_chain,
        "peptide_chain": peptide_chain,
        "receptor_residues": len(receptor_residues),
        "peptide_residues": len(peptide_residues),
        "site_distance": site_distance,
        "contact_residues": contact_residues,
        "n_contact_residues": len(contact_residues),
        "peptide_centroid_x": centroid[0],
        "peptide_centroid_y": centroid[1],
        "peptide_centroid_z": centroid[2],
        "chain_counts": chain_counts,
    }


def jaccard(a, b):
    if not a and not b:
        return 1.0

    union = a | b

    if not union:
        return 1.0

    return len(a & b) / len(union)


all_pose_rows = []
candidate_rows = []
contact_frequency_rows = []

for candidate, run_dir in RUNS.items():

    ranked_path = run_dir / "ranked_poses.csv"

    if not ranked_path.exists():
        print(f"WARNING: missing {ranked_path}")
        continue

    ranked = pd.read_csv(ranked_path)

    pose_results = []

    # Only analyze poses that appear in ranked_poses.csv,
    # i.e. poses that survived scoring.
    for _, row in ranked.iterrows():

        pose_filename = str(row["pose_filename"])
        pose_path = run_dir / "poses" / pose_filename

        if not pose_path.exists():
            # Some installations may place poses directly
            # in the run directory.
            alternate = run_dir / pose_filename

            if alternate.exists():
                pose_path = alternate
            else:
                print(
                    f"WARNING: missing pose file for "
                    f"{candidate}: {pose_filename}"
                )
                continue

        try:
            result = analyze_pose(pose_path)
        except Exception as e:
            print(
                f"WARNING: failed {candidate} {pose_filename}: {e}"
            )
            continue

        contacts = result["contact_residues"]

        # Residue-number test for your two structural anchors.
        anchor_hits = set()

        for label in contacts:
            parts = label.split(":")
            if len(parts) >= 2:
                residue_number = parts[1]
                if residue_number in ANCHOR_RESIDUES:
                    anchor_hits.add(residue_number)

        pose_record = {
            "candidate": candidate,
            "pose_filename": pose_filename,
            "rank": row["rank"],
            "cluster_id": row["cluster_id"],
            "delta_g": row["delta_g"],
            "rank_score": row["rank_score"],
            "hybrid_score": row["hybrid_score"],
            "bsa": row["bsa"],
            "n_contact_residues_reported": row["n_contact_residues"],
            "site_distance": result["site_distance"],
            "n_contact_residues_structural": (
                result["n_contact_residues"]
            ),
            "anchor_203_contact": "203" in anchor_hits,
            "anchor_205_contact": "205" in anchor_hits,
            "either_anchor_contact": bool(anchor_hits),
            "contact_residue_set": contacts,
        }

        all_pose_rows.append(pose_record)
        pose_results.append(pose_record)

    if not pose_results:
        continue

    # Contact-set consistency.
    contact_sets = [
        x["contact_residue_set"]
        for x in pose_results
    ]

    pairwise_jaccard = []

    for a, b in combinations(contact_sets, 2):
        pairwise_jaccard.append(jaccard(a, b))

    # Frequency of each receptor residue appearing at the interface.
    freq = {}

    for contacts in contact_sets:
        for residue in contacts:
            freq[residue] = freq.get(residue, 0) + 1

    n_poses = len(pose_results)

    for residue, count in sorted(
        freq.items(),
        key=lambda x: (-x[1], x[0])
    ):
        contact_frequency_rows.append({
            "candidate": candidate,
            "receptor_residue": residue,
            "poses_with_contact": count,
            "contact_fraction": count / n_poses,
        })

    # Cluster occupancy.
    clusters = pd.Series(
        [x["cluster_id"] for x in pose_results]
    )

    cluster_counts = clusters.value_counts()

    largest_cluster_fraction = (
        cluster_counts.iloc[0] / n_poses
    )

    candidate_rows.append({
        "candidate": candidate,
        "n_analyzed_poses": n_poses,
        "median_site_distance_A": np.median([
            x["site_distance"] for x in pose_results
        ]),
        "fraction_within_15A_of_site": np.mean([
            x["site_distance"] <= 15
            for x in pose_results
        ]),
        "fraction_within_20A_of_site": np.mean([
            x["site_distance"] <= 20
            for x in pose_results
        ]),
        "median_structural_contact_residues": np.median([
            x["n_contact_residues_structural"]
            for x in pose_results
        ]),
        "mean_contact_set_jaccard": (
            np.mean(pairwise_jaccard)
            if pairwise_jaccard
            else np.nan
        ),
        "median_delta_g": np.median([
            x["delta_g"] for x in pose_results
        ]),
        "median_rank_score": np.median([
            x["rank_score"] for x in pose_results
        ]),
        "n_clusters": len(cluster_counts),
        "largest_cluster_fraction": largest_cluster_fraction,
        "anchor_203_contact_fraction": np.mean([
            x["anchor_203_contact"]
            for x in pose_results
        ]),
        "anchor_205_contact_fraction": np.mean([
            x["anchor_205_contact"]
            for x in pose_results
        ]),
        "either_anchor_contact_fraction": np.mean([
            x["either_anchor_contact"]
            for x in pose_results
        ]),
    })


# Save pose-level records without the Python set object.
pose_df = pd.DataFrame(all_pose_rows)

if len(pose_df):
    pose_df["contact_residues"] = pose_df[
        "contact_residue_set"
    ].apply(
        lambda s: ";".join(sorted(s))
    )

    pose_df = pose_df.drop(
        columns=["contact_residue_set"]
    )

pose_out = Path("prototype2/p2_pose_interface_analysis.csv")
pose_df.to_csv(pose_out, index=False)

candidate_df = pd.DataFrame(candidate_rows)

# Primary ordering:
# 1. contact consistency
# 2. site proximity
# 3. clustering concentration
# This deliberately does NOT turn everything into one fake score.
candidate_df = candidate_df.sort_values(
    [
        "mean_contact_set_jaccard",
        "median_site_distance_A",
    ],
    ascending=[
        False,
        True,
    ],
)

candidate_out = Path(
    "prototype2/p2_interface_candidate_summary.csv"
)
candidate_df.to_csv(candidate_out, index=False)

freq_df = pd.DataFrame(contact_frequency_rows)

freq_out = Path(
    "prototype2/p2_contact_residue_frequency.csv"
)
freq_df.to_csv(freq_out, index=False)


print("\n" + "=" * 90)
print("P2 INTERFACE CONSISTENCY")
print("=" * 90)

if len(candidate_df):
    display_cols = [
        "candidate",
        "n_analyzed_poses",
        "median_site_distance_A",
        "fraction_within_15A_of_site",
        "fraction_within_20A_of_site",
        "median_structural_contact_residues",
        "mean_contact_set_jaccard",
        "n_clusters",
        "largest_cluster_fraction",
        "anchor_203_contact_fraction",
        "anchor_205_contact_fraction",
        "either_anchor_contact_fraction",
    ]

    print(
        candidate_df[display_cols]
        .round(3)
        .to_string(index=False)
    )

print("\nSaved:")
print(pose_out)
print(candidate_out)
print(freq_out)
