from pathlib import Path
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser

BASE = Path("prototype2/benchmark_1ycr")
POSE_DIR = BASE / "rapidock_20" / "poses"
RANKED = BASE / "rapidock_20" / "ranked_poses.csv"
EXP = BASE / "1YCR.pdb"

SEQUENCE = "ETFSDLWKLLPEN"


AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}


def chain_sequence(chain):
    residues = []
    seq = []

    for r in chain:
        if r.id[0] != " ":
            continue
        if r.resname not in AA3_TO_1:
            continue
        seq.append(AA3_TO_1[r.resname])
        residues.append(r)

    return "".join(seq), residues


def find_chain(structure, sequence):
    matches = []

    for model in structure:
        for chain in model:
            seq, residues = chain_sequence(chain)
            if seq == sequence:
                matches.append((chain.id, residues))

    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one {sequence} chain, found {len(matches)}"
        )

    return matches[0]


def ca_coords(residues):
    return np.array([r["CA"].coord for r in residues])


def rmsd(a, b):
    if a.shape != b.shape:
        raise ValueError(f"Shape mismatch: {a.shape} vs {b.shape}")
    return float(np.sqrt(np.mean(np.sum((a - b)**2, axis=1))))


parser = PDBParser(QUIET=True)
exp_struct = parser.get_structure("exp", EXP)

_, exp_res = find_chain(exp_struct, SEQUENCE)
exp_ca = ca_coords(exp_res)

ranked = pd.read_csv(RANKED)

results = []

for _, row in ranked.iterrows():

    pose_file = row["pose_filename"]
    pose_path = POSE_DIR / pose_file

    if not pose_path.exists():
        print(f"WARNING: missing {pose_path}")
        continue

    pose_struct = parser.get_structure(
        f"pose_{row['rank']}",
        pose_path
    )

    _, pose_res = find_chain(pose_struct, SEQUENCE)
    pose_ca = ca_coords(pose_res)

    if len(pose_ca) != len(exp_ca):
        print(f"WARNING: atom count mismatch for {pose_file}")
        continue

    pose_rmsd = rmsd(pose_ca, exp_ca)

    results.append({
        "rank": int(row["rank"]),
        "pose": pose_file,
        "pose_rmsd_A": pose_rmsd,
        "delta_g_kcal_mol": row["delta_g"],
        "rank_score": row["rank_score"],
        "n_contact_residues": row["n_contact_residues"],
        "is_clashed": row["is_clashed"],
        "charged_confidence": row["charged_confidence"],
    })


out = pd.DataFrame(results)
out = out.sort_values("pose_rmsd_A")

print("\n=== 1YCR RAPiDock pose recovery ===")
print(out.to_string(index=False))

best = out.iloc[0]

print("\n=== BEST POSE BY EXPERIMENTAL RMSD ===")
print(f"Pose: {best['pose']}")
print(f"RMSD: {best['pose_rmsd_A']:.3f} A")
print(f"RAPiDock rank: {int(best['rank'])}")
print(f"Delta G: {best['delta_g_kcal_mol']:.4f} kcal/mol")

rank1 = out[out["rank"] == 1].iloc[0]

print("\n=== TOP-RANKED POSE ===")
print(f"Pose: {rank1['pose']}")
print(f"RMSD: {rank1['pose_rmsd_A']:.3f} A")
print(f"Delta G: {rank1['delta_g_kcal_mol']:.4f} kcal/mol")

out_path = BASE / "rapidock_20" / "pose_rmsd_all.csv"
out.to_csv(out_path, index=False)

print(f"\nSaved: {out_path}")
