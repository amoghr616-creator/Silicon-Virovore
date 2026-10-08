from pathlib import Path
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

BASE = Path("prototype2/p2_rapidock/P2-01-clean")
POSE_DIR = BASE / "poses"
RANKED = BASE / "ranked_poses.csv"

SITE = np.array([-8.405, -5.823, 1.523], dtype=float)
SITE_RADIUS = 15.0

SEQUENCE = "MKQAVNALLVFFAGSSDAIRR"

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}

def find_peptide_chain(structure, sequence):
    matches = []

    for model in structure:
        for chain in model:
            seq = []
            ca = []

            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if residue.resname not in AA3_TO_1:
                    continue

                seq.append(AA3_TO_1[residue.resname])

                if "CA" in residue:
                    ca.append(residue["CA"].coord)

            if "".join(seq) == sequence:
                matches.append(np.array(ca, dtype=float))

    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one peptide chain matching {sequence}, "
            f"found {len(matches)}"
        )

    return matches[0]


parser = PDBParser(QUIET=True)
ranked = pd.read_csv(RANKED)

rows = []

for _, row in ranked.iterrows():

    pose_path = POSE_DIR / row["pose_filename"]

    if not pose_path.exists():
        continue

    structure = parser.get_structure("pose", pose_path)

    try:
        ca = find_peptide_chain(structure, SEQUENCE)
    except RuntimeError as e:
        print(f"{row['pose_filename']}: {e}")
        continue

    centroid = ca.mean(axis=0)
    distance = float(np.linalg.norm(centroid - SITE))

    rows.append({
        "rank": int(row["rank"]),
        "pose": row["pose_filename"],
        "site_distance_A": distance,
        "on_site_15A": distance <= SITE_RADIUS,
        "delta_g_kcal_mol": float(row["delta_g"]),
        "rank_score": float(row["rank_score"]),
        "contacts": int(row["n_contact_residues"]),
        "clashed": bool(row["is_clashed"]),
    })

out = pd.DataFrame(rows)

print("\n=== P2-01 SITE FILTER ===")
print(out.sort_values("site_distance_A").to_string(index=False))

usable = out[
    (out["on_site_15A"]) &
    (~out["clashed"])
].copy()

print("\n=== USABLE ON-SITE POSES ===")

if usable.empty:
    print("NO USABLE POSES WITHIN 15 A OF THE HYPOTHESIZED SITE")
else:
    usable = usable.sort_values(
        ["rank_score", "contacts"],
        ascending=[True, False]
    )

    print(usable.to_string(index=False))

    best = usable.iloc[0]

    print("\n=== BEST USABLE P2-01 POSE ===")
    print(f"Pose: {best['pose']}")
    print(f"Site distance: {best['site_distance_A']:.3f} A")
    print(f"Rank score: {best['rank_score']:.4f}")
    print(f"Predicted delta G: {best['delta_g_kcal_mol']:.4f} kcal/mol")
    print(f"Contacts: {best['contacts']}")

out_path = BASE / "site_filtered_poses.csv"
out.to_csv(out_path, index=False)

print(f"\nSaved: {out_path}")
