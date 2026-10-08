from pathlib import Path
import numpy as np

experimental = Path("prototype2/benchmark_1f4v/flim_chainD.pdb")
predicted = Path("prototype2/benchmark_1f4v/fliM_predicted.pdb")

def get_ca(path, chain):
    atoms = []

    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[21] != chain:
            continue
        if line[12:16].strip() != "CA":
            continue

        residue = int(line[22:26])

        xyz = np.array([
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        ])

        atoms.append((residue, xyz))

    return atoms

def kabsch_rmsd(ref, mobile):
    ref_center = ref.mean(axis=0)
    mobile_center = mobile.mean(axis=0)

    X = ref - ref_center
    Y = mobile - mobile_center

    H = Y.T @ X
    U, S, Vt = np.linalg.svd(H)

    D = np.eye(3)
    if np.linalg.det(Vt.T @ U.T) < 0:
        D[-1, -1] = -1

    R = Vt.T @ D @ U.T

    aligned = Y @ R.T

    return np.sqrt(
        np.mean(np.sum((aligned - X) ** 2, axis=1))
    )

exp = get_ca(experimental, "D")
pred = get_ca(predicted, "A")

print("Experimental CA atoms:", len(exp))
print("Predicted CA atoms:", len(pred))

if len(exp) != len(pred):
    raise SystemExit(
        f"ERROR: residue counts differ: "
        f"experimental={len(exp)}, predicted={len(pred)}"
    )

exp_xyz = np.array([x[1] for x in exp])
pred_xyz = np.array([x[1] for x in pred])

rmsd = kabsch_rmsd(exp_xyz, pred_xyz)

print()
print("=== 1F4V STRUCTURE VALIDATION ===")
print("Experimental: PDB 1F4V, chain D")
print("Predicted:    MGDSILSQAEIDALLN, ESMFold")
print("Residues compared:", len(exp))
print(f"C-alpha RMSD: {rmsd:.3f} Å")

Path(
    "prototype2/benchmark_1f4v/structure_validation.txt"
).write_text(
    f"PDB=1F4V\n"
    f"experimental_chain=D\n"
    f"predicted_sequence=MGDSILSQAEIDALLN\n"
    f"n_residues={len(exp)}\n"
    f"ca_rmsd_angstrom={rmsd:.6f}\n"
)

print(
    "Saved: "
    "prototype2/benchmark_1f4v/structure_validation.txt"
)
