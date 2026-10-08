from pathlib import Path
import urllib.request
import numpy as np

PDB_URL = "https://files.rcsb.org/download/1YCR.pdb"

experimental = Path("prototype2/benchmark_1ycr/1YCR.pdb")
predicted = Path("results/structures/ETFSDLWKLLPEN.pdb")

experimental.parent.mkdir(parents=True, exist_ok=True)

if not experimental.exists():
    print("Downloading 1YCR experimental structure...")
    urllib.request.urlretrieve(PDB_URL, experimental)

def get_ca(path, chain):
    atoms = []

    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[21] != chain:
            continue
        if line[12:16].strip() != "CA":
            continue

        try:
            residue = int(line[22:26])
            xyz = np.array([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ])
        except ValueError:
            continue

        atoms.append((residue, xyz))

    return atoms


def kabsch(ref, mobile):
    ref_center = ref.mean(axis=0)
    mobile_center = mobile.mean(axis=0)

    X = ref - ref_center
    Y = mobile - mobile_center

    H = Y.T @ X
    U, S, Vt = np.linalg.svd(H)

    d = np.linalg.det(Vt.T @ U.T)

    D = np.eye(3)
    if d < 0:
        D[-1, -1] = -1

    R = Vt.T @ D @ U.T

    aligned = Y @ R.T
    rmsd = np.sqrt(np.mean(np.sum((aligned - X) ** 2, axis=1)))

    return rmsd


exp = get_ca(experimental, "B")
pred = get_ca(predicted, "A")

print("Experimental CA atoms:", len(exp))
print("Predicted CA atoms:", len(pred))

if len(exp) != len(pred):
    raise SystemExit(
        f"ERROR: residue counts differ: experimental={len(exp)}, predicted={len(pred)}"
    )

exp_xyz = np.array([x[1] for x in exp])
pred_xyz = np.array([x[1] for x in pred])

rmsd = kabsch(exp_xyz, pred_xyz)

print()
print("=== 1YCR STRUCTURE VALIDATION ===")
print("Experimental: PDB 1YCR, chain B")
print("Predicted:    ETFSDLWKLLPEN, ESMFold")
print("Residues compared:", len(exp))
print(f"C-alpha RMSD: {rmsd:.3f} Å")
print("pLDDT: 72.96")

# Save machine-readable result
out = Path("prototype2/benchmark_1ycr/structure_validation.txt")
out.write_text(
    f"PDB=1YCR\n"
    f"experimental_chain=B\n"
    f"predicted_sequence=ETFSDLWKLLPEN\n"
    f"predicted_backend=ESMFold\n"
    f"plddt=72.96\n"
    f"n_residues={len(exp)}\n"
    f"ca_rmsd_angstrom={rmsd:.6f}\n"
)

print("Saved:", out)
