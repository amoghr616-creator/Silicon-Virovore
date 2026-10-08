from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser

BASE = Path("prototype2/p2_rapidock/P2-01-clean")

RECEPTOR = BASE / "receptor_for_rapidock.pdb"
POSE = BASE / "best_pose_vina_relaxed.pdb"
RELAXED = BASE / "best_pose_vina_relaxed.pdb"


def heavy_atoms(structure):
    atoms = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue

                for atom in residue:
                    if atom.element.upper() != "H":
                        atoms.append(atom)

    return atoms


def split_atoms(structure, peptide_sequence):
    aa = {
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
        "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
        "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
        "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
    }

    peptide = []
    receptor = []

    for model in structure:
        for chain in model:
            seq = ""
            residues = []

            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if residue.resname not in aa:
                    continue

                seq += aa[residue.resname]
                residues.append(residue)

            if seq == peptide_sequence:
                for residue in residues:
                    for atom in residue:
                        if atom.element.upper() != "H":
                            peptide.append(atom)

    # Receptor comes from separate receptor structure,
    # so this function is only used to locate the peptide.
    return peptide


def distance_stats(receptor_atoms, peptide_atoms):
    receptor_xyz = np.array([a.coord for a in receptor_atoms])
    peptide_xyz = np.array([a.coord for a in peptide_atoms])

    # Chunk to avoid making an enormous matrix.
    minimum = float("inf")
    counts = {1.5: 0, 2.0: 0, 2.5: 0, 3.0: 0}

    for p in peptide_xyz:
        d = np.sqrt(np.sum((receptor_xyz - p) ** 2, axis=1))
        minimum = min(minimum, float(d.min()))

        for cutoff in counts:
            counts[cutoff] += int(np.sum(d < cutoff))

    return minimum, counts


parser = PDBParser(QUIET=True)

receptor = parser.get_structure("receptor", RECEPTOR)
pose = parser.get_structure("pose", POSE)

receptor_atoms = heavy_atoms(receptor)
peptide_atoms = split_atoms(
    pose,
    "MKQAVNALLVFFAGSSDAIRR"
)

if not peptide_atoms:
    raise RuntimeError("Could not identify P2-01 peptide chain in best_pose.pdb")

minimum, counts = distance_stats(
    receptor_atoms,
    peptide_atoms
)

print("=== P2-01 BEST POSE GEOMETRY ===")
print("Peptide heavy atoms:", len(peptide_atoms))
print("Receptor heavy atoms:", len(receptor_atoms))
print(f"Minimum receptor-peptide heavy-atom distance: {minimum:.3f} A")

for cutoff, count in counts.items():
    print(f"Pairs < {cutoff:.1f} A: {count}")
