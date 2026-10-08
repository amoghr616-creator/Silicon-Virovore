from __future__ import annotations

import sys
import numpy as np
from Bio.PDB import PDBParser


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def get_chain_sequence(chain):
    sequence = []
    residues = []

    for residue in chain:
        if residue.id[0] != " ":
            continue
        if residue.resname not in AA3_TO_1:
            continue

        sequence.append(AA3_TO_1[residue.resname])
        residues.append(residue)

    return "".join(sequence), residues


def find_matching_chain(structure, target_sequence):
    matches = []

    for model in structure:
        for chain in model:
            sequence, residues = get_chain_sequence(chain)

            if sequence == target_sequence:
                matches.append((chain.id, residues))

    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one chain matching {target_sequence}, "
            f"but found {len(matches)}"
        )

    return matches[0]


def get_ca_coordinates(residues):
    return np.array(
        [residue["CA"].coord for residue in residues],
        dtype=float,
    )


def calculate_rmsd(docked, experimental):
    if docked.shape != experimental.shape:
        raise ValueError(
            f"Coordinate shape mismatch: "
            f"{docked.shape} vs {experimental.shape}"
        )

    differences = docked - experimental

    return float(
        np.sqrt(np.mean(np.sum(differences ** 2, axis=1)))
    )


def main():
    if len(sys.argv) != 4:
        print(
            "Usage:\n"
            "python prototype2/pose_rmsd.py "
            "<docked.pdb> <experimental.pdb> <sequence>"
        )
        sys.exit(1)

    docked_path = sys.argv[1]
    experimental_path = sys.argv[2]
    peptide_sequence = sys.argv[3].strip().upper()

    parser = PDBParser(QUIET=True)

    docked_structure = parser.get_structure("docked", docked_path)
    experimental_structure = parser.get_structure(
        "experimental",
        experimental_path,
    )

    docked_chain_id, docked_residues = find_matching_chain(
        docked_structure,
        peptide_sequence,
    )

    experimental_chain_id, experimental_residues = find_matching_chain(
        experimental_structure,
        peptide_sequence,
    )

    docked_ca = get_ca_coordinates(docked_residues)
    experimental_ca = get_ca_coordinates(experimental_residues)

    rmsd = calculate_rmsd(
        docked_ca,
        experimental_ca,
    )

    print(f"Docked peptide chain: {docked_chain_id}")
    print(f"Experimental peptide chain: {experimental_chain_id}")
    print(f"Peptide C-alpha atoms: {len(docked_ca)}")
    print(f"Direct docking-pose C-alpha RMSD: {rmsd:.3f} Å")


if __name__ == "__main__":
    main()
