from pathlib import Path
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

FILES = [
    Path("prototype4/data/haddock3/receptors/9sya_source_chain_A.pdb"),
    Path("prototype4/data/haddock3/receptors/9sya_source_chain_B.pdb"),
    Path("prototype2/9SYA_chainA.pdb"),
    Path("data/receptor/9SYA.pdb"),
]

parser = PDBParser(QUIET=True)

for path in FILES:
    print(f"\n{'=' * 72}")
    print("FILE:", path)

    if not path.exists():
        print("MISSING")
        continue

    structure = parser.get_structure(path.stem, str(path))
    model = next(structure.get_models())

    for chain in model:
        residues = [
            r for r in chain
            if is_aa(r, standard=False)
        ]

        if not residues:
            continue

        sequence = "".join(
            seq1(r.resname, custom_map={"MSE": "M"})
            for r in residues
        )

        print(
            f"\nChain {chain.id!r}: "
            f"{len(residues)} amino-acid residues"
        )
        print("First residue:", residues[0].id)
        print("Last residue:", residues[-1].id)
        print("Sequence length:", len(sequence))
        print("Sequence preview:", sequence[:35], "...", sequence[-20:])

        for target_number in (203, 205):
            matches = [
                r for r in residues
                if r.resname == "ARG" and r.id[1] == target_number
            ]

            if not matches:
                print(f"ARG {target_number}: NOT FOUND")
                continue

            for residue in matches:
                required_atoms = {
                    "NE", "NH1", "NH2", "CZ"
                }
                present = required_atoms.intersection(
                    residue.child_dict
                )

                print(
                    f"ARG {target_number}: "
                    f"chain={chain.id}, "
                    f"residue_id={residue.id}, "
                    f"guanidinium_atoms={sorted(present)}"
                )

print("\nAUDIT COMPLETE")
