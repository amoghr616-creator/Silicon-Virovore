from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np

INPUT = "prototype2/9SYA_chainA_meeko_clean.pdb"
OUTPUT = "prototype2/9SYA_pocket15.pdb"

SITE = np.array([-8.405, -5.823, 1.523], dtype=float)
RADIUS = 15.0

class PocketSelect(Select):
    def accept_residue(self, residue):

        # Keep standard amino acids only.
        if residue.id[0] != " ":
            return 0

        # Keep the whole residue if any atom is within the pocket radius.
        for atom in residue:
            d = np.linalg.norm(atom.coord - SITE)

            if d <= RADIUS:
                return 1

        return 0


parser = PDBParser(QUIET=True)
structure = parser.get_structure("9SYA", INPUT)

io = PDBIO()
io.set_structure(structure)
io.save(OUTPUT, PocketSelect())

# Report what was kept.
structure2 = parser.get_structure("pocket", OUTPUT)
chain = structure2[0]["A"]

residues = [
    r for r in chain
    if r.id[0] == " "
]

print(f"Saved: {OUTPUT}")
print(f"Pocket radius: {RADIUS:.1f} A")
print(f"Residues retained: {len(residues)}")
print(
    "Residues:",
    [(r.id[1], r.resname) for r in residues]
)
