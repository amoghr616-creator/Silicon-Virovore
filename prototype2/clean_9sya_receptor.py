from Bio.PDB import PDBParser, PDBIO, Select

INPUT = "prototype2/9SYA_chainA_final.pdb"
OUTPUT = "prototype2/9SYA_chainA_meeko_clean.pdb"

class CleanSelect(Select):
    def accept_residue(self, residue):
        # Remove incomplete GLN A:109
        if residue.id[0] == " " and residue.id[1] == 109:
            return 0
        return 1

parser = PDBParser(QUIET=True)
structure = parser.get_structure("9SYA", INPUT)

io = PDBIO()
io.set_structure(structure)
io.save(OUTPUT, CleanSelect())

print(f"Saved: {OUTPUT}")
