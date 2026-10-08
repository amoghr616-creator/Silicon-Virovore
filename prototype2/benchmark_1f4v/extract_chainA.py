from Bio.PDB import PDBParser, PDBIO, Select

INPUT = "prototype2/benchmark_1f4v/1F4V.pdb"
OUTPUT = "prototype2/benchmark_1f4v/chey_chainA.pdb"

class ChainASelect(Select):
    def accept_chain(self, chain):
        return chain.id == "A"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("1f4v", INPUT)

io = PDBIO()
io.set_structure(structure)
io.save(OUTPUT, ChainASelect())

print(f"Saved {OUTPUT}")
