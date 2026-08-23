
from openmm.app import PDBFile, Modeller
pdb = PDBFile('lead_1.pdb')
modeller = Modeller(pdb.topology, pdb.positions)
with open('lead_1_fixed_final.pdb', 'w') as out:
    PDBFile.writeFile(modeller.topology, modeller.positions, out)
