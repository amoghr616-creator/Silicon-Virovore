from openmm.app import PDBFile, Modeller

# Load your pipeline's raw output file
pdb = PDBFile('lead_1.pdb')
modeller = Modeller(pdb.topology, pdb.positions)

# Clean and write out a standard protein PDB file
with open('lead_1_pipeline_clean.pdb', 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print("Successfully created lead_1_pipeline_clean.pdb!")
