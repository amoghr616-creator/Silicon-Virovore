from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Load the original backbone/peptide
fixer = PDBFixer(filename='lead_1.pdb')

# Explicitly add missing heavy side-chain atoms
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# Write out a natively formatted PDB file that GROMACS loves
with open('lead_1_gromacs.pdb', 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

print("Saved clean GROMACS-ready file as lead_1_gromacs.pdb")
