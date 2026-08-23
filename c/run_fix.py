import subprocess

# Rebuilds a clean standard protein PDB from your original file
with open('generate_clean.py', 'w') as f:
    f.write('''
from openmm.app import PDBFile, Modeller
pdb = PDBFile('lead_1.pdb')
modeller = Modeller(pdb.topology, pdb.positions)
with open('lead_1_fixed_final.pdb', 'w') as out:
    PDBFile.writeFile(modeller.topology, modeller.positions, out)
''')

subprocess.run(['python', 'generate_clean.py'])

# Executes GROMACS with properly separated arguments
subprocess.run([
    'gmx', 'pdb2gmx',
    '-f', 'lead_1_fixed_final.pdb',
    '-o', 'lead_1_processed.gro',
    '-water', 'tip3p',
    '-ff', 'charmm36-feb2026',
    '-ignh', '-ter'
])
