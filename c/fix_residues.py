from Bio import PDB

# Mapping 1-letter codes to 3-letter codes
d3 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
      'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
      'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
      'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}

parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure('lead', 'lead_1.pdb')

for model in structure:
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip()
            if len(resname) == 1 and resname in d3:
                residue.resname = d3[resname]

io = PDB.PDBIO()
io.set_structure(structure)
io.save('lead_1_3letter.pdb')
print("Converted residues to standard 3-letter codes in lead_1_3letter.pdb!")
