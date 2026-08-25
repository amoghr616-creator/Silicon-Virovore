from Bio.PDB import PDBParser
import numpy as np

def calculate_hydrophobic_score(pdb_file_path):
    """
    Parses a PDB docking file and calculates a Hydrophobic Contact Score (HCS)
    based on close-range carbon-carbon interactions between the peptide and receptor.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file_path)
    
    # In PatchDock outputs, the receptor is typically Chain A and the ligand/peptide is Chain B.
    # Adjust chain IDs if your specific PDB uses different naming.
    model = structure[0]
    
    receptor_carbons = []
    peptide_carbons = []
    
    # Step 1: Extract all Carbon atoms (representing hydrophobic sidechains)
    for chain in model:
        for residue in chain:
            for atom in residue:
                # We target Carbon atoms ('C', 'CA', 'CB', 'CG', 'CD', etc.) for hydrophobic mapping
                if atom.element == 'C':
                    if chain.id == 'A': # Receptor
                        receptor_carbons.append(atom)
                    elif chain.id == 'B' or chain.id == ' ': # Peptide (often blank or 'B')
                        peptide_carbons.append(atom)

    if not receptor_carbons or not peptide_carbons:
        # Fallback: if chains aren't explicitly split, let's look at first vs second molecule
        chains = list(model.get_chains())
        if len(chains) >= 2:
            receptor_carbons = [a for r in chains[0].get_residues() for a in r if a.element == 'C']
            peptide_carbons = [a for r in chains[1].get_residues() for a in r if a.element == 'C']
        else:
            return 0.0, 0

    # Step 2: Calculate 3D distances between all carbon pairs
    hydrophobic_contacts = 0
    total_score = 0.0
    threshold = 4.5 # Angstroms (standard limit for hydrophobic interactions)

    for r_atom in receptor_carbons:
        r_coord = r_atom.get_coord()
        for p_atom in peptide_carbons:
            p_coord = p_atom.get_coord()
            
            # Calculate Euclidean distance in 3D space: sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
            distance = np.linalg.norm(r_coord - p_coord)
            
            if distance <= threshold:
                hydrophobic_contacts += 1
                # Closer contacts contribute exponentially more to thermodynamic stability
                total_score += (threshold - distance)

    return round(total_score, 2), hydrophobic_contacts

if __name__ == "__main__":
    # Test the script on your PDB file
    pdb_path = "docking.res.1.pdb"
    try:
        score, contacts = calculate_hydrophobic_score(pdb_path)
        print(f"[*] Analysis of {pdb_path}:")
        print(f"[-] Total Hydrophobic Contacts: {contacts}")
        print(f"[-] Calculated Hydrophobic Score: {score}")
    except Exception as e:
        print(f"[X] Error reading PDB: {e}")