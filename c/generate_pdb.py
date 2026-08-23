# generate_pdb.py
import math

# Standard back-bone dihedral angles for an ideal alpha-helix
PHI = math.radians(-57.0)
PSI = math.radians(-47.0)
BOND_N_CA = 1.458   # Angstroms
BOND_CA_C = 1.525
BOND_C_N  = 1.329

def generate_alpha_helix_pdb(sequence: str, output_filename: str = "lead_candidate_1.pdb"):
    atoms = []
    atom_id = 1
    
    # Helical translation geometry parameters per residue
    rise_per_residue = 1.5    # Angstroms along z-axis
    rotation_per_residue = math.radians(100.0) # 3.6 residues per turn
    radius = 2.3              # Radius of CA backbone cylinder in Angstroms

    with open(output_filename, "w") as f:
        f.write(f"HEADER    SILICON VIROVORE LEAD PEPTIDE 3D MODEL\n")
        f.write(f"TITLE     ALPHA-HELICAL MODEL FOR {sequence}\n")
        f.write(f"REMARK    GENERATED VIA IDEALIZED BACKBONE DIHEDRAL GEOMETRY\n")
        
        for i, aa in enumerate(sequence, 1):
            angle = (i - 1) * rotation_per_residue
            z_ca = (i - 1) * rise_per_residue
            
            # Alpha Carbon (CA)
            ca_x = radius * math.cos(angle)
            ca_y = radius * math.sin(angle)
            ca_z = z_ca
            
            # Back-bone Nitrogen (N)
            n_x = (radius - 0.5) * math.cos(angle - 0.3)
            n_y = (radius - 0.5) * math.sin(angle - 0.3)
            n_z = ca_z - 0.7
            
            # Back-bone Carbonyl (C)
            c_x = (radius + 0.4) * math.cos(angle + 0.3)
            c_y = (radius + 0.4) * math.sin(angle + 0.3)
            c_z = ca_z + 0.6
            
            # Carbonyl Oxygen (O)
            o_x = c_x + 0.8 * math.cos(angle + 0.5)
            o_y = c_y + 0.8 * math.sin(angle + 0.5)
            o_z = c_z + 0.4

            # Write PDB ATOM lines
            f.write(f"ATOM  {atom_id:5d}  N   {aa:>3s} A{i:4d}    {n_x:8.3f}{n_y:8.3f}{n_z:8.3f}  1.00 20.00           N\n")
            atom_id += 1
            f.write(f"ATOM  {atom_id:5d}  CA  {aa:>3s} A{i:4d}    {ca_x:8.3f}{ca_y:8.3f}{ca_z:8.3f}  1.00 20.00           C\n")
            atom_id += 1
            f.write(f"ATOM  {atom_id:5d}  C   {aa:>3s} A{i:4d}    {c_x:8.3f}{c_y:8.3f}{c_z:8.3f}  1.00 20.00           C\n")
            atom_id += 1
            f.write(f"ATOM  {atom_id:5d}  O   {aa:>3s} A{i:4d}    {o_x:8.3f}{o_y:8.3f}{o_z:8.3f}  1.00 20.00           O\n")
            atom_id += 1

        f.write("END\n")

    print(f"[EXPORT] PDB structure successfully written to '{output_filename}'")

if __name__ == "__main__":
    lead_seq = "MKLAVFALLVFFAGSSDLIRR"
    generate_alpha_helix_pdb(lead_seq, "MKLAVFALLVFFAGSSDLIRR.pdb")