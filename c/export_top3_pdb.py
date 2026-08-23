# export_top3_pdbs.py
import math

LEADS = {
    "lead_1": "MKLAVFALLVFFAGSSDLIRR",
    "lead_2": "MFLAVFALLVTFAGSSDKKLR",
    "lead_3": "MFRAVDALLVFFAFFFDKVLR"
}

def build_pdb(sequence: str, filename: str):
    rise_per_residue = 1.5
    rotation_per_residue = math.radians(100.0)
    radius = 2.3
    atom_id = 1

    with open(filename, "w") as f:
        f.write(f"HEADER    SILICON VIROVORE LEAD PEPTIDE MODEL\n")
        f.write(f"TITLE     ALPHA-HELICAL BACKBONE FOR {sequence}\n")
        
        for i, aa in enumerate(sequence, 1):
            angle = (i - 1) * rotation_per_residue
            ca_z = (i - 1) * rise_per_residue
            ca_x = radius * math.cos(angle)
            ca_y = radius * math.sin(angle)
            
            n_x = (radius - 0.5) * math.cos(angle - 0.3)
            n_y = (radius - 0.5) * math.sin(angle - 0.3)
            n_z = ca_z - 0.7
            
            c_x = (radius + 0.4) * math.cos(angle + 0.3)
            c_y = (radius + 0.4) * math.sin(angle + 0.3)
            c_z = ca_z + 0.6
            
            o_x = c_x + 0.8 * math.cos(angle + 0.5)
            o_y = c_y + 0.8 * math.sin(angle + 0.5)
            o_z = c_z + 0.4

            f.write(f"ATOM  {atom_id:5d}  N   {aa:>3s} A{i:4d}    {n_x:8.3f}{n_y:8.3f}{n_z:8.3f}  1.00 20.00           N\n")
            atom_id += 1
            f.write(f"ATOM  {atom_id:5d}  CA  {aa:>3s} A{i:4d}    {ca_x:8.3f}{ca_y:8.3f}{ca_z:8.3f}  1.00 20.00           C\n")
            atom_id += 1
            f.write(f"ATOM  {atom_id:5d}  C   {aa:>3s} A{i:4d}    {c_x:8.3f}{c_y:8.3f}{c_z:8.3f}  1.00 20.00           C\n")
            atom_id += 1
            f.write(f"ATOM  {atom_id:5d}  O   {aa:>3s} A{i:4d}    {o_x:8.3f}{o_y:8.3f}{o_z:8.3f}  1.00 20.00           O\n")
            atom_id += 1

        f.write("END\n")
    print(f"[EXPORT] Created {filename} ({sequence})")

if __name__ == "__main__":
    for name, seq in LEADS.items():
        build_pdb(seq, f"{name}.pdb")