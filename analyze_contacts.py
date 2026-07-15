import math

def calculate_distance(coord1, coord2):
    return math.sqrt(sum((c1 - c2) ** 2 for c1, c2 in zip(coord1, coord2)))

receptor_atoms = []
peptide_atoms = []

docking_file_path = "/Users/centurion616/Desktop/Silicon-Virovore/docking.res.1.pdb"

try:
    with open(docking_file_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith("ATOM"):
                # Robust parsing: first try standard PDB columns
                try:
                    chain = line[21].strip()
                    res_seq = int(line[22:26])
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except (ValueError, IndexError):
                    # Fallback to whitespace splitting if the format is shifted
                    parts = line.split()
                    if len(parts) >= 11:
                        atom_name = parts[2]
                        res_name = parts[3]
                        # If chain is merged or missing
                        if len(parts[4]) == 1 and parts[4].isalpha():
                            chain = parts[4]
                            res_seq = int(parts[5])
                        else:
                            chain = ''
                            res_seq = int(parts[4])
                        x, y, z = float(parts[-6]), float(parts[-5]), float(parts[-4])
                    else:
                        continue
                
                atom_data = {
                    "res_seq": res_seq,
                    "res_name": res_name,
                    "atom_name": atom_name,
                    "coords": (x, y, z),
                    "line_num": line_num
                }
                
                if chain == 'A':
                    receptor_atoms.append(atom_data)
                elif chain == '' or chain == ' ':
                    peptide_atoms.append(atom_data)
                    
except FileNotFoundError:
    print(f"Error: Could not find the file at {docking_file_path}")
    exit(1)

# --- DIAGNOSTIC PRINT STATEMENTS ---
print("="*60)
print("DIAGNOSTIC REPORT")
print("="*60)
print(f"Total Receptor (Chain A) atoms loaded: {len(receptor_atoms)}")
if receptor_atoms:
    print(f"  First receptor atom: Residue {receptor_atoms[0]['res_name']}-{receptor_atoms[0]['res_seq']}, Coords: {receptor_atoms[0]['coords']}")
    print(f"  Last receptor atom:  Residue {receptor_atoms[-1]['res_name']}-{receptor_atoms[-1]['res_seq']}, Coords: {receptor_atoms[-1]['coords']}")

print(f"\nTotal Peptide (Blank Chain) atoms loaded: {len(peptide_atoms)}")
if peptide_atoms:
    print(f"  First peptide atom:  Residue {peptide_atoms[0]['res_name']}-{peptide_atoms[0]['res_seq']}, Coords: {peptide_atoms[0]['coords']}")
    print(f"  Last peptide atom:   Residue {peptide_atoms[-1]['res_name']}-{peptide_atoms[-1]['res_seq']}, Coords: {peptide_atoms[-1]['coords']}")
print("="*60 + "\n")

# Find interactions within 4.0 Angstroms
contacts = {}
for p_atom in peptide_atoms:
    for r_atom in receptor_atoms:
        dist = calculate_distance(p_atom["coords"], r_atom["coords"])
        if dist < 4.0:
            p_res = f"{p_atom['res_name']}-{p_atom['res_seq']}"
            r_res = f"{r_atom['res_name']}-{r_atom['res_seq']}"
            if p_res not in contacts:
                contacts[p_res] = set()
            contacts[p_res].add((r_res, round(dist, 2)))

if contacts:
    print("   CRITICAL INTERACTION MAP (Distance < 4.0 Å)")
    print("="*60)
    for p_res, r_list in sorted(contacts.items(), key=lambda x: int(x[0].split('-')[1])):
        targets = ", ".join([f"{r[0]} ({r[1]}Å)" for r in sorted(r_list)])
        print(f"Peptide Residue {p_res:<8} ---> binds Receptor: {targets}")
    print("="*60)
else:
    print("No contacts found within 4.0 Å.")