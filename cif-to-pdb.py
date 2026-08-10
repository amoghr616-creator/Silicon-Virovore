import os

def convert_cif_text_to_pdb(cif_lines, output_pdb_path):
    pdb_lines = []
    atom_started = False
    
    for line in cif_lines:
        line = line.strip()
        
        if line.startswith("_atom_site."):
            atom_started = True
            continue
            
        if atom_started and (line.startswith("ATOM") or line.startswith("HETATM")):
            parts = line.split()
            if len(parts) < 13:
                continue
                
            try:
                group_pdb = parts[0][:6]     # "ATOM" or "HETATM"
                atom_id = int(parts[1])      # Atom serial number
                element = parts[2][:2]       # Chemical element (e.g., N, C, O, S)
                atom_name = parts[3][:4]     # Atom name (e.g., CA, CB, CG)
                res_name = parts[4][:3]      # Residue name (e.g., MET, LYS, ALA)
                seq_id = int(parts[5])       # Sequence ID
                chain_id = parts[8][:1]      # Chain ID (must be exactly 1 character)
                
                x = float(parts[10])
                y = float(parts[11])
                z = float(parts[12])
                
                occupancy = float(parts[13]) if len(parts) > 13 and parts[13] != '?' else 1.00
                b_factor = float(parts[14]) if len(parts) > 14 and parts[14] != '?' else 90.00
                
                # Standard PDB formatting rules for atom names:
                # 1-2 character elements (like C, N, O, S) start at column 14.
                # 4-letter atom names (like 1H5') start at column 13.
                if len(atom_name) < 4:
                    formatted_atom_name = f" {atom_name:<3}"
                else:
                    formatted_atom_name = f"{atom_name:<4}"

                # STRICT COLUMNS (PDB Official Format):
                # Col 1-6   : Record name (ATOM)
                # Col 7-11  : Atom Serial Number
                # Col 13-16 : Atom Name
                # Col 17    : Alternate Location Indicator (blank)
                # Col 18-20 : Residue Name (e.g., MET)
                # Col 21    : Chain Identifier (e.g., A)
                # Col 22-26 : Residue Sequence Number
                # Col 31-38 : X coordinate
                # Col 39-46 : Y coordinate
                # Col 47-54 : Z coordinate
                # Col 55-60 : Occupancy
                # Col 61-66 : Temperature factor (B-factor)
                # Col 77-78 : Element symbol
                
                pdb_line = (
                    f"{group_pdb:<6}"          # 1-6
                    f"{atom_id:>5}"            # 7-11
                    f" "                       # 12 (blank)
                    f"{formatted_atom_name}"   # 13-16
                    f" "                       # 17 (blank)
                    f"{res_name:>3}"           # 18-20 (Residue Name)
                    f"{chain_id:>1}"           # 21 (Chain Identifier)
                    f"{seq_id:>4}"             # 22-25 (Sequence ID)
                    f" "                       # 26 (Insertion Code)
                    f"   "                     # 27-30 (blank)
                    f"{x:>8.3f}"               # 31-38
                    f"{y:>8.3f}"               # 39-46
                    f"{z:>8.3f}"               # 47-54
                    f"{occupancy:>6.2f}"       # 55-60
                    f"{b_factor:>6.2f}"        # 61-66
                    f"          "              # 67-76 (blank)
                    f"{element:>2}"            # 77-78 (Element symbol)
                    f"  "                      # 79-80 (blank)
                )
                pdb_lines.append(pdb_line)
            except (ValueError, IndexError):
                continue
                
        elif atom_started and line.startswith("#") and len(pdb_lines) > 0:
            atom_started = False

    if not pdb_lines:
        print(f"Warning: No structural atoms could be parsed for {output_pdb_path}")
        return False

    with open(output_pdb_path, 'w') as f:
        for line in pdb_lines:
            f.write(line + "\n")
        f.write("END\n")
    
    print(f"Successfully converted structural data into: {output_pdb_path}")
    return True

def find_and_convert(cif_name, pdb_name):
    possible_paths = [
        cif_name,
        f"c/{cif_name}",
        f"../{cif_name}"
    ]
    target_path = None
    for path in possible_paths:
        if os.path.exists(path):
            target_path = path
            break
            
    if target_path:
        print(f"Found input structure file at: {target_path}")
        with open(target_path, "r") as cif_file:
            lines = cif_file.readlines()
        convert_cif_text_to_pdb(lines, pdb_name)
    else:
        print(f"Error: Could not find '{cif_name}' in your folders.")

if __name__ == "__main__":
    print("[*] Starting Batch CIF-to-PDB Conversion Process...")
    find_and_convert("structure-1-3.cif", "new_optimized.pdb")