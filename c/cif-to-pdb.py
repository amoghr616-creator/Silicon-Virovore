import sys
import os

def convert_cif_text_to_pdb(cif_lines, output_pdb_path):
    pdb_lines = []
    atom_started = False
    
    # Dynamically find header indices to prevent misalignment
    headers = {}
    header_idx = 0
    
    for line in cif_lines:
        line = line.strip()
        
        # Track loop headers
        if line.startswith("_atom_site."):
            atom_started = True
            headers[line] = header_idx
            header_idx += 1
            continue
            
        # Parse the actual data rows
        if atom_started and (line.startswith("ATOM") or line.startswith("HETATM")):
            parts = line.split()
            if len(parts) < 10:
                continue
                
            try:
                # Dynamically match columns if available, otherwise use strict fallbacks
                group_pdb = parts[0]
                atom_id = int(parts[1])
                element = parts[2]
                atom_name = parts[3]
                res_name = parts[4]
                seq_id = int(parts[5])
                chain_id = parts[8]
                
                # Match typical coordinates from standard ESMFold/AlphaFold CIF outputs
                x = float(parts[10])
                y = float(parts[11])
                z = float(parts[12])
                
                occupancy = float(parts[13]) if len(parts) > 13 and parts[13] != '?' else 1.00
                b_factor = float(parts[14]) if len(parts) > 14 and parts[14] != '?' else 90.00
                
                # Format exactly to standard PDB fixed-width specifications
                pdb_line = (
                    f"{group_pdb:<6}"          
                    f"{atom_id:>5} "           
                    f"{atom_name:<4}"          
                    f"{res_name:>3} "          
                    f"{chain_id}{seq_id:>4}    " 
                    f"{x:>8.3f}"               
                    f"{y:>8.3f}"               
                    f"{z:>8.3f}"               
                    f"{occupancy:>6.2f}"       
                    f"{b_factor:>6.2f}"        
                    f"          {element:>2}  " 
                )
                pdb_lines.append(pdb_line)
            except (ValueError, IndexError) as e:
                continue
                
        # Stop parsing when exiting the loop block
        elif atom_started and line.startswith("#") and len(pdb_lines) > 0:
            atom_started = False

    if not pdb_lines:
        print("Warning: No structural atoms could be parsed. Check your CIF format.")
        return False

    # Write out the completed PDB file
    with open(output_pdb_path, 'w') as f:
        for line in pdb_lines:
            f.write(line + "\n")
        f.write("END\n")
    
    print(f"Successfully converted complete structural data into: {output_pdb_path}")
    return True

if __name__ == "__main__":
    # Define possible search paths for your structural file in VS Code
    possible_paths = [
        "structure-1-2.cif",
        "c/structure-1-2.cif",
        "../structure-1-2.cif"
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
        convert_cif_text_to_pdb(lines, "validated_helix.pdb")
    else:
        print("Error: Could not find 'structure-1-2.cif' in your project folders.")
        print("Please ensure the file is named exactly structure-1-2.cif and sits in your directory.")