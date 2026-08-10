def fix_pdb_spacing(input_path, output_path):
    fixed_lines = []
    
    with open(input_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # We split the line by whitespace to extract the clean values safely
                parts = line.split()
                if len(parts) < 11:
                    fixed_lines.append(line)
                    continue
                
                try:
                    record = parts[0]       # ATOM or HETATM
                    serial = int(parts[1])   # Atom ID
                    name = parts[2]         # Atom Name (e.g., N, CA)
                    res_name = parts[3]     # Residue Name (e.g., MET)
                    chain = parts[4]        # Chain ID (e.g., A)
                    res_seq = int(parts[5])  # Residue Sequence Number
                    
                    x = float(parts[6])     # X Coord
                    y = float(parts[7])     # Y Coord
                    z = float(parts[8])     # Z Coord
                    occ = float(parts[9])    # Occupancy
                    temp = float(parts[10])  # B-factor
                    
                    element = parts[11] if len(parts) > 11 else name[0]
                    
                    # Correct spacing for 1-2 letter atom names
                    if len(name) < 4:
                        formatted_name = f" {name:<3}"
                    else:
                        formatted_name = f"{name:<4}"
                    
                    # Force strict, standard PDB format column alignments
                    fixed_line = (
                        f"{record:<6}"          # 1-6
                        f"{serial:>5d}"          # 7-11
                        f" "                     # 12
                        f"{formatted_name}"      # 13-16
                        f" "                     # 17 (Alternate location indicator)
                        f"{res_name:>3}"         # 18-20 (Residue Name)
                        f"{chain:1}"             # 21 (Chain ID)
                        f"{res_seq:>4d}"         # 22-25 (Residue Sequence)
                        f" "                     # 26 (Insertion code)
                        f"   "                   # 27-30
                        f"{x:>8.3f}"             # 31-38
                        f"{y:>8.3f}"             # 39-46
                        f"{z:>8.3f}"             # 47-54
                        f"{occ:>6.2f}"           # 55-60
                        f"{temp:>6.2f}"          # 61-66
                        f"          "            # 67-76
                        f"{element:>2}"          # 77-78
                        f"  "                    # 79-80
                    )
                    fixed_lines.append(fixed_line + "\n")
                except ValueError:
                    fixed_lines.append(line)
            else:
                fixed_lines.append(line)
                
    with open(output_path, 'w') as f:
        f.writelines(fixed_lines)
    print(f"File successfully fixed and saved to: {output_path}")

if __name__ == "__main__":
    fix_pdb_spacing("new_optimized.pdb", "new_optimized_fixed.pdb")