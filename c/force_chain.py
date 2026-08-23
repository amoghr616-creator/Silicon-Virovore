with open("lead_1_gromacs.pdb", "r") as f:
    lines = f.readlines()

fixed_lines = []
for line in lines:
    if line.startswith("ATOM") or line.startswith("HETATM"):
        # Force record type to ATOM (never HETATM for standard protein)
        record = "ATOM  "
        # Keep atom serial, name, altloc, resName, chainID, resSeq, coords, etc.
        # Standard PDB columns:
        # 0-6: ATOM
        # 6-11: serial
        # 12-16: atom name
        # 17-19: res name
        # 21: chain ID (force to 'A')
        # 22-26: res seq
        atom_part = line[6:21]
        res_seq_and_coords = line[22:]
        
        new_line = record + atom_part + "A" + res_seq_and_coords
        fixed_lines.append(new_line)
    else:
        fixed_lines.append(line)

with open("lead_1_fixed_chain.pdb", "w") as f:
    f.writelines(fixed_lines)

print("Created lead_1_fixed_chain.pdb with forced chain A assignments.")