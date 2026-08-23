cat << 'EOF' > fix_names.py
three_letter_map = {
    'M': 'MET', 'K': 'LYS', 'A': 'ALA', 'V': 'VAL', 'F': 'PHE',
    'L': 'LEU', 'I': 'ILE', 'G': 'GLY', 'S': 'SER', 'D': 'ASP',
    'R': 'ARG', 'N': 'ASN', 'C': 'CYS', 'E': 'GLU', 'Q': 'GLN',
    'H': 'HIS', 'P': 'PRO', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR'
}

with open("lead_1_fixed_chain.pdb", "r") as f:
    lines = f.readlines()

new_lines = []
for line in lines:
    if line.startswith("ATOM") or line.startswith("HETATM"):
        res_name = line[17:20].strip()
        if res_name in three_letter_map:
            correct_res = three_letter_map[res_name].ljust(3)
            line = line[:17] + correct_res + line[20:]
    new_lines.append(line)

with open("lead_1_3letter.pdb", "w") as f:
    f.writelines(new_lines)

print("Created lead_1_3letter.pdb successfully!")
EOF