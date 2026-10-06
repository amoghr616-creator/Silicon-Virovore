import math
from collections import defaultdict

PDB = "../data/receptor/9SYA.pdb"
CHAIN = "A"
ANCHORS = {203, 205}
CUTOFF = 6.0

atoms = []
anchors = []

with open(PDB) as f:
    for line in f:
        if not line.startswith("ATOM"):
            continue

        chain = line[21].strip()
        if chain != CHAIN:
            continue

        try:
            resid = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue

        atom = {
            "resid": resid,
            "resname": line[17:20].strip(),
            "atom": line[12:16].strip(),
            "xyz": (x, y, z),
        }

        atoms.append(atom)

        if resid in ANCHORS:
            anchors.append(atom)

def dist(a, b):
    return math.sqrt(sum(
        (a[i] - b[i]) ** 2 for i in range(3)
    ))

nearby = defaultdict(lambda: {"resname": None, "min_dist": 999})

for atom in atoms:
    if atom["resid"] in ANCHORS:
        continue

    d = min(dist(atom["xyz"], x["xyz"]) for x in anchors)

    if d <= CUTOFF:
        r = atom["resid"]
        nearby[r]["resname"] = atom["resname"]
        nearby[r]["min_dist"] = min(nearby[r]["min_dist"], d)

print("Chain:", CHAIN)
print("Anchor residues: ARG203 and ARG205")
print("Cutoff:", CUTOFF, "A")
print()

for resid in sorted(nearby):
    print(
        f"{resid:4d} {nearby[resid]['resname']:3s} "
        f"{nearby[resid]['min_dist']:6.2f} A"
    )
