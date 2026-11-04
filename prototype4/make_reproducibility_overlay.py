from pathlib import Path
import tarfile

ROOT = Path.cwd()
OUT = ROOT / "prototype4/results/p4_visualization_20261103"
OUT.mkdir(parents=True, exist_ok=True)

runs = {
    "seed1917": ROOT / (
        "prototype4/results/p4_reliability_repeats_20261103/"
        "HADDOCK3_chain_A_seed1917"
    ),
    "seed2917": ROOT / (
        "prototype4/results/p4_reliability_repeats_20261103/"
        "HADDOCK3_chain_A_seed2917"
    ),
}

pdb_paths = {}

for label, run in runs.items():
    archive = (
        run / "analysis/6_seletopclusts_analysis/summary.tgz"
    )
    if not archive.is_file():
        raise FileNotFoundError(archive)

    with tarfile.open(archive, "r:*") as tar:
        matches = [
            m for m in tar.getmembers()
            if m.isfile()
            and Path(m.name).name == "cluster_1_model_1.pdb"
        ]

        if len(matches) != 1:
            raise RuntimeError(
                f"{label}: expected one rank-1 pose, found {len(matches)}"
            )

        handle = tar.extractfile(matches[0])
        if handle is None:
            raise RuntimeError(f"Could not extract {label}")

        output = OUT / f"{label}_rank1.pdb"
        output.write_bytes(handle.read())
        pdb_paths[label] = output

pml_path = OUT / "reproducibility_overlay.pml"
png_path = OUT / "reproducibility_overlay.png"
pse_path = OUT / "reproducibility_overlay.pse"

pml = f'''
reinitialize

load "{pdb_paths['seed1917']}", run1917
load "{pdb_paths['seed2917']}", run2917

# Align both receptor structures, retaining their peptide coordinates.
align run2917 and chain A and name CA, run1917 and chain A and name CA

hide everything

# Use one aligned receptor as the background.
show surface, run1917 and chain A
color gray80, run1917 and chain A
set transparency, 0.30, run1917 and chain A

# Contrasting peptide predictions.
show cartoon, run1917 and chain B
color cyan, run1917 and chain B

show cartoon, run2917 and chain B
color magenta, run2917 and chain B

# Target residues on the reference receptor.
select arg203, run1917 and chain A and resi 203 and resn ARG
select arg205, run1917 and chain A and resi 205 and resn ARG

show sticks, arg203 or arg205
color red, arg203
color orange, arg205
set stick_radius, 0.25, arg203 or arg205

label arg203 and name CZ, "Arg203"
label arg205 and name CZ, "Arg205"
set label_size, 25
set label_color, black
set label_outline_color, white

bg_color white
set orthoscopic, 1
set ray_shadows, 0
set antialias, 2
set two_sided_lighting, 1

orient (run1917 and chain A) or (run1917 and chain B) or (run2917 and chain B)
zoom (run1917 and chain A) or (run1917 and chain B) or (run2917 and chain B), 8

png "{png_path}", 2600, 1900, dpi=300, ray=1
save "{pse_path}"
'''

pml_path.write_text(pml, encoding="utf-8")

print("Extracted structures:")
for name, path in pdb_paths.items():
    print(name, path)
print("PyMOL script:", pml_path)
print("PNG:", png_path)
print("Session:", pse_path)
