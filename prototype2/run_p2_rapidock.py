from pathlib import Path
import subprocess
import sys

CANDIDATES = {
    "P2-01": "MKQAVNALLVFFAGSSDAIRR",
    "P2-02": "DKLAVTALLVVFAESSDKLRR",
    "P2-03": "MKLAPFALLVVFAGASDWIRR",
    "P2-04": "MKLAVFALLVCFAESSDLVRR",
    "P2-05": "MKQAVFALLQVFAGSSDWIRR",
    "P2-06": "MKLAVFALRVFFAESSDAIRR",
}

ROOT = Path.cwd()
RECEPTOR = ROOT / "prototype2" / "9SYA_chainA_final.pdb"
OUTROOT = ROOT / "prototype2" / "p2_rapidock"

SITE = (-8.405, -5.823, 1.523)

if not RECEPTOR.exists():
    raise FileNotFoundError(RECEPTOR)

for name, sequence in CANDIDATES.items():

    outdir = OUTROOT / name
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "hybridock-pep", "dock",
        "--peptide", sequence,
        "--receptor", str(RECEPTOR),
        "--site", str(SITE[0]), str(SITE[1]), str(SITE[2]),
        "--box", "52",
        "--n-samples", "100",
        "--output-dir", str(outdir),
    ]

    print("\n" + "=" * 70)
    print(name, sequence)
    print("=" * 70)
    print(" ".join(cmd))

    result = subprocess.run(cmd)

    if result.returncode != 0:
        print(f"FAILED: {name}")
        continue

    print(f"COMPLETED: {name}")

print("\nFinished batch.")
