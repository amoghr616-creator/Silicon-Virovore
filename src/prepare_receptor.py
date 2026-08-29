"""
prepare_receptor.py

Stage 1 of the Silicon Virovore pipeline.

Loads the HERV-K Env receptor structure,
performs basic cleaning,
and writes a docking-ready structure.

Future versions will integrate
PDBFixer, OpenMM, and AutoDockTools.
"""

from pathlib import Path
import logging

logger = logging.getLogger(__name__)

# ------------------------------------------------------------------
# Project Paths
# ------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent

INPUT = ROOT / "data" / "receptor" / "herv_k_env.pdb"

OUTPUT_DIR = ROOT / "results" / "receptor"

OUTPUT = OUTPUT_DIR / "clean_env.pdb"


# ------------------------------------------------------------------
# Cleaning Function
# ------------------------------------------------------------------

def clean_pdb(input_file: Path, output_file: Path):
    """
    Basic receptor cleanup.

    Current version:
    - Removes water molecules (HOH)

    Future versions:
    - Remove unwanted ligands
    - Add hydrogens
    - Repair missing atoms
    - Assign charges
    """

    with open(input_file) as fin, open(output_file, "w") as fout:
        for line in fin:
            if line.startswith("HETATM") and "HOH" in line:
                continue

            fout.write(line)


# ------------------------------------------------------------------
# Main Pipeline Function
# ------------------------------------------------------------------

def prepare_receptor():

    logger.info("Preparing receptor...")

    print(f"Input: {INPUT}")
    print(f"Exists: {INPUT.exists()}")

    if not INPUT.exists():
        raise FileNotFoundError(
            f"Could not find receptor:\n{INPUT}"
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    clean_pdb(INPUT, OUTPUT)

    logger.info(f"Saved cleaned receptor -> {OUTPUT}")

    return OUTPUT


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    prepare_receptor()