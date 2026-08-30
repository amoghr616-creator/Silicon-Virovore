"""
config.py

Central configuration for the Silicon Virovore pipeline.

All global parameters should be defined here instead of
being hardcoded throughout the codebase.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


# ============================================================
# Project Directories
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = ROOT / "data"

RESULTS_DIR = ROOT / "results"

FIGURES_DIR = RESULTS_DIR / "figures"

REPORT_DIR = RESULTS_DIR / "report"

STRUCTURE_DIR = RESULTS_DIR / "structures"

DOCKING_DIR = RESULTS_DIR / "docking"

RECEPTOR_DIR = DATA_DIR / "receptor"

PEPTIDE_DIR = DATA_DIR / "peptides"


# ============================================================
# Receptor
# ============================================================

RECEPTOR_PDB = RECEPTOR_DIR / "herv_k_env.pdb"

CLEAN_RECEPTOR = RESULTS_DIR / "receptor" / "clean_env.pdb"


# ============================================================
# Evolutionary Algorithm
# ============================================================

POPULATION_SIZE = 100

GENERATIONS = 25

MUTATION_RATE = 0.05

ELITE_COUNT = 10

RANDOM_SEED = 42


# ============================================================
# Docking
# ============================================================

PEPTIDE_FRAGMENT_SIZE = 9

TIER2_DOCKING_THRESHOLD = -7.0

USE_ML_SURROGATE = True


# ============================================================
# Ranking Weights
# ============================================================

FITNESS_WEIGHT = 0.35

DOCKING_WEIGHT = 0.35

HELIX_WEIGHT = 0.10

SOLVATION_WEIGHT = 0.10

CONFIDENCE_WEIGHT = 0.10


# ============================================================
# Bootstrap
# ============================================================

BOOTSTRAP_ITERATIONS = 100


# ============================================================
# Plotting
# ============================================================

FIGURE_DPI = 300

SAVE_PDF = True

SAVE_PNG = True


# ============================================================
# Logging
# ============================================================

LOG_LEVEL = "INFO"

LOG_FILE = RESULTS_DIR / "pipeline.log"


# ============================================================
# Pipeline
# ============================================================

TOP_K = 10

SAVE_INTERMEDIATE_FILES = True

VALIDATE_RESULTS = True

GENERATE_REPORT = True

GENERATE_PLOTS = True


# ============================================================
# Optional Future Features
# ============================================================

ENABLE_PARETO = True

ENABLE_BOOTSTRAP = True

ENABLE_DIVERSITY = True

ENABLE_CORRELATION = True

ENABLE_MD = False

ENABLE_EXPERIMENTAL_DOCKING = False


# ============================================================
# Runtime Configuration Object
# ============================================================

@dataclass(slots=True)
class PipelineSettings:

    population_size: int = POPULATION_SIZE

    generations: int = GENERATIONS

    mutation_rate: float = MUTATION_RATE

    docking_threshold: float = TIER2_DOCKING_THRESHOLD

    bootstrap_iterations: int = BOOTSTRAP_ITERATIONS

    random_seed: int = RANDOM_SEED

    top_k: int = TOP_K