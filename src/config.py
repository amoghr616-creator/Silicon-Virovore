"""
config.py

Central configuration for the Silicon Virovore pipeline.

All global parameters should be defined here instead of
being hardcoded throughout the codebase.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
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

RECEPTOR_PDBQT = RECEPTOR_DIR / "receptor.pdbqt"

CLEAN_RECEPTOR = RESULTS_DIR / "receptor" / "clean_env.pdb"

# ============================================================
# Initial Peptide
# ============================================================

DEFAULT_SEED_SEQUENCE = "MKLAVFALLVFFAGSSDLIRR"

PEPTIDE_LENGTH = 21

AMINO_ACIDS = (
    "ACDEFGHIKLMNPQRSTVWY"
)

# ============================================================
# Evolutionary Algorithm
# ============================================================

POPULATION_SIZE = 15

GENERATIONS = 5

MUTATION_RATE = 0.05

ELITE_COUNT = 10

ARISE_ENABLED = True

RANDOM_SEED = 616


# ============================================================
# Docking
# ============================================================

ESMFOLD_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

ESMFOLD_TIMEOUT_SECONDS = 45

ESMFOLD_MAX_ATTEMPTS = 3

ESMFOLD_RETRY_BACKOFF_SECONDS = 2

ESMFOLD_MAX_RETRY_DELAY_SECONDS = 8

PEPTIDE_FRAGMENT_SIZE = 9

TIER2_DOCKING_THRESHOLD = -6.0

VINA_EXECUTABLE = ROOT / "vina_1.2.7_mac_aarch64"

VINA_EXHAUSTIVENESS = 1

VINA_NUM_MODES = 1

VINA_MAX_EVALS = 5000

VINA_CPU = 1

VINA_BOX_PADDING = 4.0

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

REMOVE_DUPLICATE_CANDIDATES = True

DIVERSITY_REFILL_ATTEMPTS = 5

ARISE_MIN_OBSERVATIONS = 5

SAVE_INTERMEDIATE_FILES = True

VALIDATE_RESULTS = True

GENERATE_REPORT = True

GENERATE_PLOTS = True

seed_sequence: str = DEFAULT_SEED_SEQUENCE

peptide_length: int = PEPTIDE_LENGTH

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

    elite_count: int = ELITE_COUNT

    docking_threshold: float = TIER2_DOCKING_THRESHOLD

    bootstrap_iterations: int = BOOTSTRAP_ITERATIONS

    random_seed: int = RANDOM_SEED

    top_k: int = TOP_K

    def __post_init__(self) -> None:
        validate_configuration(self)


def validate_configuration(settings: PipelineSettings | None = None) -> None:
    """Validate runtime configuration before an experiment starts."""

    settings = settings or PipelineSettings.__new__(PipelineSettings)

    population_size = getattr(settings, "population_size", POPULATION_SIZE)
    generations = getattr(settings, "generations", GENERATIONS)
    mutation_rate = getattr(settings, "mutation_rate", MUTATION_RATE)
    elite_count = getattr(settings, "elite_count", ELITE_COUNT)
    bootstrap_iterations = getattr(
        settings,
        "bootstrap_iterations",
        BOOTSTRAP_ITERATIONS,
    )
    tier2_threshold = getattr(
        settings,
        "docking_threshold",
        TIER2_DOCKING_THRESHOLD,
    )
    random_seed = getattr(settings, "random_seed", RANDOM_SEED)

    if population_size <= 0:
        raise ValueError("population_size must be greater than zero.")
    if generations <= 0:
        raise ValueError("generations must be greater than zero.")
    if not 0.0 <= mutation_rate <= 1.0:
        raise ValueError("mutation_rate must be between 0 and 1.")
    if elite_count < 0 or elite_count > population_size:
        raise ValueError(
            "elite_count must be between zero and population_size."
        )
    if bootstrap_iterations < 0:
        raise ValueError("bootstrap_iterations cannot be negative.")
    if not math.isfinite(float(tier2_threshold)):
        raise ValueError("tier2_threshold must be finite.")
    if PEPTIDE_LENGTH <= 0:
        raise ValueError("peptide_length must be greater than zero.")
    if len(AMINO_ACIDS) != 20 or len(set(AMINO_ACIDS)) != 20:
        raise ValueError("AMINO_ACIDS must contain 20 unique residues.")
    if any(residue not in AMINO_ACIDS for residue in DEFAULT_SEED_SEQUENCE):
        raise ValueError("DEFAULT_SEED_SEQUENCE contains an invalid residue.")
    if len(DEFAULT_SEED_SEQUENCE) != PEPTIDE_LENGTH:
        raise ValueError("DEFAULT_SEED_SEQUENCE has the wrong length.")
    if random_seed is not None and not isinstance(random_seed, int):
        raise ValueError("random_seed must be an integer or None.")
