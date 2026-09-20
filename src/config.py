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

# This input is user-editable, but must contain only the 20 canonical residues.
DEFAULT_SEED_SEQUENCE = "MKLAVFALLVFFAGSSDLIRR"

PILOT_SEED_SEQUENCE = "MKLAVFALLVFFAGSSDLIRR"

PEPTIDE_LENGTH = 21

AMINO_ACIDS = (
    "ACDEFGHIKLMNPQRSTVWY"
)

# ============================================================
# Evolutionary Algorithm
# ============================================================

POPULATION_SIZE = 10

GENERATIONS = 10

MUTATION_RATE = 0.05

ELITE_COUNT = 10

TOURNAMENT_SIZE = 7

ARISE_ENABLED = True

ALE_ENABLED = True

ADAPTIVE_ENABLED = True

RANDOM_SEED = 616

EXPERIMENT_CONDITION = "adaptive"

RUN_ID = ""

STRUCTURE_WORKERS = 1

DOCKING_WORKERS = 1

CACHE_ENABLED = True

# ============================================================
# ARISE / ALE hotspot experiments
# ============================================================

HOTSPOT_ENABLED = True

HOTSPOT_START = None

HOTSPOT_END = None

HOTSPOT_DISCOVERY_WINDOW_MIN = 2

HOTSPOT_DISCOVERY_WINDOW_MAX = 6

HOTSPOT_UPDATE_INTERVAL = 1

HOTSPOT_TOP_QUANTILE = 0.25

HOTSPOT_BOOTSTRAP_ITERATIONS = 100

MUTATION_GUIDANCE_MODE = "hotspot_biased"

PRESERVE_BEST = False

SCORE_THRESHOLD = None


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

    tournament_size: int = TOURNAMENT_SIZE

    docking_threshold: float = TIER2_DOCKING_THRESHOLD

    bootstrap_iterations: int = BOOTSTRAP_ITERATIONS

    random_seed: int = RANDOM_SEED

    top_k: int = TOP_K

    seed_sequence: str = DEFAULT_SEED_SEQUENCE

    adaptive_enabled: bool = ADAPTIVE_ENABLED

    arise_enabled: bool = ARISE_ENABLED

    ale_enabled: bool = ALE_ENABLED

    output_directory: Path = RESULTS_DIR

    run_id: str = RUN_ID

    condition: str = EXPERIMENT_CONDITION

    structure_workers: int = STRUCTURE_WORKERS

    docking_workers: int = DOCKING_WORKERS

    cache_enabled: bool = CACHE_ENABLED

    generate_plots: bool = GENERATE_PLOTS

    hotspot_enabled: bool = HOTSPOT_ENABLED

    hotspot_start: int | None = HOTSPOT_START

    hotspot_end: int | None = HOTSPOT_END

    hotspot_discovery_window_min: int = HOTSPOT_DISCOVERY_WINDOW_MIN

    hotspot_discovery_window_max: int = HOTSPOT_DISCOVERY_WINDOW_MAX

    hotspot_update_interval: int = HOTSPOT_UPDATE_INTERVAL

    hotspot_top_quantile: float = HOTSPOT_TOP_QUANTILE

    hotspot_bootstrap_iterations: int = HOTSPOT_BOOTSTRAP_ITERATIONS

    mutation_guidance_mode: str = MUTATION_GUIDANCE_MODE

    preserve_best: bool = PRESERVE_BEST

    score_threshold: float | None = SCORE_THRESHOLD

    def __post_init__(self) -> None:
        validate_configuration(self)


def validate_configuration(settings: PipelineSettings | None = None) -> None:
    """Validate runtime configuration before an experiment starts."""

    settings = settings or PipelineSettings.__new__(PipelineSettings)

    population_size = getattr(settings, "population_size", POPULATION_SIZE)
    generations = getattr(settings, "generations", GENERATIONS)
    mutation_rate = getattr(settings, "mutation_rate", MUTATION_RATE)
    elite_count = getattr(settings, "elite_count", ELITE_COUNT)
    tournament_size = getattr(settings, "tournament_size", TOURNAMENT_SIZE)
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
    seed_sequence = getattr(settings, "seed_sequence", DEFAULT_SEED_SEQUENCE)
    condition = getattr(settings, "condition", EXPERIMENT_CONDITION)
    structure_workers = getattr(settings, "structure_workers", STRUCTURE_WORKERS)
    docking_workers = getattr(settings, "docking_workers", DOCKING_WORKERS)
    hotspot_start = getattr(settings, "hotspot_start", HOTSPOT_START)
    hotspot_end = getattr(settings, "hotspot_end", HOTSPOT_END)
    window_min = getattr(
        settings,
        "hotspot_discovery_window_min",
        HOTSPOT_DISCOVERY_WINDOW_MIN,
    )
    window_max = getattr(
        settings,
        "hotspot_discovery_window_max",
        HOTSPOT_DISCOVERY_WINDOW_MAX,
    )
    hotspot_update_interval = getattr(
        settings,
        "hotspot_update_interval",
        HOTSPOT_UPDATE_INTERVAL,
    )
    hotspot_top_quantile = getattr(
        settings,
        "hotspot_top_quantile",
        HOTSPOT_TOP_QUANTILE,
    )
    hotspot_bootstrap_iterations = getattr(
        settings,
        "hotspot_bootstrap_iterations",
        HOTSPOT_BOOTSTRAP_ITERATIONS,
    )
    mutation_guidance_mode = getattr(
        settings,
        "mutation_guidance_mode",
        MUTATION_GUIDANCE_MODE,
    )
    score_threshold = getattr(settings, "score_threshold", SCORE_THRESHOLD)

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
    if tournament_size <= 0:
        raise ValueError("tournament_size must be greater than zero.")
    if bootstrap_iterations < 0:
        raise ValueError("bootstrap_iterations cannot be negative.")
    if not math.isfinite(float(tier2_threshold)):
        raise ValueError("tier2_threshold must be finite.")
    if PEPTIDE_LENGTH <= 0:
        raise ValueError("peptide_length must be greater than zero.")
    if len(AMINO_ACIDS) != 20 or len(set(AMINO_ACIDS)) != 20:
        raise ValueError("AMINO_ACIDS must contain 20 unique residues.")
    if len(seed_sequence) != PEPTIDE_LENGTH:
        raise ValueError("seed_sequence has the wrong length.")
    if any(residue not in AMINO_ACIDS for residue in seed_sequence):
        raise ValueError("seed_sequence contains an invalid residue.")
    if random_seed is not None and not isinstance(random_seed, int):
        raise ValueError("random_seed must be an integer or None.")
    if condition not in {"baseline", "adaptive", "pilot"}:
        raise ValueError(
            "condition must be one of: baseline, adaptive, pilot."
        )
    if structure_workers <= 0:
        raise ValueError("structure_workers must be greater than zero.")
    if docking_workers <= 0:
        raise ValueError("docking_workers must be greater than zero.")
    if hotspot_start is not None and not 0 <= hotspot_start < PEPTIDE_LENGTH:
        raise ValueError("hotspot_start must be within the sequence.")
    if hotspot_end is not None and not 0 < hotspot_end <= PEPTIDE_LENGTH:
        raise ValueError("hotspot_end must be within the sequence.")
    if hotspot_start is not None and hotspot_end is not None:
        if hotspot_end <= hotspot_start:
            raise ValueError("hotspot_end must be greater than hotspot_start.")
    if window_min < 2 or window_max < window_min:
        raise ValueError("Invalid hotspot discovery window range.")
    if window_max > PEPTIDE_LENGTH:
        raise ValueError("hotspot discovery window exceeds peptide length.")
    if hotspot_update_interval <= 0:
        raise ValueError("hotspot_update_interval must be greater than zero.")
    if not 0.0 < hotspot_top_quantile < 1.0:
        raise ValueError("hotspot_top_quantile must be between 0 and 1.")
    if hotspot_bootstrap_iterations < 0:
        raise ValueError("hotspot_bootstrap_iterations cannot be negative.")
    if mutation_guidance_mode not in {
        "uniform",
        "hotspot_only",
        "hotspot_biased",
    }:
        raise ValueError(
            "mutation_guidance_mode must be uniform, hotspot_only, "
            "or hotspot_biased."
        )
    if score_threshold is not None and not math.isfinite(float(score_threshold)):
        raise ValueError("score_threshold must be finite when configured.")
