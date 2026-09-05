"""
bridge.py

Python interface to the Silicon Virovore native C backend.
"""

from pathlib import Path
import ctypes
import platform

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.models import Candidate
from src.config import (
    PEPTIDE_FRAGMENT_SIZE,
    POPULATION_SIZE,
    MUTATION_RATE,
)

# ==========================================================
# Locate Native Library
# ==========================================================

ROOT = Path(__file__).resolve().parent

SYSTEM = platform.system()

if SYSTEM == "Darwin":
    library_candidates = [
        ROOT / "libsafety.dylib",
        ROOT / "libsafety.so",
    ]
elif SYSTEM == "Linux":
    library_candidates = [
        ROOT / "libsafety.so",
    ]
else:
    raise RuntimeError(f"Unsupported operating system: {SYSTEM}")

LIB_PATH = next((p for p in library_candidates if p.exists()), None)

if LIB_PATH is None:
    raise FileNotFoundError(
        "Could not locate the Silicon Virovore native backend.\n"
        "Run `make` before executing the pipeline."
    )

lib = ctypes.CDLL(str(LIB_PATH))

print(f"[Bridge] Loaded native backend: {LIB_PATH.resolve()}")
# ==========================================================
# Constants
# ==========================================================

SEQ_LEN = 21
FRAGMENT_SIZE = PEPTIDE_FRAGMENT_SIZE

# ==========================================================
# C Function Signatures
# ==========================================================

lib.c_check_sequence_fitness.argtypes = [
    ctypes.c_char_p,
]
lib.c_check_sequence_fitness.restype = ctypes.c_double

PopulationArray = (
    (ctypes.c_char * (SEQ_LEN + 1))
)

lib.c_generate_mutated_population.restype = None

ImportanceArray = ctypes.c_double * SEQ_LEN

PopulationType = (
    (ctypes.c_char * (SEQ_LEN + 1))
    * POPULATION_SIZE
)

lib.c_generate_adaptive_population.argtypes = [
    ctypes.c_char_p,
    ImportanceArray,
    PopulationType,
    ctypes.c_int,
]

lib.c_generate_adaptive_population.restype = None
# ==========================================================
# Python Wrappers
# ==========================================================

def process_candidate_peptide(sequence: str) -> Candidate:
    """
    Evaluate a peptide using the native C backend.

    Returns
    -------
    Candidate
        Shared pipeline data model.
    """

    sequence = sequence.upper()

    if len(sequence) != SEQ_LEN:
        raise ValueError(
            f"Sequence must be exactly {SEQ_LEN} amino acids."
        )

    score = lib.c_check_sequence_fitness(
        sequence.encode("utf-8")
    )

    fragments = [
        sequence[i:i + FRAGMENT_SIZE]
        for i in range(len(sequence) - FRAGMENT_SIZE + 1)
    ]

    return Candidate(
        sequence=sequence,
        c_score=score,
        fragments=fragments,
    )

def generate_c_population(
    seed_sequence: str,
    pop_size: int = POPULATION_SIZE,
    mutation_rate: float = MUTATION_RATE,
) -> list[str]:
    """
    Generate a mutated peptide population using the native C engine.
    """

    seed_sequence = seed_sequence.upper()

    if len(seed_sequence) != SEQ_LEN:
        raise ValueError(
            f"Seed sequence must be {SEQ_LEN} amino acids."
        )

    PopulationType = (
        (ctypes.c_char * (SEQ_LEN + 1))
        * pop_size
    )

    population = PopulationType()

    lib.c_generate_mutated_population.argtypes = [
        ctypes.c_char_p,
        PopulationType,
        ctypes.c_int,
        ctypes.c_double,
    ]

    lib.c_generate_mutated_population(
        seed_sequence.encode("utf-8"),
        population,
        pop_size,
        mutation_rate,
    )

    return [
        population[i].value.decode("utf-8")
        for i in range(pop_size)
    ]

def generate_adaptive_population(
    seed_sequence: str,
    importance_map: list[float],
    pop_size: int = POPULATION_SIZE,
) -> list[str]:

    seed_sequence = seed_sequence.upper()

    PopulationType = (
        (ctypes.c_char * (SEQ_LEN + 1))
        * pop_size
    )

    population = PopulationType()

    ImportanceType = ctypes.c_double * SEQ_LEN

    importance = ImportanceType(*importance_map)

    lib.c_generate_adaptive_population(
        seed_sequence.encode("utf-8"),
        importance,
        population,
        pop_size,
    )

    return [
        population[i].value.decode("utf-8")
        for i in range(pop_size)
    ]