"""
bridge.py

Python interface to the Silicon Virovore native C backend.
"""

from pathlib import Path
import ctypes
import threading
import platform
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
_C_ENGINE_CACHE_VERSION = "engine-v1"
_fitness_cache: dict[tuple[str, int, str], float] = {}
_fitness_cache_hits = 0
_fitness_cache_misses = 0
_fitness_cache_lock = threading.Lock()

# ==========================================================
# C Function Signatures
# ==========================================================

lib.c_check_sequence_fitness.argtypes = [
    ctypes.c_char_p,
]
lib.c_check_sequence_fitness.restype = ctypes.c_double

lib.c_seed_random.argtypes = [ctypes.c_uint]
lib.c_seed_random.restype = None

PopulationRow = ctypes.c_char * (SEQ_LEN + 1)
PopulationPointer = ctypes.POINTER(PopulationRow)

lib.c_generate_mutated_population.restype = None

ImportanceArray = ctypes.c_double * SEQ_LEN

lib.c_generate_adaptive_population.argtypes = [
    ctypes.c_char_p,
    ImportanceArray,
    PopulationPointer,
    ctypes.c_int,
]

lib.c_generate_adaptive_population.restype = None

lib.c_generate_policy_population.argtypes = [
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    PopulationPointer,
    ctypes.c_int,
]
lib.c_generate_policy_population.restype = None
# ==========================================================
# Python Wrappers
# ==========================================================

def seed_native_random(seed: int) -> None:
    """Seed the native C RNG used by both population generators."""

    lib.c_seed_random(ctypes.c_uint(seed))


def clear_evaluation_cache() -> None:
    """Clear deterministic native-evaluation cache and counters."""

    global _fitness_cache_hits, _fitness_cache_misses
    with _fitness_cache_lock:
        _fitness_cache.clear()
        _fitness_cache_hits = 0
        _fitness_cache_misses = 0


def evaluation_cache_stats() -> dict[str, int]:
    with _fitness_cache_lock:
        return {
            "hits": _fitness_cache_hits,
            "misses": _fitness_cache_misses,
            "entries": len(_fitness_cache),
        }


def process_candidate_peptide(
    sequence: str,
    cache_enabled: bool = True,
) -> Candidate:
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

    global _fitness_cache_hits, _fitness_cache_misses
    cache_key = (sequence, FRAGMENT_SIZE, _C_ENGINE_CACHE_VERSION)
    cache_hit = False
    if cache_enabled:
        with _fitness_cache_lock:
            score = _fitness_cache.get(cache_key)
        if score is not None:
            cache_hit = True
            _fitness_cache_hits += 1
        else:
            _fitness_cache_misses += 1
            score = lib.c_check_sequence_fitness(
                sequence.encode("utf-8")
            )
            with _fitness_cache_lock:
                _fitness_cache[cache_key] = score
    else:
        _fitness_cache_misses += 1
        score = lib.c_check_sequence_fitness(
            sequence.encode("utf-8")
        )
    '''print(
        f"[C FITNESS] {sequence} -> {score:.6f}"
    )'''
    fragments = [
        sequence[i:i + FRAGMENT_SIZE]
        for i in range(len(sequence) - FRAGMENT_SIZE + 1)
    ]

    candidate = Candidate(
        sequence=sequence,
        c_score=score,
        fragments=fragments,
    )
    candidate.metadata["c_fitness_cache_hit"] = cache_hit
    return candidate

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
        PopulationRow
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
        PopulationRow
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


def generate_policy_population(
    seed_sequence: str,
    mutation_rate: float,
    guidance_mode: str = "uniform",
    hotspot_start: int | None = None,
    hotspot_end: int | None = None,
    pop_size: int = POPULATION_SIZE,
) -> list[str]:
    """Generate baseline or frozen-hotspot-guided candidates.

    The native implementation gives every candidate the same rounded number
    of mutation events. This is separate from the legacy adaptive generator.
    """

    seed_sequence = seed_sequence.upper()
    if len(seed_sequence) != SEQ_LEN:
        raise ValueError(f"Seed sequence must be {SEQ_LEN} amino acids.")

    mode_values = {
        "uniform": 0,
        "hotspot_only": 1,
        "hotspot_biased": 2,
    }
    if guidance_mode not in mode_values:
        raise ValueError("Unknown mutation guidance mode.")

    PopulationType = PopulationRow * pop_size
    population = PopulationType()
    lib.c_generate_policy_population(
        seed_sequence.encode("utf-8"),
        mutation_rate,
        -1 if hotspot_start is None else int(hotspot_start),
        -1 if hotspot_end is None else int(hotspot_end),
        mode_values[guidance_mode],
        population,
        pop_size,
    )
    return [population[i].value.decode("utf-8") for i in range(pop_size)]
