import ctypes
import os

# Load shared library
lib_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "libsafety.so"))
c_lib = ctypes.CDLL(lib_path)

# Configure argument & return types for safety check
c_lib.c_check_sequence_fitness.argtypes = [ctypes.c_char_p]
c_lib.c_check_sequence_fitness.restype = ctypes.c_double

# Configure types for dynamic C population generator
SEQ_LEN = 21
POP_SIZE = 10
PopulationArray = (ctypes.c_char * (SEQ_LEN + 1)) * POP_SIZE

c_lib.c_generate_mutated_population.argtypes = [
    ctypes.c_char_p,
    PopulationArray,
    ctypes.c_int,
    ctypes.c_double
]
c_lib.c_generate_mutated_population.restype = None

def generate_c_population(seed_seq: str, pop_size: int = 10, mutation_rate: float = 0.15) -> list[str]:
    """
    Calls the native C PRNG engine to mutate a seed candidate sequence directly in memory.
    """
    pop_buffer = PopulationArray()
    c_lib.c_generate_mutated_population(
        seed_seq.encode('utf-8'),
        pop_buffer,
        pop_size,
        mutation_rate
    )
    return [pop_buffer[i].value.decode('utf-8') for i in range(pop_size)]

def run_c_safety_scan(sequence: str, threshold: float = 0.0) -> tuple[bool, float]:
    score = c_lib.c_check_sequence_fitness(sequence.encode('utf-8'))
    return (score >= threshold), score

def slice_into_9mers(sequence: str, window_size: int = 9) -> list[str]:
    if len(sequence) < window_size:
        return []
    return [sequence[i:i + window_size] for i in range(len(sequence) - window_size + 1)]

def process_candidate_peptide(sequence: str) -> dict | None:
    is_safe, score = run_c_safety_scan(sequence)
    if not is_safe:
        return None
    return {
        "sequence": sequence,
        "c_score": score,
        "fragments": slice_into_9mers(sequence)
    }