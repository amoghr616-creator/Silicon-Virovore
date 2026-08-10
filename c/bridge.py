import ctypes
import os
import sys

# ==============================================================================
# 1. SHARED LIBRARY LOADING & C SYMBOL BINDING
# ==============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_NAME = "libsafety.so" if sys.platform != "win32" else "libsafety.dll"
LIB_PATH = os.path.join(SCRIPT_DIR, LIB_NAME)

if not os.path.exists(LIB_PATH):
    raise FileNotFoundError(
        f"Compiled library '{LIB_NAME}' not found in {SCRIPT_DIR}.\n"
        f"Run 'gcc -O3 -shared -fPIC -o libsafety.so *.c' in terminal first."
    )

c_lib = ctypes.CDLL(LIB_PATH)

# Bind to the safe C wrapper function
if hasattr(c_lib, "c_check_sequence_fitness"):
    c_eval_fitness = c_lib.c_check_sequence_fitness
    print("[C HANDSHAKE SUCCESS] Bound Python ctypes to safe C wrapper: 'c_check_sequence_fitness'")
elif hasattr(c_lib, "evaluate_fitness"):
    c_eval_fitness = c_lib.evaluate_fitness
    print("[C HANDSHAKE SUCCESS] Bound Python ctypes directly to C function: 'evaluate_fitness'")
else:
    raise AttributeError(
        f"Could not find 'c_check_sequence_fitness' in {LIB_NAME}.\n"
        f"Ensure you added the C wrapper function to one of your .c files and recompiled."
    )

# Set signature: double c_check_sequence_fitness(const char* sequence)
c_eval_fitness.argtypes = [ctypes.c_char_p]
c_eval_fitness.restype = ctypes.c_double


# ==============================================================================
# 2. C FITNESS & SAFETY SCAN WRAPPER
# ==============================================================================
def evaluate_c_fitness(sequence: str) -> float:
    """
    Passes sequence string directly to C memory buffer.
    Returns calculated fitness score without disk I/O.
    """
    seq_bytes = sequence.encode('utf-8')
    score = c_eval_fitness(seq_bytes)
    return float(score)


def run_c_safety_scan(sequence: str, threshold: float = 0.0) -> tuple[bool, float]:
    """
    Evaluates sequence against C safety/fitness threshold.
    """
    fitness_score = evaluate_c_fitness(sequence)
    is_safe = fitness_score >= threshold
    return is_safe, fitness_score


# ==============================================================================
# 3. TAB 4: 9-MER FRAGMENT SLICER (Resolves 21-mer Sampling Bottleneck)
# ==============================================================================
def slice_into_9mers(sequence: str, window_size: int = 9) -> list[str]:
    """
    Splits 21-mer sequence into overlapping rigid 9-mer anchor fragments
    to guide Rosetta FlexPepDock / AutoDock Vina reconstruction.
    """
    if len(sequence) < window_size:
        raise ValueError(f"Sequence length ({len(sequence)}) is shorter than window size ({window_size})")
    
    fragments = [sequence[i : i + window_size] for i in range(len(sequence) - window_size + 1)]
    return fragments


# ==============================================================================
# 4. MASTER ORCHESTRATOR LOOP
# ==============================================================================
def process_candidate_peptide(sequence_21mer: str) -> dict | None:
    print(f"\n==================================================")
    print(f" Processing Candidate 21-mer: {sequence_21mer}")
    print(f"==================================================")
    
    # Phase A: C Memory Fitness & Safety Evaluation
    is_safe, score = run_c_safety_scan(sequence_21mer, threshold=0.0)
    print(f" [C EVALUATION] Raw C Fitness Score: {score:.4f}")
    
    if not is_safe:
        print(" [REJECTED] Candidate failed C fitness/safety threshold.")
        return None
    
    print(" [PASSED] C Fitness & Safety Scan.")
    
    # Phase B: Tab 4 Fragment Generation for 3D Docking Reconstruction
    fragments = slice_into_9mers(sequence_21mer, window_size=9)
    print(f" [GENERATED] {len(fragments)} Overlapping 9-mer Anchor Fragments:")
    for idx, frag in enumerate(fragments):
        print(f"   Anchor Fragment {idx+1:02d}: {frag}")
        
    return {
        "sequence": sequence_21mer,
        "c_score": score,
        "fragments": fragments
    }


if __name__ == "__main__":
    # Standard 20 canonical amino acids + 1 valid amino acid = valid 21-mer
    test_sequence = "ACDEFGHIKLMNPQRSTVWYA"  # All valid single-letter AA codes
    result = process_candidate_peptide(test_sequence)