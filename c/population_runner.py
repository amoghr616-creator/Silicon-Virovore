import time
import json
from bridge import process_candidate_peptide
from docking_vina import PeptideDockingScorer

def run_population_screening(candidate_pool: list[str], output_file: str = "lead_candidates.json") -> list[dict]:
    """
    Streams a population of candidate 21-mers through the full pipeline:
      1. C High-Speed Safety Gatekeeper
      2. 9-mer Anchor Fragment Slicing
      3. Tier 1 ML Surrogate Docking Evaluation
      4. Tier 2 Vina Activation (if threshold hit)
    """
    scorer = PeptideDockingScorer(use_ml_surrogate=True)
    lead_candidates = []

    print(f"\n==================================================")
    print(f" STARTING POPULATION SCREENING ({len(candidate_pool)} Candidates)")
    print(f"==================================================")

    start_time = time.time()

    for idx, seq in enumerate(candidate_pool):
        print(f"\n--- Processing Candidate [{idx+1}/{len(candidate_pool)}]: {seq} ---")
        
        # Run C scan & 9-mer slicing via bridge
        candidate_data = process_candidate_peptide(seq)
        
        if not candidate_data:
            print(" [REJECTED] Failed C safety/validation scan.")
            continue

        # Evaluate docking scores
        results = scorer.evaluate_candidate(candidate_data)
        lead_candidates.append(results)

    elapsed = time.time() - start_time
    print(f"\n==================================================")
    print(f" SCREENING COMPLETE")
    print(f" Processed: {len(candidate_pool)} candidates in {elapsed:.2f}s")
    print(f" Qualified Leads: {len(lead_candidates)}")
    print(f"==================================================")

    # Save output log
    with open(output_file, "w") as f:
        json.dump(lead_candidates, f, indent=4)
    print(f"[EXPORTED] Saved results to '{output_file}'")

    return lead_candidates


if __name__ == "__main__":
    # Test batch of diverse candidate 21-mers
    sample_population = [
        "ACDEFGHIKLMNPQRSTVWYA",
        "KLMNPQRSTVWACDEFGHIKL",
        "YVWTSRQPONMLKJIHGFEDC",
        "ACDEFGHIKLMNPQRSTVWY1"   # Triggers invalid AA rejection safely
    ]

    run_population_screening(sample_population)