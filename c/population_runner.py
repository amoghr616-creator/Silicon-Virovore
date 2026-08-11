import time
import json
from bridge import process_candidate_peptide, generate_c_population
from docking_vina import PeptideDockingScorer

def run_evolutionary_loop(seed_seq: str, generations: int = 3, pop_size: int = 5, mutation_rate: float = 0.20):
    scorer = PeptideDockingScorer(use_ml_surrogate=True)
    current_seed = seed_seq
    all_time_best = None

    print(f"\n==================================================")
    print(f" LAUNCHING GA EVOLUTION RUN ({generations} Generations, Pop Size {pop_size})")
    print(f" Initial Seed: {current_seed}")
    print(f"==================================================")

    start_time = time.time()

    for gen in range(1, generations + 1):
        print(f"\n>>> GENERATION {gen}/{generations} <<<")
        
        # 1. Generate mutated population in C memory
        candidates = generate_c_population(current_seed, pop_size=pop_size, mutation_rate=mutation_rate)
        
        gen_results = []
        for seq in candidates:
            c_data = process_candidate_peptide(seq)
            if not c_data:
                continue
            
            dock_res = scorer.evaluate_candidate(c_data)
            gen_results.append(dock_res)

        if not gen_results:
            print(" [WARNING] All candidates failed safety scan. Reverting to original seed.")
            continue

        # 2. Rank population by strongest (lowest) anchor Delta G
        gen_results.sort(key=lambda x: x["strongest_anchor_delta_g"])
        top_lead = gen_results[0]
        
        # 3. Update top seed for next generation
        current_seed = top_lead["sequence"]
        
        if all_time_best is None or top_lead["strongest_anchor_delta_g"] < all_time_best["strongest_anchor_delta_g"]:
            all_time_best = top_lead

        print(f"\n [GEN {gen} LEADER]")
        print(f"   Sequence: {top_lead['sequence']}")
        print(f"   Strongest Anchor ΔG: {top_lead['strongest_anchor_delta_g']} kcal/mol")

    elapsed = time.time() - start_time
    
    print(f"\n==================================================")
    print(f" EVOLUTION RUN COMPLETE ({elapsed:.2f}s)")
    print(f" All-Time Strongest Lead: {all_time_best['sequence']}")
    print(f" Best Anchor ΔG: {all_time_best['strongest_anchor_delta_g']} kcal/mol")
    print(f"==================================================")

    with open("evolution_leads.json", "w") as f:
        json.dump(all_time_best, f, indent=4)

if __name__ == "__main__":
    initial_seed = "ACDEFGHIKLMNPQRSTVWYA"
    run_evolutionary_loop(seed_seq=initial_seed, generations=3, pop_size=5, mutation_rate=0.20)