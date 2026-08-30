import json
import time
from .bridge import generate_c_population, process_candidate_peptide
from src.docking_vina import PeptideDockingScorer

def calculate_composite_fitness(c_score: float, anchor_dg: float) -> float:
    if c_score <= 0.01:
        return 0.0
    return round(c_score * (1.0 + (abs(anchor_dg) * 0.2)), 4)

def run_seed_campaign(
    seed_name: str,
    seed_seq: str,
    pop_size: int = 50,
    generations: int = 100,
    initial_mutation_rate: float = 0.30
) -> dict:
    print(f"\n" + "=" * 66)
    print(f" 🧬 LAUNCHING CAMPAIGN: {seed_name}")
    print(f" Seed Sequence : {seed_seq}")
    print("=" * 66)

    scorer = PeptideDockingScorer()
    current_seed = seed_seq
    mutation_rate = initial_mutation_rate
    min_mutation_rate = 0.05
    mutation_decay = (initial_mutation_rate - min_mutation_rate) / generations

    all_time_champ = None
    all_time_score = -1.0
    history = []

    start_time = time.time()

    for gen in range(1, generations + 1):
        population_strings = generate_c_population(
            seed_seq=current_seed,
            pop_size=pop_size,
            mutation_rate=mutation_rate
        )

        evaluated_population = []
        for seq in population_strings:
            c_data = process_candidate_peptide(seq)
            if c_data is None:
                continue

            dock_res = scorer.evaluate_candidate(c_data)
            c_score = dock_res["c_score"]
            anchor_dg = dock_res["strongest_anchor_delta_g"]
            comp_score = calculate_composite_fitness(c_score, anchor_dg)

            dock_res["composite_fitness"] = comp_score
            evaluated_population.append(dock_res)

        if not evaluated_population:
            continue

        evaluated_population.sort(key=lambda x: x["composite_fitness"], reverse=True)
        gen_best = evaluated_population[0]

        if gen_best["composite_fitness"] > all_time_score:
            all_time_score = gen_best["composite_fitness"]
            all_time_champ = gen_best
            current_seed = gen_best["sequence"]

        history.append({
            "gen": gen,
            "best_seq": gen_best["sequence"],
            "comp_fit": gen_best["composite_fitness"],
            "anchor_dg": gen_best["strongest_anchor_delta_g"]
        })

        mutation_rate = max(min_mutation_rate, mutation_rate - mutation_decay)

    elapsed = round(time.time() - start_time, 2)
    print(f" [DONE] Campaign '{seed_name}' completed in {elapsed}s")
    print(f" Top Sequence  : {all_time_champ['sequence']}")
    print(f" Anchor Motif  : {all_time_champ['strongest_anchor']} ({all_time_champ['strongest_anchor_delta_g']} kcal/mol)")
    print(f" Composite Fit : {all_time_champ['composite_fitness']:.4f}")

    return {
        "seed_name": seed_name,
        "initial_seed": seed_seq,
        "champion": all_time_champ,
        "execution_time_s": elapsed
    }

def main():
    # 5 Diverse Seed Archetypes targeting different structural start points
    seeds = {
        "Native Env Target": "MKLAVDALLVTFAGSSDKKRR",
        "High-Hydrophobic Seed": "QCDPLGHDKLDNFQASTVWYA",
        "Alpha-Helical Seed": "EALEKALKEALEKALKEALEK",
        "Aromatic Stacking Seed": "YWWFARMWYAFWYAQASTVWY",
        "Charged Amphipathic Seed": "KRLRDALLVTFAKRLRDKKRR"
    }

    results = []
    for name, seq in seeds.items():
        res = run_seed_campaign(seed_name=name, seed_seq=seq, pop_size=50, generations=100)
        results.append(res)

    # Sort campaigns by final champion composite fitness
    results.sort(key=lambda x: x["champion"]["composite_fitness"], reverse=True)

    print("\n" + "=" * 80)
    print(" 🏆 MULTI-SEED SWEEP SUMMARY RANKINGS")
    print("=" * 80)
    print(f"{'Rank':<5} | {'Campaign Name':<24} | {'Champion Sequence':<22} | {'Anchor ΔG':<10} | {'Comp Fit':<8}")
    print("-" * 80)

    for idx, r in enumerate(results, 1):
        champ = r["champion"]
        print(
            f"{idx:<5} | {r['seed_name']:<24} | {champ['sequence']:<22} | "
            f"{champ['strongest_anchor_delta_g']:<10} | {champ['composite_fitness']:<8.3f}"
        )
    print("=" * 80)

    # Export combined results
    with open("multi_seed_sweep_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\n[EXPORT] Multi-seed sweep telemetry saved to multi_seed_sweep_results.json")

if __name__ == "__main__":
    main()