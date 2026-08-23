# verify_top3.py
import json
import os

def check_top_candidates(json_file="multi_seed_sweep_results.json"):
    if not os.path.exists(json_file):
        print(f"❌ Error: Could not find '{json_file}'. Make sure multi_seed_sweep.py has been run!")
        return

    with open(json_file, "r") as f:
        results = json.load(f)

    # Sort strictly by composite fitness descending
    sorted_by_fitness = sorted(results, key=lambda x: x["champion"]["composite_fitness"], reverse=True)

    print("=" * 80)
    print(" 📊 VERIFIED MULTI-SEED SWEEP RANKINGS (Parsed from multi_seed_sweep_results.json)")
    print("=" * 80)
    print(f"{'Rank':<5} | {'Campaign Name':<24} | {'Champion Sequence':<22} | {'Anchor ΔG':<10} | {'Comp Fit':<8}")
    print("-" * 80)

    for idx, r in enumerate(sorted_by_fitness, 1):
        champ = r["champion"]
        print(
            f"{idx:<5} | {r['seed_name']:<24} | {champ['sequence']:<22} | "
            f"{champ['strongest_anchor_delta_g']:<10} | {champ['composite_fitness']:<8.4f}"
        )
    print("=" * 80)

    print("\n🏆 True Top 3 Lead Panel:")
    for i in range(min(3, len(sorted_by_fitness))):
        champ = sorted_by_fitness[i]["champion"]
        print(f"  Lead #{i+1}: {champ['sequence']} (Fitness: {champ['composite_fitness']:.4f}, Anchor: {champ['strongest_anchor']} at {champ['strongest_anchor_delta_g']} kcal/mol)")

if __name__ == "__main__":
    check_top_candidates()