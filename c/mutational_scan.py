# mutational_scan.py
from .bridge import process_candidate_peptide
from src.docking_vina import PeptideDockingScorer

AA_POOL = "ACDEFGHIKLMNPQRSTVWY"
CHAMPION = "MKLAVFALLVFFAGSSDLIRR"

def calculate_composite_fitness(c_score: float, anchor_dg: float) -> float:
    if c_score <= 0.01:
        return 0.0
    return round(c_score * (1.0 + (abs(anchor_dg) * 0.2)), 4)

def run_mutational_scan(base_seq: str):
    print(f"[*] Running Exhaustive Mutational Scan for: {base_seq}")
    scorer = PeptideDockingScorer()
    
    # Evaluate baseline champion
    base_c = process_candidate_peptide(base_seq)
    base_dock = scorer.evaluate_candidate(base_c)
    base_fitness = calculate_composite_fitness(base_dock["c_score"], base_dock["strongest_anchor_delta_g"])
    
    print(f"[*] Champion Baseline Fitness: {base_fitness:.4f}\n")
    
    better_mutants = []
    total_scanned = 0

    for i in range(len(base_seq)):
        original_aa = base_seq[i]
        for aa in AA_POOL:
            if aa == original_aa:
                continue
            
            mutant_seq = base_seq[:i] + aa + base_seq[i+1:]
            c_data = process_candidate_peptide(mutant_seq)
            if not c_data:
                continue
                
            dock_res = scorer.evaluate_candidate(c_data)
            fitness = calculate_composite_fitness(dock_res["c_score"], dock_res["strongest_anchor_delta_g"])
            total_scanned += 1

            if fitness > base_fitness:
                better_mutants.append({
                    "pos": i + 1,
                    "swap": f"{original_aa}{i+1}{aa}",
                    "sequence": mutant_seq,
                    "fitness": fitness,
                    "delta": round(fitness - base_fitness, 4)
                })

    print(f"[DONE] Scanned {total_scanned} single-point mutants.")
    
    if not better_mutants:
        print("🏆 VERIFIED: Champion sequence is at a TRUE LOCAL OPTIMUM! No single mutation improved fitness.")
    else:
        print(f"⚠️ FOUND {len(better_mutants)} MUTANTS WITH HIGHER FITNESS:")
        better_mutants.sort(key=lambda x: x["fitness"], reverse=True)
        for m in better_mutants[:5]:
            print(f"   - Variant {m['swap']}: Fitness {m['fitness']} (+{m['delta']}) -> {m['sequence']}")

if __name__ == "__main__":
    run_mutational_scan(CHAMPION)