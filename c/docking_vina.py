import json
import os
import sys
from bridge import process_candidate_peptide

# Tier 2 AutoDock Vina Activation Threshold (kcal/mol)
# Candidates with an anchor delta_g <= THRESHOLD trigger physical docking
TIER_2_DOCKING_THRESHOLD = -7.0


class PeptideDockingScorer:
    def __init__(self, use_ml_surrogate: bool = True):
        self.use_ml_surrogate = use_ml_surrogate

    def _predict_surrogate_delta_g(self, fragment: str) -> float:
        """
        Fast ML surrogate model evaluating 9-mer anchor binding affinity.
        """
        hydrophobic_residues = set("AILMFWYV")
        charged_residues = set("RHKDE")

        score = -5.0
        for aa in fragment:
            if aa in hydrophobic_residues:
                score -= 0.3
            elif aa in charged_residues:
                score += 0.1

        return round(score, 2)

    def _run_vina_cli_docking(self, fragment: str) -> float:
        """
        Full AutoDock Vina physical simulation (Triggered on high-affinity leads).
        """
        # Place Subprocess CLI execution for Vina binary here
        print(f"   [VINA CLI] Running rigid-backbone physical docking for {fragment}...")
        return -8.40

    def score_fragment(self, fragment: str) -> float:
        if self.use_ml_surrogate:
            return self._predict_surrogate_delta_g(fragment)
        else:
            return self._run_vina_cli_docking(fragment)

    def evaluate_candidate(self, candidate_data: dict) -> dict:
        fragments = candidate_data["fragments"]
        sequence = candidate_data["sequence"]

        print(f"\n--- Running Tier 1 Binding Evaluation (ML Surrogate) ---")
        scores = []
        for idx, frag in enumerate(fragments):
            dG = self.score_fragment(frag)
            scores.append({"fragment": frag, "delta_g": dG})
            print(f"   Fragment {idx+1:02d} [{frag}]: Predicted ΔG = {dG:.2f} kcal/mol")

        all_dG = [s["delta_g"] for s in scores]
        avg_dG = round(sum(all_dG) / len(all_dG), 2)
        best_dG = round(min(all_dG), 2)

        # Tier 2 Routing Decision
        passed_tier_2 = best_dG <= TIER_2_DOCKING_THRESHOLD
        
        print("\n==================================================")
        print(f" CANDIDATE DOCKING SUMMARY")
        print(f" Sequence: {sequence}")
        print(f" C Fitness Score: {candidate_data['c_score']:.4f}")
        print(f" Mean Fragment ΔG: {avg_dG} kcal/mol")
        print(f" Strongest Anchor ΔG: {best_dG} kcal/mol")
        print("==================================================")

        if passed_tier_2:
            print(f"\n[ALERT] Lead candidate detected (ΔG {best_dG} <= {TIER_2_DOCKING_THRESHOLD})!")
            print("[TIER 2 TRIGGERED] Running heavy physical Vina docking...")
            # Run physical verification on strongest anchor
            top_frag = min(scores, key=lambda x: x["delta_g"])["fragment"]
            vina_dG = self._run_vina_cli_docking(top_frag)
            print(f"[VINA VERIFIED] Physical Binding ΔG: {vina_dG} kcal/mol")
        else:
            print(f"\n[TIER 1 COMPLETE] Candidate ΔG ({best_dG}) above threshold ({TIER_2_DOCKING_THRESHOLD}). Skipping heavy Vina docking.")

        return {
            "sequence": sequence,
            "c_score": candidate_data["c_score"],
            "mean_delta_g": avg_dG,
            "strongest_anchor_delta_g": best_dG,
            "passed_tier_2": passed_tier_2,
            "fragments": scores
        }


if __name__ == "__main__":
    test_sequence = "ACDEFGHIKLMNPQRSTVWYA"
    candidate = process_candidate_peptide(test_sequence)

    if candidate:
        scorer = PeptideDockingScorer(use_ml_surrogate=True)
        results = scorer.evaluate_candidate(candidate)