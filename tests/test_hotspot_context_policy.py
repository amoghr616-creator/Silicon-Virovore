import unittest

from src.ale_policy import resolve_mutation_policy
from src.context_policy import (
    SearchContext,
    decide_contextual_policy,
    run_contextual_step,
)
from src.hotspot import (
    _bootstrap_stability,
    _split_history,
    discover_hotspot_model,
)
from src.epistasis_entropy import EpistasisEntropyController
from src.population_runner import generate_candidates


def synthetic_records(run_count=4):
    records = []
    for run_index in range(run_count):
        for generation in range(2):
            for candidate_index in range(8):
                sequence = list("A" * 10)
                if candidate_index < 4:
                    sequence[3:6] = "KLM"
                    score = 1.0 + candidate_index * 0.01
                else:
                    sequence[3:6] = "CCC"
                    score = 0.1 + candidate_index * 0.01
                records.append({
                    "sequence": "".join(sequence),
                    "score": score,
                    "generation": generation,
                    "run_id": f"run-{run_index}",
                    "arise_score_comparable_across_generations": True,
                })
    return records


class HotspotContextPolicyTests(unittest.TestCase):
    def test_split_keeps_runs_whole_and_single_run_has_no_holdout(self):
        records = synthetic_records()
        discovery, holdout = _split_history(records, seed=5)
        discovery_runs = {record["run_id"] for record in discovery}
        holdout_runs = {record["run_id"] for record in holdout}
        self.assertTrue(discovery_runs.isdisjoint(holdout_runs))
        self.assertTrue(discovery)
        self.assertTrue(holdout)

        one_run = [record for record in records if record["run_id"] == "run-0"]
        single_discovery, single_holdout = _split_history(one_run, seed=5)
        self.assertEqual(single_discovery, one_run)
        self.assertEqual(single_holdout, [])
        self.assertEqual(
            _bootstrap_stability(
                one_run,
                range(0, 8),
                range(2, 5),
                0.25,
                20,
                5,
            ),
            {},
        )

    def test_permutation_gate_disables_ale_hotspot_on_failure(self):
        model = discover_hotspot_model(
            synthetic_records(),
            discovery_generation=1,
            run_id="run-0",
            bootstrap_iterations=5,
            permutation_iterations=1,
            require_significant_hotspot=True,
            random_seed=4,
        )
        self.assertEqual(model["status"], "significance_unvalidated")
        self.assertFalse(model["hotspot_significance_passed"])
        self.assertFalse(model["hotspot_significance_validated"])
        self.assertIsNone(model["selected_hotspot"])
        policy = resolve_mutation_policy(
            "hotspot_biased",
            model,
            enabled=True,
        )
        self.assertFalse(policy["guided"])
        self.assertEqual(policy["mode"], "uniform")

        unvalidated_candidate = resolve_mutation_policy(
            "hotspot_biased",
            {
                "selected_hotspot": {"start_position": 2, "end_position": 5},
                "hotspot_permutation_p_value": 0.001,
            },
            enabled=True,
        )
        self.assertFalse(unvalidated_candidate["guided"])
        self.assertEqual(
            unvalidated_candidate["guidance_disabled_reason"],
            "hotspot_significance_unvalidated",
        )

        unvalidated_context = SearchContext(
            hotspot_p_value=0.001,
            hotspot_stability=1.0,
            stagnation_generations=0,
            diversity=0.9,
        )
        decision = decide_contextual_policy(
            unvalidated_context,
            remaining_evaluations=100,
        )
        self.assertEqual(decision.action, "baseline")
        self.assertEqual(decision.reason, "hotspot_significance_unvalidated")

    def test_frozen_hotspot_passes_independent_holdout_validation(self):
        records = []
        for run_index in range(16):
            for generation in range(2):
                for candidate_index in range(8):
                    is_high = candidate_index in (2, 3)
                    sequence = list("A" * 10)
                    sequence[3:5] = "KL" if is_high else "CC"
                    score = (
                        1.0 + candidate_index * 0.01
                        if is_high
                        else 0.1 + candidate_index * 0.01
                    )
                    records.append({
                        "sequence": "".join(sequence),
                        "score": score,
                        "generation": generation,
                        "run_id": f"independent-{run_index}",
                        "arise_score_comparable_across_generations": True,
                    })

        model = discover_hotspot_model(
            records,
            discovery_generation=1,
            run_id="independent-0",
            top_quantile=0.25,
            window_min=2,
            window_max=4,
            bootstrap_iterations=50,
            permutation_iterations=199,
            permutation_alpha=0.05,
            require_significant_hotspot=True,
            random_seed=9,
        )

        self.assertTrue(model["hotspot_significance_validated"])
        self.assertTrue(model["hotspot_significance_passed"])
        self.assertIsNotNone(model["selected_hotspot"])
        self.assertGreaterEqual(model["holdout_run_count"], 3)
        self.assertGreaterEqual(model["holdout_sample_size"], 20)
        self.assertLess(model["holdout_permutation_p_value"], 0.05)
        self.assertGreaterEqual(
            model["selected_hotspot"]["bootstrap_stability"], 0.70
        )
        self.assertFalse(model["selection_uses_heldout"])


    def test_exploratory_pair_mode_uses_unvalidated_hotspot_explicitly(self):
        parent = "MKLAVFALLVFFAGSSDLIRR"
        model = {
            "selected_hotspot": {
                "start_position": 4,
                "end_position": 6,
                "window_length": 2,
                "hotspot_score": 0.14,
                "bootstrap_stability": 0.53,
            },
            "hotspot_significance_validated": False,
            "discovery_generation": 1,
        }
        candidates = generate_candidates(
            parent,
            population_size=16,
            mutation_rate=0.05,
            arise_enabled=True,
            cache_enabled=True,
            hotspot_model=model,
            hotspot_enabled=True,
            ale_enabled=True,
            mutation_guidance_mode=(
                "hotspot_entropy_substitution_epistasis_exploratory"
            ),
            generation=2,
            run_id="test-explicit-exploratory-hotspot",
            preserve_sequence=parent,
            total_generations=10,
            random_seed=123,
            operator_controller=EpistasisEntropyController(),
        )
        assert len(candidates) == 16
        assert all(
            c.metadata["hotspot_positions_used"] == [4, 5]
            for c in candidates
        )
        assert all(
            c.metadata["hotspot_guidance_evidence_status"]
            == "exploratory_unvalidated"
            for c in candidates
        )
        assert model["hotspot_significance_validated"] is False

    def test_context_selector_uses_evidence_and_budget(self):
        reliable = SearchContext(
            hotspot_p_value=0.01,
            hotspot_stability=0.9,
            stagnation_generations=0,
            diversity=0.7,
            hotspot_significance_validated=True,
        )
        adaptive = decide_contextual_policy(
            reliable,
            remaining_evaluations=100,
        )
        self.assertEqual(adaptive.action, "adaptive")

        weak = SearchContext(
            hotspot_p_value=None,
            hotspot_stability=0.9,
            stagnation_generations=10,
            diversity=0.1,
        )
        self.assertEqual(
            decide_contextual_policy(
                weak,
                remaining_evaluations=100,
            ).action,
            "baseline",
        )

        audit_context = SearchContext(
            hotspot_p_value=0.01,
            hotspot_stability=0.9,
            stagnation_generations=4,
            diversity=0.1,
            hotspot_significance_validated=True,
        )
        self.assertEqual(
            decide_contextual_policy(
                audit_context,
                remaining_evaluations=100,
                branch_budget=2,
                replications=2,
            ).action,
            "audit",
        )
        self.assertEqual(
            decide_contextual_policy(
                audit_context,
                remaining_evaluations=20,
                branch_budget=2,
                replications=2,
            ).action,
            "baseline",
        )

    def test_contextual_audit_charges_exact_cost(self):
        result = run_contextual_step(
            context=SearchContext(0.01, 0.9, 4, 0.1, True),
            start_state=(0, 0, 0),
            objective=lambda state, rng: float(sum(state)),
            preferred_position=0,
            remaining_evaluations=100,
            branch_budget=2,
            replications=2,
            seed=9,
        )
        self.assertEqual(result.action, "audit")
        self.assertEqual(result.evaluations_spent, 8)
        self.assertEqual(result.remaining_evaluations, 92)
        self.assertIsNotNone(result.audit_result)
        self.assertEqual(len(result.next_state), 3)


if __name__ == "__main__":
    unittest.main()