import unittest

from scripts.benchmark_combined_adaptive_search import run_benchmark ; 
from scripts.benchmark_combined_adaptive_search import make_landscape


from src.p2_adaptive_search import (
    CombinedAdaptiveController,
    IndependentSignalValidator,
    MatchedPolicyTrialRunner,
    SignalObservation,
)


def signal_records(split, groups, values):
    records = []
    order = 0
    for group in groups:
        for signal, outcome in values:
            order += 1
            records.append(SignalObservation(
                observation_id=f"{split}-{group}-{signal}",
                independent_group_id=group,
                run_id=f"run-{group}",
                lineage_id=f"lineage-{group}",
                context_id="synthetic",
                signal=signal,
                outcome=float(outcome),
                split=split,
                observation_order=order,
            ))
    return records


class CombinedAdaptiveSearchTests(unittest.TestCase):
    def test_related_groups_cannot_cross_discovery_and_validation(self):
        discovery = signal_records("discovery", ["g1", "g2", "g3", "g4"], [("position:0", 1), ("position:1", 0)])
        validation = signal_records("validation", ["g4", "g5", "g6", "g7"], [("position:0", 1), ("position:1", 0)])
        validator = IndependentSignalValidator(permutations=20)
        with self.assertRaises(ValueError):
            validator.fit_and_validate(discovery, validation)
    def test_deceptive_landscapes_vary_and_have_known_optima(self):
        landscapes = [
            make_landscape("deceptive", instance_id)
            for instance_id in range(20)
        ]

        # Different instances should contain different trap structures
        # and different global target configurations.
        self.assertGreater(
            len({landscape.deceptive_blocks for landscape in landscapes}),
            1,
        )
        self.assertGreater(
            len({landscape.deceptive_target for landscape in landscapes}),
            1,
        )

        for landscape in landscapes:
            target = landscape.deceptive_target
            optimum = landscape.optimum()

            self.assertAlmostEqual(
                landscape.noiseless_value(target),
                optimum,
            )

            trap = tuple(1 - bit for bit in target)
            trap_score = landscape.noiseless_value(trap)

            self.assertLess(trap_score, optimum)

            # The complement is a local trap: any one-bit move away
            # from it initially reduces the objective.
            for position in range(landscape.feature_count):
                neighbor = list(trap)
                neighbor[position] = 1 - neighbor[position]

                self.assertLess(
                    landscape.noiseless_value(tuple(neighbor)),
                    trap_score,
                )
    def test_max_statistic_correction_is_deterministic_and_can_pass_known_signal(self):
        values = [("position:0", 1), ("position:1", 0), ("position:2", 0)]
        discovery = signal_records("discovery", [f"d{i}" for i in range(6)], values)
        validation = signal_records("validation", [f"v{i}" for i in range(6)], values)
        first = IndependentSignalValidator(permutations=100, seed=9).fit_and_validate(discovery, validation)
        second = IndependentSignalValidator(permutations=100, seed=9).fit_and_validate(discovery, validation)
        self.assertEqual(first, second)
        self.assertEqual(first.selected_signal, "position:0")
        self.assertTrue(first.passes)
        self.assertGreaterEqual(first.adjusted_p_values["position:0"], 0.0)

    def test_matched_trials_allow_reversible_mutations(self):
        runner = MatchedPolicyTrialRunner(branch_budget=10)
        observed_states = []

        def objective(state, rng):
            observed_states.append(tuple(state))
            return float(sum(state))

        result = runner.run_trial(
            trial_id="reversible-test",
            independent_group_id="independent-instance-1",
            seed=19,
            start_state=(0, 0, 0),
            objective=objective,
            preferred_position=0,
        )

        self.assertEqual(result.baseline_evaluations, 10)
        self.assertEqual(result.adaptive_evaluations, 10)

        # The first 10 observations are the baseline branch.
        baseline_states = observed_states[:10]
        self.assertEqual(len(baseline_states), 10)

        # At least one transition must reduce the number of active bits.
        active_counts = [sum(state) for state in baseline_states]
        self.assertTrue(
            any(
                after < before
                for before, after in zip(
                    active_counts, active_counts[1:]
                )
            ),
            "baseline search should be able to reverse a bit flip",
        )

    def test_insufficient_signal_evidence_fails_closed(self):
        discovery = signal_records("discovery", ["d1", "d2"], [("position:0", 1), ("position:1", 0)])
        validation = signal_records("validation", ["v1", "v2"], [("position:0", 1), ("position:1", 0)])
        result = IndependentSignalValidator(permutations=10).fit_and_validate(discovery, validation)
        self.assertFalse(result.passes)
        self.assertIn("insufficient_discovery_groups", result.reason)

    def test_matched_trials_start_identically_and_use_equal_budgets(self):
        runner = MatchedPolicyTrialRunner(branch_budget=7)
        result = runner.run_trial(
            trial_id="t1",
            independent_group_id="instance-1",
            seed=4,
            start_state=(0, 0, 0, 0),
            objective=lambda state, rng: float(sum(state)),
            preferred_position=0,
        )
        self.assertEqual(result.start_state, (0, 0, 0, 0))
        self.assertEqual(result.baseline_evaluations, result.adaptive_evaluations)
        self.assertEqual(result.baseline_evaluations, 7)

    def test_signal_success_does_not_imply_policy_success(self):
        values = [("position:0", 1), ("position:1", 0)]
        signal = IndependentSignalValidator(permutations=20).fit_and_validate(
            signal_records("discovery", [f"d{i}" for i in range(20)], values),
            signal_records("validation", [f"v{i}" for i in range(20)], values),
        )
        runner = MatchedPolicyTrialRunner(branch_budget=4)
        failed_policy = runner.compare([], minimum_trials=2)
        controller = CombinedAdaptiveController(signal_result=signal, policy_result=failed_policy)
        decision = controller.decide(checkpoint="evaluation:1", mode="combined")
        self.assertTrue(signal.passes)
        self.assertFalse(failed_policy.passes)
        self.assertFalse(decision.adaptive_selected)

    def test_validation_failure_blocks_combined_signal_use(self):
        values = [("position:0", 1), ("position:1", 0)]
        validator = IndependentSignalValidator(minimum_effect=2.0, permutations=20)
        signal = validator.fit_and_validate(
            signal_records("discovery", [f"d{i}" for i in range(5)], values),
            signal_records("validation", [f"v{i}" for i in range(5)], values),
        )
        runner = MatchedPolicyTrialRunner(branch_budget=4)
        trials = [runner.run_trial(
            trial_id=f"t{i}", independent_group_id=f"g{i}", seed=i,
            start_state=(0, 0), objective=lambda state, rng: float(sum(state)), preferred_position=0,
        ) for i in range(4)]
        policy = runner.compare(trials, minimum_trials=2, minimum_improvement=-1.0, required_win_rate=0.0)
        decision = CombinedAdaptiveController(signal_result=signal, policy_result=policy).decide(checkpoint="evaluation:1", mode="combined")
        self.assertFalse(signal.passes)
        self.assertTrue(policy.passes)
        self.assertFalse(decision.adaptive_selected)

    def test_smoke_benchmark_is_equal_budget_and_reproducible(self):
        first = run_benchmark(development_instances=12, final_instances=2, seeds=[3, 5], total_budget=200, branch_budget=4)
        second = run_benchmark(development_instances=12, final_instances=2, seeds=[3, 5], total_budget=200, branch_budget=4)
        self.assertEqual({run["total_evaluations"] for run in first["runs"]}, {200})
        self.assertEqual(first["audits"], second["audits"])
        comparable = [(run["family"], run["instance_id"], run["seed"], run["mode"], run["best_objective"], run["regret"]) for run in first["runs"]]
        comparable_second = [(run["family"], run["instance_id"], run["seed"], run["mode"], run["best_objective"], run["regret"]) for run in second["runs"]]
        self.assertEqual(comparable, comparable_second)
        self.assertTrue(all(row["total_evaluations_equal"] for row in first["summary"]))
                # Verify signal discovery, signal validation, and matched policy
        # trials use disjoint development groups.
        for family, audit in first["audits"].items():
            groups = audit["development_group_ids"]

            discovery = set(groups["discovery"])
            validation = set(groups["signal_validation"])
            policy_trials = set(groups["policy_trials"])

            self.assertTrue(discovery.isdisjoint(validation))
            self.assertTrue(discovery.isdisjoint(policy_trials))
            self.assertTrue(validation.isdisjoint(policy_trials))

            costs = audit["audit_cost_by_mode"]
            self.assertEqual(costs["baseline"], 0)
            self.assertGreaterEqual(
                costs["combined"], costs["signal_validation_only"]
            )

        for family, audit in first["audits"].items():
            for mode, expected_cost in audit["audit_cost_by_mode"].items():
                matching = [
                    run for run in first["runs"]
                    if run["family"] == family and run["mode"] == mode
                ]
                self.assertEqual(
                    sum(run["audit_evaluations"] for run in matching),
                    expected_cost,
                )
                self.assertTrue(
                    all(
                        run["one_time_audit_evaluations"] == expected_cost
                        for run in matching
                    )
                )
                self.assertTrue(
                    all(
                        run["total_evaluations"] == 200
                        for run in matching
                    )
                )

if __name__ == "__main__":
    unittest.main()
