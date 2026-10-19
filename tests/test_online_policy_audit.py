import random
import unittest

from src.online_policy_audit import run_online_audit


class OnlinePolicyAuditTests(unittest.TestCase):
    def test_equal_branch_budgets_and_total_cost(self):
        result = run_online_audit(
            start_state=(0, 0, 0, 0),
            objective=lambda state, rng: float(sum(state)),
            preferred_position=0,
            branch_budget=8,
            replications=4,
            seed=11,
        )

        self.assertEqual(result.total_objective_evaluations, 64)
        self.assertEqual(len(result.trials), 4)

        for trial in result.trials:
            self.assertEqual(
                trial.baseline.start_state,
                trial.adaptive.start_state,
            )
            self.assertEqual(
                trial.baseline.evaluations,
                trial.adaptive.evaluations,
            )
            self.assertEqual(trial.baseline.evaluations, 8)
            self.assertEqual(trial.adaptive.evaluations, 8)
            self.assertEqual(
                len(trial.baseline.evaluated_states),
                trial.baseline.evaluations,
            )
            self.assertEqual(
                len(trial.adaptive.evaluated_states),
                trial.adaptive.evaluations,
            )
            self.assertEqual(
                trial.baseline.evaluated_states[0],
                trial.baseline.start_state,
            )

    def test_result_is_reproducible_for_fixed_seed(self):
        kwargs = dict(
            start_state=(0, 0, 0, 0),
            objective=lambda state, rng: float(sum(state)),
            preferred_position=0,
            branch_budget=8,
            replications=4,
            seed=23,
        )

        first = run_online_audit(**kwargs)
        second = run_online_audit(**kwargs)

        self.assertEqual(
            first.mean_paired_advantage,
            second.mean_paired_advantage,
        )
        self.assertEqual(first.selected_policy, second.selected_policy)

    def test_no_preferred_signal_does_not_select_adaptive(self):
        result = run_online_audit(
            start_state=(0, 0, 0, 0),
            objective=lambda state, rng: float(sum(state)),
            preferred_position=None,
            branch_budget=8,
            replications=4,
            seed=7,
        )

        self.assertEqual(result.selected_policy, "baseline")

    def test_invalid_inputs_fail(self):
        with self.assertRaises(ValueError):
            run_online_audit(
                start_state=(),
                objective=lambda state, rng: 0.0,
                preferred_position=None,
            )

        with self.assertRaises(ValueError):
            run_online_audit(
                start_state=(0, 0),
                objective=lambda state, rng: 0.0,
                preferred_position=3,
            )


if __name__ == "__main__":
    unittest.main()
