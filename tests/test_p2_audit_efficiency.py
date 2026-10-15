import unittest

from scripts.benchmark_p2_audit_efficiency import (
    PLAN_DEFINITIONS,
    _sign_permutation_p,
)


class AuditEfficiencyTests(unittest.TestCase):
    def test_predeclared_audit_costs_are_bounded_and_monotone(self):
        costs = {
            allocation: plan["discovery_groups"] * 8
            + plan["validation_groups"] * 8
            + 2 * plan["branch_budget"] * plan["matched_trials"]
            for allocation, plan in PLAN_DEFINITIONS.items()
        }
        self.assertEqual(costs, {0: 0, 5: 16, 10: 40, 20: 80, 80: 320})
        self.assertEqual(sorted(costs), [0, 5, 10, 20, 80])

    def test_sequential_sign_test_is_deterministic(self):
        values = [0.2, 0.4, 0.3, 0.5, 0.1]
        self.assertEqual(_sign_permutation_p(values, seed=13, permutations=200), _sign_permutation_p(values, seed=13, permutations=200))
        self.assertLess(_sign_permutation_p(values, seed=13, permutations=200), 0.5)


if __name__ == "__main__":
    unittest.main()
