import unittest

from scripts.benchmark_online_policy_audit import run_benchmark


class OnlinePolicyBenchmarkTests(unittest.TestCase):
    def test_all_methods_respect_total_budget(self):
        payload = run_benchmark(
            development_instances=12,
            final_instances=2,
            seeds=[3],
            total_budget=100,
            branch_budget=2,
            audit_replications=2,
            audit_interval=8,
        )
        self.assertEqual(
            set(payload["methods"]),
            {"baseline", "fixed_adaptive", "online_audit"},
        )
        for run in payload["runs"]:
            self.assertEqual(run["total_evaluations"], 100)
            self.assertEqual(
                run["search_evaluations"]
                + run["online_audit_evaluations"]
                + run["allocated_discovery_evaluations"],
                100,
            )

    def test_online_audits_have_exact_cost(self):
        payload = run_benchmark(
            development_instances=12,
            final_instances=2,
            seeds=[3],
            total_budget=100,
            branch_budget=2,
            audit_replications=2,
            audit_interval=8,
        )
        online_runs = [
            run for run in payload["runs"]
            if run["mode"] == "online_audit"
        ]
        for run in online_runs:
            self.assertEqual(
                run["online_audit_evaluations"],
                run["audit_count"] * 2 * 2 * 2,
            )

    def test_fixed_seed_reproduces_quality_metrics(self):
        kwargs = dict(
            development_instances=12,
            final_instances=2,
            seeds=[3],
            total_budget=100,
            branch_budget=2,
            audit_replications=2,
            audit_interval=8,
        )
        first = run_benchmark(**kwargs)
        second = run_benchmark(**kwargs)
        first_values = [
            (
                run["family"],
                run["instance_id"],
                run["mode"],
                run["best_expected_objective"],
                run["regret"],
            )
            for run in first["runs"]
        ]
        second_values = [
            (
                run["family"],
                run["instance_id"],
                run["mode"],
                run["best_expected_objective"],
                run["regret"],
            )
            for run in second["runs"]
        ]
        self.assertEqual(first_values, second_values)


if __name__ == "__main__":
    unittest.main()
