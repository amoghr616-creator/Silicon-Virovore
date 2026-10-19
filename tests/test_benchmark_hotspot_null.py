import unittest
import csv
from pathlib import Path
import tempfile

from scripts.benchmark_hotspot_null import (
    _run_trial,
    _blocked_max_statistic_p_value,
    _permute_scores_within_blocks,
    LINEAGE_BLOCK_FIELDS,
    interval_overlap,
    generate_synthetic_history,
    run_benchmark,
    wilson_interval,
    write_outputs,
)


class HotspotNullBenchmarkTests(unittest.TestCase):
    def test_generators_produce_reproducible_detector_schema(self):
        for scenario in (
            "exchangeable_null",
            "structured_null",
            "injected_hotspot",
        ):
            first, target = generate_synthetic_history(
                scenario,
                seed=113,
                effect_size=0.8,
            )
            second, second_target = generate_synthetic_history(
                scenario,
                seed=113,
                effect_size=0.8,
            )
            self.assertEqual(first, second)
            self.assertEqual(target, second_target)
            self.assertEqual(len(first), 8 * 4 * 12)
            for record in first:
                self.assertIn("candidate_id", record)
                self.assertIn("parent_sequence", record)
                self.assertTrue(record["sequence"])
                self.assertIsInstance(record["score"], float)
                self.assertIn("generation", record)
                self.assertIn("run_id", record)
                self.assertTrue(
                    record["arise_score_comparable_across_generations"]
                )
        self.assertEqual(target, (4, 7))

    def test_metrics_return_expected_binomial_and_window_values(self):
        lower, upper = wilson_interval(0, 10)
        self.assertEqual(lower, 0.0)
        self.assertAlmostEqual(upper, 0.2775, places=3)
        self.assertEqual(interval_overlap((4, 7), (4, 7)), 1.0)
        self.assertAlmostEqual(interval_overlap((3, 7), (4, 8)), 0.6)
        self.assertEqual(interval_overlap((0, 2), (4, 7)), 0.0)
        self.assertIsNone(interval_overlap(None, (4, 7)))

    def test_actual_detector_trial_is_reproducible(self):
        first = _run_trial(
            "injected_hotspot",
            trial_seed=991,
            permutation_iterations=9,
            effect_size=1.0,
        )
        second = _run_trial(
            "injected_hotspot",
            trial_seed=991,
            permutation_iterations=9,
            effect_size=1.0,
        )
        self.assertEqual(first, second)
        self.assertIn(
            first["candidate_shuffle_p_value"],
            (None, *[i / 10 for i in range(1, 11)]),
        )
        self.assertIn(
            first["lineage_generation_block_shuffle_p_value"],
            (None, *[i / 10 for i in range(1, 11)]),
        )

    def test_block_permutation_preserves_hierarchy_and_score_multisets(self):
        records, _ = generate_synthetic_history("structured_null", seed=91)
        permuted = _permute_scores_within_blocks(
            records,
            seed=22,
            block_fields=LINEAGE_BLOCK_FIELDS,
        )
        self.assertEqual(
            [record["sequence"] for record in records],
            [record["sequence"] for record in permuted],
        )
        self.assertEqual(
            [tuple(record[field] for field in LINEAGE_BLOCK_FIELDS) for record in records],
            [tuple(record[field] for field in LINEAGE_BLOCK_FIELDS) for record in permuted],
        )
        original_blocks = {}
        permuted_blocks = {}
        for record in records:
            key = tuple(record[field] for field in LINEAGE_BLOCK_FIELDS)
            original_blocks.setdefault(key, []).append(record["score"])
        for record in permuted:
            key = tuple(record[field] for field in LINEAGE_BLOCK_FIELDS)
            permuted_blocks.setdefault(key, []).append(record["score"])
        self.assertEqual(
            {key: sorted(scores) for key, scores in original_blocks.items()},
            {key: sorted(scores) for key, scores in permuted_blocks.items()},
        )

    def test_block_test_is_reproducible_and_rejects_missing_block_keys(self):
        records, _ = generate_synthetic_history("structured_null", seed=92)
        discovery = [record for record in records if record["run_id"] != "synthetic-run-000"]
        kwargs = dict(
            iterations=19,
            seed=36,
            block_fields=LINEAGE_BLOCK_FIELDS,
        )
        self.assertEqual(
            _blocked_max_statistic_p_value(discovery, **kwargs),
            _blocked_max_statistic_p_value(discovery, **kwargs),
        )
        self.assertIsNone(
            _blocked_max_statistic_p_value(
                [{"sequence": "AAAA", "score": 1.0}],
                iterations=5,
                seed=1,
                block_fields=LINEAGE_BLOCK_FIELDS,
            )
        )
        with self.assertRaises(ValueError):
            _blocked_max_statistic_p_value(
                records,
                iterations=0,
                seed=1,
                block_fields=LINEAGE_BLOCK_FIELDS,
            )

    def test_seeded_structured_null_exposes_candidate_shuffle_failure(self):
        trial = _run_trial(
            "structured_null",
            trial_seed=21261016,
            permutation_iterations=199,
        )
        self.assertEqual(trial["candidate_shuffle_p_value"], 0.005)
        self.assertTrue(trial["candidate_shuffle_significant"])
        self.assertEqual(
            trial["lineage_generation_block_shuffle_p_value"],
            0.075,
        )
        self.assertFalse(
            trial["lineage_generation_block_shuffle_significant"]
        )

    def test_small_benchmark_runs_all_scenarios_and_reports_intervals(self):
        result = run_benchmark(
            mode="quick",
            trials=1,
            permutations=2,
            seed=77,
        )
        self.assertEqual(result["trial_count_per_scenario"], 1)
        self.assertEqual(result["permutation_iterations_per_trial"], 2)
        self.assertEqual(len(result["summaries"]), 10)
        self.assertEqual(len(result["trials"]), 5)
        for row in result["summaries"]:
            self.assertEqual(row["trials"], 1)
            self.assertEqual(len(row["rate_ci95"]), 2)
            self.assertTrue(0.0 <= row["rate_ci95"][0] <= row["rate_ci95"][1] <= 1.0)
        alternatives = [
            row for row in result["summaries"]
            if row["scenario"] == "injected_hotspot"
        ]
        self.assertEqual(len(alternatives), 6)
        self.assertTrue(all("detection_power_ci95" in row for row in alternatives))
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory) / "benchmark"
            write_outputs(result, output_dir)
            with (output_dir / "summary.csv").open(newline="") as handle:
                written_rows = list(csv.DictReader(handle))
            self.assertEqual(len(written_rows), 10)
            self.assertTrue((output_dir / "report.md").is_file())


if __name__ == "__main__":
    unittest.main()
