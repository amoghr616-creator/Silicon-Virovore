import unittest

from prototype2.prepare_validation_panels import select_nonredundant_panel


class ValidationPanelSelectionTests(unittest.TestCase):
    def test_selects_by_native_fitness_and_deduplicates_sequence(self):
        records = [
            {
                "candidate_id": "r:g1:c1",
                "generation": 1,
                "sequence": "AAAA",
                "native_fitness": 0.8,
                "score": 0.7,
                "vina_score": -99.0,
            },
            {
                "candidate_id": "r:g2:c2",
                "generation": 2,
                "sequence": "AAAA",
                "native_fitness": 0.9,
                "score": 0.6,
                "vina_score": 100.0,
            },
            {
                "candidate_id": "r:g1:c3",
                "generation": 1,
                "sequence": "CCCC",
                "native_fitness": 0.7,
                "score": 0.5,
                "vina_score": -100.0,
            },
        ]
        panel = select_nonredundant_panel(
            records,
            run_id="r",
            arm="baseline",
            panel_size=2,
        )
        self.assertEqual([row["candidate_id"] for row in panel], ["r:g2:c2", "r:g1:c3"])
        self.assertEqual(panel[0]["original_fitness_score"], 0.9)
        self.assertEqual(panel[0]["sequence_occurrences_in_source_run"], 2)
        self.assertFalse(panel[0]["unique_in_source_run"])
        self.assertTrue(all(row["unique_in_selected_panel"] for row in panel))
        self.assertTrue(all(row["selection_score_source"] == "candidate_history.native_fitness" for row in panel))

    def test_insufficient_unique_sequences_fails_closed(self):
        records = [
            {
                "candidate_id": "one",
                "generation": 1,
                "sequence": "AAAA",
                "native_fitness": 0.8,
                "score": 0.7,
            },
            {
                "candidate_id": "two",
                "generation": 2,
                "sequence": "AAAA",
                "native_fitness": 0.9,
                "score": 0.6,
            },
        ]
        with self.assertRaises(ValueError):
            select_nonredundant_panel(
                records,
                run_id="r",
                arm="baseline",
                panel_size=2,
            )


if __name__ == "__main__":
    unittest.main()
