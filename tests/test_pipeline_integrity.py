import json
import tempfile
import unittest
from pathlib import Path
import sys
from unittest.mock import patch

from c.bridge import (
    generate_adaptive_population,
    generate_c_population,
    seed_native_random,
)
from src.config import PipelineSettings
from src.docking_vina import MLSurrogateBackend, PeptideDockingScorer
from src.models import Candidate, candidate_has_valid_structure
from src.population_runner import ARISEEngine
from src.predict_structure import predict_population_structures
from src.ranking import CandidateRanker
from src.validation import Validator
from run_pipeline import _candidate_provenance, evidence_summary


def pdb_text(plddt: float = 70.0) -> str:
    return (
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00"
        f" {plddt:5.2f}           C\nEND\n"
    )


class FakeRequests:
    class RequestException(Exception):
        pass

    class ConnectionError(RequestException):
        pass

    @staticmethod
    def post(*_args, **_kwargs):
        raise AssertionError("FakeRequests.post must be configured by the test")


class TestPipelineIntegrity(unittest.TestCase):
    def test_consensus_preserves_negative_sign(self):
        scorer = object.__new__(PeptideDockingScorer)
        self.assertEqual(scorer._consensus_score(-7.0, -6.5), -6.8)

    def test_surrogate_is_deterministic_and_local(self):
        backend = MLSurrogateBackend()
        with patch("src.docking_vina.subprocess.run") as run:
            first = backend.score_fragment("AILMFWYVA")
            second = backend.score_fragment("AILMFWYVA")
            run.assert_not_called()
        self.assertEqual(first, second)
        self.assertNotEqual(first, backend.score_fragment("DEKRSTN"))

    def test_structure_failure_does_not_stop_population(self):
        candidates = [Candidate("A" * 21), Candidate("C" * 21)]

        def post(_url, *, data, **_kwargs):
            if data == candidates[0].sequence:
                raise FakeRequests.ConnectionError("temporary failure")
            return type("Response", (), {
                "status_code": 200,
                "text": pdb_text(75.0),
            })()

        with tempfile.TemporaryDirectory() as directory:
            with patch("src.predict_structure.STRUCTURE_DIR", Path(directory)):
                with patch.dict(sys.modules, {"requests": FakeRequests}):
                    with patch.object(FakeRequests, "post", side_effect=post):
                        with patch("src.predict_structure.time.sleep"):
                            result = predict_population_structures(candidates)
                        self.assertTrue(candidate_has_valid_structure(result[1]))

        self.assertEqual(len(result), 2)
        self.assertFalse(result[0].metadata["structure_available"])
        self.assertIsNone(result[0].structure_confidence)

    def test_all_structure_failures_have_no_fake_confidence(self):
        candidates = [Candidate("A" * 21), Candidate("C" * 21)]
        with tempfile.TemporaryDirectory() as directory:
            with patch("src.predict_structure.STRUCTURE_DIR", Path(directory)):
                with patch.dict(sys.modules, {"requests": FakeRequests}):
                    with patch.object(
                        FakeRequests,
                        "post",
                        side_effect=FakeRequests.ConnectionError("offline"),
                    ):
                        with patch("src.predict_structure.time.sleep"):
                            result = predict_population_structures(candidates)

        self.assertEqual(sum(c.structure_path is not None for c in result), 0)
        self.assertTrue(all(c.structure_confidence is None for c in result))
        self.assertTrue(all(
            c.metadata["structure_status"] == "unavailable"
            for c in result
        ))

    def test_cached_structure_requires_valid_file_and_confidence(self):
        candidate = Candidate("A" * 21)
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            cached = directory / f"{candidate.sequence}.pdb"
            cached.write_text(pdb_text(81.0))
            with patch("src.predict_structure.STRUCTURE_DIR", directory):
                with patch.dict(sys.modules, {"requests": FakeRequests}):
                    with patch.object(FakeRequests, "post") as post:
                        result = predict_population_structures([candidate])[0]
                        post.assert_not_called()
                        self.assertTrue(candidate_has_valid_structure(result))
        self.assertEqual(result.metadata["structure_status"], "cached")

    def test_tier2_eligible_without_structure_is_not_validated(self):
        scorer = PeptideDockingScorer()
        scorer.vina_backend = object()
        candidate = Candidate("AILMFWYVA")
        result = scorer.evaluate_candidate(candidate)
        self.assertTrue(result.metadata["tier2_eligible"])
        self.assertEqual(result.metadata["tier2_status"], "eligible_no_structure")
        self.assertFalse(result.metadata["tier2_validated"])
        self.assertFalse(result.passed_tier_2)
        self.assertIsNone(result.vina_delta_g)
        self.assertEqual(result.metadata["vina_status"], "not_attempted")

    def test_tier2_with_structure_and_backend_is_validated(self):
        scorer = PeptideDockingScorer()
        scorer.vina_backend = type(
            "FakeVina",
            (),
            {"score_fragment": lambda self, fragment, structure: -8.1},
        )()
        candidate = Candidate("AILMFWYVA")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "candidate.pdb"
            path.write_text(pdb_text(75.0))
            candidate.structure_path = path
            candidate.structure_confidence = 75.0
            candidate.metadata["structure_available"] = True
            result = scorer.evaluate_candidate(candidate)
        self.assertTrue(result.metadata["tier2_validated"])
        self.assertEqual(result.metadata["tier2_status"], "validated")
        self.assertTrue(result.passed_tier_2)
        self.assertEqual(result.vina_delta_g, -8.1)
        self.assertEqual(result.metadata["vina_status"], "success")

    def test_zero_vina_is_retained_only_as_a_real_success(self):
        scorer = PeptideDockingScorer()
        scorer.vina_backend = type(
            "ZeroVina",
            (),
            {"score_fragment": lambda self, fragment, structure: 0.0},
        )()
        candidate = Candidate("AILMFWYVA")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "candidate.pdb"
            path.write_text(pdb_text(75.0))
            candidate.structure_path = path
            candidate.structure_confidence = 75.0
            candidate.metadata["structure_available"] = True
            result = scorer.evaluate_candidate(candidate)
        self.assertTrue(result.metadata["tier2_validated"])
        self.assertEqual(result.metadata["vina_status"], "success")
        self.assertEqual(result.vina_delta_g, 0.0)

    def test_failed_fragment_is_not_a_fake_score(self):
        class PartialBackend:
            def score_fragment(self, fragment):
                if fragment.startswith("C"):
                    raise RuntimeError("fragment failed")
                return -6.0

        scorer = PeptideDockingScorer(backend=PartialBackend())
        scorer.vina_backend = None
        result = scorer.evaluate_candidate(Candidate("ACDEFGHIKL"))
        self.assertEqual(result.docking_scores, [-6.0])
        self.assertEqual(result.metadata["docking_status"], "incomplete")
        self.assertIn("CDEFGHIKL", result.metadata["docking_failed_fragments"])

    def test_native_population_is_reproducible_when_reseeded(self):
        seed_native_random(616)
        first = generate_c_population("A" * 21, pop_size=6, mutation_rate=0.2)
        seed_native_random(616)
        second = generate_c_population("A" * 21, pop_size=6, mutation_rate=0.2)
        self.assertEqual(first, second)

    def test_raw_and_normalized_docking_ordering(self):
        raw = [-7.0, -6.5, -5.5]
        normalized = CandidateRanker().normalize(raw, reverse=True)
        self.assertGreater(normalized[0], normalized[1])
        self.assertGreater(normalized[1], normalized[2])

    def test_arise_receives_overall_score_and_excludes_missing_docking(self):
        engine = ARISEEngine()
        valid = Candidate("A" * 21, overall_score=0.8)
        valid.metadata["ranking_observation_valid"] = True
        invalid = Candidate("C" * 21, overall_score=0.9)
        invalid.metadata["ranking_observation_valid"] = False
        engine.observe_generation([valid, invalid])
        self.assertEqual(valid.metadata["arise_observation_source"], "overall_score")
        self.assertEqual(valid.metadata["arise_observation_score"], 0.8)
        self.assertNotIn("arise_observation_score", invalid.metadata)
        self.assertEqual(engine.importance_observations()[0], 1)

    def test_final_candidate_provenance_uses_candidate_fields(self):
        candidate = Candidate("A" * 21, overall_score=0.73, confidence=0.61)
        candidate.structure_confidence = 82.05
        candidate.metadata.update({
            "generation": 4,
            "generation_seed": "C" * 21,
            "structure_available": True,
            "structure_status": "success",
            "consensus_docking": -6.3,
            "tier2_eligible": True,
            "tier2_attempted": True,
            "tier2_validated": True,
            "tier2_status": "validated",
            "vina_status": "success",
            "ranking_components": {"docking": 0.2},
            "generation_overall_score": 0.74,
            "generation_confidence": 0.59,
            "generation_ranking_components": {"docking": 0.21},
        })
        candidate.strongest_anchor_delta_g = -6.5
        candidate.mean_delta_g = -6.0
        candidate.vina_delta_g = 0.0
        provenance = _candidate_provenance(candidate)
        self.assertEqual(provenance["generation_first_observed"], 4)
        self.assertEqual(provenance["pLDDT"], 82.05)
        self.assertEqual(provenance["vina_score"], 0.0)
        self.assertEqual(provenance["overall_score"], 0.73)
        self.assertEqual(provenance["generation_overall_score"], 0.74)

    def test_evidence_summary_matches_candidate_states(self):
        first = Candidate("A" * 21)
        first.metadata.update({
            "docking_status": "complete",
            "tier2_eligible": True,
            "tier2_attempted": True,
            "tier2_validated": True,
            "vina_status": "success",
        })
        first.vina_delta_g = 0.0
        second = Candidate("C" * 21)
        second.metadata.update({
            "docking_status": "complete",
            "tier2_eligible": False,
            "tier2_attempted": False,
            "tier2_validated": False,
            "vina_status": "not_attempted",
        })
        counts = evidence_summary([first, second])
        self.assertEqual(counts["tier2_eligible"], 1)
        self.assertEqual(counts["tier2_attempted"], 1)
        self.assertEqual(counts["tier2_validated"], 1)
        self.assertEqual(counts["vina_success"], 1)
        self.assertEqual(counts["vina_missing"], 1)

    def test_latest_audit_preserves_one_summary_per_generation(self):
        audit_path = Path("results/experiment_audit.json")
        if not audit_path.exists():
            self.skipTest("latest experiment audit is not present")
        audit = json.loads(audit_path.read_text())
        summaries = audit["optimization"]["generation_summaries"]
        generations = [summary["generation"] for summary in summaries]
        self.assertEqual(generations, list(range(1, len(generations) + 1)))

    def test_adaptive_population_accepts_partial_refill_size(self):
        seed_native_random(616)
        result = generate_adaptive_population(
            "A" * 21,
            [0.5] * 21,
            pop_size=1,
        )
        self.assertEqual(len(result), 1)
        self.assertEqual(len(result[0]), 21)

    def test_configuration_rejects_invalid_values(self):
        with self.assertRaises(ValueError):
            PipelineSettings(population_size=0)
        with self.assertRaises(ValueError):
            PipelineSettings(mutation_rate=1.5)

    def test_validator_rejects_inconsistent_structure_state(self):
        candidate = Candidate("A" * 21)
        candidate.metadata["structure_available"] = False
        candidate.structure_confidence = 73.0
        result = Validator().validate([candidate], tempfile.gettempdir())
        self.assertFalse(result.passed)
        self.assertTrue(any("pLDDT" in error for error in result.errors))


if __name__ == "__main__":
    unittest.main()
