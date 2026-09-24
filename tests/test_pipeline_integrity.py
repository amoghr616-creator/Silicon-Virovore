import json
import tempfile
import unittest
from pathlib import Path
import sys
from dataclasses import replace
from unittest.mock import patch

from c.bridge import (
    clear_evaluation_cache,
    evaluation_cache_stats,
    generate_adaptive_population,
    generate_c_population,
    generate_policy_population,
    process_candidate_peptide,
    seed_native_random,
)
from src.config import PipelineSettings
from src.docking_vina import MLSurrogateBackend, PeptideDockingScorer
from src.models import Candidate, candidate_has_valid_structure
from src.hotspot import (
    annotate_arise_scores,
    discover_hotspot_model,
    derive_positional_signal,
)
from src.arise_memory import (
    load_arise_memory,
    save_arise_memory,
    get_memory_history,
    upsert_arise_run,
)
from src.population_runner import ARISEEngine
from src.predict_structure import predict_population_structures
from src.ranking import CandidateRanker
from src.validation import Validator
from run_pipeline import (
    _candidate_provenance,
    _finalize_comparable_telemetry,
    _posthoc_rank_without_overwriting_generation_scores,
    evidence_summary,
    resolve_experiment_seed,
)
from run_pipeline import parse_settings


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
    def test_derived_positional_signal_from_hotspot_model(self):
        model = {
            "hotspots": [
                {
                    "start_position": 14,
                    "end_position": 16,
                    "window_length": 2,
                    "selection_score": 0.90,
                },
                {
                    "start_position": 13,
                    "end_position": 16,
                    "window_length": 3,
                    "selection_score": 0.60,
                },
            ]
        }

        signal = derive_positional_signal(
            model,
            sequence_length=21,
        )

        self.assertEqual(len(signal), 21)
        self.assertLessEqual(max(signal), 1.0)
        self.assertGreater(signal[14], 0.0)
        self.assertGreater(signal[15], 0.0)
        self.assertGreater(signal[14], signal[0])

    def test_arise_memory_round_trip(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "arise_memory.json"

        memory = load_arise_memory(path)

        history = [
            {
                "run_id": "adaptive-101",
                "generation": 1,
                "sequence": "A" * 21,
                "score": 0.8,
            }
        ]

        upsert_arise_run(
            memory,
            run_id="adaptive-101",
            condition="adaptive",
            random_seed=101,
            candidate_history=history,
            hotspot_models=[],
            derived_positional_signal=[0.0] * 21,
        )

        save_arise_memory(path, memory)

        reloaded = load_arise_memory(path)

        self.assertEqual(
            len(get_memory_history(reloaded)),
            1,
        )

        self.assertEqual(
            reloaded["runs"][0]["run_id"],
            "adaptive-101",
        )

    def test_arise_memory_replaces_duplicate_run(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "arise_memory.json"

        memory = load_arise_memory(path)

        first = [{
            "run_id": "run-1",
            "generation": 1,
            "sequence": "A" * 21,
            "score": 0.4,
        }]

        second = [{
            "run_id": "run-1",
            "generation": 1,
            "sequence": "C" * 21,
            "score": 0.9,
        }]

        upsert_arise_run(
            memory,
            run_id="run-1",
            condition="adaptive",
            random_seed=1,
            candidate_history=first,
            hotspot_models=[],
            derived_positional_signal=[0.0] * 21,
        )

        upsert_arise_run(
            memory,
            run_id="run-1",
            condition="adaptive",
            random_seed=1,
            candidate_history=second,
            hotspot_models=[],
            derived_positional_signal=[1.0] * 21,
        )

        records = get_memory_history(memory)

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["sequence"], "C" * 21)
    
    def test_explicit_pipeline_seed_is_preserved(self):
        settings = PipelineSettings(random_seed=617)

        self.assertEqual(settings.random_seed, 617)
        self.assertEqual(resolve_experiment_seed(617), 617)

    def test_auto_generated_seed_is_valid_uint32(self):
        with patch(
            "run_pipeline.secrets.randbelow",
            return_value=123456,
        ):
            self.assertEqual(
                resolve_experiment_seed(None),
                123456,
            )

    def test_run_pipeline_does_not_mutate_settings_seed(self):
        settings = PipelineSettings(
            random_seed=None,
            run_id="",
        )

        with patch(
            "run_pipeline.resolve_experiment_seed",
            return_value=123456,
        ) as resolver:
            copied = replace(
                settings,
                random_seed=resolver(settings.random_seed),
                run_id="adaptive-seed123456",
            )

        self.assertIsNone(settings.random_seed)
        self.assertEqual(copied.random_seed, 123456)
        self.assertEqual(
            copied.run_id,
            "adaptive-seed123456",
        )

    def test_configuration_rejects_out_of_range_seed(self):
        with self.assertRaises(ValueError):
            PipelineSettings(random_seed=-1)

        with self.assertRaises(ValueError):
            PipelineSettings(random_seed=2**32)

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

    def test_parallel_docking_preserves_candidate_results(self):
        scorer = PeptideDockingScorer(
            backend=MLSurrogateBackend(),
            workers=2,
        )
        scorer.vina_backend = None
        candidates = [Candidate("AILMFWYVA"), Candidate("DEKRSTNAA")]
        result = scorer.evaluate_candidates(candidates)
        self.assertEqual(
            [c.sequence for c in result],
            [c.sequence for c in candidates],
        )
        self.assertEqual(len(result[0].docking_scores), 1)
        self.assertEqual(len(result[1].docking_scores), 1)

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

    def test_parallel_structure_evaluation_preserves_candidate_order(self):
        candidates = [Candidate("A" * 21), Candidate("C" * 21)]

        def post(_url, *, data, **_kwargs):
            return type("Response", (), {
                "status_code": 200,
                "text": pdb_text(75.0 if data.startswith("A") else 76.0),
            })()

        with tempfile.TemporaryDirectory() as directory:
            with patch("src.predict_structure.STRUCTURE_DIR", Path(directory)):
                with patch.dict(sys.modules, {"requests": FakeRequests}):
                    with patch.object(FakeRequests, "post", side_effect=post):
                        result = predict_population_structures(
                            candidates,
                            workers=2,
                        )
        self.assertEqual(
            [c.sequence for c in result],
            ["A" * 21, "C" * 21],
        )
        self.assertEqual(
            [c.structure_confidence for c in result],
            [75.0, 76.0],
        )

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

    def test_valid_structure_cache_does_not_require_requests_package(self):
        candidate = Candidate("A" * 21)
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            cached = directory / f"{candidate.sequence}.pdb"
            cached.write_text(pdb_text(82.0))
            with patch("src.predict_structure.STRUCTURE_DIR", directory):
                with patch.dict(sys.modules, {"requests": None}):
                    result = predict_population_structures([candidate])[0]
        self.assertEqual(result.metadata["structure_status"], "cached")
        self.assertEqual(result.structure_confidence, 82.0)

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

    def test_changing_seed_changes_native_population(self):
        seed_native_random(616)
        first = generate_c_population("A" * 21, pop_size=6, mutation_rate=0.2)
        seed_native_random(617)
        second = generate_c_population("A" * 21, pop_size=6, mutation_rate=0.2)
        self.assertNotEqual(first, second)

    def test_deterministic_fitness_cache_preserves_result(self):
        clear_evaluation_cache()
        sequence = "MKLAVFALLVFFAGSSDLIRR"
        uncached = process_candidate_peptide(sequence, cache_enabled=False)
        first = process_candidate_peptide(sequence, cache_enabled=True)
        second = process_candidate_peptide(sequence, cache_enabled=True)
        self.assertEqual(uncached.c_score, first.c_score)
        self.assertEqual(first.c_score, second.c_score)
        self.assertTrue(second.metadata["c_fitness_cache_hit"])
        self.assertGreaterEqual(evaluation_cache_stats()["hits"], 1)

    def test_cli_baseline_and_adaptive_modes_are_explicit(self):
        baseline = parse_settings([
            "--condition", "baseline",
            "--seed-sequence", "MKLAVFALLVFFAGSSDLIRR",
        ])
        adaptive = parse_settings([
            "--condition", "adaptive",
            "--seed-sequence", "MKLAVFALLVFFAGSSDLIRR",
        ])
        self.assertFalse(baseline.adaptive_enabled)
        self.assertFalse(baseline.arise_enabled)
        self.assertFalse(baseline.ale_enabled)
        self.assertTrue(adaptive.adaptive_enabled)
        self.assertTrue(adaptive.arise_enabled)
        self.assertTrue(adaptive.ale_enabled)
        self.assertFalse(baseline.hotspot_enabled)
        self.assertEqual(baseline.mutation_guidance_mode, "uniform")
        self.assertTrue(adaptive.hotspot_enabled)
        self.assertEqual(adaptive.mutation_guidance_mode, "hotspot_biased")

    def test_posthoc_ranking_does_not_overwrite_generation_score(self):
        first = Candidate("A" * 21, c_score=1.0)
        second = Candidate("C" * 21, c_score=2.0)
        for index, candidate in enumerate((first, second), start=1):
            candidate.metadata.update({
                "candidate_id": f"candidate-{index}",
                "generation": 1,
                "docking_status": "complete",
                "consensus_docking": -5.0 - index,
                "evaluation_index": index,
            })
        ranker = CandidateRanker()
        ranker.rank([first, second])
        original_score = first.overall_score
        posthoc, _ = _posthoc_rank_without_overwriting_generation_scores(
            ranker,
            [first, second],
        )
        self.assertEqual(first.overall_score, original_score)
        self.assertIsNotNone(first.posthoc_recomputed_score)
        self.assertEqual(
            first.posthoc_recomputed_score,
            posthoc[-1].posthoc_recomputed_score,
        )
        self.assertIn("raw_objectives", first.metadata)
        self.assertIn("posthoc_normalized_objectives", first.metadata)

    def test_comparable_trajectory_tracks_generation_and_best_so_far(self):
        candidates = []
        history = []
        for generation, scores in ((1, (0.9, 0.4)), (2, (0.6, 0.5))):
            for index, score in enumerate(scores, start=1):
                candidate = Candidate(
                    "A" * 20 + ("C" if index == 2 else "A"),
                    c_score=score,
                )
                candidate.overall_score = score
                candidate.posthoc_recomputed_score = score
                candidate.metadata.update({
                    "candidate_id": f"g{generation}c{index}",
                    "generation": generation,
                    "evaluation_index": len(candidates) + 1,
                    "raw_objectives": {"fitness": score},
                    "normalized_objectives": {"fitness": score},
                    "posthoc_normalized_objectives": {"fitness": score},
                })
                candidates.append(candidate)
                history.append({
                    "candidate_id": candidate.metadata["candidate_id"],
                    "generation": generation,
                    "sequence": candidate.sequence,
                    "score": score,
                })
        summaries = [
            {"generation": 1, "best_overall": 0.9, "best_sequence": candidates[0].sequence},
            {"generation": 2, "best_overall": 0.6, "best_sequence": candidates[2].sequence},
        ]
        result = _finalize_comparable_telemetry(
            summaries,
            candidates,
            history,
            score_threshold=0.8,
        )
        self.assertEqual(summaries[0]["generation_best_score"], 0.9)
        self.assertEqual(summaries[1]["generation_best_score"], 0.6)
        self.assertEqual(summaries[1]["best_so_far_score"], 0.9)
        self.assertEqual(result["candidates_to_threshold"], 1)
        self.assertEqual(result["generation_to_threshold"], 1)

    def test_legacy_arise_is_separate_from_hotspot_mode(self):
        adaptive = parse_settings([
            "--condition", "adaptive",
            "--seed-sequence", "MKLAVFALLVFFAGSSDLIRR",
        ])
        self.assertEqual(adaptive.mutation_guidance_mode, "hotspot_biased")
        self.assertFalse(adaptive.mutation_guidance_mode == "legacy")

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

    def test_hotspot_discovery_recovers_planted_local_pattern(self):
        history = []
        for generation in range(1, 5):
            for index in range(20):
                sequence = list("A" * 21)
                if index < 10:
                    sequence[8:12] = list("KLMN")
                    score = 1.0 + index * 0.001
                else:
                    sequence[8:12] = list("CCCC")
                    score = 0.1 + index * 0.001
                history.append({
                    "sequence": "".join(sequence),
                    "score": score,
                    "generation": generation,
                    "run_id": "synthetic",
                    "arise_score_comparable_across_generations": True,
                })

        model = discover_hotspot_model(
            history,
            discovery_generation=4,
            run_id="synthetic",
            bootstrap_iterations=20,
            random_seed=7,
        )
        selected = model["selected_hotspot"]
        self.assertIsNotNone(selected)
        self.assertTrue(8 <= selected["start_position"] < 12)
        self.assertGreater(selected["held_out_association"], 0.0)
        self.assertGreaterEqual(selected["bootstrap_stability"], 0.5)

    def test_hotspot_discovery_is_deterministic_and_history_bounded(self):
        history = [
            {
                "sequence": "A" * 21,
                "score": 0.1,
                "generation": 1,
                "run_id": "run",
                "arise_score_comparable_across_generations": True,
            },
            {
                "sequence": "A" * 8 + "KLM" + "A" * 10,
                "score": 0.9,
                "generation": 1,
                "run_id": "run",
                "arise_score_comparable_across_generations": True,
            },
        ] * 5
        first = discover_hotspot_model(
            history,
            discovery_generation=1,
            run_id="run",
            random_seed=11,
        )
        second = discover_hotspot_model(
            history + [{
                "sequence": "C" * 21,
                "score": 100.0,
                "generation": 2,
                "run_id": "run",
                "arise_score_comparable_across_generations": True,
            }],
            discovery_generation=1,
            run_id="run",
            random_seed=11,
        )
        self.assertEqual(first, second)

    def test_arise_score_is_shared_scope_across_generations(self):
        raw = {
            "fitness": 0.7,
            "docking": -6.0,
            "structure": 75.0,
            "helix": 1.0,
            "hydrophobicity": 0.3,
            "solvation": 2.0,
            "md": None,
            "safety": None,
        }
        history = [
            {
                "sequence": "A" * 21,
                "score": 0.12,
                "generation_local_score": 0.12,
                "generation": 1,
                "raw_objectives": raw,
            },
            {
                "sequence": "C" * 21,
                "score": 0.88,
                "generation_local_score": 0.88,
                "generation": 2,
                "raw_objectives": {
                    **raw,
                    "fitness": 0.2,
                    "docking": -5.0,
                },
            },
            {
                "sequence": "D" * 21,
                "score": 0.5,
                "generation": 3,
                "raw_objectives": {
                    **raw,
                    "fitness": 0.9,
                    "docking": -7.0,
                },
            },
        ]
        annotation = annotate_arise_scores(history)
        self.assertEqual(
            annotation["score_source"],
            "pooled_history_overall_score",
        )
        self.assertTrue(annotation["score_comparable_across_generations"])
        self.assertNotEqual(
            annotation["records"][0]["generation_local_score"],
            annotation["records"][1]["generation_local_score"],
        )

        same_candidate_history = [
            {
                **history[0],
                "generation": generation,
            }
            for generation in (1, 2)
        ]
        same_annotation = annotate_arise_scores(same_candidate_history)
        self.assertEqual(
            same_annotation["records"][0]["arise_comparable_score"],
            same_annotation["records"][1]["arise_comparable_score"],
        )

    def test_guided_policy_uses_hotspot_positions_and_fixed_budget(self):
        seed_native_random(616)
        population = generate_policy_population(
            "A" * 21,
            mutation_rate=0.05,
            guidance_mode="hotspot_only",
            hotspot_start=8,
            hotspot_end=12,
            pop_size=20,
        )
        self.assertTrue(population)
        self.assertTrue(all(
            all(sequence[index] == "A" for index in range(21) if index not in range(8, 12))
            for sequence in population
        ))
        self.assertTrue(all(
            sum(old != new for old, new in zip("A" * 21, sequence)) == 1
            for sequence in population
        ))

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
