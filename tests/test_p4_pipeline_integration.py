import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import src.p4_pipeline_integration as p4


def make_candidate(sequence):
    return SimpleNamespace(
        sequence=sequence,
        metadata={},
        c_score=0.5,
        overall_score=0.5,
        rank=0,
    )


def write_manifest(tmp_path, sequence, *, expected_hash=None):
    pose = tmp_path / "pose.pdb"
    pose.write_text("TEST STRUCTURE\n")

    actual_hash = hashlib.sha256(pose.read_bytes()).hexdigest()

    manifest = {
        "schema_version": 1,
        "target": {
            "target_structure_id": "TEST_TARGET",
            "receptor_chain": "R",
            "peptide_chain": "P",
            "hotspots": [
                {
                    "chain": "R",
                    "residue_number": 42,
                    "residue_name": "GLU",
                }
            ],
        },
        "complexes": [
            {
                "sequence": sequence,
                "path": pose.name,
                "sha256": expected_hash or actual_hash,
                "provenance": "test_fixture",
            }
        ],
    }

    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest))
    return path, pose, actual_hash



def test_selection_comparison_uses_same_cohort_and_reports_target_effect():
    sequence_winner = make_candidate("SEQUENCE_ONLY_WINNER")
    target_winner = make_candidate("TARGET_ORDERED_WINNER")
    sequence_winner.metadata["p4_target_evidence"] = {
        "status": "success",
        "target_geometry_score_mean": 0.35,
        "hotspot_proximity_score": 0.40,
    }
    target_winner.metadata["p4_target_evidence"] = {
        "status": "success",
        "target_geometry_score_mean": 0.85,
        "hotspot_proximity_score": 0.90,
    }

    result = p4.build_selection_comparison(
        [sequence_winner, target_winner],
        [target_winner, sequence_winner],
        {
            "pareto_reordered_candidates": True,
            "selection_method": "pareto_front_then_target_geometry_then_hotspot_proximity",
        },
    )

    assert result["same_candidate_cohort_verified"] is True
    assert result["selection_changed"] is True
    assert result["selection_effect"] == "target_ordering_changed_generation_winner"
    assert result["target_ordered_winner_sequence_only_rank"] == 2
    assert result["winner_metric_deltas"][
        "target_winner_minus_sequence_only_winner_target_geometry_score_mean"
    ] == pytest.approx(0.50)
    assert "does not validate binding affinity" in result["interpretation_limit"]


def test_selection_comparison_preserves_duplicate_cohort_multiplicity():
    candidate_a = make_candidate("A")
    candidate_b = make_candidate("B")
    result = p4.build_selection_comparison(
        [candidate_a, candidate_a, candidate_b],
        [candidate_a, candidate_b, candidate_b],
        {"pareto_reordered_candidates": True},
    )
    assert result["same_candidate_cohort_verified"] is False
    assert result["selection_effect"] == "invalid_comparison_cohort_mismatch"
    assert result["winner_metric_deltas"] == {}


def test_selection_comparison_rejects_different_candidate_cohorts():
    sequence_winner = make_candidate("BASELINE")
    target_winner = make_candidate("UNRELATED")
    result = p4.build_selection_comparison(
        [sequence_winner],
        [target_winner],
        {"pareto_reordered_candidates": True},
    )

    assert result["same_candidate_cohort_verified"] is False
    assert result["selection_effect"] == "invalid_comparison_cohort_mismatch"
    assert result["winner_metric_deltas"] == {}

def test_generic_hotspot_configuration():
    hotspots = p4._normalise_hotspots(
        [
            {
                "chain": "R",
                "residue_number": 42,
                "residue_name": "GLU",
            },
            {
                "chain": "R",
                "residue_number": 51,
                "residue_name": "ASP",
            },
        ],
        "R",
    )

    assert len(hotspots) == 2
    assert hotspots[0].label() == "R:42:GLU"
    assert hotspots[1].label() == "R:51:ASP"


def test_duplicate_hotspots_are_rejected():
    with pytest.raises(ValueError, match="Duplicate hotspot"):
        p4._normalise_hotspots(
            [
                {"chain": "R", "residue_number": 42,
                 "residue_name": "GLU"},
                {"chain": "R", "residue_number": 42,
                 "residue_name": "GLU"},
            ],
            "R",
        )


def test_manifest_matches_exact_candidate_sequence(tmp_path):
    sequence = "ACDEFGHIK"
    manifest, pose, digest = write_manifest(tmp_path, sequence)

    matching = make_candidate(sequence)
    unmatched = make_candidate("ACDEFGHIL")

    report = p4._register_manifest(
        [matching, unmatched],
        manifest,
    )

    assert report["matched_candidate_count"] == 1
    assert report["matched_complex_count"] == 1
    assert matching.metadata["target_complexes"][0]["sha256"] == digest
    assert Path(matching.metadata["target_complexes"][0]["path"]) == pose
    assert unmatched.metadata["target_complexes"] == []


def test_manifest_rejects_incorrect_sha256(tmp_path):
    sequence = "ACDEFGHIK"
    manifest, _, _ = write_manifest(
        tmp_path,
        sequence,
        expected_hash="0" * 64,
    )

    candidate = make_candidate(sequence)
    report = p4._register_manifest([candidate], manifest)

    assert report["matched_candidate_count"] == 0
    assert candidate.metadata["target_complexes"] == []
    assert len(report["rejected_manifest_records"]) == 1
    assert "SHA-256 mismatch" in (
        report["rejected_manifest_records"][0]["error"]
    )


def test_reliability_requires_independent_runs():
    entries = [
        {
            "independent_run_id": "run_1",
            "random_seed": 101,
            "protocol_sha256": "a" * 64,
            "_interface_fingerprint_5A": ["A"],
        },
        {
            "independent_run_id": "run_1",
            "random_seed": 101,
            "protocol_sha256": "a" * 64,
            "_interface_fingerprint_5A": ["A"],
        },
    ]

    result = p4._reliability_diagnostic(entries)

    assert result["status"] == "insufficient_independent_runs"
    assert result["calibrated"] is False
    assert result["feedback_active"] is False


def test_reliability_rejects_pose_without_run_id():
    entries = [
        {
            "independent_run_id": "run_1",
            "random_seed": 101,
            "protocol_sha256": "a" * 64,
            "_interface_fingerprint_5A": ["A"],
        },
        {
            "independent_run_id": "run_2",
            "random_seed": 202,
            "protocol_sha256": "a" * 64,
            "_interface_fingerprint_5A": ["A"],
        },
        {
            "independent_run_id": None,
            "random_seed": 303,
            "protocol_sha256": "a" * 64,
            "_interface_fingerprint_5A": ["A"],
        },
    ]
    result = p4._reliability_diagnostic(entries)
    assert result["status"] == "missing_independent_run_metadata"
    assert result["calibrated"] is False
    assert result["feedback_active"] is False


def test_reliability_rejects_none_or_malformed_protocol_hash():
    entries = []
    for run_id, seed in (("run_1", 101), ("run_2", 202)):
        for _ in range(2):
            entries.append({
                "independent_run_id": run_id,
                "random_seed": seed,
                "protocol_sha256": None,
                "_interface_fingerprint_5A": ["A"],
            })
    result = p4._reliability_diagnostic(entries)
    assert result["status"] == "protocol_mismatch_or_missing"
    assert result["calibrated"] is False


def test_reliability_rejects_missing_seed_on_any_pose():
    entries = []
    for run_id, seed in (("run_1", 101), ("run_2", 202)):
        for pose_index in range(2):
            entries.append({
                "independent_run_id": run_id,
                "random_seed": seed if pose_index == 0 else None,
                "protocol_sha256": "a" * 64,
                "_interface_fingerprint_5A": ["A"],
            })
    result = p4._reliability_diagnostic(entries)
    assert result["status"] == "invalid_run_seed_metadata"
    assert result["calibrated"] is False


def test_reliability_diagnostic_remains_uncalibrated(monkeypatch):
    def fake_cross_run_reliability(runs):
        assert len(runs) == 2
        assert all(len(run) == 2 for run in runs)
        return {
            "reliability_proxy": 0.8,
            "mean_cross_run_jaccard_similarity": 0.9,
        }

    monkeypatch.setattr(
        "prototype4.p4_closed_loop.cross_run_reliability",
        fake_cross_run_reliability,
    )

    entries = []

    for run_id, seed in (("run_1", 101), ("run_2", 202)):
        for pose_index in range(2):
            entries.append({
                "independent_run_id": run_id,
                "random_seed": seed,
                "protocol_sha256": "a" * 64,
                "_interface_fingerprint_5A": [
                    f"R:{pose_index}:GLU|P:1:ALA"
                ],
            })

    result = p4._reliability_diagnostic(entries)

    assert result["status"] == "diagnostic_only"
    assert result["calibrated"] is False
    assert result["feedback_active"] is False
    assert result["n_independent_runs"] == 2


def test_missing_target_evidence_is_reported_and_not_ranked(tmp_path, monkeypatch):
    monkeypatch.delenv(
        "SILICON_VIROVORE_P4_TARGET_MANIFEST",
        raising=False,
    )
    monkeypatch.delenv(
        "SILICON_VIROVORE_P4_TARGET_PARETO",
        raising=False,
    )

    candidate = make_candidate("ACDEFGHIK")

    ranked, report = p4.apply_p4_target_objective(
        [candidate],
        generation=1,
        output_directory=tmp_path,
        pareto_enabled=True,
    )

    assert report["target_evidence_coverage"] == 0.0
    assert report["pareto_reordered_candidates"] is False
    assert candidate.metadata["p4_target_evidence"]["status"] == (
        "missing_complex"
    )
    assert Path(report["report_path"]).is_file()
    assert ranked == [candidate]


def _write_pipeline_test_complex(tmp_path):
    """Create a valid, controlled ALA-receptor / GLY-peptide complex."""
    def atom(serial, residue, chain, residue_number, x, y, z):
        return (
            f"ATOM  {serial:5d} {'CA':4s} {residue:>3s} "
            f"{chain}{residue_number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            f"{1.00:6.2f}{20.00:6.2f}          {'C':>2s}\n"
        )

    lines = []
    positions = (
        (0, 0, 0),
        (0, 10, 0),
        (0, 0, 10),
    )

    serial = 1

    for residue_number, (x, y, z) in enumerate(positions, 1):
        lines.append(
            atom(serial, "ALA", "A", residue_number, x, y, z)
        )
        serial += 1

    for residue_number, (x, y, z) in enumerate(positions, 1):
        lines.append(
            atom(serial, "GLY", "B", residue_number, x + 3, y, z)
        )
        serial += 1

    lines.append("END\n")

    path = tmp_path / "candidate_complex.pdb"
    path.write_text("".join(lines), encoding="utf-8")
    return path


def _write_pipeline_manifest(tmp_path, pose, hotspot, target_id):
    digest = hashlib.sha256(pose.read_bytes()).hexdigest()

    manifest = {
        "schema_version": 1,
        "target": {
            "target_structure_id": target_id,
            "receptor_chain": "A",
            "peptide_chain": "B",
            "hotspots": [hotspot],
        },
        "complexes": [
            {
                "sequence": "GGG",
                "path": pose.name,
                "sha256": digest,
                "provenance": "controlled_integration_test",
            }
        ],
    }

    path = tmp_path / "target_manifest.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    return path


def test_end_to_end_custom_target_scoring_and_pareto(tmp_path, monkeypatch):
    pose = _write_pipeline_test_complex(tmp_path)

    manifest = _write_pipeline_manifest(
        tmp_path,
        pose,
        {
            "chain": "A",
            "residue_number": 1,
            "residue_name": "ALA",
        },
        "SYNTHETIC_TARGET_1",
    )

    monkeypatch.setenv(
        "SILICON_VIROVORE_P4_TARGET_MANIFEST",
        str(manifest),
    )

    candidate = make_candidate("GGG")

    ranked, report = p4.apply_p4_target_objective(
        [candidate],
        generation=1,
        output_directory=tmp_path / "run",
        pareto_enabled=True,
    )

    evidence = candidate.metadata["p4_target_evidence"]

    assert report["manifest"]["matched_candidate_count"] == 1
    assert report["target_evidence_coverage"] == 1.0
    assert report["pareto_reordered_candidates"] is True

    assert evidence["status"] == "success"
    assert evidence["target_structure_id"] == "SYNTHETIC_TARGET_1"
    assert evidence["pose_count"] == 1

    assert evidence["target_geometry_score_mean"] >= 0.0
    assert evidence["target_geometry_score_summary"]["count"] == 1
    assert evidence["hotspot_distance_summaries_A"]["A:1:ALA"]["count"] == 1
    assert evidence["hotspot_distances_A_mean"]["A:1:ALA"] == pytest.approx(
        3.0, abs=0.01
    )

    assert candidate.metadata["p4_reliability"]["status"] == (
        "missing_independent_run_metadata"
    )
    assert candidate.metadata["p4_reliability"]["calibrated"] is False
    assert ranked == [candidate]

    report_path = Path(report["report_path"])
    assert report_path.is_file()



def test_live_candidate_complexes_take_precedence_over_historical_manifest(
    tmp_path, monkeypatch
):
    """A stale manifest must not erase poses generated in the current run."""
    sequence = "GGG"
    (tmp_path / "historical").mkdir()
    (tmp_path / "live").mkdir()
    historical_pose = _write_pipeline_test_complex(tmp_path / "historical")
    live_pose = _write_pipeline_test_complex(tmp_path / "live")

    historical_manifest = _write_pipeline_manifest(
        tmp_path / "historical",
        historical_pose,
        {
            "chain": "A",
            "residue_number": 1,
            "residue_name": "ALA",
        },
        "SYNTHETIC_HISTORICAL_TARGET",
    )
    monkeypatch.setenv(
        "SILICON_VIROVORE_P4_TARGET_MANIFEST",
        str(historical_manifest),
    )

    candidate = make_candidate(sequence)
    candidate.metadata.update({
        "candidate_id": "live-candidate",
        "target_structure_id": "SYNTHETIC_LIVE_TARGET",
        "target_receptor_chain": "A",
        "target_peptide_chain": "B",
        "target_hotspots": [
            {"chain": "A", "residue_number": 1, "residue_name": "ALA"}
        ],
        "target_complexes": [{
            "path": str(live_pose.resolve()),
            "sha256": hashlib.sha256(live_pose.read_bytes()).hexdigest(),
            "target_structure_id": "SYNTHETIC_LIVE_TARGET",
            "receptor_chain": "A",
            "peptide_chain": "B",
            "provenance": "candidate_specific_HADDOCK3_prediction",
            "independent_run_id": "live-seed-1917",
            "random_seed": 1917,
            "protocol_sha256": "live-protocol-sha256",
        }],
    })

    scored_paths = []

    def fake_score(received_candidate, entry):
        assert received_candidate is candidate
        scored_paths.append(Path(entry["path"]).resolve())
        assert entry["target_structure_id"] == "SYNTHETIC_LIVE_TARGET"
        return {
            "path": str(Path(entry["path"]).resolve()),
            "target_structure_id": entry["target_structure_id"],
            "score_version": "handoff-regression-test",
            "target_geometry_score": 0.7,
            "hotspot_distances_A": {"A:1:ALA": 3.0},
            "hotspot_proximity_score": 0.8,
            "minimum_interchain_distance_A": 2.8,
            "severe_clash_count": 0,
            "interface_residue_pairs_5A": 5,
            "independent_run_id": entry.get("independent_run_id"),
            "random_seed": entry.get("random_seed"),
            "protocol_sha256": entry.get("protocol_sha256"),
            "_interface_fingerprint_5A": ["A:1:ALA|B:1:GLY"],
        }

    monkeypatch.setattr(p4, "_score_complex", fake_score)

    ranked, report = p4.apply_p4_target_objective(
        [candidate],
        generation=1,
        output_directory=tmp_path / "run",
        pareto_enabled=False,
    )

    assert scored_paths == [live_pose.resolve()]
    assert scored_paths[0] != historical_pose.resolve()
    assert report["manifest"]["live_metadata_candidate_count"] == 1
    assert report["manifest"]["live_metadata_complex_count"] == 1
    assert report["manifest"]["manifest_fallback_candidate_count"] == 0
    assert report["target_evidence_coverage"] == 1.0
    candidate_record = report["candidates"][0]
    assert candidate_record["target_complex_source"] == "live_candidate_metadata"
    assert candidate_record["registered_complex_count"] == 1
    assert candidate_record["complex_generation"] is None
    assert candidate.metadata["p4_target_evidence"]["status"] == "success"
    assert (
        candidate.metadata["p4_target_evidence"]["target_structure_id"]
        == "SYNTHETIC_LIVE_TARGET"
    )
    assert ranked == [candidate]


def test_target_switch_uses_new_hotspot_definition(tmp_path, monkeypatch):
    pose = _write_pipeline_test_complex(tmp_path)

    manifest = _write_pipeline_manifest(
        tmp_path,
        pose,
        {
            "chain": "A",
            "residue_number": 2,
            "residue_name": "ALA",
        },
        "SYNTHETIC_TARGET_2",
    )

    monkeypatch.setenv(
        "SILICON_VIROVORE_P4_TARGET_MANIFEST",
        str(manifest),
    )

    candidate = make_candidate("GGG")

    _, report = p4.apply_p4_target_objective(
        [candidate],
        generation=1,
        output_directory=tmp_path / "run",
        pareto_enabled=False,
    )

    evidence = candidate.metadata["p4_target_evidence"]

    assert report["target_evidence_coverage"] == 1.0
    assert evidence["status"] == "success"
    assert evidence["target_structure_id"] == "SYNTHETIC_TARGET_2"

    distances = evidence["hotspot_distances_A_mean"]

    assert set(distances) == {"A:2:ALA"}
    assert distances["A:2:ALA"] == pytest.approx(3.0, abs=0.01)



def test_strict_target_optimization_requires_repeatability_and_targets_parent_order(
    tmp_path, monkeypatch
):
    """Target objectives order the next parent only with full repeatability evidence."""
    scores = {
        "HIGH_FITNESS_WEAK_TARGET": (0.40, 0.40),
        "LOWER_FITNESS_STRONG_TARGET": (0.85, 0.90),
    }

    def fake_score(candidate, entry):
        geometry, proximity = scores[candidate.sequence]
        return {
            "path": entry["path"],
            "target_structure_id": "9SYA",
            "score_version": "test",
            "target_geometry_score": geometry,
            "hotspot_distances_A": {
                "A:203:ARG": 3.0,
                "A:205:ARG": 3.5,
            },
            "hotspot_proximity_score": proximity,
            "minimum_interchain_distance_A": 2.8,
            "severe_clash_count": 0,
            "interface_residue_pairs_5A": 12,
            "independent_run_id": entry["independent_run_id"],
            "random_seed": entry["random_seed"],
            "protocol_sha256": entry["protocol_sha256"],
            "_interface_fingerprint_5A": ["A:203:ARG|B:1:ALA"],
        }

    def fake_reliability(entries):
        return {
            "status": "diagnostic_only",
            "n_independent_runs": 2,
            "reliability_proxy": 0.60,
            "calibrated": False,
            "feedback_active": False,
        }

    monkeypatch.setattr(p4, "_score_complex", fake_score)
    monkeypatch.setattr(p4, "_reliability_diagnostic", fake_reliability)
    monkeypatch.delenv("SILICON_VIROVORE_P4_TARGET_MANIFEST", raising=False)

    high_fitness = make_candidate("HIGH_FITNESS_WEAK_TARGET")
    high_fitness.c_score = 0.95
    high_fitness.overall_score = 0.95
    low_fitness = make_candidate("LOWER_FITNESS_STRONG_TARGET")
    low_fitness.c_score = 0.80
    low_fitness.overall_score = 0.80

    for candidate in (high_fitness, low_fitness):
        candidate.metadata["target_structure_id"] = "9SYA"
        candidate.metadata["target_receptor_chain"] = "A"
        candidate.metadata["target_peptide_chain"] = "B"
        candidate.metadata["target_hotspots"] = [
            {"chain": "A", "residue_number": 203, "residue_name": "ARG"},
            {"chain": "A", "residue_number": 205, "residue_name": "ARG"},
        ]
        candidate.metadata["target_complexes"] = [
            {
                "path": str(tmp_path / f"{candidate.sequence}_seed1917" / "pose1.pdb"),
                "target_structure_id": "9SYA",
                "receptor_chain": "A",
                "peptide_chain": "B",
                "independent_run_id": "seed1917",
                "random_seed": 1917,
                "protocol_sha256": "a" * 64,
            },
            {
                "path": str(tmp_path / f"{candidate.sequence}_seed2917" / "pose1.pdb"),
                "target_structure_id": "9SYA",
                "receptor_chain": "A",
                "peptide_chain": "B",
                "independent_run_id": "seed2917",
                "random_seed": 2917,
                "protocol_sha256": "a" * 64,
            },
        ]

    ranked, report = p4.apply_p4_target_objective(
        [high_fitness, low_fitness],
        generation=1,
        output_directory=tmp_path / "run",
        pareto_enabled=True,
        require_reproducibility=True,
    )

    assert report["target_evidence_coverage"] == 1.0
    assert report["target_selection_coverage"] == 1.0
    assert report["reliability_used_for_ranking"] is False
    assert report["pareto_reordered_candidates"] is True
    assert report["selection_method"] == (
        "pareto_front_then_target_geometry_then_hotspot_proximity"
    )
    assert ranked[0] is low_fitness
    assert ranked[0].metadata["p4_target_selection_applied"] is True
    assert ranked[0].metadata["p4_reliability"]["calibrated"] is False


def test_strict_target_optimization_withholds_pareto_if_any_repeatability_is_missing(
    tmp_path, monkeypatch
):
    def fake_score(candidate, entry):
        return {
            "path": entry["path"],
            "target_structure_id": "9SYA",
            "score_version": "test",
            "target_geometry_score": 0.7,
            "hotspot_distances_A": {
                "A:203:ARG": 3.0,
                "A:205:ARG": 3.5,
            },
            "hotspot_proximity_score": 0.8,
            "minimum_interchain_distance_A": 2.8,
            "severe_clash_count": 0,
            "interface_residue_pairs_5A": 12,
            "independent_run_id": entry["independent_run_id"],
            "random_seed": entry["random_seed"],
            "protocol_sha256": entry["protocol_sha256"],
            "_interface_fingerprint_5A": ["A:203:ARG|B:1:ALA"],
        }

    def fake_reliability(entries):
        if "CANDIDATE_TWO" in entries[0]["path"]:
            return {
                "status": "insufficient_independent_runs",
                "n_independent_runs": 1,
                "calibrated": False,
                "feedback_active": False,
            }
        return {
            "status": "diagnostic_only",
            "n_independent_runs": 2,
            "reliability_proxy": 0.7,
            "calibrated": False,
            "feedback_active": False,
        }

    monkeypatch.setattr(p4, "_score_complex", fake_score)
    monkeypatch.setattr(p4, "_reliability_diagnostic", fake_reliability)
    monkeypatch.delenv("SILICON_VIROVORE_P4_TARGET_MANIFEST", raising=False)

    candidates = [make_candidate("CANDIDATE_ONE"), make_candidate("CANDIDATE_TWO")]
    for candidate in candidates:
        candidate.metadata["target_structure_id"] = "9SYA"
        candidate.metadata["target_receptor_chain"] = "A"
        candidate.metadata["target_peptide_chain"] = "B"
        candidate.metadata["target_hotspots"] = [
            {"chain": "A", "residue_number": 203, "residue_name": "ARG"},
            {"chain": "A", "residue_number": 205, "residue_name": "ARG"},
        ]
        candidate.metadata["target_complexes"] = [
            {
                "path": str(tmp_path / f"{candidate.sequence}_seed1917" / "pose1.pdb"),
                "target_structure_id": "9SYA",
                "receptor_chain": "A",
                "peptide_chain": "B",
                "independent_run_id": "seed1917",
                "random_seed": 1917,
                "protocol_sha256": "a" * 64,
            },
            {
                "path": str(tmp_path / f"{candidate.sequence}_seed2917" / "pose1.pdb"),
                "target_structure_id": "9SYA",
                "receptor_chain": "A",
                "peptide_chain": "B",
                "independent_run_id": "seed2917",
                "random_seed": 2917,
                "protocol_sha256": "a" * 64,
            },
        ]

    ranked, report = p4.apply_p4_target_objective(
        candidates,
        generation=1,
        output_directory=tmp_path / "run",
        pareto_enabled=True,
        require_reproducibility=True,
    )

    assert report["target_evidence_coverage"] == 1.0
    assert report["target_selection_coverage"] == 0.5
    assert report["pareto_reordered_candidates"] is False
    assert "repeatability diagnostics" in report["reason"]
    assert ranked == candidates
