import csv
import json
from pathlib import Path
from types import SimpleNamespace

from src.p4_pipeline_integration import (
    _summarize_pose_source_groups,
)
from src.p4_generation_audit import (
    write_p4_generation_assessment,
)


def test_pose_source_groups_are_descriptive_and_not_independent():
    poses = [
        {
            "path": "/tmp/chain_A/pose1.pdb",
            "target_geometry_score": 0.7,
            "hotspot_proximity_score": 0.2,
            "interface_residue_pairs_5A": 10,
            "minimum_interchain_distance_A": 3.0,
            "hotspot_distances_A": {
                "A:203:ARG": 4.0,
                "A:205:ARG": 8.0,
            },
        },
        {
            "path": "/tmp/chain_A/pose2.pdb",
            "target_geometry_score": 0.9,
            "hotspot_proximity_score": 0.6,
            "interface_residue_pairs_5A": 14,
            "minimum_interchain_distance_A": 5.0,
            "hotspot_distances_A": {
                "A:203:ARG": 6.0,
                "A:205:ARG": 10.0,
            },
        },
        {
            "path": "/tmp/chain_B/pose1.pdb",
            "target_geometry_score": 0.5,
            "hotspot_proximity_score": 0.3,
            "interface_residue_pairs_5A": 4,
            "minimum_interchain_distance_A": 2.8,
            "hotspot_distances_A": {
                "A:203:ARG": 7.0,
                "A:205:ARG": 9.0,
            },
        },
    ]

    result = _summarize_pose_source_groups(poses)

    assert result["grouping_basis"] == (
        "pose_file_parent_directory_name"
    )
    assert result["independence_claimed"] is False
    assert result["group_count"] == 2

    groups = {
        item["source_group"]: item
        for item in result["groups"]
    }

    assert groups["chain_A"]["pose_count"] == 2
    assert groups["chain_A"]["target_geometry_score"]["mean"] == 0.8
    assert groups["chain_A"]["hotspot_proximity_score"]["mean"] == 0.4
    assert (
        groups["chain_A"]["hotspot_distance_summaries_A"]
        ["A:203:ARG"]["mean"] == 5.0
    )

    assert groups["chain_B"]["pose_count"] == 1
    assert groups["chain_B"]["target_geometry_score"]["mean"] == 0.5


def test_integrated_json_and_csv_preserve_pose_source_groups(tmp_path):
    source_groups = _summarize_pose_source_groups([
        {
            "path": "/tmp/chain_A/pose1.pdb",
            "target_geometry_score": 0.8,
            "hotspot_proximity_score": 0.4,
            "interface_residue_pairs_5A": 12,
            "minimum_interchain_distance_A": 3.1,
            "hotspot_distances_A": {
                "A:203:ARG": 4.0,
                "A:205:ARG": 7.0,
            },
        }
    ])

    candidate = SimpleNamespace(
        rank=1,
        metadata={
            "candidate_id": "test-candidate",
            "p4_target_evidence": {
                "status": "success",
                "target_geometry_score_mean": 0.8,
                "hotspot_proximity_score": 0.4,
                "pose_source_groups": source_groups,
            },
            "p4_reliability": {
                "status": "insufficient_independent_runs",
                "calibrated": False,
                "feedback_active": False,
            },
        },
    )

    report = write_p4_generation_assessment(
        [candidate],
        {
            "generation": 1,
            "manifest": {"status": "loaded"},
            "pareto_reordered_candidates": False,
            "legacy_fitness_scores_modified": False,
            "reliability_diagnostic_candidate_count": 0,
        },
        tmp_path,
    )

    assert (
        report["candidates"][0]["pose_source_groups"]
        == source_groups
    )

    csv_path = (
        tmp_path
        / "p4_integrated_assessment"
        / "generation_0001.csv"
    )

    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 1
    csv_groups = json.loads(rows[0]["pose_source_groups"])
    assert csv_groups == source_groups
