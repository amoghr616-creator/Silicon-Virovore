import json
from pathlib import Path

from scripts.p4_matched_experiment import make_pool, parse_seeds, plan
from scripts.p4_target_coverage_audit import audit_run


def test_make_pool_is_deterministic_unique_and_seed_specific():
    first = make_pool(631, 30, "MKLAVFALLVFFAGSSDLIRR")
    repeat = make_pool(631, 30, "MKLAVFALLVFFAGSSDLIRR")
    other = make_pool(632, 30, "MKLAVFALLVFFAGSSDLIRR")
    assert first == repeat
    assert len(first) == len(set(first)) == 30
    assert first != other
    assert first[0] == "MKLAVFALLVFFAGSSDLIRR"


def test_parse_seeds_supports_ranges_and_lists():
    assert parse_seeds("631-633") == [631, 632, 633]
    assert parse_seeds("631,635,640") == [631, 635, 640]


def test_plan_creates_shared_pools_and_syntax_checked_launcher(tmp_path):
    import argparse

    out = tmp_path / "matched"
    args = argparse.Namespace(
        seeds="631,632",
        experiment_dir=out,
        population_size=20,
        generations=2,
        shortlist_size=4,
        dock_seeds=(1917, 2917),
        mutation_mode="hotspot_entropy_substitution_epistasis_exploratory_v2_matched",
        ncores=2,
        seed_sequence="MKLAVFALLVFFAGSSDLIRR",
    )
    assert plan(args) == 0
    manifest = json.loads((out / "experiment_plan.json").read_text())
    assert manifest["seeds"] == [631, 632]
    assert manifest["sequence_evaluations_per_run_expected"] == 40
    assert manifest["arms"] == ["P3", "P4-TA"]
    assert len(manifest["pool_records"]) == 2
    launcher = (out / "launch_matched_p3_p4ta.sh").read_text()
    assert "--p4-exploratory-target-optimization" in launcher
    assert "SILICON_VIROVORE_P4_DOCK_SEEDS" in launcher
    assert "SILICON_VIROVORE_INITIAL_POOL_FASTA" in launcher
    assert "P3_seed631" in launcher and "P4_TA_seed631" in launcher
    assert (out / "run_ledger.tsv").is_file()


def test_target_coverage_audit_preserves_candidate_level_failures(tmp_path):
    run_dir = tmp_path / "run"
    evidence_dir = run_dir / "p4_target_evidence"
    evidence_dir.mkdir(parents=True)
    report = {
        "generation": 1,
        "pareto_reordered_candidates": False,
        "candidates": [
            {
                "candidate_id": "ok",
                "sequence": "MKLAVFALLVFFAGSSDLIRR",
                "structure_shortlist": {"selected": True, "rank": 1},
                "complex_generation": {"status": "success"},
                "registered_complex_count": 10,
                "status": "success",
                "target_geometry_score_mean": 0.7,
                "hotspot_proximity_score": 0.8,
                "reliability": {"status": "diagnostic_only", "n_independent_runs": 2},
                "target_selection_eligible": True,
            },
            {
                "candidate_id": "miss",
                "sequence": "MKLAVFALLVFFAGSSDLIRW",
                "structure_shortlist": {"selected": True, "rank": 2},
                "complex_generation": {"status": "no_valid_complexes"},
                "registered_complex_count": 0,
                "status": "missing_complex",
                "reliability": {"status": "unavailable_without_valid_target_evidence"},
                "target_selection_eligible": False,
                "target_selection_reason": "no valid structures",
            },
        ],
    }
    (evidence_dir / "generation_0001.json").write_text(json.dumps(report))
    summary = audit_run(run_dir, tmp_path / "out")
    assert summary["shortlist_candidate_records"] == 2
    assert summary["valid_target_evidence_records"] == 1
    assert summary["repeatability_eligible_records"] == 1
    candidates = (tmp_path / "out" / "p4_target_coverage_candidates.csv").read_text()
    assert "missing_complex" in candidates
