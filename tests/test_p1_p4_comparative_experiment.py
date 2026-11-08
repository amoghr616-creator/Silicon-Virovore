import argparse
import hashlib
import json
from pathlib import Path

from scripts import p1_p4_comparative_experiment as experiment


def _historical_manifest(tmp_path: Path) -> Path:
    complex_path = tmp_path / "historical_complex.pdb"
    complex_path.write_text("REMARK synthetic test fixture only\n", encoding="utf-8")
    digest = hashlib.sha256(complex_path.read_bytes()).hexdigest()
    payload = {
        "schema_version": 1,
        "target": {
            "target_structure_id": "9SYA",
            "receptor_chain": "A",
            "peptide_chain": "B",
            "hotspots": [
                {"chain": "A", "residue_number": 203, "residue_name": "ARG"},
                {"chain": "A", "residue_number": 205, "residue_name": "ARG"},
            ],
        },
        "complexes": [{
            "sequence": "MKLAVEALLVFFAGSSDLIRR",
            "path": str(complex_path),
            "sha256": digest,
        }],
    }
    path = tmp_path / "historical_pose_manifest.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_matched_pool_is_seed_deterministic_and_contains_registered_sequences():
    seed = "MKLAVFALLVFFAGSSDLIRR"
    historical = "MKLAVEALLVFFAGSSDLIRR"
    first = experiment.generate_matched_pool(631, 100, seed, historical)
    repeat = experiment.generate_matched_pool(631, 100, seed, historical)
    other = experiment.generate_matched_pool(632, 100, seed, historical)
    assert first == repeat
    assert first != other
    assert len(first) == len(set(first)) == 100
    assert first[:2] == [seed, historical]


def test_plan_freezes_shared_pools_and_registers_eight_runs(tmp_path, monkeypatch):
    manifest = _historical_manifest(tmp_path)
    receptor = tmp_path / "receptor.pdb"
    receptor.write_text("REMARK synthetic receptor fixture\n", encoding="utf-8")
    haddock_config = tmp_path / "haddock.cfg"
    haddock_config.write_text("# synthetic config fixture\n", encoding="utf-8")
    monkeypatch.setattr(experiment, "DEFAULT_RECEPTOR", receptor)
    monkeypatch.setattr(experiment, "DEFAULT_HADDOCK_CONFIG", haddock_config)
    out = tmp_path / "experiment"
    args = argparse.Namespace(
        seeds="631,632",
        dock_seeds="1917,2917",
        experiment_dir=out,
        population_size=100,
        generations=2,
        shortlist_size=2,
        refill_cap=2,
        haddock_ncores=4,
        seed_sequence="MKLAVFALLVFFAGSSDLIRR",
        historical_manifest=manifest,
        initial_pool_dir=tmp_path / "no_existing_pools",
        external_validation_summary=tmp_path / "no_external_summary.json",
    )
    assert experiment.plan_experiment(args) == 0
    plan = json.loads((out / "experiment_plan.json").read_text())
    assert plan["planned_sequence_evaluations_total"] == 1600
    assert plan["planned_haddock_execution_slots_total"] == 16
    assert plan["planned_esmfold_candidate_attempts_max"] == 16
    assert len(plan["pool_records"]) == 2
    target_snapshot = plan["historical_target_manifest"]["snapshot"]
    assert target_snapshot["status"] == "snapshotted_and_hash_verified"
    assert Path(target_snapshot["snapshot_manifest_path"]).is_file()
    assert len(target_snapshot["complex_records"]) == 1
    assert Path(target_snapshot["complex_records"][0]["snapshot_path"]).is_file()
    assert len(list((out / "run_manifests").glob("*.json"))) == 8
    hashes = {row["sha256"] for row in plan["pool_records"]}
    assert len(hashes) == 2  # Seed-specific pools differ; all four arms share each one.
    for seed in plan["seeds"]:
        rows = [row for row in json.loads((out / "experiment_plan.json").read_text())["pool_records"] if row["seed"] == seed]
        assert len(rows) == 1
        assert experiment.sha256_file(Path(rows[0]["path"])) == rows[0]["sha256"]
    assert not (out / "LAUNCH_STARTED").exists()


def test_arm_definitions_encode_valid_matched_proposal_modes():
    from dataclasses import replace
    from src.config import PipelineSettings, validate_configuration

    p1 = experiment.ARMS["P1"]["mutation_guidance_mode"]
    p2 = experiment.ARMS["P2"]["mutation_guidance_mode"]
    p3 = experiment.ARMS["P3"]["mutation_guidance_mode"]
    base_settings = PipelineSettings()
    for mode in (p1, p2, p3):
        validate_configuration(replace(base_settings, mutation_guidance_mode=mode))
    assert experiment.ARMS["P1"]["adaptive_enabled"] is False
    assert experiment.ARMS["P2"]["adaptive_enabled"] is True
    assert experiment.ARMS["P3"]["structural_assessment"] is True
    assert experiment.ARMS["P3"]["target_ordering"] is False
    assert experiment.ARMS["P4-TA"]["target_ordering"] is True


def test_p3_and_p4_get_same_structure_budget_but_only_p4_orders(tmp_path):
    manifest = _historical_manifest(tmp_path)
    plan_data = {
        "seed_sequence": "MKLAVFALLVFFAGSSDLIRR",
        "population_size": 100,
        "generations": 2,
        "shortlist_size": 2,
        "refill_cap": 2,
        "haddock_seeds": [1917, 2917],
        "haddock_ncores": 4,
        "historical_target_manifest": {"path": str(manifest), "sha256": experiment.sha256_file(manifest)},
        "p4_gate_summary_path": str(tmp_path / "NO_REGISTERED_SUMMARY.json"),
        "pool_records": [{"seed": 631, "path": str(tmp_path / "seed_631.fasta"), "sha256": "a" * 64}],
    }
    pool = Path(plan_data["pool_records"][0]["path"])
    pool.write_text(">a\nMKLAVFALLVFFAGSSDLIRR\n", encoding="utf-8")
    plan_data["pool_records"][0]["sha256"] = experiment.sha256_file(pool)
    p3_command, p3_env, _, _, _ = experiment._build_command(plan_data, tmp_path, 631, "P3")
    p4_command, p4_env, _, _, _ = experiment._build_command(plan_data, tmp_path, 631, "P4-TA")
    assert p3_command[p3_command.index("--seed-sequence") + 1] == plan_data["seed_sequence"]
    assert p4_command[p4_command.index("--seed-sequence") + 1] == plan_data["seed_sequence"]
    assert p4_env["SILICON_VIROVORE_P4_EXTERNAL_VALIDATION_SUMMARY"] == plan_data["p4_gate_summary_path"]
    for env in (p3_env, p4_env):
        assert env["SILICON_VIROVORE_P4_GENERATE_COMPLEXES"] == "1"
        assert env["SILICON_VIROVORE_P4_STRUCTURAL_GENERATION"] == "1"
        assert env["SILICON_VIROVORE_P4_STRUCTURE_LIMIT"] == "2"
        assert env["SILICON_VIROVORE_P4_MAX_DOCK_CANDIDATES"] == "2"
        assert env["SILICON_VIROVORE_P4_DOCK_SEEDS"] == "1917,2917"
        assert env["SILICON_VIROVORE_P4_TARGET_MANIFEST"] == str(manifest)
    assert p3_env["SILICON_VIROVORE_P4_TARGET_PARETO"] == "0"
    assert p4_env["SILICON_VIROVORE_P4_TARGET_PARETO"] == "1"
    assert "--no-arise-memory-enabled" in p4_command
    assert "--no-arise-memory-use-history" in p4_command
    memory_path = Path(p4_command[p4_command.index("--arise-memory-path") + 1])
    assert memory_path.name == "P4_TA_seed631.json"
    assert p4_env["SILICON_VIROVORE_P4_STRUCTURAL_GENERATION"] == "1"


def test_exact_sign_flip_is_descriptive_for_two_matched_seeds():
    assert experiment._exact_two_sided_signflip([0.2, 0.3]) == 0.5
    assert experiment._exact_two_sided_signflip([]) is None


def test_pipeline_uses_shared_pool_only_for_generation_one():
    source = (experiment.ROOT / "run_pipeline.py").read_text(encoding="utf-8")
    assert "initial_pool_sequences if generation == 1 else None" in source
    assert "proposal_operator_active" in source


def test_pipeline_patch_declares_generation_scoping_and_nonfatal_abstention():
    source = (experiment.ROOT / "run_pipeline.py").read_text(encoding="utf-8")
    assert "SILICON_VIROVORE_P4_STRUCTURAL_GENERATION" in source
    assert "p4_structural_generation_active = p4_generation_requested" in source
    assert "P4 target-selection abstention" in source
    assert "--arise-memory-path" in source


def test_external_validation_snapshot_freezes_summary_folds_and_artifacts(tmp_path):
    source_dir = tmp_path / "external_source"
    source_dir.mkdir()
    source_csv = source_dir / "data.csv"
    source_csv.write_text("group,measurement\nG1,1.0\n", encoding="utf-8")
    import hashlib
    source_hash = hashlib.sha256(source_csv.read_bytes()).hexdigest()
    predictions = source_dir / "predictions.csv"
    predictions.write_text("group,p3,p4\nG1,0.1,0.2\n", encoding="utf-8")
    summary = {
        "source": str(source_csv),
        "source_sha256": source_hash,
        "rows": 1,
        "unique_peptide_sequences": 1,
        "independent_binding_groups": 1,
        "binding_groups": ["G1"],
        "validation": "leave-one-Binding-Group-out",
        "folds": [{"held_out_binding_group": "G1"}],
        "artifacts": {"predictions": "predictions.csv"},
    }
    summary_path = source_dir / "summary.json"
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    snapshot = experiment.snapshot_external_validation(summary_path, tmp_path / "frozen")
    assert snapshot["status"] == "summary_snapshotted"
    assert snapshot["fold_count"] == 1
    assert snapshot["source_dataset_sha256"] == source_hash
    assert len(snapshot["copied_artifacts"]) == 1
    assert Path(snapshot["summary_path"]).is_file()
    assert Path(snapshot["source_dataset_path"]).is_file()


def test_missing_external_summary_is_recorded_as_unavailable_without_fake_metrics(tmp_path):
    snapshot = experiment.snapshot_external_validation(
        tmp_path / "missing_summary.json", tmp_path / "external_snapshot"
    )
    assert snapshot["status"] == "not_registered_at_plan_time"
    sentinel = Path(snapshot["summary_path"])
    assert sentinel.is_file()
    payload = json.loads(sentinel.read_text())
    assert payload["status"] == "not_registered_at_plan_time"
    assert "paired_primary_effect" not in payload


def test_target_coverage_requires_exact_sequence_and_manifest_verification(tmp_path):
    pool = tmp_path / "seed_631.fasta"
    pool.write_text(">seed\nMKLAVFALLVFFAGSSDLIRR\n", encoding="utf-8")
    pool_hash = hashlib.sha256(pool.read_bytes()).hexdigest()
    experiment_dir = tmp_path / "experiment_exact"
    run_id = "P3_seed631"
    run_dir = experiment_dir / "runs" / run_id
    target_dir = run_dir / "p4_target_evidence"
    target_dir.mkdir(parents=True)
    (target_dir / "generation_0001.json").write_text(json.dumps({
        "candidates": [
            {"candidate_id": "verified", "sequence": "MKLAVEALLVFFAGSSDLIRR",
             "status": "success", "structure_shortlist": {"selected": True},
             "complex_generation": {"sequence_verified": True,
                                    "execution_manifests_verified": True,
                                    "status": "success", "independent_run_count": 2},
             "target_selection_eligible": True,
             "target_selection_reason": "diagnostic present", "reliability": {}},
            {"candidate_id": "unverified", "sequence": "MKLAVFALLVFFAGSSDLIRR",
             "status": "success", "structure_shortlist": {"selected": True},
             "complex_generation": {"sequence_verified": True,
                                    "execution_manifests_verified": False,
                                    "status": "success", "independent_run_count": 2},
             "target_selection_eligible": False,
             "target_selection_reason": "manifest verification failed", "reliability": {}},
        ],
        "pareto_reordered_candidates": False,
    }), encoding="utf-8")
    record = experiment._load_run_data(
        experiment_dir,
        {"pool_records": [{"seed": 631, "path": str(pool), "sha256": pool_hash}]},
        631, "P3",
    )
    evidence = record["p4_generation1"]
    assert evidence["shortlist_count"] == 2
    assert evidence["valid_target_evidence_count"] == 1
    assert evidence["exact_sequence_evidence_coverage"] == 0.5
    assert evidence["abstention_count"] == 1
    assert record["target_candidates"][1]["execution_manifests_verified"] is False


def test_analysis_preserves_structural_failures_even_without_final_run_audit(tmp_path):
    pool = tmp_path / "seed_631.fasta"
    pool.write_text(">seed\nMKLAVFALLVFFAGSSDLIRR\n", encoding="utf-8")
    import hashlib
    pool_hash = hashlib.sha256(pool.read_bytes()).hexdigest()
    experiment_dir = tmp_path / "experiment"
    run_id = "P4_TA_seed631"
    run_dir = experiment_dir / "runs" / run_id
    run_dir.mkdir(parents=True)
    manifests = experiment_dir / "run_manifests"
    manifests.mkdir(parents=True)
    (manifests / f"{run_id}.json").write_text(
        json.dumps({"status": "FAILED", "exit_code": 1}), encoding="utf-8"
    )
    refill_dir = run_dir / "p4_structure_refill"
    refill_dir.mkdir()
    (refill_dir / "generation_0001.json").write_text(json.dumps({
        "refill_attempted_count": 1,
        "unattempted_queue_count": 1,
        "candidates": [
            {"queue_rank": 1, "role": "primary", "candidate_id": "c1",
             "attempted": True, "status": "http_504_failed", "structure_available": False,
             "error": "HTTP 504", "selected": True, "selection_rank": 1},
            {"queue_rank": 2, "role": "refill", "candidate_id": "c2",
             "attempted": False, "status": "not_attempted_quota_filled",
             "structure_available": False, "selected": False},
        ],
    }), encoding="utf-8")
    target_dir = run_dir / "p4_target_evidence"
    target_dir.mkdir()
    (target_dir / "generation_0001.json").write_text(json.dumps({
        "candidates": [{
            "candidate_id": "c1", "sequence": "MKLAVEALLVFFAGSSDLIRR",
            "status": "missing_complex", "structure_shortlist": {"selected": True},
            "complex_generation": {"status": "failed", "sequence_verified": False,
                                   "execution_manifests_verified": False},
            "target_selection_eligible": False,
            "target_selection_reason": "Target geometry evidence unavailable or invalid",
            "reliability": {},
        }],
        "pareto_reordered_candidates": False,
        "reason": "abstained",
    }), encoding="utf-8")

    plan = {
        "pool_records": [{"seed": 631, "path": str(pool), "sha256": pool_hash}],
    }
    record = experiment._load_run_data(experiment_dir, plan, 631, "P4-TA")
    assert record["run_status"] == "FAILED"
    assert len(record["structure_attempts"]) == 2
    assert record["structure_attempts"][0]["status"] == "http_504_failed"
    assert record["p4_generation1"]["abstention_count"] == 1
    assert record["target_candidates"][0]["abstention_reason"] == "Target geometry evidence unavailable or invalid"


def test_scored_evaluation_count_uses_finite_native_fitness_values():
    history = [
        {"native_fitness": 0.25},
        {"native_fitness": "0.75"},
        {"native_fitness": None},
        {"native_fitness": float("nan")},
        {"native_fitness": float("inf")},
        {"native_fitness": "not-a-score"},
        "malformed-row",
    ]
    assert experiment._count_scored_native_fitness(history) == 2
