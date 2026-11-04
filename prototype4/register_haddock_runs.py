#!/usr/bin/env python3
import hashlib
import json
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / "prototype4/data/p4_closed_loop_registry.sqlite"
AUDIT = ROOT / "prototype4/results/haddock_pose_audit_20261031"

STAGES = [
    "0_topoaa", "1_rigidbody", "2_seletop",
    "3_flexref", "4_emref", "5_clustfcc", "6_seletopclusts",
]

def file_hash(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def manifest_hash(items):
    payload = "\n".join(
        f"{item['path']}\t{item['sha256']}" for item in sorted(
            items, key=lambda x: x["path"]
        )
    )
    return hashlib.sha256(payload.encode()).hexdigest()

def artifact_list(paths):
    items = []
    for path in sorted(set(paths)):
        if not path.is_file():
            raise FileNotFoundError(path)
        items.append({
            "path": str(path.relative_to(ROOT)),
            "sha256": file_hash(path),
            "size_bytes": path.stat().st_size,
        })
    if not items:
        raise RuntimeError("Empty artifact list")
    return items

def register(con, chain):
    run_root = ROOT / f"prototype4/results/haddock3_chain_{chain}_pilot_20261031"
    receptor = run_root / "data/0_topoaa" / f"9sya_source_chain_{chain}.pdb"
    peptide = run_root / "data/0_topoaa/clean_peptide_ensemble.pdb"

    restraint_root = ROOT / "prototype4/data/haddock3/restraints"
    restraints = [
        restraint_root / f"chain_{chain}_receptor.act-pass",
        restraint_root / f"chain_{chain}_peptide.act-pass",
        restraint_root / f"chain_{chain}_arg203_arg205.tbl",
    ]

    source_receptor = (
        ROOT / f"prototype4/data/haddock3/receptors/9sya_source_chain_{chain}.pdb"
    )

    inputs = artifact_list([source_receptor, receptor, peptide, *restraints])

    config_paths = [
        run_root / "data/configurations/enhanced_haddock_params.json"
    ]
    for stage in STAGES:
        stage_root = run_root / stage
        config_paths += list(stage_root.glob("params.cfg"))
        config_paths += list(stage_root.glob("io.json"))

    configs = artifact_list(config_paths)

    output_paths = []
    for stage in ("1_rigidbody", "3_flexref", "4_emref", "6_seletopclusts"):
        output_paths += list((run_root / stage).glob("*.pdb.gz"))

    audit_poses = sorted(
        (AUDIT / f"chain_{chain}").glob("cluster_*_model_1.pdb")
    )
    outputs = artifact_list(output_paths + audit_poses)

    counts = {
        stage: len(list((run_root / stage).glob("*.pdb.gz")))
        for stage in ("1_rigidbody", "3_flexref", "4_emref", "6_seletopclusts")
    }
    counts["audited_pose_files"] = len(audit_poses)

    if counts["1_rigidbody"] != 60:
        raise RuntimeError(f"Expected 60 rigidbody outputs for chain {chain}: {counts}")
    if counts["3_flexref"] != 5 or counts["4_emref"] != 5:
        raise RuntimeError(f"Unexpected refinement output counts: {counts}")
    if counts["6_seletopclusts"] != 5 or counts["audited_pose_files"] != 5:
        raise RuntimeError(f"Unexpected selected pose counts: {counts}")

    config_sha = manifest_hash(configs)
    input_sha = manifest_hash(inputs)
    output_sha = manifest_hash(outputs)

    run_id = f"haddock3_9sya_chain_{chain}_20261031_{config_sha[:12]}"

    metadata = {
        "target_structure_id": "9SYA",
        "receptor_chain": chain,
        "peptide_chain": "B",
        "hotspots": ["ARG:203", "ARG:205"],
        "software_version_source": "HADDOCK3 environment version previously verified",
        "run_root": str(run_root.relative_to(ROOT)),
        "workflow_stages": STAGES,
        "stage_output_counts": counts,
        "independent_execution_count": 1,
        "selected_pose_count": len(audit_poses),
        "pose_independence_warning":
            "Five selected poses belong to this single HADDOCK execution; "
            "they are not independent runs or independent trajectories.",
        "random_seed": None,
        "random_seed_note": "Not recorded in the inspected provenance artifacts.",
        "start_time": None,
        "completion_time": None,
        "timing_note": "Exact timestamps not recovered from preserved artifacts.",
        "rmsd_note": "No reference-based L-RMSD recorded; database l_rmsd is NULL.",
        "input_artifacts": inputs,
        "configuration_artifacts": configs,
        "output_artifacts": outputs,
        "config_sha256": config_sha,
        "input_sha256": input_sha,
        "output_sha256": output_sha,
    }

    values = (
        run_id,
        None,
        f"9SYA:chain_{chain}:ARG203_ARG205",
        f"single_HADDOCK3_execution_chain_{chain}",
        "HADDOCK3",
        "2026.10.1",
        "HADDOCK3 workflow; original shell invocation not preserved",
        config_sha,
        None,
        input_sha,
        output_sha,
        None,
        None,
        "completed",
        None,
        None,
        json.dumps(metadata, sort_keys=True, separators=(",", ":")),
    )

    existing = con.execute(
        "SELECT * FROM structure_runs WHERE run_id=?", (run_id,)
    ).fetchone()

    if existing:
        columns = [r[1] for r in con.execute(
            "PRAGMA table_info(structure_runs)"
        )]
        observed = tuple(existing[c] for c in columns)
        if observed != values:
            raise RuntimeError(f"Existing run record conflicts: {run_id}")
        print(f"UNCHANGED {run_id}")
    else:
        con.execute(
            """
            INSERT INTO structure_runs (
                run_id, cycle_uid, target_key, run_role, software,
                software_version, command_string, config_sha256,
                random_seed, input_sha256, output_sha256, trajectory_path,
                l_rmsd, status, started_at, completed_at, metadata_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            values,
        )
        print(f"REGISTERED {run_id}")

    print(f"  Input artifacts: {len(inputs)}")
    print(f"  Config artifacts: {len(configs)}")
    print(f"  Output artifacts: {len(outputs)}")
    print(f"  Output counts: {counts}")
    print(f"  Config SHA-256: {config_sha}")
    print(f"  Input SHA-256:  {input_sha}")
    print(f"  Output SHA-256: {output_sha}")

def main():
    if not DB.is_file():
        raise FileNotFoundError(DB)

    con = sqlite3.connect(DB)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA foreign_keys=ON")
    try:
        con.execute("BEGIN IMMEDIATE")
        register(con, "A")
        register(con, "B")

        errors = con.execute("PRAGMA foreign_key_check").fetchall()
        if errors:
            raise RuntimeError(f"Foreign-key violations: {errors[:5]}")

        con.commit()
        print("\nSTRUCTURE RUN REGISTRATION PASSED")
        print("Registered HADDOCK executions:",
              con.execute("SELECT COUNT(*) FROM structure_runs").fetchone()[0])
    except Exception:
        con.rollback()
        raise
    finally:
        con.close()

if __name__ == "__main__":
    main()
