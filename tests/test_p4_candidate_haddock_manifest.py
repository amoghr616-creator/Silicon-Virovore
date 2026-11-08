import hashlib
import json
from pathlib import Path

from src.p4_candidate_haddock import _valid_execution_manifest


def test_execution_manifest_is_reused_only_when_archive_and_protocol_match(tmp_path):
    archive = tmp_path / "summary.tgz"
    archive.write_bytes(b"fake archive bytes")
    payload = {
        "schema_version": 1,
        "status": "success",
        "execution_id": "exec-unique-1",
        "independent_run_id": "candidate:HADDOCK3:seed1917",
        "random_seed": 1917,
        "protocol_sha256": "a" * 64,
        "receptor_sha256": "b" * 64,
        "peptide_input_sha256": "c" * 64,
        "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
        "return_code": 0,
        "verified_pose_sha256": {"cluster_1/model_1.pdb": "d" * 64},
    }
    manifest = tmp_path / "execution_manifest.json"
    manifest.write_text(json.dumps(payload))

    recovered = _valid_execution_manifest(
        manifest,
        seed=1917,
        protocol_sha256="a" * 64,
        receptor_sha256="b" * 64,
        peptide_sha256="c" * 64,
        archive_path=archive,
    )
    assert recovered["execution_id"] == "exec-unique-1"

    assert _valid_execution_manifest(
        manifest,
        seed=2917,
        protocol_sha256="a" * 64,
        receptor_sha256="b" * 64,
        peptide_sha256="c" * 64,
        archive_path=archive,
    ) is None

    invalid_poses = dict(payload)
    invalid_poses["verified_pose_sha256"] = {}
    manifest.write_text(json.dumps(invalid_poses))
    assert _valid_execution_manifest(
        manifest,
        seed=1917,
        protocol_sha256="a" * 64,
        receptor_sha256="b" * 64,
        peptide_sha256="c" * 64,
        archive_path=archive,
    ) is None

    invalid_return = dict(payload)
    invalid_return["return_code"] = 1
    manifest.write_text(json.dumps(invalid_return))
    assert _valid_execution_manifest(
        manifest,
        seed=1917,
        protocol_sha256="a" * 64,
        receptor_sha256="b" * 64,
        peptide_sha256="c" * 64,
        archive_path=archive,
    ) is None

    manifest.write_text(json.dumps(payload))
    archive.write_bytes(b"changed archive")
    assert _valid_execution_manifest(
        manifest,
        seed=1917,
        protocol_sha256="a" * 64,
        receptor_sha256="b" * 64,
        peptide_sha256="c" * 64,
        archive_path=archive,
    ) is None
