"""Run Prototype 2 docking only after receptor identity preflight.

Outputs are written to ``p2_rapidock_validated`` so historical runs under
``p2_rapidock`` remain untouched. Each completed candidate receives an
independent structural-validation JSON file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

try:
    from .clean_9sya_receptor import prepare_receptor
    from .target_audit import (
        TargetMappingError,
        audit_residue_mapping,
        coordinate_frame_audit,
        diagnose_pose,
    )
except ImportError:  # direct script execution
    from clean_9sya_receptor import prepare_receptor
    from target_audit import TargetMappingError, audit_residue_mapping, coordinate_frame_audit, diagnose_pose


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data" / "receptor" / "9SYA.pdb"
DEFAULT_RECEPTOR = ROOT / "prototype2" / "target_aware_audit" / "9SYA_chainA_repaired.pdb"
DEFAULT_OUTROOT = ROOT / "prototype2" / "p2_rapidock_validated"
SITE = (-8.405, -5.823, 1.523)

CANDIDATES = {
    "P2-01": "MKQAVNALLVFFAGSSDAIRR",
    "P2-02": "DKLAVTALLVVFAESSDKLRR",
    "P2-03": "MKLAPFALLVVFAGASDWIRR",
    "P2-04": "MKLAVFALLVCFAESSDLVRR",
    "P2-05": "MKQAVFALLQVFAGSSDWIRR",
    "P2-06": "MKLAVFALRVFFAESSDAIRR",
}


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _source_version(source_root: Path) -> str | None:
    pyproject = source_root / "pyproject.toml"
    if not pyproject.exists():
        return None
    match = re.search(r'^version\s*=\s*["\']([^"\']+)["\']', pyproject.read_text(errors="replace"), re.MULTILINE)
    return match.group(1) if match else None


def _backend_preflight(
    sequence: str,
    box: float,
    n_samples: int,
    *,
    receptor: Path,
    output_root: Path,
    seed: int | None = None,
) -> dict[str, object]:
    """Check every dependency and input needed before launching docking.

    This is intentionally stricter than the backend's own length router: a
    21-mer may not silently fall back to a different checkpoint in a
    reproducibility run.
    """

    source_root = ROOT / "hybridock-pep"
    model_dir = source_root / "third_party" / "RAPiDock" / "train_models" / "CGTensorProductEquivariantModel"
    threshold = 13
    longer = model_dir / "longer_local.pt"
    local = model_dir / "rapidock_local.pt"
    executable = shutil.which("hybridock-pep")
    errors: list[str] = []
    warnings: list[str] = []
    help_text = ""
    version_probe: dict[str, object] = {"status": "not_run"}
    if executable is None:
        errors.append("hybridock-pep is not available on PATH; install or activate the intended backend environment")
        warnings.append("backend executable/version/argument support are unresolved because no executable was found")
    else:
        try:
            version = subprocess.run(
                [executable, "--version"], capture_output=True, text=True, timeout=10, check=False,
            )
            version_probe = {
                "status": "completed",
                "returncode": version.returncode,
                "stdout": version.stdout.strip(),
                "stderr": version.stderr.strip(),
            }
            if version.returncode != 0:
                warnings.append(
                    "backend does not expose a successful --version response; "
                    "source package version is recorded instead"
                )
        except (OSError, subprocess.TimeoutExpired) as exc:
            version_probe = {"status": "failed", "error": str(exc)}
            warnings.append(f"backend version probe failed: {exc}")
        try:
            help_result = subprocess.run(
                [executable, "dock", "--help"], capture_output=True, text=True, timeout=10, check=False,
            )
            help_text = (help_result.stdout or "") + "\n" + (help_result.stderr or "")
            required_flags = (
                "--peptide", "--receptor", "--site", "--box", "--output-dir",
                "--n-samples", "--seed",
            )
            missing_flags = [flag for flag in required_flags if flag not in help_text]
            if help_result.returncode != 0 or missing_flags:
                errors.append(
                    "backend dock help does not confirm the required arguments: "
                    + (", ".join(missing_flags) or f"exit {help_result.returncode}")
                )
        except (OSError, subprocess.TimeoutExpired) as exc:
            errors.append(f"backend argument-support probe failed: {exc}")

    if not source_root.exists():
        errors.append(f"backend source directory is missing: {source_root}")
    if not model_dir.exists():
        errors.append(f"RAPiDock runtime model directory is missing: {model_dir}")
    if not local.exists():
        errors.append(f"required RAPiDock checkpoint is missing from runtime model directory: {local}")
    long_required = len(sequence) >= threshold
    if long_required and not longer.exists():
        errors.append(
            f"{len(sequence)}-residue peptide requires {longer.name} at threshold {threshold}; "
            "the checkpoint is absent, so no fallback will be used"
        )

    mapping = audit_residue_mapping(receptor, reference_path=REFERENCE, chain="A")
    if not mapping["mapping_valid"]:
        errors.append("ARG203/ARG205 mapping is unresolved: " + "; ".join(mapping["failure_reasons"]))
    try:
        coordinate_frame = coordinate_frame_audit(
            [receptor], reference_path=REFERENCE, chain="A", site_center=SITE,
        )["representations"][0]
    except (OSError, ValueError, TargetMappingError) as exc:
        coordinate_frame = {"coordinate_frame_comparable": False, "error": str(exc)}
    if not coordinate_frame.get("coordinate_frame_comparable", False):
        errors.append("receptor target coordinates are not comparable to the verified reference frame")

    site_valid = len(SITE) == 3 and all(isinstance(value, (int, float)) for value in SITE)
    if not site_valid or not all(float(value) == float(value) for value in SITE):
        errors.append(f"target site coordinates are invalid: {SITE}")
    if box <= 0:
        errors.append(f"box size must be positive: {box}")
    if n_samples <= 0:
        errors.append(f"n-samples must be positive: {n_samples}")

    output_root = output_root.resolve()
    if output_root.exists() and any(output_root.iterdir()):
        errors.append(f"output directory is non-empty and would risk overwriting results: {output_root}")
    elif not output_root.parent.exists() or not output_root.parent.is_dir():
        errors.append(f"output parent directory does not exist: {output_root.parent}")

    selected = longer if long_required else local
    return {
        "ready": not errors,
        "errors": errors,
        "warnings": warnings,
        "backend_command_on_path": executable,
        "backend_source": str(source_root),
        "backend_source_present": source_root.exists(),
        "backend_source_version": _source_version(source_root),
        "runtime_model_directory": str(model_dir),
        "runtime_model_directory_present": model_dir.exists(),
        "backend_version_probe": version_probe,
        "backend_dock_help_confirmed": bool(help_text),
        "required_dock_arguments": ["--peptide", "--receptor", "--site", "--box", "--output-dir", "--n-samples", "--seed"],
        "peptide_length": len(sequence),
        "long_checkpoint_threshold": threshold,
        "long_checkpoint_expected": long_required,
        "long_checkpoint_present": longer.exists(),
        "selected_checkpoint_if_run": str(selected) if selected.exists() else None,
        "fallback_used": False,
        "checkpoint_provenance": {
            "rapidock_base_source": "RAPiDock authors, Zenodo 10.5281/zenodo.14193621",
            "rapidock_base_license": "CC-BY-4.0",
            "longer_local_source_status": "no published artifact or legitimate download source identified",
            "longer_local_policy": "required for 13+ residue strict preflight; fallback is disabled",
            "runtime_resolution": "third_party/RAPiDock/train_models/CGTensorProductEquivariantModel",
        },
        "rapidock_local_checkpoint_present": local.exists(),
        "checkpoint_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (local, longer) if path.exists()
        },
        "receptor": str(receptor),
        "receptor_mapping_valid": mapping["mapping_valid"],
        "coordinate_frame": coordinate_frame,
        "site_center_A": list(SITE),
        "declared_box_size_A": box,
        "n_samples": n_samples,
        "seed": seed,
        "site_constrains_stage_1_sampling": False,
        "stage_1_sampling_scope": "receptor-wide RAPiDock sampling",
        "off_pocket_filter": {"metric": "C-alpha centroid distance to site", "cutoff_A": 35.0},
        "auto_box_behavior": {"user_box_is_minimum": True, "maximum_edge_A": 100.0, "safety_margin_A": 4.0},
        "effective_box_size_recorded_by_archived_metadata": False,
        "output_root": str(output_root),
    }


def _preflight_report_path(output_root: Path) -> Path:
    """Keep a preflight report out of the directory reserved for run outputs."""
    return output_root.parent / f"{output_root.name}_preflight.json"


def _build_dock_command(
    name: str,
    sequence: str,
    receptor: Path,
    output_dir: Path,
    box: float,
    n_samples: int,
    seed: int | None,
) -> list[str]:
    command = [
        "hybridock-pep", "dock",
        "--peptide", sequence,
        "--receptor", str(receptor),
        "--site", str(SITE[0]), str(SITE[1]), str(SITE[2]),
        "--box", str(box),
        "--n-samples", str(n_samples),
        "--output-dir", str(output_dir),
    ]
    if seed is not None:
        command.extend(["--seed", str(seed)])
    return command


def run(args: argparse.Namespace) -> int:
    receptor = args.receptor
    if receptor == DEFAULT_RECEPTOR and not receptor.exists():
        prepare_receptor(output_path=receptor)
    preflight = _backend_preflight(
        CANDIDATES["P2-01"], args.box, args.n_samples,
        receptor=receptor, output_root=args.output_root, seed=args.seed,
    )
    report_path = _preflight_report_path(args.output_root)
    _write_json(report_path, preflight)
    if not preflight["ready"]:
        raise RuntimeError("Docking preflight failed:\n- " + "\n- ".join(preflight["errors"]) + f"\nReport: {report_path}")

    args.output_root.mkdir(parents=True, exist_ok=True)
    summary = {
        "started_at": datetime.now(timezone.utc).isoformat(),
        "reference": str(REFERENCE),
        "receptor": str(receptor),
        "receptor_mapping_valid": True,
        "output_root": str(args.output_root),
        "backend_preflight": preflight,
        "candidates": [],
    }
    for name, sequence in CANDIDATES.items():
        outdir = args.output_root / name
        outdir.mkdir(parents=True, exist_ok=True)
        candidate_seed = None if args.seed is None else args.seed + len(summary["candidates"])
        cmd = _build_dock_command(
            name, sequence, receptor, outdir, args.box, args.n_samples, candidate_seed,
        )
        print("\n" + "=" * 70)
        print(name, sequence)
        print(" ".join(cmd))
        result = subprocess.run(cmd, check=False)
        row = {
            "candidate": name,
            "sequence": sequence,
            "returncode": result.returncode,
            "output_dir": str(outdir),
            "seed": candidate_seed,
            "best_pose_exists": (outdir / "best_pose.pdb").exists(),
            "structural_validation": None,
        }
        if result.returncode == 0 and (outdir / "best_pose.pdb").exists():
            pose_receptor = outdir / "receptor_for_rapidock_full.pdb"
            if not pose_receptor.exists():
                pose_receptor = outdir / "receptor_for_rapidock.pdb"
            if not pose_receptor.exists():
                pose_receptor = receptor
            try:
                validation = diagnose_pose(
                    pose_receptor,
                    outdir / "best_pose.pdb",
                    sequence,
                    reference_path=REFERENCE,
                    site_center=SITE,
                    box_size_A=args.box,
                )
                _write_json(outdir / "structural_validation.json", validation)
                row["structural_validation"] = {
                    "status": validation.get("pose_geometry_status"),
                    "target_contact_atom_count": validation.get("target_contact_atom_count", 0),
                    "severe_clash_count": validation.get("severe_clash_count"),
                    "failure_categories": validation.get("failure_categories", []),
                }
            except Exception as exc:  # retain run output, but fail closed in metadata
                failure = {"status": "unavailable", "error": str(exc)}
                _write_json(outdir / "structural_validation.json", failure)
                row["structural_validation"] = failure
        summary["candidates"].append(row)
        print("COMPLETED:" if result.returncode == 0 else "FAILED:", name)
    summary["finished_at"] = datetime.now(timezone.utc).isoformat()
    _write_json(args.output_root / "batch_summary.json", summary)
    print(f"\nBatch summary: {args.output_root / 'batch_summary.json'}")
    return 0 if all(item["returncode"] == 0 for item in summary["candidates"]) else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor", type=Path, default=DEFAULT_RECEPTOR)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTROOT)
    parser.add_argument("--n-samples", type=int, default=100)
    parser.add_argument("--box", type=float, default=52.0)
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Base RNG seed; each candidate receives base seed plus its candidate index",
    )
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="write backend/site/checkpoint provenance without launching docking",
    )
    args = parser.parse_args()
    if args.preflight_only:
        receptor = args.receptor
        if receptor == DEFAULT_RECEPTOR and not receptor.exists():
            prepare_receptor(output_path=receptor)
        preflight = _backend_preflight(
            CANDIDATES["P2-01"], args.box, args.n_samples,
            receptor=receptor, output_root=args.output_root, seed=args.seed,
        )
        report_path = _preflight_report_path(args.output_root)
        _write_json(report_path, preflight)
        print(json.dumps(preflight, indent=2, sort_keys=True))
        print(f"Preflight report: {report_path}")
        return 0 if preflight["ready"] else 2
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
