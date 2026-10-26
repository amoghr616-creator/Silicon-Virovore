#!/usr/bin/env python3
"""Small geometry-control pilot for target_aware_structural_score_v1.1_clash_gated.

This pilot scores a supplied receptor/peptide complex and creates translated
pose controls. These are computational geometric decoys, not experimentally
verified nonbinders and not an independent-complex held-out benchmark.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.target_aware_scoring import (
    SCORE_VERSION, StructuralScoreConfig, aligned_peptide_ca_rmsd,
    parse_hotspot_spec, score_peptide_receptor_complex,
    sha256_file, translate_peptide_chain_pdb,
)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                keys.append(key); seen.add(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list)) else v for k, v in row.items()})


def _try_native_c_score(sequence: str, requested: bool) -> tuple[float | None, str]:
    if not requested:
        return None, "not_requested"
    if len(sequence) != 21:
        return None, f"not_computed_native_bridge_requires_exactly_21_aa_received_{len(sequence)}"
    if any(aa not in "ACDEFGHIKLMNPQRSTVWY" for aa in sequence):
        return None, "not_computed_noncanonical_sequence"
    try:
        from c.bridge import evaluate_sequence_metrics
        result = evaluate_sequence_metrics(sequence)
        return float(result["fitness_score"]), "success_native_sequence_prior_only"
    except Exception as exc:
        return None, f"native_bridge_unavailable:{type(exc).__name__}:{exc}"


def _pairwise_wins(reference: float, decoys: list[float]) -> dict[str, Any]:
    vals = [1.0 if reference > d else 0.5 if reference == d else 0.0 for d in decoys]
    return {"decoys_compared": len(decoys),
            "pairwise_ranking_accuracy_with_ties_half": (sum(vals) / len(vals) if vals else None),
            "wins": sum(v == 1.0 for v in vals), "ties": sum(v == 0.5 for v in vals),
            "losses": sum(v == 0.0 for v in vals)}



def _chain_heavy_atom_coordinates(pdb_path: Path, chain_id: str):
    """Read finite heavy-atom coordinates for one chain, first MODEL only."""
    import math
    lines = pdb_path.read_text(encoding="utf-8", errors="strict").splitlines()
    has_models = any(line.startswith("MODEL ") for line in lines)
    in_first_model = not has_models
    coords = []
    for line in lines:
        if has_models and line.startswith("MODEL "):
            if in_first_model:
                break
            in_first_model = True
            continue
        if has_models and line.startswith("ENDMDL") and in_first_model:
            break
        if has_models and not in_first_model:
            continue
        if line[:6].strip().upper() not in {"ATOM", "HETATM"} or len(line) < 54:
            continue
        if line[21:22].strip() != chain_id:
            continue
        if line[16:17] not in {" ", "A", ""}:
            continue
        atom_name = line[12:16].strip().upper()
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if element in {"H", "D"} or (not element and atom_name.startswith(("H", "D"))):
            continue
        try:
            xyz = tuple(float(line[a:b].strip()) for a, b in ((30, 38), (38, 46), (46, 54)))
        except ValueError as exc:
            raise ValueError(f"Malformed coordinates in chain {chain_id}: {line!r}") from exc
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError(f"Non-finite coordinates in chain {chain_id}: {line!r}")
        coords.append(xyz)
    if not coords:
        raise ValueError(f"No heavy-atom coordinates found for chain {chain_id!r} in {pdb_path}")
    return coords


def _far_separating_translations(source: Path, receptor_chain: str,
                                 peptide_chain: str, margin_A: float = 30.0):
    """Return deterministic translations putting the whole peptide outside receptor AABB.

    These are intentionally trivial non-contact controls, not hard decoys and
    not experimentally verified nonbinders. They test if no-interface poses
    receive low structural scores; they do not test native-pose recognition.
    """
    import math
    if not math.isfinite(margin_A) or margin_A < 8.0:
        raise ValueError("margin_A must be finite and at least 8.0 Å")
    receptor = _chain_heavy_atom_coordinates(source, receptor_chain)
    peptide = _chain_heavy_atom_coordinates(source, peptide_chain)
    rmin = [min(x[i] for x in receptor) for i in range(3)]
    rmax = [max(x[i] for x in receptor) for i in range(3)]
    pmin = [min(x[i] for x in peptide) for i in range(3)]
    pmax = [max(x[i] for x in peptide) for i in range(3)]
    # Move the complete peptide bounding box to the +X, -X, and +Y exterior.
    return [
        (rmax[0] + margin_A - pmin[0], 0.0, 0.0),
        (rmin[0] - margin_A - pmax[0], 0.0, 0.0),
        (0.0, rmax[1] + margin_A - pmin[1], 0.0),
    ]

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complex-pdb", required=True, help="One PDB with receptor and posed peptide in one coordinate frame")
    parser.add_argument("--receptor-chain", required=True)
    parser.add_argument("--peptide-chain", required=True)
    parser.add_argument("--target-id", default="unspecified")
    parser.add_argument("--hotspot", action="append", default=[], help="Verified CHAIN:RESNUM[:ICODE]:RESNAME; repeat as needed")
    parser.add_argument("--rmsd-moving-pdb", help="Optional independently generated pose for receptor-aligned peptide CA RMSD")
    parser.add_argument("--out-root", default="prototype3/results/target_aware_structural_v1")
    parser.add_argument("--seed", type=int, default=20261026)
    parser.add_argument("--skip-c-score", action="store_true", help="Do not try native baseline c_score, even for a 21-aa peptide")
    args = parser.parse_args()

    source = Path(args.complex_pdb).resolve()
    if not source.is_file():
        parser.error(f"complex PDB does not exist: {source}")
    try:
        hotspots = tuple(parse_hotspot_spec(value) for value in args.hotspot)
    except ValueError as exc:
        parser.error(str(exc))

    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out_root = Path(args.out_root)
    if not out_root.is_absolute():
        out_root = ROOT / out_root
    out_root.mkdir(parents=True, exist_ok=True)
    run_dir = out_root / f"pilot_{timestamp}_{os.getpid()}"
    run_dir.mkdir(parents=False, exist_ok=False)  # Never overwrite an earlier run.
    controls_dir = run_dir / "pose_controls"
    controls_dir.mkdir()

    config = StructuralScoreConfig(hotspot_residues=hotspots)
    # Prespecified controls are deliberately outside the receptor bounding box.
    # These trivial no-contact controls test geometry handling, not hard-decoy ranking.
    control_margin_A = 30.0
    translations = _far_separating_translations(
        source, args.receptor_chain, args.peptide_chain, margin_A=control_margin_A)
    created_controls = []
    for idx, (dx, dy, dz) in enumerate(translations, start=1):
        output = controls_dir / f"translated_decoy_{idx}.pdb"
        control = translate_peptide_chain_pdb(source, output, args.peptide_chain, dx, dy, dz)
        control["control_type"] = "computational_far_separated_pose_control_not_experimental_nonbinder"
        control["minimum_axis_separation_A"] = control_margin_A
        created_controls.append(control)

    pose_specs = [{"pose_id": "supplied_reference_pose", "path": source, "label": "input_reference_pose",
                   "control_type": "user_supplied_structure_not_independently_validated"}]
    pose_specs += [{"pose_id": f"translated_decoy_{i+1}", "path": Path(item["output"]),
                    "label": "computational_far_separated_control", "control_type": item["control_type"]}
                   for i, item in enumerate(created_controls)]

    raw_records = []
    for pose in pose_specs:
        rec = score_peptide_receptor_complex(
            pose["path"], config=config, receptor_chain=args.receptor_chain,
            peptide_chain=args.peptide_chain, target_structure_id=args.target_id,
            structure_provenance=pose["control_type"])
        c_score, c_status = _try_native_c_score(rec.get("peptide_sequence") or "", not args.skip_c_score)
        rec["c_score"] = c_score
        rec["c_score_status"] = c_status
        rec["pose_id"] = pose["pose_id"]
        rec["pose_label"] = pose["label"]
        rec["complex_group_id"] = f"{args.target_id}:{sha256_file(source)[:12]}"
        rec["input_pdb_sha256"] = sha256_file(pose["path"])
        rec["relative_pose_path"] = str(Path(pose["path"]).relative_to(run_dir) if Path(pose["path"]).is_relative_to(run_dir) else Path(pose["path"]))
        raw_records.append(rec)

    score_rows = []
    for rec in raw_records:
        m = rec.get("structural_metrics", {})
        score_rows.append({
            "pose_id": rec["pose_id"], "pose_label": rec["pose_label"],
            "complex_group_id": rec["complex_group_id"], "target_structure_id": rec["target_structure_id"],
            "peptide_sequence": rec["peptide_sequence"], "c_score": rec["c_score"],
            "c_score_status": rec["c_score_status"],
            "target_aware_structural_score": rec["target_aware_structural_score"],
            "structural_score_status": rec["structural_score_status"], "score_valid": rec["score_valid"],
            "receptor_chain": rec["receptor_chain"], "peptide_chain": rec["peptide_chain"],
            "structure_provenance": rec["structure_provenance"], "structural_metrics": m,
            "eligibility_reasons": rec["eligibility_reasons"], "input_pdb_sha256": rec["input_pdb_sha256"],
            "score_version": rec["score_version"],
        })
    _write_csv(run_dir / "structural_scores.csv", score_rows)

    manifest_rows = []
    for pose in pose_specs:
        manifest_rows.append({
            "pose_id": pose["pose_id"], "complex_group_id": f"{args.target_id}:{sha256_file(source)[:12]}",
            "target_structure_id": args.target_id, "complex_pdb": str(pose["path"]),
            "complex_pdb_sha256": sha256_file(pose["path"]), "receptor_chain": args.receptor_chain,
            "peptide_chain": args.peptide_chain, "label": pose["label"], "control_type": pose["control_type"],
            "independent_biological_unit": f"{args.target_id}:{sha256_file(source)[:12]}",
            "selection_rule": "one supplied pose plus three deterministic translations placing the peptide bounding box at least 30 A outside the receptor bounding box; trivial non-contact controls, no score-based selection",
            "independent_experimental_nonbinder": False,
        })
    _write_csv(run_dir / "benchmark_manifest.csv", manifest_rows)
    _write_csv(run_dir / "split_manifest.csv", [{
        "complex_group_id": manifest_rows[0]["complex_group_id"], "split": "development_pilot_not_locked_test",
        "seed": args.seed, "independent_complex_groups_available": 1,
        "reason": "single input complex and its geometric controls are insufficient for an independent held-out generalization test",
    }])

    ref = next(r for r in raw_records if r["pose_id"] == "supplied_reference_pose")
    decoys = [r for r in raw_records if r["pose_label"] in {"computational_decoy", "computational_far_separated_control"}]
    valid_decoys = [r for r in decoys if r["score_valid"] and r["target_aware_structural_score"] is not None]
    valid_reference = ref["score_valid"] and ref["target_aware_structural_score"] is not None
    metrics = {
        "score_version": SCORE_VERSION,
        "experiment_type": "small_development_pilot_geometry_sensitivity_not_confirmatory_validation",
        "valid_input_reference_pose_count": int(bool(valid_reference)),
        "independently_verified_experimental_reference_complexes": None,
        "computational_translated_pose_decoys": len(decoys),
        "control_construction": "three deterministic far-separated translations; geometry sanity controls, not hard decoys",
        "control_margin_A": control_margin_A,
        "valid_decoy_scores": len(valid_decoys),
        "invalid_decoy_scores": len(decoys) - len(valid_decoys),
        "invalid_decoy_pose_ids": [r["pose_id"] for r in decoys
                                   if not r["score_valid"] or r["target_aware_structural_score"] is None],
        "independent_complex_groups": 1,
        "locked_held_out_test_performed": False,
        "input_reference_provenance_independently_verified": False,
        "structural_pairwise_reference_vs_translated_decoys": _pairwise_wins(
            float(ref["target_aware_structural_score"]), [float(r["target_aware_structural_score"]) for r in valid_decoys]
        ) if valid_reference else None,
        "contact_only_pairwise_reference_vs_translated_decoys": _pairwise_wins(
            float(ref["structural_metrics"]["simple_baseline_objectives_0_1"]["contact_only"]),
            [float(r["structural_metrics"]["simple_baseline_objectives_0_1"]["contact_only"]) for r in valid_decoys]
        ) if valid_reference else None,
        "geometry_only_pairwise_reference_vs_translated_decoys": _pairwise_wins(
            float(ref["structural_metrics"]["simple_baseline_objectives_0_1"]["geometry_only"]),
            [float(r["structural_metrics"]["simple_baseline_objectives_0_1"]["geometry_only"]) for r in valid_decoys]
        ) if valid_reference else None,
        "clash_only_pairwise_reference_vs_translated_decoys": _pairwise_wins(
            float(ref["structural_metrics"]["simple_baseline_objectives_0_1"]["clash_only_no_clash_is_best"]),
            [float(r["structural_metrics"]["simple_baseline_objectives_0_1"]["clash_only_no_clash_is_best"]) for r in valid_decoys]
        ) if valid_reference else None,
        "c_score_pairwise_reference_vs_translated_decoys": _pairwise_wins(
            float(ref["c_score"]), [float(r["c_score"]) for r in decoys if r["c_score"] is not None]
        ) if ref["c_score"] is not None and all(r["c_score"] is not None for r in decoys) else None,
        "c_score_comparison_status": ref["c_score_status"],
        "interpretation": [
            "Translated poses are computational geometric controls, not experimental nonbinders.",
            "The native c_score is sequence-only; for the same peptide its value should not change under pose translation.",
            "A one-complex pilot cannot establish generalization or affinity prediction.",
        ],
    }
    if args.rmsd_moving_pdb:
        try:
            metrics["reference_pose_alignment_audit"] = aligned_peptide_ca_rmsd(
                source, Path(args.rmsd_moving_pdb), receptor_chain=args.receptor_chain, peptide_chain=args.peptide_chain)
        except Exception as exc:
            metrics["reference_pose_alignment_audit"] = {"status": "failed", "reason": f"{type(exc).__name__}:{exc}"}
    else:
        metrics["reference_pose_alignment_audit"] = {"status": "not_run", "reason": "--rmsd-moving-pdb not supplied; do not infer a new RMSD from translated controls"}
    (run_dir / "benchmark_metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    per_complex = [{
        "complex_group_id": ref["complex_group_id"],
        "target_structure_id": args.target_id,
        "input_reference_pose_valid": valid_reference,
        "number_of_translated_pose_decoys": len(decoys),
        "number_of_valid_decoys": len(valid_decoys),
        "number_of_invalid_decoys": len(decoys) - len(valid_decoys),
        "structural_pairwise_accuracy": (metrics["structural_pairwise_reference_vs_translated_decoys"] or {}).get("pairwise_ranking_accuracy_with_ties_half"),
        "contact_only_pairwise_accuracy": (metrics["contact_only_pairwise_reference_vs_translated_decoys"] or {}).get("pairwise_ranking_accuracy_with_ties_half"),
        "geometry_only_pairwise_accuracy": (metrics["geometry_only_pairwise_reference_vs_translated_decoys"] or {}).get("pairwise_ranking_accuracy_with_ties_half"),
        "clash_only_pairwise_accuracy": (metrics["clash_only_pairwise_reference_vs_translated_decoys"] or {}).get("pairwise_ranking_accuracy_with_ties_half"),
        "c_score_pairwise_accuracy": (metrics["c_score_pairwise_reference_vs_translated_decoys"] or {}).get("pairwise_ranking_accuracy_with_ties_half"),
        "independent_unit_count": 1,
        "interpretation": "single-complex geometric sensitivity pilot, not held-out validation",
    }]
    _write_csv(run_dir / "per_complex_metrics.csv", per_complex)

    try:
        from importlib.metadata import version as pkg_version
        numpy_version = pkg_version("numpy")
    except Exception:
        numpy_version = "unknown"
    (run_dir / "run_config.json").write_text(json.dumps({
        "score_version": SCORE_VERSION, "formula": ref.get("structural_metrics", {}).get("score_formula"),
        "config": {**config.__dict__, "hotspot_residues": [h.__dict__ for h in hotspots]},
        "seed": args.seed, "pose_control_translations_A": [list(t) for t in translations],
        "pose_control_margin_A": control_margin_A,
        "pose_control_definition": "axis-aligned peptide bounding box moved at least the configured margin beyond receptor bounding box; controls deliberately lack an interface and are not hard decoys",
        "input_complex_pdb": str(source), "input_complex_sha256": sha256_file(source),
        "receptor_chain": args.receptor_chain, "peptide_chain": args.peptide_chain,
        "target_structure_id": args.target_id, "c_score_attempted": not args.skip_c_score,
        "pose_generation": "not performed: only supplied pose and deterministic translated controls scored",
        "timestamp_utc": timestamp,
    }, indent=2, sort_keys=True, default=list) + "\n", encoding="utf-8")
    try:
        commit_result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                       check=True, capture_output=True, text=True)
        git_commit = commit_result.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        git_commit = None
    (run_dir / "environment.json").write_text(json.dumps({
        "python": sys.version, "platform": platform.platform(), "numpy": numpy_version,
        "score_module": SCORE_VERSION, "git_commit": git_commit,
        "note": "Working-tree state is preserved and separately recorded by the patch preflight installer.",
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    report_lines = [
        "# P3 Target-Aware Structural Score Pilot", "",
        f"- Score version: `{SCORE_VERSION}`",
        f"- Target ID label supplied by caller: `{args.target_id}`",
        f"- Input complex SHA-256: `{sha256_file(source)}`",
        f"- Structurally valid input reference-pose score: `{valid_reference}`",
        "- Independently verified experimental reference complexes: not determined by this pilot",
        f"- Control construction: deterministic far-separated geometric controls (minimum axis separation `{control_margin_A:g} Å`); not hard decoys",
        f"- Valid translated decoy scores: `{len(valid_decoys)}/{len(decoys)}`",
        f"- Invalid translated decoys rejected: `{len(decoys) - len(valid_decoys)}/{len(decoys)}`",
        f"- Independent complex groups: `1` (pilot only; no locked test)",
        f"- Structural pairwise reference-over-decoy accuracy: `{metrics['structural_pairwise_reference_vs_translated_decoys']}`",
        f"- Native `c_score` comparison: `{metrics['c_score_comparison_status']}`",
        f"- Contact-only pairwise accuracy: `{metrics['contact_only_pairwise_reference_vs_translated_decoys']}`",
        f"- Geometry-only pairwise accuracy: `{metrics['geometry_only_pairwise_reference_vs_translated_decoys']}`",
        f"- Clash-only pairwise accuracy: `{metrics['clash_only_pairwise_reference_vs_translated_decoys']}`",
        "", "## Interpretation", "",
        "This run tests deterministic sensitivity to artificial peptide translations, not recognition of experimentally validated binders. The three controls are deliberately placed outside the receptor bounding box with at least 30 A axis separation; they are trivial no-contact controls, not hard decoys or experimentally verified nonbinders. No affinity prediction, biological efficacy, or cross-complex generalization claim is supported.",
        "The score uses coordinate geometry only. Solvent-accessible surface area and electrostatics are not calculated. Any inter-chain heavy-atom pair below the configured severe-clash cutoff (2.0 Å by default) invalidates the official structural score: the score is null, diagnostic metrics are retained, and the pose is excluded from structural ranking metrics.",
        "A valid reference pose is whatever complex PDB was supplied; inspect its provenance and chain mapping before calling it an experimental reference.",
        "", "## Reproduction", "",
        "```bash", f"python prototype3/run_target_aware_pilot.py --complex-pdb '{source}' --receptor-chain '{args.receptor_chain}' --peptide-chain '{args.peptide_chain}' --target-id '{args.target_id}' --seed {args.seed}", "```", "",
        "Generated files: `benchmark_manifest.csv`, `split_manifest.csv`, `structural_scores.csv`, `benchmark_metrics.json`, `run_config.json`, `environment.json`, and translated control PDBs.",
    ]
    (run_dir / "REPORT.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(json.dumps({"run_dir": str(run_dir), "benchmark_metrics": metrics, "files": sorted(p.name for p in run_dir.iterdir())}, indent=2, sort_keys=True))
    return 0 if valid_reference else 2


if __name__ == "__main__":
    raise SystemExit(main())
