"""Production target-aware scoring and reliability reporting.

This module does not generate docking poses or predict experimental affinity.
It evaluates sequence-matched, pre-generated receptor-peptide complexes.

Reliability remains diagnostic-only until externally calibrated.
"""

from __future__ import annotations

import hashlib
import json
import json
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Any

from src.target_aware_scoring import (
    HotspotSpec,
    StructuralScoreConfig,
    score_peptide_receptor_complex,
)


TARGET_VERSION = "target_config_v3_generic"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _normalise_hotspots(raw: Any, receptor_chain: str) -> tuple[HotspotSpec, ...]:
    if not isinstance(raw, list) or not raw:
        raise ValueError("At least one target hotspot must be configured")

    result = []
    seen = set()

    for item in raw:
        if not isinstance(item, dict):
            raise ValueError("Every hotspot must be a JSON object")

        chain = str(item.get("chain", receptor_chain)).strip()
        residue_number = int(item["residue_number"])
        residue_name = str(item["residue_name"]).strip().upper()
        insertion_code = str(item.get("insertion_code", "")).strip()

        if not chain or not residue_name or residue_number <= 0:
            raise ValueError(f"Invalid hotspot definition: {item}")

        if chain != receptor_chain:
            raise ValueError(
                f"Hotspot chain {chain!r} does not match receptor chain "
                f"{receptor_chain!r}"
            )

        spec = HotspotSpec(
            chain=chain,
            residue_number=residue_number,
            residue_name=residue_name,
            insertion_code=insertion_code,
        )

        if spec.label() in seen:
            raise ValueError(f"Duplicate hotspot: {spec.label()}")

        seen.add(spec.label())
        result.append(spec)

    return tuple(result)


def _default_legacy_hotspots(target_id: str, receptor_chain: str) -> list[dict]:
    """Preserve the historical HERV-K configuration for old callers."""
    if target_id != "9SYA":
        return []

    return [
        {
            "chain": receptor_chain,
            "residue_number": 203,
            "residue_name": "ARG",
        },
        {
            "chain": receptor_chain,
            "residue_number": 205,
            "residue_name": "ARG",
        },
    ]


def _register_manifest(
    candidates: list[Any],
    manifest_path: str | Path,
) -> dict[str, Any]:
    path = Path(manifest_path).expanduser().resolve()

    if not path.is_file():
        raise FileNotFoundError(f"Target manifest not found: {path}")

    manifest = json.loads(path.read_text(encoding="utf-8"))

    if manifest.get("schema_version") != 1:
        raise ValueError("Unsupported target manifest schema_version")

    target = manifest.get("target")
    complexes = manifest.get("complexes")

    if not isinstance(target, dict) or not isinstance(complexes, list):
        raise ValueError("Manifest must contain target and complexes")

    target_id = str(target.get("target_structure_id", "")).strip()
    receptor_chain = str(target.get("receptor_chain", "")).strip()
    peptide_chain = str(target.get("peptide_chain", "")).strip()

    if not target_id or not receptor_chain or not peptide_chain:
        raise ValueError("Target ID and both chain IDs are required")

    if receptor_chain == peptide_chain:
        raise ValueError("Receptor and peptide chains must be different")

    hotspots = _normalise_hotspots(
        target.get("hotspots"),
        receptor_chain,
    )

    by_sequence: dict[str, list[dict[str, Any]]] = defaultdict(list)
    rejected = []

    for index, raw in enumerate(complexes):
        try:
            if not isinstance(raw, dict):
                raise ValueError("Complex record must be an object")

            sequence = str(raw.get("sequence", "")).strip().upper()
            raw_path = raw.get("path")
            expected_hash = str(raw.get("sha256", "")).lower()

            if not sequence or not raw_path:
                raise ValueError("Sequence and path are required")

            if len(sequence) < 2:
                raise ValueError("Invalid peptide sequence")

            if not expected_hash:
                raise ValueError("A SHA-256 hash is required")

            complex_path = Path(raw_path).expanduser()

            if not complex_path.is_absolute():
                complex_path = path.parent / complex_path

            complex_path = complex_path.resolve()

            if not complex_path.is_file():
                raise FileNotFoundError(complex_path)

            actual_hash = _sha256(complex_path)

            if actual_hash != expected_hash:
                raise ValueError(
                    f"SHA-256 mismatch for {complex_path.name}"
                )

            by_sequence[sequence].append({
                "path": str(complex_path),
                "sha256": actual_hash,
                "sequence": sequence,
                "target_structure_id": target_id,
                "receptor_chain": receptor_chain,
                "peptide_chain": peptide_chain,
                "hotspots": [
                    {
                        "chain": h.chain,
                        "residue_number": h.residue_number,
                        "residue_name": h.residue_name,
                        "insertion_code": h.insertion_code,
                    }
                    for h in hotspots
                ],
                "provenance": str(
                    raw.get("provenance", "manifest_registered_complex")
                ),
                "independent_run_id": raw.get("independent_run_id"),
                "random_seed": raw.get("random_seed"),
                "protocol_sha256": raw.get("protocol_sha256"),
                "execution_id": raw.get("execution_id"),
                "execution_manifest_path": raw.get("execution_manifest_path"),
                "execution_manifest_sha256": raw.get("execution_manifest_sha256"),
            })

        except Exception as exc:
            rejected.append({
                "record_index": index,
                "error": f"{type(exc).__name__}: {exc}",
            })

    matched_records = 0
    matched_candidates = 0
    live_metadata_candidates = 0
    live_metadata_complexes = 0
    manifest_fallback_candidates = 0
    manifest_fallback_complexes = 0
    manifest_sequence_match_candidates = 0

    for candidate in candidates:
        meta = candidate.metadata
        sequence = str(candidate.sequence).upper()
        manifest_entries = by_sequence.get(sequence, [])
        if manifest_entries:
            manifest_sequence_match_candidates += 1

        # Live structures generated for this candidate in the current run
        # take precedence over a configured historical manifest. The old
        # behavior unconditionally assigned manifest_entries here, silently
        # replacing freshly generated HADDOCK poses with an empty list when
        # the historical manifest did not contain the current sequence.
        live_entries = meta.get("target_complexes")
        if isinstance(live_entries, dict):
            live_entries = [live_entries]
        elif live_entries:
            live_entries = list(live_entries)
        else:
            live_entries = []

        if live_entries:
            entries = live_entries
            live_metadata_candidates += 1
            live_metadata_complexes += len(entries)
            source = "live_candidate_metadata"
        else:
            entries = manifest_entries
            if entries:
                manifest_fallback_candidates += 1
                manifest_fallback_complexes += len(entries)
            source = "historical_manifest_fallback"

        # Fill missing target metadata from the manifest, but never replace
        # target identity or hotspot configuration already attached to the
        # live candidate.
        defaults = {
            "target_structure_id": target_id,
            "target_receptor_chain": receptor_chain,
            "target_peptide_chain": peptide_chain,
            "target_hotspots": target["hotspots"],
        }
        for key, value in defaults.items():
            if meta.get(key) in (None, "", []):
                meta[key] = value

        meta["target_complexes"] = entries
        meta["target_complex_source"] = source

        if entries:
            matched_candidates += 1
            matched_records += len(entries)

    return {
        "status": "loaded",
        "manifest_path": str(path),
        "target_structure_id": target_id,
        "candidate_count": len(candidates),
        # Effective records available to the evaluator, regardless of source.
        "matched_candidate_count": matched_candidates,
        "matched_complex_count": matched_records,
        "unmatched_candidate_count": len(candidates) - matched_candidates,
        "source_precedence": "live_candidate_metadata_then_manifest_fallback",
        "live_metadata_candidate_count": live_metadata_candidates,
        "live_metadata_complex_count": live_metadata_complexes,
        "manifest_fallback_candidate_count": manifest_fallback_candidates,
        "manifest_fallback_complex_count": manifest_fallback_complexes,
        "manifest_sequence_match_candidate_count": (
            manifest_sequence_match_candidates
        ),
        "rejected_manifest_records": rejected,
    }



def build_selection_comparison(
    sequence_only_candidates: list[Any],
    target_ordered_candidates: list[Any],
    target_report: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Compare generation winners before and after target-aware ordering.

    Both rankings must contain the same candidate objects from the same
    generation. This isolates the selection effect from candidate generation
    and avoids comparing scores across independently sampled populations.
    """
    report = target_report or {}
    # Compare the full identity multiset, not just a set: duplicate entries
    # must not let [A, A, B] masquerade as the same cohort as [A, B, B].
    same_cohort = (
        len(sequence_only_candidates) == len(target_ordered_candidates)
        and sorted(id(item) for item in sequence_only_candidates)
        == sorted(id(item) for item in target_ordered_candidates)
    )

    def number(value: Any) -> float | None:
        try:
            result = float(value)
        except (TypeError, ValueError):
            return None
        return result if math.isfinite(result) else None

    def winner_record(candidate: Any | None) -> dict[str, Any] | None:
        if candidate is None:
            return None
        meta = getattr(candidate, "metadata", {}) or {}
        evidence = meta.get("p4_target_evidence", {}) or {}
        return {
            "candidate_id": meta.get("candidate_id"),
            "sequence": getattr(candidate, "sequence", None),
            "overall_score": number(getattr(candidate, "overall_score", None)),
            "c_score": number(getattr(candidate, "c_score", None)),
            "target_geometry_score_mean": number(
                evidence.get("target_geometry_score_mean")
            ),
            "hotspot_proximity_score": number(
                evidence.get("hotspot_proximity_score")
            ),
            "target_evidence_status": evidence.get("status", "unavailable"),
        }

    sequence_winner = (
        sequence_only_candidates[0] if sequence_only_candidates else None
    )
    target_winner = (
        target_ordered_candidates[0] if target_ordered_candidates else None
    )
    target_ordering_applied = bool(
        report.get("pareto_reordered_candidates", False)
    )
    changed = (
        sequence_winner is not None
        and target_winner is not None
        and sequence_winner is not target_winner
    )

    seq_record = winner_record(sequence_winner)
    target_record = winner_record(target_winner)
    deltas = {}
    if same_cohort and seq_record is not None and target_record is not None:
        for key in (
            "overall_score",
            "c_score",
            "target_geometry_score_mean",
            "hotspot_proximity_score",
        ):
            left = seq_record.get(key)
            right = target_record.get(key)
            deltas[f"target_winner_minus_sequence_only_winner_{key}"] = (
                right - left
                if left is not None and right is not None
                else None
            )

    if not same_cohort:
        effect = "invalid_comparison_cohort_mismatch"
    elif not target_ordering_applied:
        effect = "target_ordering_not_applied_sequence_only_order_retained"
    elif changed:
        effect = "target_ordering_changed_generation_winner"
    else:
        effect = "target_ordering_applied_but_generation_winner_unchanged"

    baseline_rank_of_target_winner = None
    if same_cohort and target_winner is not None:
        baseline_rank_of_target_winner = next(
            index + 1
            for index, item in enumerate(sequence_only_candidates)
            if item is target_winner
        )

    return {
        "comparison_scope": "same_generation_same_candidate_cohort",
        "candidate_count": len(sequence_only_candidates),
        "same_candidate_cohort_verified": same_cohort,
        "sequence_only_order_source": "CandidateRanker_before_P4_target_ordering",
        "target_ordering_applied": target_ordering_applied,
        "target_selection_method": report.get("selection_method"),
        "selection_changed": bool(changed),
        "selection_effect": effect,
        "sequence_only_winner": seq_record,
        "target_ordered_winner": target_record,
        "target_ordered_winner_sequence_only_rank": (
            baseline_rank_of_target_winner
        ),
        "winner_metric_deltas": deltas,
        "interpretation_limit": (
            "This is a within-generation selection comparison using the same "
            "candidate cohort. It does not validate binding affinity, target "
            "accuracy, or improvement on independent experimental data."
        ),
    }

def _interface_fingerprint(
    path: Path,
    receptor_chain: str,
    peptide_chain: str,
    cutoff_A: float = 5.0,
) -> set[str]:
    """Return unique receptor-peptide residue pairs within the cutoff."""
    import numpy as np

    from src.target_aware_scoring import _parse_pdb_atoms

    atoms = _parse_pdb_atoms(path, {receptor_chain, peptide_chain})

    receptor = [
        atom for atom in atoms
        if atom.chain == receptor_chain
        and atom.element.upper() not in {"H", "D"}
    ]
    peptide = [
        atom for atom in atoms
        if atom.chain == peptide_chain
        and atom.element.upper() not in {"H", "D"}
    ]

    if not receptor or not peptide:
        raise ValueError("Cannot build interface fingerprint from empty chains")

    rxyz = np.asarray([atom.xyz for atom in receptor], dtype=float)
    pxyz = np.asarray([atom.xyz for atom in peptide], dtype=float)

    distances_squared = np.sum(
        (rxyz[:, None, :] - pxyz[None, :, :]) ** 2,
        axis=2,
    )

    rows, cols = np.where(distances_squared <= cutoff_A ** 2)

    def residue_token(atom: Any) -> str:
        return (
            f"{atom.chain}:{atom.residue_number}"
            f"{atom.insertion_code or ''}:{atom.residue_name}"
        )

    return {
        f"{residue_token(receptor[i])}|{residue_token(peptide[j])}"
        for i, j in zip(rows, cols)
    }


def _reliability_diagnostic(
    complex_results: list[dict[str, Any]],
) -> dict[str, Any]:
    """Calculate an uncalibrated repeatability diagnostic.

    Every scored pose must carry a nonempty run ID, an integer random seed,
    and a valid shared SHA-256 protocol hash. At least two documented runs
    and two poses per run are required. Missing metadata is a hard abstention:
    it is never silently excluded from the diagnostic while remaining in the
    candidate's structural score.
    """
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)

    if not complex_results:
        return {
            "status": "no_complex_results",
            "independent_run_count": 0,
            "calibrated": False,
            "feedback_active": False,
            "reason": "No scored complexes are available for repeatability analysis.",
        }

    # Fail closed if even one scored pose is missing its run identity.
    for index, result in enumerate(complex_results):
        raw_run_id = result.get("independent_run_id")
        if raw_run_id is None or not str(raw_run_id).strip():
            return {
                "status": "missing_independent_run_metadata",
                "independent_run_count": len(groups),
                "calibrated": False,
                "feedback_active": False,
                "reason": f"Pose {index} has no independent_run_id.",
            }
        groups[str(raw_run_id).strip()].append(result)

    if len(groups) < 2:
        return {
            "status": "insufficient_independent_runs",
            "independent_run_count": len(groups),
            "calibrated": False,
            "feedback_active": False,
            "reason": "At least two documented independent runs are required.",
        }

    for run_id, entries in groups.items():
        if len(entries) < 2:
            return {
                "status": "insufficient_poses_per_run",
                "independent_run_count": len(groups),
                "calibrated": False,
                "feedback_active": False,
                "reason": f"Run {run_id} has fewer than two poses.",
            }

    # Protocol identity must be an actual SHA-256 digest on every pose.
    protocol_values = []
    for index, result in enumerate(complex_results):
        value = result.get("protocol_sha256")
        digest = str(value).strip().lower() if value is not None else ""
        if len(digest) != 64 or any(ch not in "0123456789abcdef" for ch in digest):
            return {
                "status": "protocol_mismatch_or_missing",
                "independent_run_count": len(groups),
                "calibrated": False,
                "feedback_active": False,
                "reason": f"Pose {index} lacks a valid 64-character SHA-256 protocol hash.",
            }
        protocol_values.append(digest)

    protocols = set(protocol_values)
    if len(protocols) != 1:
        return {
            "status": "protocol_mismatch_or_missing",
            "independent_run_count": len(groups),
            "calibrated": False,
            "feedback_active": False,
            "reason": "Independent runs must share one exact protocol SHA-256.",
        }

    seeds: dict[str, int] = {}

    for run_id, entries in groups.items():
        run_seeds = set()
        for entry in entries:
            raw_seed = entry.get("random_seed")
            if isinstance(raw_seed, bool) or raw_seed is None:
                return {
                    "status": "invalid_run_seed_metadata",
                    "independent_run_count": len(groups),
                    "calibrated": False,
                    "feedback_active": False,
                    "reason": f"Run {run_id} contains a pose without a valid random seed.",
                }
            try:
                seed_value = int(raw_seed)
                if float(raw_seed) != seed_value or seed_value < 0:
                    raise ValueError("seed is not a nonnegative integer")
            except (TypeError, ValueError, OverflowError):
                return {
                    "status": "invalid_run_seed_metadata",
                    "independent_run_count": len(groups),
                    "calibrated": False,
                    "feedback_active": False,
                    "reason": f"Run {run_id} contains a non-integer random seed.",
                }
            run_seeds.add(seed_value)

        if len(run_seeds) != 1:
            return {
                "status": "invalid_run_seed_metadata",
                "independent_run_count": len(groups),
                "calibrated": False,
                "feedback_active": False,
                "reason": f"Run {run_id} lacks one consistent random seed.",
            }

        seeds[run_id] = next(iter(run_seeds))

    if len(set(seeds.values())) != len(seeds):
        return {
            "status": "nonindependent_seed_metadata",
            "independent_run_count": len(groups),
            "calibrated": False,
            "feedback_active": False,
            "reason": "Independent run IDs must use distinct random seeds.",
        }

    execution_ids = set()
    for run_id, entries in groups.items():
        run_execution_ids = {str(entry.get("execution_id") or "").strip() for entry in entries}
        run_manifest_hashes = {str(entry.get("execution_manifest_sha256") or "").strip().lower() for entry in entries}
        if len(run_execution_ids) != 1 or not next(iter(run_execution_ids), ""):
            return {
                "status": "missing_or_conflicting_execution_identity",
                "independent_run_count": len(groups),
                "calibrated": False,
                "feedback_active": False,
                "reason": f"Run {run_id} lacks one shared execution_id across its poses.",
            }
        if len(run_manifest_hashes) != 1:
            return {
                "status": "missing_or_conflicting_execution_manifest",
                "independent_run_count": len(groups),
                "calibrated": False,
                "feedback_active": False,
                "reason": f"Run {run_id} lacks one shared execution-manifest hash.",
            }
        digest = next(iter(run_manifest_hashes))
        if len(digest) != 64 or any(ch not in "0123456789abcdef" for ch in digest):
            return {
                "status": "missing_or_invalid_execution_manifest",
                "independent_run_count": len(groups),
                "calibrated": False,
                "feedback_active": False,
                "reason": f"Run {run_id} has no valid execution-manifest SHA-256.",
            }
        execution_ids.add(next(iter(run_execution_ids)))

    if len(execution_ids) != len(groups):
        return {
            "status": "reused_execution_identity",
            "independent_run_count": len(groups),
            "calibrated": False,
            "feedback_active": False,
            "reason": "Each independent run must have a distinct recorded execution_id.",
        }

    runs = []
    pose_fingerprint_summary = {}

    for run_id in sorted(groups):
        fingerprints = []
        fingerprint_errors = []

        for entry in groups[run_id]:
            fingerprint = entry.get("_interface_fingerprint_5A")
            fingerprint_error = entry.get("fingerprint_error")

            # A missing fingerprint or a parser/calculation error is not a
            # measured zero-contact pose. Keep those cases fail-closed.
            if fingerprint_error:
                fingerprint_errors.append(str(fingerprint_error))
                continue
            if fingerprint is None:
                return {
                    "status": "missing_interface_fingerprint",
                    "independent_run_count": len(groups),
                    "calibrated": False,
                    "feedback_active": False,
                    "reason": f"Run {run_id} has a pose with no fingerprint field.",
                }

            # An empty collection is a valid observed result: no inter-chain
            # residue pairs were within the 5 A cutoff. Preserve it as an
            # empty set in the entropy/consensus calculation; do not treat it
            # as missing data or silently discard it.
            fingerprints.append(set(fingerprint))

        if fingerprint_errors:
            return {
                "status": "fingerprint_calculation_failed",
                "independent_run_count": len(groups),
                "calibrated": False,
                "feedback_active": False,
                "reason": (
                    f"Run {run_id} has {len(fingerprint_errors)} fingerprint "
                    "calculation error(s): " + "; ".join(fingerprint_errors[:3])
                ),
            }

        nonempty_count = sum(bool(fingerprint) for fingerprint in fingerprints)
        empty_count = len(fingerprints) - nonempty_count
        if len(fingerprints) < 2 or nonempty_count < 2:
            return {
                "status": "insufficient_nonempty_interface_fingerprints",
                "independent_run_count": len(groups),
                "calibrated": False,
                "feedback_active": False,
                "reason": (
                    f"Run {run_id} has {nonempty_count} nonempty interface "
                    "fingerprint(s); at least two are required per run."
                ),
                "pose_fingerprint_summary": {
                    run_id: {
                        "pose_count": len(fingerprints),
                        "nonempty_fingerprint_count": nonempty_count,
                        "empty_fingerprint_count": empty_count,
                    }
                },
            }

        pose_fingerprint_summary[run_id] = {
            "pose_count": len(fingerprints),
            "nonempty_fingerprint_count": nonempty_count,
            "empty_fingerprint_count": empty_count,
            "nonempty_fraction": nonempty_count / len(fingerprints),
        }
        runs.append(fingerprints)

    try:
        from prototype4.p4_closed_loop import cross_run_reliability

        result = cross_run_reliability(runs)

    except Exception as exc:
        return {
            "status": "calculation_failed",
            "independent_run_count": len(groups),
            "calibrated": False,
            "feedback_active": False,
            "reason": f"{type(exc).__name__}: {exc}",
        }

    result.update({
        "status": "diagnostic_only",
        "n_independent_runs": len(groups),
        "random_seeds": seeds,
        "protocol_sha256": next(iter(protocols)),
        "pose_fingerprint_summary": pose_fingerprint_summary,
        "empty_fingerprint_policy": (
            "Measured empty 5 A contact sets are retained as zero-contact "
            "observations; every independent run must contain at least two "
            "nonempty pose fingerprints. Missing fields and calculation "
            "errors remain ineligible."
        ),
        "calibrated": False,
        "feedback_active": False,
        "interpretation": (
            "Interface repeatability diagnostic only. "
            "Not experimental binding probability or structural accuracy."
        ),
    })

    return result


def _score_complex(
    candidate: Any,
    entry: dict[str, Any],
    *,
    require_execution_manifest: bool = False,
) -> dict[str, Any]:
    path = Path(entry["path"]).expanduser().resolve()

    if not path.is_file():
        raise FileNotFoundError(path)

    expected_hash = entry.get("sha256")

    if expected_hash and _sha256(path) != str(expected_hash).lower():
        raise ValueError("Complex SHA-256 changed after manifest registration")

    require_execution_manifest = bool(
        require_execution_manifest or entry.get("_require_execution_manifest")
    )
    execution_id = entry.get("execution_id")
    execution_manifest_path = entry.get("execution_manifest_path")
    execution_manifest_sha256 = entry.get("execution_manifest_sha256")
    if require_execution_manifest:
        if not execution_id or not execution_manifest_path or not execution_manifest_sha256:
            raise ValueError("Independent execution manifest metadata is required for P4-TA.")
        manifest_path = Path(execution_manifest_path).expanduser().resolve()
        if not manifest_path.is_file():
            raise FileNotFoundError(f"Execution manifest missing: {manifest_path}")
        if _sha256(manifest_path) != str(execution_manifest_sha256).lower():
            raise ValueError("Execution-manifest SHA-256 mismatch.")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        checks = {
            "status": "success",
            "return_code": 0,
            "execution_id": execution_id,
            "independent_run_id": entry.get("independent_run_id"),
            "random_seed": entry.get("random_seed"),
            "protocol_sha256": entry.get("protocol_sha256"),
            "sequence_sha256": hashlib.sha256(
                str(candidate.sequence).upper().encode("utf-8")
            ).hexdigest(),
        }
        mismatches = [key for key, expected in checks.items() if manifest.get(key) != expected]
        if mismatches:
            raise ValueError("Execution-manifest provenance mismatch: " + ", ".join(mismatches))
        registered_pose_hashes = manifest.get("verified_pose_sha256")
        if not isinstance(registered_pose_hashes, dict) or _sha256(path) not in registered_pose_hashes.values():
            raise ValueError("Complex file hash is absent from the execution manifest's verified poses.")
        archive_value = manifest.get("archive_path")
        archive_hash = manifest.get("archive_sha256")
        if not archive_value or not archive_hash:
            raise ValueError("Execution manifest lacks archived-output provenance.")
        archive_path = Path(archive_value).expanduser().resolve()
        if not archive_path.is_file() or _sha256(archive_path) != str(archive_hash).lower():
            raise ValueError("HADDOCK output archive is missing or its hash changed.")

    target_id = str(entry.get("target_structure_id", "")).strip()
    receptor_chain = str(entry.get("receptor_chain", "A")).strip()
    peptide_chain = str(entry.get("peptide_chain", "B")).strip()

    if not target_id:
        raise ValueError("Target structure ID is missing")

    hotspots_raw = entry.get("hotspots")

    if not hotspots_raw:
        hotspots_raw = _default_legacy_hotspots(
            target_id,
            receptor_chain,
        )

    hotspots = _normalise_hotspots(hotspots_raw, receptor_chain)

    config = StructuralScoreConfig(
        hotspot_residues=hotspots,
    )

    result = score_peptide_receptor_complex(
        path,
        config=config,
        receptor_chain=receptor_chain,
        peptide_chain=peptide_chain,
        target_structure_id=target_id,
        structure_provenance=str(
            entry.get("provenance", "manifest_complex")
        ),
    )

    actual_sequence = str(result.get("peptide_sequence", "")).upper()

    if actual_sequence != str(candidate.sequence).upper():
        raise ValueError(
            "Candidate/complex sequence mismatch: "
            f"{candidate.sequence} != {actual_sequence}"
        )

    if result.get("score_valid") is not True:
        raise ValueError(
            "Structural validation failed: "
            + ";".join(result.get("eligibility_reasons", []))
        )

    metrics = result.get("structural_metrics", {})
    hotspot_rows = metrics.get("hotspot_metrics", [])

    observed = {
        row.get("hotspot"): float(row["minimum_distance_A"])
        for row in hotspot_rows
        if row.get("minimum_distance_A") is not None
    }

    required_labels = {hotspot.label() for hotspot in hotspots}
    missing = required_labels - set(observed)

    if missing:
        raise ValueError(
            "Missing configured hotspot measurements: "
            + ", ".join(sorted(missing))
        )

    hotspot_distances = {
        label: observed[label]
        for label in sorted(required_labels)
    }

    proximity = sum(
        math.exp(-max(distance - 4.0, 0.0) / 2.0)
        for distance in hotspot_distances.values()
    ) / len(hotspot_distances)

    score = float(result["target_aware_structural_score"])

    if not math.isfinite(score):
        raise ValueError("Structural score is not finite")

    try:
        fingerprint = _interface_fingerprint(
            path,
            receptor_chain,
            peptide_chain,
            cutoff_A=5.0,
        )
        fingerprint_error = None
    except Exception as exc:
        fingerprint = set()
        fingerprint_error = f"{type(exc).__name__}: {exc}"

    output = {
        "path": str(path),
        "target_structure_id": target_id,
        "score_version": result.get("score_version"),
        "complex_sha256": result.get(
            "structure_provenance", {}
        ).get("sha256"),
        "target_geometry_score": score,
        "hotspot_distances_A": hotspot_distances,
        "hotspot_proximity_score": proximity,
        "minimum_interchain_distance_A": metrics.get(
            "minimum_interchain_heavy_atom_distance_A"
        ),
        "severe_clash_count": metrics.get(
            "severe_clash_heavy_atom_pairs"
        ),
        "interface_residue_pairs_5A": metrics.get(
            "distinct_interchain_residue_pairs", {}
        ).get("5A"),
        "independent_run_id": entry.get("independent_run_id"),
        "random_seed": entry.get("random_seed"),
        "protocol_sha256": entry.get("protocol_sha256"),
        "execution_id": execution_id,
        "execution_manifest_sha256": execution_manifest_sha256,
        "_interface_fingerprint_5A": sorted(fingerprint),
        "fingerprint_error": fingerprint_error,
    }

    # Preserve compatibility with historical HERV-K reports.
    for label, distance in hotspot_distances.items():
        if label.endswith(":203:ARG"):
            output["arg203_min_A"] = distance
        elif label.endswith(":205:ARG"):
            output["arg205_min_A"] = distance

    return output


def _pareto_fronts(
    objectives: list[tuple[float, ...]],
) -> list[int]:
    remaining = set(range(len(objectives)))
    ranks = [-1] * len(objectives)
    rank = 0

    while remaining:
        front = []

        for i in sorted(remaining):
            dominated = False

            for j in remaining:
                if i == j:
                    continue

                a, b = objectives[j], objectives[i]

                if all(x >= y for x, y in zip(a, b)) and any(
                    x > y for x, y in zip(a, b)
                ):
                    dominated = True
                    break

            if not dominated:
                front.append(i)

        if not front:
            raise RuntimeError("Pareto front construction failed")

        for i in front:
            ranks[i] = rank
            remaining.remove(i)

        rank += 1

    return ranks



def _summarize_distribution(values: list[float]) -> dict[str, Any]:
    """Summarize variation across scored poses without calling them replicates."""
    ordered = sorted(float(value) for value in values)

    if not ordered:
        raise ValueError("Cannot summarize an empty pose distribution")

    n = len(ordered)
    middle = n // 2

    if n % 2:
        median = ordered[middle]
    else:
        median = (ordered[middle - 1] + ordered[middle]) / 2.0

    return {
        "count": n,
        "mean": sum(ordered) / n,
        "median": median,
        "minimum": ordered[0],
        "maximum": ordered[-1],
    }


def _summarize_pose_source_groups(
    complex_results: list[dict[str, Any]],
) -> dict[str, Any]:
    """Summarize pose metrics by source directory.

    Directory labels are descriptive file-origin strata only.
    They are not independent runs, experimental replicates,
    or evidence of structural accuracy.
    """
    grouped: dict[str, list[dict[str, Any]]] = {}

    for item in complex_results:
        raw_path = item.get("path")
        group = (
            Path(str(raw_path)).parent.name
            if raw_path else "unknown_source"
        )
        group = group or "unknown_source"
        grouped.setdefault(group, []).append(item)

    def metric_summary(
        poses: list[dict[str, Any]],
        key: str,
    ) -> dict[str, Any] | None:
        values = []

        for pose in poses:
            try:
                value = float(pose.get(key))
            except (TypeError, ValueError):
                continue

            if math.isfinite(value):
                values.append(value)

        return _summarize_distribution(values) if values else None

    summaries = []

    for group, poses in sorted(grouped.items()):
        hotspot_labels = sorted({
            label
            for pose in poses
            for label in (pose.get("hotspot_distances_A") or {})
        })

        hotspot_summaries = {}

        for label in hotspot_labels:
            values = []

            for pose in poses:
                distances = pose.get("hotspot_distances_A") or {}

                if label not in distances:
                    continue

                try:
                    value = float(distances[label])
                except (TypeError, ValueError):
                    continue

                if math.isfinite(value):
                    values.append(value)

            hotspot_summaries[label] = (
                _summarize_distribution(values) if values else None
            )

        summaries.append({
            "source_group": group,
            "pose_count": len(poses),
            "target_geometry_score": metric_summary(
                poses, "target_geometry_score"
            ),
            "hotspot_proximity_score": metric_summary(
                poses, "hotspot_proximity_score"
            ),
            "interface_residue_pairs_5A": metric_summary(
                poses, "interface_residue_pairs_5A"
            ),
            "minimum_interchain_distance_A": metric_summary(
                poses, "minimum_interchain_distance_A"
            ),
            "hotspot_distance_summaries_A": hotspot_summaries,
        })

    return {
        "grouping_basis": "pose_file_parent_directory_name",
        "independence_claimed": False,
        "interpretation": (
            "Descriptive file-origin grouping only. "
            "These groups are not independent runs, experimental "
            "replicates, or structural-accuracy validation."
        ),
        "group_count": len(summaries),
        "groups": summaries,
    }


def apply_p4_target_objective(
    candidates: list[Any],
    *,
    generation: int,
    output_directory: str | Path,
    pareto_enabled: bool | None = None,
    require_reproducibility: bool = False,
    selection_scope: str = "full_population",
    external_validation_gate: dict[str, Any] | None = None,
    allow_unvalidated_exploratory_selection: bool = False,
) -> tuple[list[Any], dict[str, Any]]:
    """Score registered complexes and optionally apply target-aware ordering.

    ``require_reproducibility`` is the strict live-optimization mode. Every
    candidate in the declared ``selection_scope`` must have valid target
    geometry plus a repeatability diagnostic from at least two independent
    docking runs (>=2 scored poses per run, unique seeds, and a shared protocol
    hash). The repeatability diagnostic remains uncalibrated and is never used
    as a ranking weight.

    ``selection_scope="shortlist"`` is a budget-limited treatment: only the
    predeclared ``p4_structure_shortlist`` is eligible for target-aware
    reordering. Candidates outside that shortlist keep their original P3
    positions. This avoids claiming complete-population optimization when only
    a subset has structural evidence, while still permitting a genuine,
    controlled target-aware selection effect within the evaluated subset.
    """
    if selection_scope not in {"full_population", "shortlist"}:
        raise ValueError("selection_scope must be 'full_population' or 'shortlist'.")

    if external_validation_gate is None and require_reproducibility:
        from src.p4_external_validation_gate import load_external_validation_gate
        external_validation_gate = load_external_validation_gate()
    external_validation_gate = external_validation_gate or {
        "status": "not_required_for_diagnostic_mode",
        "allow_target_selection": not require_reproducibility,
        "reason": "External transfer gate applies to target-aware selection mode.",
    }
    external_validation_passed = bool(
        external_validation_gate.get("allow_target_selection")
    ) and external_validation_gate.get("status") != "not_required_for_diagnostic_mode"
    exploratory_override_used = bool(
        require_reproducibility
        and allow_unvalidated_exploratory_selection
        and not external_validation_passed
    )
    external_validation_allows_selection = (
        external_validation_passed or exploratory_override_used
        if require_reproducibility else True
    )

    manifest_report = {
        "status": "not_configured",
        "matched_candidate_count": 0,
        "matched_complex_count": 0,
    }

    manifest_path = os.environ.get(
        "SILICON_VIROVORE_P4_TARGET_MANIFEST"
    )

    if manifest_path:
        manifest_report = _register_manifest(candidates, manifest_path)

    if pareto_enabled is None:
        pareto_enabled = (
            os.environ.get(
                "SILICON_VIROVORE_P4_TARGET_PARETO", "1"
            ) == "1"
        )

    records = []
    eligible = []
    selection_eligible = []
    reliability_count = 0

    for candidate in candidates:
        meta = candidate.metadata
        raw_entries = meta.get("target_complexes")

        if raw_entries is None:
            raw_path = meta.get("target_complex_path")
            raw_entries = [{"path": raw_path}] if raw_path else []

        if isinstance(raw_entries, dict):
            raw_entries = [raw_entries]

        registered_complex_count = (
            len(raw_entries)
            if isinstance(raw_entries, (list, tuple))
            else int(bool(raw_entries))
        )
        complex_source = meta.get(
            "target_complex_source",
            "candidate_metadata" if registered_complex_count else "none",
        )

        entry_results = []
        errors = []

        for raw in raw_entries:
            try:
                if isinstance(raw, str):
                    raw = {"path": raw}

                target_id = raw.get(
                    "target_structure_id",
                    meta.get("target_structure_id", "9SYA"),
                )

                receptor_chain = raw.get(
                    "receptor_chain",
                    meta.get("target_receptor_chain", "A"),
                )

                peptide_chain = raw.get(
                    "peptide_chain",
                    meta.get("target_peptide_chain", "B"),
                )

                hotspots = raw.get(
                    "hotspots",
                    meta.get("target_hotspots"),
                )

                if not hotspots:
                    hotspots = _default_legacy_hotspots(
                        str(target_id),
                        str(receptor_chain),
                    )

                entry = {
                    "path": raw["path"],
                    "sha256": raw.get("sha256"),
                    "target_structure_id": target_id,
                    "receptor_chain": receptor_chain,
                    "peptide_chain": peptide_chain,
                    "hotspots": hotspots,
                    "provenance": raw.get(
                        "provenance",
                        meta.get(
                            "target_complex_provenance",
                            "external_complex",
                        ),
                    ),
                    "independent_run_id": raw.get("independent_run_id"),
                    "random_seed": raw.get("random_seed"),
                    "protocol_sha256": raw.get("protocol_sha256"),
                    "execution_id": raw.get("execution_id"),
                    "execution_manifest_path": raw.get("execution_manifest_path"),
                    "execution_manifest_sha256": raw.get("execution_manifest_sha256"),
                    "_require_execution_manifest": bool(require_reproducibility),
                }

                entry_results.append(_score_complex(candidate, entry))

            except Exception as exc:
                errors.append(f"{type(exc).__name__}: {exc}")

        if entry_results and not errors:
            hotspot_labels = sorted(
                set.intersection(*[
                    set(item["hotspot_distances_A"])
                    for item in entry_results
                ])
            )

            if not hotspot_labels:
                errors.append("No common configured hotspots across complexes")

        if entry_results and not errors:
            geometry = sum(
                item["target_geometry_score"] for item in entry_results
            ) / len(entry_results)

            mean_hotspots = {
                label: sum(
                    item["hotspot_distances_A"][label]
                    for item in entry_results
                ) / len(entry_results)
                for label in hotspot_labels
            }

            proximity = sum(
                item["hotspot_proximity_score"]
                for item in entry_results
            ) / len(entry_results)

            reliability = _reliability_diagnostic(entry_results)

            meta["p4_target_evidence"] = {
                "version": TARGET_VERSION,
                "status": "success",
                "target_structure_id": entry_results[0]["target_structure_id"],
                "pose_count": len(entry_results),
                "target_geometry_score_mean": geometry,
                "target_geometry_score_summary": _summarize_distribution([
                    item["target_geometry_score"]
                    for item in entry_results
                ]),
                "hotspot_distances_A_mean": mean_hotspots,
                "hotspot_distance_summaries_A": {
                    label: _summarize_distribution([
                        item["hotspot_distances_A"][label]
                        for item in entry_results
                    ])
                    for label in hotspot_labels
                },
                "hotspot_proximity_score": proximity,
                "hotspot_proximity_score_summary": _summarize_distribution([
                    item["hotspot_proximity_score"]
                    for item in entry_results
                ]),
                "arg203_min_A_mean": mean_hotspots.get(
                    "A:203:ARG"
                ),
                "arg205_min_A_mean": mean_hotspots.get(
                    "A:205:ARG"
                ),
                "reliability_calibrated": False,
                "pose_source_groups": _summarize_pose_source_groups(
                    entry_results
                ),
                "complex_results": [
                    {
                        key: value
                        for key, value in item.items()
                        if key != "_interface_fingerprint_5A"
                    }
                    for item in entry_results
                ],
            }

            meta["p4_reliability"] = reliability

            reliability_proxy = None
            try:
                value = float(reliability.get("reliability_proxy"))
                # Remove floating-point noise at the normalized [0, 1]
                # boundary without changing the measured Jaccard statistic.
                tolerance = 1e-12
                if math.isfinite(value) and -tolerance <= value <= 1.0 + tolerance:
                    reliability_proxy = min(1.0, max(0.0, value))
            except (TypeError, ValueError):
                pass

            reproducibility_valid = (
                reliability.get("status") == "diagnostic_only"
                and reliability_proxy is not None
                and int(reliability.get("n_independent_runs", 0)) >= 2
            )
            if reproducibility_valid:
                reliability_count += 1

            selection_valid = (
                not require_reproducibility or reproducibility_valid
            )
            meta["p4_target_selection_eligible"] = bool(selection_valid)
            meta["p4_target_selection_reason"] = (
                "Target geometry and independent-run repeatability diagnostic present"
                if selection_valid and require_reproducibility
                else "Target geometry valid; repeatability not required in diagnostic mode"
                if selection_valid
                else "Target selection withheld: independent-run repeatability evidence incomplete"
            )
            meta["p4_reliability_proxy_for_reporting_only"] = reliability_proxy

            eligible.append(candidate)
            if selection_valid:
                selection_eligible.append(candidate)

        else:
            status = (
                "missing_complex"
                if not raw_entries
                else "invalid_complex"
            )

            meta["p4_target_evidence"] = {
                "version": TARGET_VERSION,
                "status": status,
                "errors": errors,
                "reliability_calibrated": False,
            }

            meta["p4_reliability"] = {
                "status": "unavailable_without_valid_target_evidence",
                "calibrated": False,
                "feedback_active": False,
            }
            meta["p4_target_selection_eligible"] = False
            meta["p4_target_selection_reason"] = (
                "Target geometry evidence unavailable or invalid"
            )
            meta["p4_reliability_proxy_for_reporting_only"] = None

        records.append({
            "sequence": candidate.sequence,
            "candidate_id": meta.get("candidate_id"),
            "target_complex_source": complex_source,
            "registered_complex_count": registered_complex_count,
            "complex_generation": meta.get("p4_complex_generation"),
            "structure_shortlist": meta.get("p4_structure_shortlist", {}),
            **meta["p4_target_evidence"],
            "reliability": meta["p4_reliability"],
            "target_selection_eligible": meta.get(
                "p4_target_selection_eligible", False
            ),
            "target_selection_reason": meta.get(
                "p4_target_selection_reason"
            ),
            "reliability_proxy_for_reporting_only": meta.get(
                "p4_reliability_proxy_for_reporting_only"
            ),
        })

    coverage = len(eligible) / len(candidates) if candidates else 0.0
    selection_coverage = (
        len(selection_eligible) / len(candidates) if candidates else 0.0
    )

    reordered = False
    reason = "Target evidence recorded; legacy ranking retained."
    shortlist_candidates = [
        candidate for candidate in candidates
        if (candidate.metadata if isinstance(candidate.metadata, dict) else {})
        .get("p4_structure_shortlist", {}).get("selected") is True
    ]
    scope_candidates = (
        shortlist_candidates if selection_scope == "shortlist" else list(candidates)
    )
    scope_ids = {id(candidate) for candidate in scope_candidates}
    eligible_scope = [candidate for candidate in eligible if id(candidate) in scope_ids]
    selection_eligible_scope = [
        candidate for candidate in selection_eligible if id(candidate) in scope_ids
    ]

    if pareto_enabled and require_reproducibility and not external_validation_allows_selection:
        reason = (
            "Target-aware Pareto selection withheld by the independent external-validation gate: "
            + str(external_validation_gate.get("reason", "no qualifying validation evidence"))
        )
    elif pareto_enabled:
        scores_valid = bool(scope_candidates) and all(
            getattr(candidate, "c_score", None) is not None
            and math.isfinite(float(candidate.c_score))
            for candidate in scope_candidates
        )

        required_candidates = (
            selection_eligible_scope if require_reproducibility else eligible_scope
        )
        if scope_candidates and len(required_candidates) == len(scope_candidates) and scores_valid:
            objectives = []
            for candidate in scope_candidates:
                target = candidate.metadata["p4_target_evidence"]
                fitness = max(0.0, min(1.0, float(candidate.c_score)))
                objectives.append((
                    fitness,
                    float(target["target_geometry_score_mean"]),
                    float(target["hotspot_proximity_score"]),
                ))

            fronts = _pareto_fronts(objectives)
            scope_order = {id(candidate): i for i, candidate in enumerate(scope_candidates)}
            for candidate, front in zip(scope_candidates, fronts):
                candidate.metadata["p4_target_pareto_front"] = front

            if require_reproducibility:
                ordered_scope = sorted(
                    scope_candidates,
                    key=lambda candidate: (
                        candidate.metadata["p4_target_pareto_front"],
                        -float(candidate.metadata["p4_target_evidence"]["target_geometry_score_mean"]),
                        -float(candidate.metadata["p4_target_evidence"]["hotspot_proximity_score"]),
                        -max(0.0, min(1.0, float(candidate.c_score))),
                        -float(candidate.overall_score),
                        scope_order[id(candidate)],
                    ),
                )
                selection_method = "pareto_front_then_target_geometry_then_hotspot_proximity"
            else:
                ordered_scope = sorted(
                    scope_candidates,
                    key=lambda candidate: (
                        candidate.metadata["p4_target_pareto_front"],
                        -float(candidate.overall_score),
                        scope_order[id(candidate)],
                    ),
                )
                selection_method = "pareto_front_then_legacy_overall_score"

            if selection_scope == "shortlist":
                # Replace only the original shortlist slots. No candidate with
                # missing structural evidence can be promoted or demoted.
                positions = [i for i, c in enumerate(candidates) if id(c) in scope_ids]
                reordered_candidates = list(candidates)
                for position, selected_candidate in zip(positions, ordered_scope):
                    reordered_candidates[position] = selected_candidate
                candidates = reordered_candidates
            else:
                candidates = ordered_scope

            for rank, candidate in enumerate(ordered_scope, start=1):
                candidate.metadata["p4_target_scope_rank"] = rank
                candidate.metadata["p4_target_selection_method"] = selection_method
                candidate.metadata["p4_target_selection_applied"] = True
            reordered = True
            reason = (
                f"All {len(scope_candidates)} candidates in the {selection_scope} scope have valid "
                "target geometry and independent-run repeatability diagnostics; target-aware "
                "Pareto ordering applied within the evaluated scope."
                if require_reproducibility
                else f"All {len(scope_candidates)} candidates in the {selection_scope} scope have valid "
                "target geometry; Pareto ordering applied within the evaluated scope."
            )
        elif not scope_candidates:
            reason = (
                "Target optimization withheld: the declared selection scope is empty."
            )
        elif require_reproducibility and len(selection_eligible_scope) != len(scope_candidates):
            reason = (
                "Target optimization withheld: every candidate in the declared "
                "selection scope requires valid target geometry and repeatability "
                "diagnostics from at least two independent runs with at least two poses per run."
            )
        elif len(eligible_scope) != len(scope_candidates):
            reason = (
                "Pareto selection withheld: complete valid target evidence is required "
                "for every candidate in the declared selection scope."
            )
        else:
            reason = (
                "Pareto selection withheld: every candidate in the declared scope "
                "requires a finite sequence-level c_score."
            )

    report = {
        "generation": generation,
        "version": TARGET_VERSION,
        "candidate_count": len(candidates),
        "selection_scope": selection_scope,
        "selection_scope_candidate_count": len(scope_candidates),
        "population_target_evidence_coverage": coverage,
        "population_target_selection_coverage": selection_coverage,
        "target_evidence_count": len(eligible_scope),
        "target_evidence_coverage": (
            len(eligible_scope) / len(scope_candidates) if scope_candidates else 0.0
        ),
        "target_selection_eligible_count": len(selection_eligible_scope),
        "target_selection_coverage": (
            len(selection_eligible_scope) / len(scope_candidates) if scope_candidates else 0.0
        ),
        "require_reproducibility_for_selection": bool(require_reproducibility),
        "external_validation_gate": external_validation_gate,
        "external_validation_allows_selection": bool(external_validation_passed),
        "exploratory_override_used": bool(exploratory_override_used),
        "selection_validation_status": (
            "EXPLORATORY_NOT_VALIDATED" if exploratory_override_used
            else "EXTERNAL_TRANSFER_GATE_PASSED" if external_validation_passed
            else "NOT_REQUIRED_DIAGNOSTIC_MODE" if not require_reproducibility
            else "VALIDATION_NOT_PASSED"
        ),
        "reliability_diagnostic_candidate_count": reliability_count,
        "reliability_used_for_ranking": False,
        "pareto_requested": bool(pareto_enabled),
        "pareto_reordered_candidates": reordered,
        "selection_method": (
            "pareto_front_then_target_geometry_then_hotspot_proximity"
            if reordered and require_reproducibility
            else "pareto_front_then_legacy_overall_score"
            if reordered
            else "withheld"
        ),
        "reason": reason,
        "manifest": manifest_report,
        "reliability_diagnostic_candidate_count": reliability_count,
        "reliability_feedback_active": False,
        "reliability_calibrated": False,
        "legacy_fitness_scores_modified": False,
        "candidates": records,
    }

    out = Path(output_directory) / "p4_target_evidence"
    out.mkdir(parents=True, exist_ok=True)

    report_path = out / f"generation_{generation:04d}.json"
    report_path.write_text(
        json.dumps(report, indent=2, default=str) + "\n",
        encoding="utf-8",
    )

    report["report_path"] = str(report_path)

    return candidates, report
