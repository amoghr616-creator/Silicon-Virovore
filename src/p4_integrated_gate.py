"""P4 candidate-level integration gate.

Combines target-specific structural evidence, calibrated reliability
evidence, and family-validated epistasis evidence.

No synthetic measurements, inferred affinity, or arbitrary weighted
scores are generated. Missing evidence causes abstention.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections import defaultdict
from typing import Any


SEQUENCE_KEYS = (
    "candidate_sequence",
    "peptide_sequence",
    "sequence",
)

TARGET_KEYS = (
    "target_specific_score",
    "target_score",
    "docking_score",
)

UNCERTAINTY_KEYS = (
    "calibrated_expected_abs_error_kcal_mol",
    "predicted_error_upper_kcal_mol",
    "prediction_uncertainty",
    "uncertainty",
)

EPISTASIS_KEYS = (
    "predicted_delta_epsilon",
    "delta_epsilon_pred",
    "epistasis_prediction",
)


def canonical_sequence(record: dict) -> str | None:
    for key in SEQUENCE_KEYS:
        value = record.get(key)
        if isinstance(value, str):
            sequence = re.sub(r"\s+", "", value).upper()
            if sequence and re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", sequence):
                return sequence
    return None


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def number(record: dict, keys: tuple) -> float | None:
    for key in keys:
        try:
            value = float(record.get(key))
            if math.isfinite(value):
                return value
        except (TypeError, ValueError):
            pass
    return None


def flag(value: Any) -> bool:
    return value is True or (
        isinstance(value, str)
        and value.strip().lower() in {"true", "yes", "1"}
    )


def family_validation(record: dict) -> bool:
    split = str(
        record.get("validation_split_unit", record.get("split_unit", ""))
    ).lower().replace("-", "_")

    status = str(record.get("validation_status", "")).lower()

    return (
        split in {"family", "family_lofo", "lofo"}
        or status in {
            "family_held_out",
            "family_held_out_validated",
            "lofo_validated",
        }
    )


def index_evidence(rows: list[dict]) -> dict[str, dict]:
    indexed = {}
    duplicates = defaultdict(int)

    for row in rows:
        sequence = canonical_sequence(row)
        if not sequence:
            continue

        digest = sequence_hash(sequence)
        duplicates[digest] += 1

        # Duplicate evidence is not silently overwritten.
        if digest not in indexed:
            indexed[digest] = row
        else:
            indexed[digest] = {
                "__duplicate_evidence__": True,
                "__records__": [indexed[digest], row],
            }

    return indexed


def assess_target(record: dict | None) -> dict:
    if not record:
        return {"valid": False, "reason": "No candidate-matched target evidence"}

    if record.get("__duplicate_evidence__"):
        return {"valid": False, "reason": "Ambiguous duplicate target records"}

    def residue_numbers(value: Any) -> set[str]:
        text = str(value)
        found = set()
        token = ""
        for char in text:
            if char.isdigit():
                token += char
            elif token:
                found.add(token)
                token = ""
        if token:
            found.add(token)
        return found

    # Native P4 schema: geometry and hotspot evidence attached to the
    # exact candidate sequence by apply_p4_target_objective().
    native = record.get("p4_target_evidence")
    if isinstance(native, dict):
        pdb = str(
            native.get("target_structure_id", "")
        ).upper().strip()
        geometry = number(native, ("target_geometry_score_mean",))
        distances = native.get("hotspot_distances_A_mean", {})
        distances = distances if isinstance(distances, dict) else {}

        observed = set()
        for label, distance in distances.items():
            if number({"value": distance}, ("value",)) is not None:
                observed.update(residue_numbers(label))

        required = {"203", "205"}
        missing = sorted(required - observed)
        valid = (
            native.get("status") == "success"
            and pdb == "9SYA"
            and geometry is not None
            and not missing
        )

        reasons = []
        if native.get("status") != "success":
            reasons.append("Native target scoring did not succeed")
        if pdb != "9SYA":
            reasons.append("Target is not 9SYA")
        if geometry is None:
            reasons.append("Finite target geometry score missing")
        if missing:
            reasons.append(
                "Missing Arg203/Arg205 geometry: " + ",".join(missing)
            )

        return {
            "valid": valid,
            "reason": (
                "Exact-sequence 9SYA structural geometry validated"
                if valid else "; ".join(reasons)
            ),
            "pdb_id": pdb,
            "target_score": geometry,
            "evidence_type": "target_specific_structural_geometry",
            "hotspots_verified": sorted(required - set(missing)),
            "affinity_inferred": False,
        }

    # Original evidence schema retained for backward compatibility.
    pdb = str(
        record.get("target_pdb_id", record.get("pdb_id", ""))
    ).upper().replace(" ", "")

    target_specific = flag(record.get("docking_score_is_target_specific"))
    anchor = number(record, ("strongest_anchor_delta_g",))
    score = number(record, TARGET_KEYS)

    hotspot_value = record.get(
        "hotspot_residues",
        record.get("target_hotspots", record.get("hotspots", [])),
    )
    hotspot_numbers = residue_numbers(hotspot_value)

    required = {"203", "205"}
    missing = sorted(required - hotspot_numbers)
    valid = (
        pdb == "9SYA"
        and target_specific
        and anchor is not None
        and not missing
    )

    reasons = []
    if pdb != "9SYA":
        reasons.append("Target is not explicitly identified as 9SYA")
    if not target_specific:
        reasons.append("Target-specific scoring flag absent")
    if anchor is None:
        reasons.append("Strongest anchor delta missing")
    if missing:
        reasons.append("Missing hotspot evidence: " + ",".join(missing))

    return {
        "valid": valid,
        "reason": "Target evidence validated" if valid else "; ".join(reasons),
        "pdb_id": pdb,
        "target_score": score,
        "strongest_anchor_delta_g": anchor,
        "hotspots_verified": sorted(required - set(missing)),
        "affinity_inferred": False,
    }


def assess_reliability(record: dict | None) -> dict:
    if not record:
        return {"valid": False, "reason": "No candidate-matched reliability prediction"}

    if record.get("__duplicate_evidence__"):
        return {"valid": False, "reason": "Ambiguous duplicate reliability records"}

    dimension = record.get("feature_dimension")
    uncertainty = number(record, UNCERTAINTY_KEYS)

    calibration = str(record.get("calibration_status", "")).lower()
    applicability = str(
        record.get("applicability_domain_status", "")
    ).lower()

    calibrated = (
        calibration in {"calibrated", "validated", "family_calibrated"}
        or flag(record.get("uncertainty_calibrated"))
    )

    in_domain = applicability in {"in_domain", "accepted", "within_domain"}

    valid = (
        dimension == 410
        and uncertainty is not None
        and calibrated
        and in_domain
        and family_validation(record)
    )

    reasons = []
    if dimension != 410:
        reasons.append("410-dimensional feature provenance missing")
    if uncertainty is None:
        reasons.append("Candidate uncertainty estimate missing")
    if not calibrated:
        reasons.append("Calibration not verified")
    if not in_domain:
        reasons.append("Candidate applicability domain not verified")
    if not family_validation(record):
        reasons.append("Family-held-out validation provenance missing")

    return {
        "valid": valid,
        "reason": "Reliability evidence validated" if valid else "; ".join(reasons),
        "feature_dimension": dimension,
        "uncertainty": uncertainty,
        "applicability_domain": applicability,
    }


def assess_epistasis(record: dict | None) -> dict:
    if not record:
        return {"valid": False, "reason": "No candidate-matched epistasis prediction"}

    if record.get("__duplicate_evidence__"):
        return {"valid": False, "reason": "Ambiguous duplicate epistasis records"}

    delta = number(record, EPISTASIS_KEYS)
    model_hash = str(
        record.get("epistasis_model_sha256", record.get("model_sha256", ""))
    )

    family_count = record.get(
        "validation_family_count",
        record.get("n_validation_families", 0),
    )

    try:
        family_count = int(family_count)
    except (TypeError, ValueError):
        family_count = 0

    valid = (
        delta is not None
        and bool(re.fullmatch(r"[0-9a-fA-F]{64}", model_hash))
        and family_count >= 10
        and family_validation(record)
    )

    reasons = []
    if delta is None:
        reasons.append("Candidate epistasis prediction missing")
    if not re.fullmatch(r"[0-9a-fA-F]{64}", model_hash):
        reasons.append("Epistasis model provenance hash missing")
    if family_count < 10:
        reasons.append("Insufficient documented validation families")
    if not family_validation(record):
        reasons.append("Family-held-out validation not documented")

    return {
        "valid": valid,
        "reason": "Epistasis evidence validated" if valid else "; ".join(reasons),
        "predicted_delta_epsilon": delta,
        "validation_family_count": family_count,
        "model_sha256": model_hash or None,
    }


def adjudicate_candidate(
    candidate: dict,
    target_index: dict[str, dict],
    reliability_index: dict[str, dict],
    epistasis_index: dict[str, dict],
) -> dict:
    sequence = canonical_sequence(candidate)

    if not sequence:
        return {
            "candidate_id": candidate.get("candidate_id", "unknown"),
            "decision": "ABSTAIN",
            "reason": "Candidate sequence missing or invalid",
        }

    digest = sequence_hash(sequence)

    target = assess_target(target_index.get(digest))
    reliability = assess_reliability(reliability_index.get(digest))
    epistasis = assess_epistasis(epistasis_index.get(digest))

    complete = (
        target["valid"]
        and reliability["valid"]
        and epistasis["valid"]
    )

    return {
        "candidate_id": candidate.get(
            "candidate_id", candidate.get("uid", digest[:12])
        ),
        "sequence_sha256": digest,
        "sequence_length": len(sequence),
        "decision": (
            "ELIGIBLE_FOR_MULTI_OBJECTIVE_REVIEW"
            if complete else "ABSTAIN"
        ),
        "target_optimization": target,
        "reliability": reliability,
        "epistasis": epistasis,
        "scalar_fitness_generated": False,
        "pruned": False,
    }


def adjudicate_population(
    candidates: list[dict],
    target_rows: list[dict],
    reliability_rows: list[dict],
    epistasis_rows: list[dict],
) -> list[dict]:
    target_index = index_evidence(target_rows)
    reliability_index = index_evidence(reliability_rows)
    epistasis_index = index_evidence(epistasis_rows)

    return [
        adjudicate_candidate(
            candidate,
            target_index,
            reliability_index,
            epistasis_index,
        )
        for candidate in candidates
    ]

def adjudicate_candidate_objects(candidates: list[Any]) -> list[dict]:
    """Adjudicate live Candidate objects using their attached P4 evidence.

    Missing or uncalibrated channels abstain. This gate does not generate
    scalar fitness and does not remove candidates from the legacy optimizer.
    """
    candidate_rows = []
    target_by_sequence = {}
    reliability_by_sequence = {}
    epistasis_by_sequence = {}

    def add_evidence(index: dict, sequence: str, record: dict) -> None:
        previous = index.get(sequence)
        if previous is None:
            index[sequence] = record
        elif previous != record:
            if isinstance(previous, list):
                previous.append(record)
            else:
                index[sequence] = [previous, record]

    for candidate in candidates:
        sequence = str(getattr(candidate, "sequence", "")).upper()
        candidate_id = getattr(candidate, "metadata", {}).get(
            "candidate_id", "unknown"
        )
        candidate_rows.append({
            "sequence": sequence,
            "candidate_id": candidate_id,
        })

        metadata = getattr(candidate, "metadata", {})
        if not isinstance(metadata, dict) or not sequence:
            continue

        target = metadata.get("p4_target_evidence")
        if isinstance(target, dict) and target:
            add_evidence(
                target_by_sequence,
                sequence,
                {"sequence": sequence, "p4_target_evidence": target},
            )

        reliability = metadata.get("p4_reliability")
        if isinstance(reliability, dict) and reliability:
            add_evidence(
                reliability_by_sequence,
                sequence,
                {"sequence": sequence, **reliability},
            )

        epistasis = metadata.get("p4_epistasis_evidence")
        if isinstance(epistasis, dict) and epistasis:
            add_evidence(
                epistasis_by_sequence,
                sequence,
                {"sequence": sequence, **epistasis},
            )

    def flatten(index: dict) -> list[dict]:
        rows = []
        for value in index.values():
            rows.extend(value if isinstance(value, list) else [value])
        return rows

    return adjudicate_population(
        candidate_rows,
        flatten(target_by_sequence),
        flatten(reliability_by_sequence),
        flatten(epistasis_by_sequence),
    )



def recommendation_eligible_candidates(candidates: list[Any]) -> list[Any]:
    """Return evidence-qualified candidates in their original ranking order.

    This creates a separate review shortlist. It does not modify optimizer
    fitness, candidate ranking, or evolutionary parent selection.
    """
    eligible = []
    for candidate in candidates:
        metadata = getattr(candidate, "metadata", None)
        if not isinstance(metadata, dict):
            continue
        gate = metadata.get("p4_integrated_gate")
        if (
            isinstance(gate, dict)
            and gate.get("decision") == "ELIGIBLE_FOR_MULTI_OBJECTIVE_REVIEW"
        ):
            eligible.append(candidate)
    return eligible
