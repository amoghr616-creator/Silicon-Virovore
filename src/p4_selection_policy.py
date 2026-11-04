"""Cross-generation comparison using the active P4 objective.

All three P4 objectives are maximized:
    1. bounded sequence-level C fitness
    2. target geometry score
    3. hotspot proximity score

Pareto dominance is evaluated first. Incomparable candidates use the
same tie-break order as the active P4 integration.
"""

from collections.abc import Mapping
import math


_MISSING = object()

_NESTED_FIELDS = (
    "target_evidence",
    "p4_target_evidence",
    "p4_evidence",
    "selection_evidence",
)

_FIELD_NAMES = {
    "fitness": (
        "bounded_c_score",
        "c_score_fitness",
        "fitness",
        "c_score",
        "sequence_fitness",
    ),
    "geometry": (
        "target_geometry_score_mean",
        "target_geometry_score",
    ),
    "hotspot": (
        "hotspot_proximity_score",
        "hotspot_score",
    ),
}


def _read(obj, key):
    if isinstance(obj, Mapping):
        return obj.get(key, _MISSING)
    return getattr(obj, key, _MISSING)


def _find_value(candidate, names):
    containers = [candidate]

    for nested_name in _NESTED_FIELDS:
        nested = _read(candidate, nested_name)
        if nested is not _MISSING and nested is not None:
            containers.append(nested)

    for container in containers:
        for name in names:
            value = _read(container, name)
            if value is not _MISSING and value is not None:
                return value

    return _MISSING


def _finite_number(value):
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None

    return number if math.isfinite(number) else None


def p4_objective_vector(candidate):
    """Return (fitness, geometry, hotspot), or None if evidence is incomplete."""

    eligible = _find_value(candidate, ("selection_eligible",))

    if eligible is not _MISSING and eligible is not None:
        if not bool(eligible):
            return None

    fitness = _finite_number(
        _find_value(candidate, _FIELD_NAMES["fitness"])
    )
    geometry = _finite_number(
        _find_value(candidate, _FIELD_NAMES["geometry"])
    )
    hotspot = _finite_number(
        _find_value(candidate, _FIELD_NAMES["hotspot"])
    )

    if fitness is None or geometry is None or hotspot is None:
        return None

    # P4's sequence-level fitness is bounded to [0, 1].
    # Reject malformed values rather than silently clipping them.
    if not 0.0 <= fitness <= 1.0:
        return None

    return fitness, geometry, hotspot


def _dominates(a, b):
    """True when a Pareto-dominates b, assuming higher is better."""
    return (
        all(x >= y for x, y in zip(a, b))
        and any(x > y for x, y in zip(a, b))
    )


def p4_candidate_preference(candidate, incumbent):
    """Return True/False for P4 preference, or None if neither has P4 data."""

    candidate_vector = p4_objective_vector(candidate)
    incumbent_vector = p4_objective_vector(incumbent)

    if candidate_vector is None and incumbent_vector is None:
        return None

    if candidate_vector is None:
        return False

    if incumbent_vector is None:
        return True

    if _dominates(candidate_vector, incumbent_vector):
        return True

    if _dominates(incumbent_vector, candidate_vector):
        return False

    # Incomparable candidates: reproduce the active P4 tie-break order.
    # Geometry, hotspot proximity, then bounded sequence fitness.
    candidate_key = (
        candidate_vector[1],
        candidate_vector[2],
        candidate_vector[0],
    )
    incumbent_key = (
        incumbent_vector[1],
        incumbent_vector[2],
        incumbent_vector[0],
    )

    if candidate_key != incumbent_key:
        return candidate_key > incumbent_key

    # Legacy score is only the final tie-breaker.
    candidate_legacy = _finite_number(
        _read(candidate, "overall_score")
    )
    incumbent_legacy = _finite_number(
        _read(incumbent, "overall_score")
    )

    candidate_legacy = (
        candidate_legacy if candidate_legacy is not None else float("-inf")
    )
    incumbent_legacy = (
        incumbent_legacy if incumbent_legacy is not None else float("-inf")
    )

    return candidate_legacy > incumbent_legacy