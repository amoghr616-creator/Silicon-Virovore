"""ALE mutation-policy interface.

The native generator executes the policy so baseline and guided conditions use
the same fixed mutation-event budget. This module resolves the frozen ARISE
model into a small, auditable policy description.
"""

from __future__ import annotations

from typing import Any


GUIDANCE_MODES = frozenset({
    "uniform",
    "hotspot_only",
    "hotspot_biased",
    "hotspot_biased_annealed",
})

HOTSPOT_PROBABILITY_MIN = 0.25
HOTSPOT_PROBABILITY_MAX = 0.70
FIXED_HOTSPOT_PROBABILITY = 0.75


def annealed_hotspot_probability(
    generation: int,
    total_generations: int,
) -> float:
    """Linearly increase hotspot-biased mutation from 0.25 to 0.70.

    The schedule is declared before an experiment and depends only on the
    generation index, not on observed results. A one-generation run uses the
    lower bound; a run's last generation uses the upper bound.
    """
    if total_generations <= 1:
        return HOTSPOT_PROBABILITY_MIN
    progress = (generation - 1) / (total_generations - 1)
    progress = max(0.0, min(1.0, progress))
    return HOTSPOT_PROBABILITY_MIN + progress * (
        HOTSPOT_PROBABILITY_MAX - HOTSPOT_PROBABILITY_MIN
    )


def resolve_mutation_policy(
    mode: str,
    hotspot_model: dict[str, Any] | None,
    *,
    enabled: bool,
    generation: int = 1,
    total_generations: int = 10,
) -> dict[str, Any]:
    """Resolve a frozen hotspot model into an executable ALE policy."""

    if mode not in GUIDANCE_MODES:
        raise ValueError(
            "mutation_guidance_mode must be one of: "
            "uniform, hotspot_only, hotspot_biased, "
            "hotspot_biased_annealed"
        )

    selected = (hotspot_model or {}).get("selected_hotspot")
    significance_validated = (
        (hotspot_model or {}).get("hotspot_significance_validated") is True
    )
    if (
        not enabled
        or mode == "uniform"
        or not selected
        or not significance_validated
    ):
        return {
            "mode": "uniform",
            "requested_mode": mode,
            "guided": False,
            "guidance_disabled_reason": (
                "hotspot_significance_unvalidated"
                if selected and not significance_validated
                else "hotspot_unavailable_or_guidance_disabled"
            ),
            "hotspot_start": None,
            "hotspot_end": None,
            "hotspot_discovery_generation": None,
            "hotspot_probability": 0.0,
        }

    if mode == "hotspot_only":
        hotspot_probability = 1.0
    elif mode == "hotspot_biased_annealed":
        hotspot_probability = annealed_hotspot_probability(
            generation,
            total_generations,
        )
    else:
        hotspot_probability = FIXED_HOTSPOT_PROBABILITY

    return {
        "mode": mode,
        "requested_mode": mode,
        "guided": True,
        "hotspot_start": int(selected["start_position"]),
        "hotspot_end": int(selected["end_position"]),
        "hotspot_discovery_generation": hotspot_model.get(
            "discovery_generation"
        ),
        "hotspot_score": selected.get("hotspot_score"),
        "held_out_association": selected.get("held_out_association"),
        "bootstrap_stability": selected.get("bootstrap_stability"),
        "correlation_length": hotspot_model.get("correlation_length"),
        "hotspot_probability": hotspot_probability,
    }


def mutation_positions(parent: str, sequence: str) -> list[int]:
    """Return zero-based positions changed from the recorded parent."""

    return [
        index
        for index, (old, new) in enumerate(zip(parent, sequence))
        if old != new
    ]
