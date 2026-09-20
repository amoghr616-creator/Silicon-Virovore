"""ALE mutation-policy interface.

The native generator executes the policy so baseline and guided conditions use
the same fixed mutation-event budget. This module resolves the frozen ARISE
model into a small, auditable policy description.
"""

from __future__ import annotations

from typing import Any


GUIDANCE_MODES = frozenset({"uniform", "hotspot_only", "hotspot_biased"})


def resolve_mutation_policy(
    mode: str,
    hotspot_model: dict[str, Any] | None,
    *,
    enabled: bool,
) -> dict[str, Any]:
    """Resolve a frozen hotspot model into an executable ALE policy."""

    if mode not in GUIDANCE_MODES:
        raise ValueError(
            "mutation_guidance_mode must be one of: "
            "uniform, hotspot_only, hotspot_biased"
        )

    selected = (hotspot_model or {}).get("selected_hotspot")
    if not enabled or mode == "uniform" or not selected:
        return {
            "mode": "uniform",
            "requested_mode": mode,
            "guided": False,
            "hotspot_start": None,
            "hotspot_end": None,
            "hotspot_discovery_generation": None,
        }

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
    }


def mutation_positions(parent: str, sequence: str) -> list[int]:
    """Return zero-based positions changed from the recorded parent."""

    return [
        index
        for index, (old, new) in enumerate(zip(parent, sequence))
        if old != new
    ]
