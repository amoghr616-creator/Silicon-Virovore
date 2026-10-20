import pytest

from src.ale_policy import resolve_mutation_policy


@pytest.fixture
def unvalidated_hotspot():
    return {
        "selected_hotspot": {
            "start_position": 4,
            "end_position": 6,
            "window_length": 2,
            "hotspot_score": 0.17,
            "held_out_association": None,
            "bootstrap_stability": 0.0,
        },
        "hotspot_significance_validated": False,
        "discovery_generation": 1,
        "correlation_length": 2,
    }


def test_existing_annealed_mode_rejects_unvalidated_hotspot(
    unvalidated_hotspot,
):
    policy = resolve_mutation_policy(
        "hotspot_biased_annealed",
        unvalidated_hotspot,
        enabled=True,
        generation=2,
        total_generations=3,
    )

    assert policy["mode"] == "uniform"
    assert policy["guided"] is False
    assert policy["hotspot_probability"] == 0.0


def test_exploratory_mode_uses_hotspot_without_claiming_validation(
    unvalidated_hotspot,
):
    policy = resolve_mutation_policy(
        "hotspot_biased_annealed_exploratory",
        unvalidated_hotspot,
        enabled=True,
        generation=2,
        total_generations=3,
    )

    assert policy["mode"] == "hotspot_biased_annealed"
    assert policy["requested_mode"] == (
        "hotspot_biased_annealed_exploratory"
    )
    assert policy["guided"] is True
    assert policy["hotspot_start"] == 4
    assert policy["hotspot_end"] == 6
    assert policy["hotspot_probability"] == pytest.approx(0.475)
    assert policy["hotspot_significance_validated"] is False
    assert policy["guidance_evidence_status"] == (
        "exploratory_unvalidated_hotspot"
    )


def test_exploratory_mode_still_requires_enabled_guidance(
    unvalidated_hotspot,
):
    policy = resolve_mutation_policy(
        "hotspot_biased_annealed_exploratory",
        unvalidated_hotspot,
        enabled=False,
        generation=2,
        total_generations=3,
    )

    assert policy["mode"] == "uniform"
    assert policy["guided"] is False