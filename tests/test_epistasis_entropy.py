from types import SimpleNamespace

import pytest

from src.epistasis_entropy import EpistasisEntropyController


def test_shannon_entropy_per_position():
    entropy = EpistasisEntropyController.shannon_entropy(
        ["AA", "AC"]
    )

    assert entropy[0] == pytest.approx(0.0)
    assert entropy[1] == pytest.approx(1.0)


def test_pairwise_finite_difference_uses_c_score():
    controller = EpistasisEntropyController()

    # E_01 = 4 - 1 - 2 + 0 = 1.
    candidates = [
        SimpleNamespace(sequence="AAAA", c_score=0.0, overall_score=100.0),
        SimpleNamespace(sequence="CAAA", c_score=1.0, overall_score=0.0),
        SimpleNamespace(sequence="ACAA", c_score=2.0, overall_score=0.0),
        SimpleNamespace(sequence="CCAA", c_score=4.0, overall_score=0.0),
    ]

    diagnostics = controller.observe_generation(
        candidates,
        parent_sequence="AAAA",
    )

    assert diagnostics["epistasis_updates"] == 1
    assert controller.pair_count[0][1] == 1
    assert controller.pair_mean[0][1] == pytest.approx(1.0)


def test_factorial_population_contains_parent_and_unique_sequences():
    controller = EpistasisEntropyController()

    population = controller.generate_population(
        "MKLAVFALLVFFAGSSDLIRR",
        pop_size=16,
        seed=1234,
        mode="epistasis_entropy",
    )

    assert len(population) == 16
    assert population[0] == "MKLAVFALLVFFAGSSDLIRR"
    assert len(set(population)) == 16


def test_substitution_guided_v2_uses_smaller_pair_budget_and_guided_singles():
    controller = EpistasisEntropyController()
    parent = "MKLAVFALLVFFAGSSDLIRR"
    population = controller.generate_population(
        parent,
        pop_size=16,
        seed=1234,
        mode="hotspot_entropy_substitution_epistasis_exploratory_v2",
    )

    assert len(population) == 16
    assert population[0] == parent
    assert len(set(population)) == 16
    mutation_counts = [
        sum(a != b for a, b in zip(parent, sequence))
        for sequence in population[1:]
    ]
    assert mutation_counts.count(1) == 13
    assert mutation_counts.count(2) == 2
    assert all(count in {1, 2} for count in mutation_counts)


def test_matched_substitution_mode_preserves_p2_pair_budget():
    from src.config import PipelineSettings
    from src.epistasis_entropy import (
        EpistasisEntropyController,
        PAIRWISE_MODES,
    )

    matched_mode = (
        "hotspot_entropy_substitution_epistasis_exploratory_v2_matched"
    )
    p2_mode = "hotspot_epistasis_entropy_exploratory"

    assert matched_mode in PAIRWISE_MODES
    PipelineSettings(mutation_guidance_mode=matched_mode)

    parent = "MKLAVFALLVFFAGSSDLIRR"
    hotspot_positions = {2, 7}

    p2 = EpistasisEntropyController()
    p3 = EpistasisEntropyController()

    p2_population = p2.generate_population(
        parent,
        pop_size=16,
        seed=1234,
        mode=p2_mode,
        hotspot_positions=hotspot_positions,
        hotspot_boost=2.0,
    )

    p3_population = p3.generate_population(
        parent,
        pop_size=16,
        seed=1234,
        mode=matched_mode,
        hotspot_positions=hotspot_positions,
        hotspot_boost=2.0,
    )

    assert len(p2_population) == len(p3_population) == 16
    assert len(set(p2_population)) == len(set(p3_population)) == 16
    assert p2_population[0] == p3_population[0] == parent

    # Both modes must allocate one factorial pair at population size 16.
    assert len(p2.last_selected_pairs) == 1
    assert len(p3.last_selected_pairs) == 1

    # P2 uses random residue proposals; matched P3 uses learned proposals.
    assert p2.last_selected_substitutions == []
    assert len(p3.last_selected_substitutions) >= 1

    # One pair produces two single mutants and one double mutant.
    # Remaining offspring are single mutants in both conditions.
    for population in (p2_population, p3_population):
        counts = [
            sum(a != b for a, b in zip(parent, sequence))
            for sequence in population[1:]
        ]
        assert counts.count(1) == 14
        assert counts.count(2) == 1
