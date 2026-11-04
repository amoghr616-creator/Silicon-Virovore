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
