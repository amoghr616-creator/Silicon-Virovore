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