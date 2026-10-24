from random import Random

from src.epistasis_entropy import EpistasisEntropyController


PARENT = "MKLAVFALLVFFAGSSDLIRR"


def test_hotspot_increases_position_weight_without_zeroing_other_positions():
    controller = EpistasisEntropyController()

    uniform = controller.position_weights(
        4,
        use_entropy=False,
    )
    guided = controller.position_weights(
        4,
        use_entropy=False,
        hotspot_positions={1, 2},
        hotspot_boost=2.0,
    )

    assert abs(sum(guided) - 1.0) < 1e-12
    assert guided[1] > uniform[1]
    assert guided[2] > uniform[2]
    assert guided[0] > 0
    assert guided[3] > 0


def test_triple_guidance_prioritizes_validated_hotspot_pair():
    controller = EpistasisEntropyController()

    population = controller.generate_population(
        PARENT,
        pop_size=16,
        seed=12345,
        mode="hotspot_epistasis_entropy",
        hotspot_positions={4, 5},
        hotspot_boost=2.0,
    )

    assert len(population) == 16
    assert len(set(population)) == 16
    assert population[0] == PARENT
    assert controller.last_pair_weights
    assert controller.last_pair_weights[0][0] == (4, 5)
    assert controller.last_selected_pairs



def test_exploratory_triple_guidance_weights_selected_hotspot_pair():
    controller = EpistasisEntropyController()
    population = controller.generate_population(
        PARENT,
        pop_size=16,
        seed=12345,
        mode="hotspot_epistasis_entropy_exploratory",
        hotspot_positions={4, 5},
        hotspot_boost=2.0,
    )

    assert len(population) == 16
    assert len(set(population)) == 16
    assert controller.last_pair_weights
    assert controller.last_pair_weights[0][0] == (4, 5)


def test_substitution_specific_pair_score_learns_amino_acid_choices():
    controller = EpistasisEntropyController()
    parent = "MKLAVFALLVFFAGSSDLIRR"
    i, j = 4, 5
    aa_i = next(aa for aa in "ACDEFGHIKLMNPQRSTVWY" if aa != parent[i])
    aa_j = next(aa for aa in "ACDEFGHIKLMNPQRSTVWY" if aa != parent[j])

    key_i = (i, parent[i], aa_i)
    key_j = (j, parent[j], aa_j)
    pair_key = (i, parent[i], aa_i, j, parent[j], aa_j)

    controller.substitution_mean[key_i] = 0.03
    controller.substitution_count[key_i] = 100
    controller.substitution_mean[key_j] = 0.03
    controller.substitution_count[key_j] = 100
    controller.substitution_pair_mean[pair_key] = 0.04
    controller.substitution_pair_count[pair_key] = 100

    estimates = dict(controller._substitution_pair_scores(parent, i, j))
    assert estimates[(aa_i, aa_j)] == max(estimates.values())


def test_new_mode_generates_unique_population_and_records_guided_choices():
    controller = EpistasisEntropyController()
    parent = "MKLAVFALLVFFAGSSDLIRR"
    observed_epistasis_flags = []
    original_pair_weights = controller._pair_weights

    def spy_pair_weights(n, positions, use_epistasis):
        observed_epistasis_flags.append(use_epistasis)
        return original_pair_weights(n, positions, use_epistasis)

    controller._pair_weights = spy_pair_weights
    population = controller.generate_population(
        parent,
        pop_size=16,
        seed=12345,
        mode="hotspot_entropy_substitution_epistasis_exploratory",
        hotspot_positions={4, 5},
        hotspot_boost=1.25,
    )

    assert len(population) == 16
    assert len(set(population)) == 16
    assert population[0] == parent
    assert observed_epistasis_flags == [True]
    assert len(controller.last_selected_pairs) == 1
    assert len(controller.last_selected_substitutions) == 1


def test_learned_substitution_gain_biases_site_and_single_mutant_choice():
    controller = EpistasisEntropyController()
    parent = PARENT
    controller.substitution_random_probability = 0.0
    key = (0, parent[0], "A")
    controller.substitution_mean[key] = 0.50
    controller.substitution_count[key] = 100

    assert controller._choose_single_substitution(
        parent, 0, Random(9)
    ) == "A"

    captured = {}
    original_pair_weights = controller._pair_weights

    def spy_pair_weights(n, positions, use_epistasis):
        captured["positions"] = list(positions)
        captured["use_epistasis"] = use_epistasis
        return original_pair_weights(n, positions, use_epistasis)

    controller._pair_weights = spy_pair_weights
    population = controller.generate_population(
        parent,
        pop_size=16,
        seed=991,
        mode="hotspot_entropy_substitution_epistasis_exploratory",
    )
    assert len(population) == 16
    assert captured["use_epistasis"] is True
    assert captured["positions"][0] > max(captured["positions"][1:])
