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
