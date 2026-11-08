import pytest

from src.config import PipelineSettings


def test_p4_target_optimization_is_explicit_and_default_off():
    assert PipelineSettings().p4_target_optimization_enabled is False
    assert PipelineSettings(p4_target_optimization_enabled=True).p4_target_optimization_enabled is True


def test_p4_target_optimization_cannot_be_combined_with_vina_tier2():
    with pytest.raises(ValueError, match="disable standalone Vina tier-2"):
        PipelineSettings(
            p4_target_optimization_enabled=True,
            enable_vina_tier2=True,
        )


def test_exploratory_mode_requires_target_optimization():
    with pytest.raises(ValueError, match="requires P4 target optimization"):
        PipelineSettings(p4_exploratory_target_optimization_enabled=True)


def test_exploratory_mode_is_explicit_and_disabled_by_default():
    assert PipelineSettings().p4_exploratory_target_optimization_enabled is False
    settings = PipelineSettings(
        p4_target_optimization_enabled=True,
        p4_exploratory_target_optimization_enabled=True,
    )
    assert settings.p4_exploratory_target_optimization_enabled is True
