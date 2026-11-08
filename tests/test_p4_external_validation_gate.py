from src.p4_external_validation_gate import evaluate_external_validation_summary


def _summary(benefit, ci, p=0.01, improved=8, worsened=2, groups=10):
    return {
        "rows": 100,
        "outcome": "direct experimental binding free energy",
        "source_sha256": "b" * 64,
        "paired_primary_effect": {
            "macro_mae_benefit_kcal_mol": benefit,
            "group_bootstrap_95ci": ci,
            "exact_one_sided_signflip_p": p,
            "groups_improved": improved,
            "groups_worsened": worsened,
            "n_independent_groups": groups,
        },
    }


def test_gate_passes_only_when_all_predeclared_conditions_pass():
    result = evaluate_external_validation_summary(_summary(0.2, [0.05, 0.3]))
    assert result["allow_target_selection"] is True
    assert result["status"] == "external_transfer_gate_passed"


def test_gate_blocks_negative_or_uncertain_external_effect():
    result = evaluate_external_validation_summary(_summary(-0.05, [-0.12, -0.01], p=0.99, improved=3, worsened=7))
    assert result["allow_target_selection"] is False
    assert result["status"] == "external_transfer_gate_failed"


def test_gate_requires_ten_independent_groups():
    result = evaluate_external_validation_summary(_summary(0.2, [0.05, 0.3], groups=9))
    assert result["allow_target_selection"] is False
    assert result["checks"]["minimum_independent_groups"] is False
