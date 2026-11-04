import numpy as np
import pandas as pd
import pytest

from prototype4.selective_abstention import analyze_selective_abstention


def make_oof_data(n_groups=10, rows_per_group=12):
    rows = []
    for g in range(n_groups):
        # The uncertainty score orders low-error rows first within each family.
        for i in range(rows_per_group):
            error = 0.05 + (i / max(rows_per_group - 1, 1)) * 1.5
            rows.append({
                "feature_uid": f"feature-{g}-{i}",
                "heldout_group": f"family-{g}",
                "split_unit": "family",
                "y_true": float(error),
                "y_pred": 0.0,
                "calibrated_abs_error": float(error),
            })
    return pd.DataFrame(rows)


def test_selective_retention_beats_family_matched_random_on_signal_data():
    table, summary = analyze_selective_abstention(
        make_oof_data(),
        random_repetitions=200,
        family_bootstrap_repetitions=1000,
        seed=7,
    )
    row = table.loc[table["target_coverage"] == 0.5].iloc[0]
    assert summary["status"] == "family_matched_selective_analysis_complete"
    assert row["random_minus_selective_mae"] > 0
    assert row["random_minus_selective_family_bootstrap_ci95_low"] > 0
    assert row["realized_coverage"] == pytest.approx(0.5)


def test_requires_family_held_out_predictions():
    frame = make_oof_data()
    frame["split_unit"] = "complex"
    with pytest.raises(ValueError, match="family-held-out"):
        analyze_selective_abstention(
            frame, random_repetitions=100, family_bootstrap_repetitions=1000
        )


def test_rejects_duplicate_oof_prediction_ids():
    frame = make_oof_data()
    frame.loc[1, "feature_uid"] = frame.loc[0, "feature_uid"]
    with pytest.raises(ValueError, match="exactly one"):
        analyze_selective_abstention(
            frame, random_repetitions=100, family_bootstrap_repetitions=1000
        )
