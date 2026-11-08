"""Regression tests for measured zero-contact poses in P4 repeatability."""

import src.p4_pipeline_integration as p4


def _entries(*, nonempty_per_run=4, fingerprint_error=None):
    entries = []
    for run_index, (run_id, seed) in enumerate((("run_1917", 1917), ("run_2917", 2917))):
        for pose_index in range(5):
            if pose_index < nonempty_per_run:
                fingerprint = [f"R:{pose_index}:ARG|P:1:GLU"]
            else:
                fingerprint = []
            entries.append({
                "independent_run_id": run_id,
                "random_seed": seed,
                "protocol_sha256": "a" * 64,
                "execution_id": f"execution_{run_index}",
                "execution_manifest_sha256": ("b" if run_index == 0 else "c") * 64,
                "_interface_fingerprint_5A": fingerprint,
                "fingerprint_error": (
                    fingerprint_error if run_id == "run_1917" and pose_index == 0 else None
                ),
            })
    return entries


def test_empty_fingerprints_are_retained_as_measured_zero_contact_poses(monkeypatch):
    def fake_cross_run_reliability(runs):
        assert len(runs) == 2
        assert [len(run) for run in runs] == [5, 5]
        assert [sum(not fingerprint for fingerprint in run) for run in runs] == [1, 1]
        return {
            "reliability_proxy": 0.61,
            "mean_cross_run_jaccard_similarity": 0.72,
        }

    monkeypatch.setattr(
        "prototype4.p4_closed_loop.cross_run_reliability",
        fake_cross_run_reliability,
    )

    result = p4._reliability_diagnostic(_entries())

    assert result["status"] == "diagnostic_only"
    assert result["calibrated"] is False
    assert result["feedback_active"] is False
    assert result["pose_fingerprint_summary"]["run_1917"] == {
        "pose_count": 5,
        "nonempty_fingerprint_count": 4,
        "empty_fingerprint_count": 1,
        "nonempty_fraction": 0.8,
    }
    assert result["pose_fingerprint_summary"]["run_2917"]["empty_fingerprint_count"] == 1


def test_each_run_still_needs_at_least_two_nonempty_fingerprints():
    result = p4._reliability_diagnostic(_entries(nonempty_per_run=1))

    assert result["status"] == "insufficient_nonempty_interface_fingerprints"
    assert result["calibrated"] is False
    assert result["feedback_active"] is False


def test_fingerprint_calculation_errors_remain_ineligible():
    result = p4._reliability_diagnostic(
        _entries(fingerprint_error="ValueError: invalid coordinate data")
    )

    assert result["status"] == "fingerprint_calculation_failed"
    assert "invalid coordinate data" in result["reason"]
    assert result["calibrated"] is False
    assert result["feedback_active"] is False


def test_missing_fingerprint_field_remains_ineligible():
    entries = _entries()
    entries[0]["_interface_fingerprint_5A"] = None
    result = p4._reliability_diagnostic(entries)

    assert result["status"] == "missing_interface_fingerprint"
    assert result["calibrated"] is False
    assert result["feedback_active"] is False
