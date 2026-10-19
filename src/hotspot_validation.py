"""Independent-run validation for hotspot evidence.

Input record format:
{
    "run_id": "independent-run-01",
    "window_scores": {
        "window_A": 0.12,
        "window_B": 0.07,
    }
}

Important:
- run_id must identify genuinely independent runs, not individual descendants.
- Scores must be calculated by a prespecified, consistent method.
- This module does not prove biological or structural relevance.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Any
import random
import statistics


@dataclass
class HotspotValidationResult:
    selected_window: str
    observed_mean: float
    corrected_p_value: float
    independent_runs: int
    windows_tested: int
    passed: bool


def validate_hotspot(
    records: list[dict[str, Any]],
    *,
    min_independent_runs: int = 8,
    permutations: int = 5000,
    alpha: float = 0.05,
    seed: int = 2026,
) -> HotspotValidationResult:
    """Test whether the strongest window exceeds a run-level null.

    Uses one aggregate score per independent run and window, then a
    max-statistic sign-flip permutation test across all tested windows.

    This assumes independent run-level deviations are reasonably
    symmetric under the null. Treat the result as exploratory if that
    assumption is not justified for the score being tested.
    """
    if permutations < 1:
        raise ValueError("permutations must be >= 1")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1")

    # Aggregate duplicate records within each run so descendants do not
    # count as independent replicates.
    by_run: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )

    for record in records:
        run_id = str(record["run_id"])
        for window, score in record["window_scores"].items():
            value = float(score)
            if not (value == value and abs(value) != float("inf")):
                raise ValueError("Scores must be finite numbers")
            by_run[run_id][str(window)].append(value)

    run_ids = sorted(by_run)
    if len(run_ids) < min_independent_runs:
        raise ValueError(
            f"Need at least {min_independent_runs} independent runs; "
            f"received {len(run_ids)}"
        )

    # Only compare windows represented in every independent run.
    common_windows = set.intersection(
        *(set(by_run[run_id]) for run_id in run_ids)
    )
    windows = sorted(common_windows)
    if not windows:
        raise ValueError("No windows are present in every independent run")

    matrix = {
        window: [
            statistics.mean(by_run[run_id][window])
            for run_id in run_ids
        ]
        for window in windows
    }

    means = {
        window: statistics.mean(values)
        for window, values in matrix.items()
    }
    selected = max(means, key=means.get)
    observed = means[selected]

    # Center each window on zero under the null. For every permutation,
    # apply the same run-level sign flip across all windows. Taking the
    # maximum preserves the correction for searching multiple windows.
    rng = random.Random(seed)
    null_maxima = []

    for _ in range(permutations):
        signs = [
            -1.0 if rng.random() < 0.5 else 1.0
            for _ in run_ids
        ]
        permuted_means = [
            statistics.mean(
                sign * value
                for sign, value in zip(signs, matrix[window])
            )
            for window in windows
        ]
        null_maxima.append(max(permuted_means))

    p_value = (
        1 + sum(value >= observed for value in null_maxima)
    ) / (permutations + 1)

    return HotspotValidationResult(
        selected_window=selected,
        observed_mean=observed,
        corrected_p_value=p_value,
        independent_runs=len(run_ids),
        windows_tested=len(windows),
        passed=p_value < alpha,
    )