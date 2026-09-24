"""Persistent ARISE memory.

ARISE memory stores prior computational observations and learned hotspot
models across runs.

Memory can be recorded automatically while remaining excluded from new
discovery by default. This prevents accidental cross-experiment leakage.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


MEMORY_VERSION = 1


def _empty_memory() -> dict[str, Any]:
    return {
        "memory_version": MEMORY_VERSION,
        "candidate_history": [],
        "runs": [],
    }


def load_arise_memory(
    path: Path,
    *,
    max_records: int = 5000,
    max_runs: int = 50,
) -> dict[str, Any]:
    if not path.exists():
        return _empty_memory()

    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return _empty_memory()

    if not isinstance(payload, dict):
        return _empty_memory()

    payload.setdefault("memory_version", MEMORY_VERSION)
    payload.setdefault("candidate_history", [])
    payload.setdefault("runs", [])

    payload["candidate_history"] = list(
        payload["candidate_history"]
    )[-max_records:]

    payload["runs"] = list(
        payload["runs"]
    )[-max_runs:]

    return payload


def save_arise_memory(
    path: Path,
    memory: dict[str, Any],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    temporary_path = path.with_suffix(
        path.suffix + ".tmp"
    )

    temporary_path.write_text(
        json.dumps(
            memory,
            indent=2,
            default=str,
        ) + "\n"
    )

    temporary_path.replace(path)


def get_memory_history(
    memory: dict[str, Any],
) -> list[dict]:
    return list(
        memory.get("candidate_history", [])
    )


def upsert_arise_run(
    memory: dict[str, Any],
    *,
    run_id: str,
    condition: str,
    random_seed: int | None,
    candidate_history: list[dict],
    hotspot_models: list[dict],
    derived_positional_signal: list[float],
    max_records: int = 5000,
    max_runs: int = 50,
) -> dict[str, Any]:
    """Insert or replace one run in persistent memory."""

    run_id = str(run_id)

    # Remove previous copy of this run.
    historical_records = [
        record
        for record in memory.get("candidate_history", [])
        if str(record.get("run_id")) != run_id
    ]

    new_records = [
        dict(record)
        for record in candidate_history
    ]

    combined_records = (
        historical_records + new_records
    )[-max_records:]

    historical_runs = [
        run
        for run in memory.get("runs", [])
        if str(run.get("run_id")) != run_id
    ]

    run_record = {
        "run_id": run_id,
        "condition": condition,
        "random_seed": random_seed,
        "candidate_count": len(candidate_history),
        "hotspot_model_count": len(hotspot_models),
        "hotspot_models": hotspot_models,
        "derived_positional_signal": derived_positional_signal,
    }

    memory["memory_version"] = MEMORY_VERSION
    memory["candidate_history"] = combined_records
    memory["runs"] = (
        historical_runs + [run_record]
    )[-max_runs:]

    return memory


def aggregate_memory_positional_signal(
    memory: dict[str, Any],
    sequence_length: int,
) -> dict[str, Any]:
    """Average derived positional signals across remembered runs."""

    signals = []

    for run in memory.get("runs", []):
        signal = run.get(
            "derived_positional_signal"
        )

        if not isinstance(signal, list):
            continue

        if len(signal) != sequence_length:
            continue

        try:
            signals.append([
                float(value)
                for value in signal
            ])
        except (TypeError, ValueError):
            continue

    if not signals:
        return {
            "mean_signal": [0.0] * sequence_length,
            "run_count": 0,
            "support_count": [0] * sequence_length,
        }

    mean_signal = []

    for position in range(sequence_length):
        values = [
            signal[position]
            for signal in signals
        ]

        mean_signal.append(
            round(
                sum(values) / len(values),
                4,
            )
        )

    maximum = max(mean_signal, default=0.0)

    if maximum > 0:
        mean_signal = [
            round(value / maximum, 4)
            for value in mean_signal
        ]

    return {
        "mean_signal": mean_signal,
        "run_count": len(signals),
        "support_count": [
            sum(
                signal[position] > 0
                for signal in signals
            )
            for position in range(sequence_length)
        ],
    }