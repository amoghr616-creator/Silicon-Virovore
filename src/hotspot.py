"""Data-driven local sequence hotspot discovery for ARISE.

The model describes reproducible computational associations in candidate
histories. It is deliberately not a residue-level causal or structural
binding attribution model.
"""

from __future__ import annotations

from collections import Counter
import math
import random
from typing import Iterable

from src.ranking import CandidateRanker


ARISE_SCORE_SOURCE = "pooled_history_overall_score"


def annotate_arise_scores(history: list[dict]) -> dict:
    """Attach one shared-scope ARISE score to each history record.

    Raw objective components are normalized across the complete history
    available at discovery time.  This keeps candidates from different
    generations on the same scale while preserving their generation-local
    scores separately.
    """

    records = list(history)
    raw_records = [
        record.get("raw_objectives")
        for record in records
    ]
    has_raw_objectives = bool(records) and all(
        isinstance(raw, dict)
        for raw in raw_records
    )

    if has_raw_objectives:
        scores = CandidateRanker().score_raw_objective_records(raw_records)
        for record, score in zip(records, scores):
            record["arise_comparable_score"] = score
            record["arise_score_source"] = ARISE_SCORE_SOURCE
            record["arise_score_comparable_across_generations"] = True
        return {
            "records": records,
            "score_source": ARISE_SCORE_SOURCE,
            "score_comparable_across_generations": True,
        }

    # Explicitly allow external/synthetic callers to provide a score that
    # they have already established as comparable.  Unmarked scores are not
    # safe for pooled discovery and are rejected by discover_hotspot_model.
    explicitly_comparable = bool(records) and all(
        record.get("arise_score_comparable_across_generations") is True
        and record.get("score") is not None
        for record in records
    )
    if explicitly_comparable:
        for record in records:
            record["arise_comparable_score"] = record.get("score")
            record["arise_score_source"] = record.get(
                "arise_score_source",
                "provided_comparable_score",
            )
        return {
            "records": records,
            "score_source": records[0]["arise_score_source"],
            "score_comparable_across_generations": True,
        }

    return {
        "records": records,
        "score_source": "unmarked_history_score",
        "score_comparable_across_generations": False,
    }


def _valid_records(
    history: Iterable[dict],
    discovery_generation: int | None = None,
    generation_limit_run_id: str | None = None,
) -> list[dict]:
    """
    Validate discovery records.

    For the current run, only generations <= discovery_generation are used.

    Historical records from previous runs are allowed to contribute their
    complete history because their generation numbers are local to those runs.
    """

    records = []

    for record in history:
        sequence = str(record.get("sequence", "")).upper()

        score = record.get(
            "arise_comparable_score",
            record.get("score"),
        )

        generation = record.get("generation")

        if not sequence or score is None:
            continue

        try:
            score = float(score)
            generation = int(generation)
        except (TypeError, ValueError):
            continue

        if not math.isfinite(score):
            continue

        # Only constrain generations belonging to the current run.
        # Older runs have their own generation numbering.
        if (
            discovery_generation is not None
            and generation_limit_run_id is not None
            and str(record.get("run_id")) == str(generation_limit_run_id)
            and generation > discovery_generation
        ):
            continue

        records.append({
            **record,
            "sequence": sequence,
            "score": score,
        })

    return records


def _split_history(records: list[dict], seed: int) -> tuple[list[dict], list[dict]]:
    """Make a deterministic history-only discovery/holdout split."""

    if len(records) < 8:
        return list(records), []

    indices = list(range(len(records)))
    random.Random(seed).shuffle(indices)
    holdout_count = max(1, len(records) // 4)
    holdout_indices = set(indices[:holdout_count])
    discovery = [record for index, record in enumerate(records) if index not in holdout_indices]
    holdout = [record for index, record in enumerate(records) if index in holdout_indices]
    return discovery, holdout


def _high_background(
    records: list[dict],
    top_quantile: float,
) -> tuple[list[dict], list[dict]]:
    if len(records) < 2:
        return [], []
    ordered = sorted(
        enumerate(records),
        key=lambda item: (-item[1]["score"], item[0]),
    )
    high_count = max(1, int(math.ceil(len(records) * top_quantile)))
    high_count = min(high_count, len(records) - 1)
    high_indices = {index for index, _record in ordered[:high_count]}
    high = [record for index, record in enumerate(records) if index in high_indices]
    background = [record for index, record in enumerate(records) if index not in high_indices]
    return high, background


def _js_divergence(left: Counter[str], right: Counter[str]) -> float:
    """Jensen-Shannon divergence with light additive smoothing."""

    keys = set(left) | set(right)
    if not keys:
        return 0.0
    alpha = 0.5
    left_total = sum(left.values()) + alpha * len(keys)
    right_total = sum(right.values()) + alpha * len(keys)
    divergence = 0.0
    for key in keys:
        p = (left.get(key, 0) + alpha) / left_total
        q = (right.get(key, 0) + alpha) / right_total
        midpoint = (p + q) / 2.0
        divergence += 0.5 * p * math.log2(p / midpoint)
        divergence += 0.5 * q * math.log2(q / midpoint)
    return divergence


def _window_features(
    records: list[dict],
    start: int,
    length: int,
) -> tuple[list[Counter[str]], list[Counter[str]]]:
    residue_features = [Counter() for _ in range(length)]
    pair_features = [Counter() for _ in range(max(0, length - 1))]
    for record in records:
        sequence = record["sequence"][start:start + length]
        if len(sequence) != length:
            continue
        for offset, residue in enumerate(sequence):
            residue_features[offset][residue] += 1
        for offset in range(length - 1):
            pair_features[offset][sequence[offset:offset + 2]] += 1
    return residue_features, pair_features


def _window_association(
    records: list[dict],
    start: int,
    length: int,
    top_quantile: float,
) -> dict | None:
    high, background = _high_background(records, top_quantile)
    if not high or not background:
        return None

    high_residues, high_pairs = _window_features(high, start, length)
    background_residues, background_pairs = _window_features(
        background,
        start,
        length,
    )
    residue_js = (
        sum(
            _js_divergence(left, right)
            for left, right in zip(high_residues, background_residues)
        )
        / max(1, len(high_residues))
    )
    pair_js = (
        sum(
            _js_divergence(left, right)
            for left, right in zip(high_pairs, background_pairs)
        )
        / max(1, len(high_pairs))
        if high_pairs
        else 0.0
    )
    score = 0.5 * residue_js + 0.5 * pair_js
    return {
        "start_position": start,
        "end_position": start + length,
        "window_length": length,
        "hotspot_score": round(score, 8),
        "residue_distribution_association": round(residue_js, 8),
        "adjacent_pair_association": round(pair_js, 8),
        "high_performing_count": len(high),
        "background_count": len(background),
    }


def _bootstrap_stability(
    records: list[dict],
    starts: range,
    window_lengths: range,
    top_quantile: float,
    iterations: int,
    seed: int,
    hotspot_start: int | None = None,
    hotspot_end: int | None = None,
) -> dict[tuple[int, int], float]:
    if not records or iterations <= 0:
        return {}
    rng = random.Random(seed)
    selections: Counter[tuple[int, int]] = Counter()
    for _ in range(iterations):
        sample = [records[rng.randrange(len(records))] for _ in records]
        candidates = []
        for length in window_lengths:
            for start in starts:
                if hotspot_start is not None and start < hotspot_start:
                    continue
                if hotspot_end is not None and start + length > hotspot_end:
                    continue
                association = _window_association(
                    sample,
                    start,
                    length,
                    top_quantile,
                )
                if association is not None:
                    candidates.append(association)
        if candidates:
            best = max(candidates, key=lambda item: item["hotspot_score"])
            selections[(best["start_position"], best["window_length"])] += 1
    return {
        key: round(value / iterations, 6)
        for key, value in selections.items()
    }


def discover_hotspot_model(
    history: Iterable[dict],
    *,
    discovery_generation: int,
    run_id: str,
    top_quantile: float = 0.25,
    window_min: int = 2,
    window_max: int = 6,
    bootstrap_iterations: int = 100,
    random_seed: int = 616,
    hotspot_start: int | None = None,
    hotspot_end: int | None = None,
) -> dict:
    """Discover computationally associated local windows from past history.

    Window selection uses only the discovery split. Held-out association is
    reported for reproducibility and is never used to select the guided window.
    Positions are zero-based and ``end_position`` is exclusive.
    """

    annotation = annotate_arise_scores(
        history if isinstance(history, list) else list(history)
    )
    score_source = annotation["score_source"]
    score_comparable = annotation["score_comparable_across_generations"]
    records = _valid_records(
    annotation["records"],
    discovery_generation,
    generation_limit_run_id=run_id,
)
    source_run_ids = sorted({
        str(record.get("run_id"))
        for record in records
        if record.get("run_id") is not None
    })
    if run_id and run_id not in source_run_ids:
        source_run_ids.append(run_id)

    empty = {
        "model_version": "arise-hotspot-v2",
        "status": "insufficient_history",
        "discovery_generation": discovery_generation,
        "source_run_ids": source_run_ids,
        "discovery_sample_size": len(records),
        "discovery_sample_size_train": 0,
        "holdout_sample_size": 0,
        "hotspots": [],
        "selected_hotspot": None,
        "correlation_length": None,
        "method": "JS divergence over residue and adjacent-pair distributions",
        "causal_interpretation": False,
        "arise_score_source": score_source,
        "arise_score_comparable_across_generations": score_comparable,
    }
    if not score_comparable:
        empty["status"] = "incomparable_score"
        return empty
    if len(records) < 6:
        return empty

    sequence_length = len(records[0]["sequence"])
    window_min = max(2, int(window_min))
    window_max = min(sequence_length, int(window_max))
    if window_min > window_max:
        return empty
    if not 0.0 < top_quantile < 1.0:
        raise ValueError("top_quantile must be between 0 and 1.")

    discovery, holdout = _split_history(records, random_seed)
    min_start = max(0, hotspot_start) if hotspot_start is not None else 0
    max_start = sequence_length - window_min + 1
    if hotspot_end is not None:
        max_start = min(max_start, hotspot_end - window_min + 1)
    starts = range(min_start, max(min_start, max_start))
    lengths = range(window_min, window_max + 1)
    stability = _bootstrap_stability(
        discovery,
        starts,
        lengths,
        top_quantile,
        bootstrap_iterations,
        random_seed,
        hotspot_start,
        hotspot_end,
    )

    candidates = []
    for length in lengths:
        for start in starts:
            if hotspot_start is not None and start < hotspot_start:
                continue
            if hotspot_end is not None and start + length > hotspot_end:
                continue
            association = _window_association(
                discovery,
                start,
                length,
                top_quantile,
            )
            if association is None:
                continue
            heldout = _window_association(
                holdout,
                start,
                length,
                top_quantile,
            )
            association["held_out_association"] = (
                heldout["hotspot_score"] if heldout is not None else None
            )
            association["bootstrap_stability"] = stability.get(
                (start, length),
                0.0,
            )
            association["selection_score"] = round(
                association["hotspot_score"]
                * (0.5 + 0.5 * association["bootstrap_stability"]),
                8,
            )
            candidates.append(association)

    if not candidates:
        return empty

    selected = max(
        candidates,
        key=lambda item: (
            item["selection_score"],
            item["hotspot_score"],
            -item["window_length"],
            -item["start_position"],
        ),
    )
    ranked = sorted(
    candidates,
    key=lambda item: (
        item["selection_score"],
        item["hotspot_score"],
        -item["window_length"],
        -item["start_position"],
    ),
    reverse=True,
)

    for rank, fragment in enumerate(ranked, start=1):
        fragment["fragment_rank"] = rank

    selected = dict(selected)
    selected = dict(selected)
    selected.update({
        "discovery_generation": discovery_generation,
        "source_run_ids": source_run_ids,
        "discovery_sample_size": len(records),
    })
    return {
        "model_version": "arise-hotspot-v2",
        "status": "selected",
        "discovery_generation": discovery_generation,
        "source_run_ids": source_run_ids,
        "discovery_sample_size": len(records),
        "discovery_sample_size_train": len(discovery),
        "holdout_sample_size": len(holdout),
        "top_quantile": top_quantile,
        "window_min": window_min,
        "window_max": window_max,
        "bootstrap_iterations": bootstrap_iterations,
        "correlation_length": selected["window_length"],
        # Preserve the complete fragment landscape for analysis.
        "hotspots": ranked[:10],
        "fragment_scores": ranked,
        "selected_hotspot": selected,   
        "method": "JS divergence over residue and adjacent-pair distributions",
        "causal_interpretation": False,
        "selection_uses_heldout": False,
        "positions_zero_based_end_exclusive": True,
        "arise_score_source": score_source,
        "arise_score_comparable_across_generations": score_comparable,
    }
def derive_positional_signal(
    hotspot_model: dict | None,
    sequence_length: int,
) -> list[float]:
    """
    Convert the ARISE hotspot model into a normalized positional signal.

    This is an interpretability/reporting signal only.
    It is NOT used by ALE for mutation.
    """

    signal = [0.0] * sequence_length

    if not hotspot_model:
        return signal

    hotspots = hotspot_model.get("hotspots", [])

    for hotspot in hotspots:
        try:
            start = int(hotspot["start_position"])
            end = int(hotspot["end_position"])
            weight = float(
                hotspot.get(
                    "selection_score",
                    hotspot.get("hotspot_score", 0.0),
                )
            )
        except (KeyError, TypeError, ValueError):
            continue

        if weight <= 0.0 or end <= start:
            continue

        start = max(0, start)
        end = min(sequence_length, end)
        length = end - start

        if length <= 0:
            continue

        # Divide by window length so longer windows do not receive
        # a larger total contribution simply because they cover
        # more positions.
        per_position = weight / length

        for position in range(start, end):
            signal[position] += per_position

    maximum = max(signal, default=0.0)

    if maximum <= 0.0:
        return signal

    return [
        round(value / maximum, 4)
        for value in signal
    ]