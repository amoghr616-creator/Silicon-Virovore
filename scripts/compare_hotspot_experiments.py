"""Compare baseline and hotspot-guided experiment audits."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics


def _load(path: Path) -> tuple[dict, list[dict]]:
    audit = json.loads(path.read_text())
    history_path = path.parent / audit["execution"]["candidate_history_file"]
    return audit, json.loads(history_path.read_text())


def _metrics(paths: list[Path], threshold: float) -> dict:
    runs = []
    for path in paths:
        audit, history = _load(path)
        reached = [item for item in history if item["score"] >= threshold]
        generations = [item["generation"] for item in reached]
        summaries = audit["optimization"]["generation_summaries"]
        guided = [
            item for item in history
            if item.get("hotspot_guided")
        ]
        inside = outside = 0
        for item in guided:
            start = item.get("hotspot_start")
            end = item.get("hotspot_end")
            if start is None or end is None:
                continue
            for position in item.get("mutation_positions", []):
                if start <= position < end:
                    inside += 1
                else:
                    outside += 1
        final_generation = max(
            (item["generation"] for item in history),
            default=None,
        )
        final_scores = [
            item["score"]
            for item in history
            if item["generation"] == final_generation
        ]
        best_scores = [summary.get("best_overall") for summary in summaries]
        runs.append({
            "audit": str(path),
            "candidates_evaluated_to_threshold": (
                min((index + 1 for index, item in enumerate(history)
                     if item["score"] >= threshold), default=None)
            ),
            "generations_to_threshold": min(generations, default=None),
            "final_generation_best_score": max(final_scores, default=None),
            "run_best_score": max(
                (item["score"] for item in history),
                default=None,
            ),
            "convergence_rate": (
                (best_scores[-1] - best_scores[0]) / max(1, len(best_scores) - 1)
                if len(best_scores) > 1 else None
            ),
            "mean_generation_diversity": statistics.mean(
                summary.get("mean_pairwise_identity", 0.0)
                for summary in summaries
            ) if summaries else None,
            "hotspot_stability": [
                model.get("selected_hotspot", {}).get("bootstrap_stability")
                for model in audit["optimization"].get("hotspot_models", [])
                if model.get("selected_hotspot")
            ],
            "guided_mutations_inside_hotspot": inside,
            "guided_mutations_outside_hotspot": outside,
        })
    final_scores = [run["run_best_score"] for run in runs if run["run_best_score"] is not None]
    return {
        "runs": runs,
        "run_to_run_consistency": {
            "mean_best_score": statistics.mean(final_scores) if final_scores else None,
            "std_best_score": statistics.pstdev(final_scores) if len(final_scores) > 1 else 0.0,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, nargs="+", required=True)
    parser.add_argument("--guided", type=Path, nargs="+", required=True)
    parser.add_argument("--threshold", type=float, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {
        "threshold": args.threshold,
        "baseline": _metrics(args.baseline, args.threshold),
        "hotspot_guided": _metrics(args.guided, args.threshold),
    }
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
