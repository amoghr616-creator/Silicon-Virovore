"""Prepare matched, nonredundant candidate panels without docking-based selection."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
RUN_PAIRS = (
    ("baseline-616-A", "adaptive-616"),
    ("baseline-618-A", "adaptive-618-A"),
)
PANEL_SIZE = 32
MATCHED_CONFIGURATION_KEYS = (
    "random_seed",
    "population_size",
    "generations",
    "mutation_rate",
)


def select_nonredundant_panel(
    records: list[dict],
    *,
    run_id: str,
    arm: str,
    panel_size: int = PANEL_SIZE,
) -> list[dict]:
    """Select unique sequences by original native fitness, never docking score."""
    if panel_size < 1:
        raise ValueError("panel_size must be >= 1")
    ordered = sorted(
        records,
        key=lambda record: (
            -float(record["native_fitness"]),
            str(record["candidate_id"]),
        ),
    )
    occurrences: dict[str, int] = {}
    for record in records:
        sequence = str(record["sequence"])
        occurrences[sequence] = occurrences.get(sequence, 0) + 1

    selected: list[dict] = []
    seen: set[str] = set()
    for record in ordered:
        sequence = str(record["sequence"])
        if sequence in seen:
            continue
        seen.add(sequence)
        selected.append({
            "arm": arm,
            "source_run": run_id,
            "candidate_id": record["candidate_id"],
            "generation": int(record["generation"]),
            "sequence": sequence,
            "original_fitness_score": float(record["native_fitness"]),
            "original_composite_score": float(record["score"]),
            "selection_score_source": "candidate_history.native_fitness",
            "sequence_occurrences_in_source_run": occurrences[sequence],
            "unique_in_source_run": occurrences[sequence] == 1,
            "unique_in_selected_panel": True,
        })
        if len(selected) == panel_size:
            break
    if len(selected) != panel_size:
        raise ValueError(
            f"{run_id} has only {len(selected)} unique selectable sequences; "
            f"need {panel_size}"
        )
    return selected


def _load_run(run_id: str) -> tuple[dict, list[dict]]:
    run_dir = ROOT / "results" / run_id
    audit_path = run_dir / "experiment_audit.json"
    history_path = run_dir / "candidate_history.json"
    audit = json.loads(audit_path.read_text())
    history = json.loads(history_path.read_text())
    return audit["configuration"], history


def build_panels(panel_size: int = PANEL_SIZE) -> tuple[list[dict], dict]:
    panels: list[dict] = []
    summaries: dict = {"selection_rule": "top native_fitness, descending; candidate_id ascending for ties; deduplicate exact sequence before taking panel", "docking_used_for_selection": False, "panel_size_per_run_arm": panel_size, "matched_pairs": []}

    for baseline_id, adaptive_id in RUN_PAIRS:
        baseline_config, baseline_records = _load_run(baseline_id)
        adaptive_config, adaptive_records = _load_run(adaptive_id)
        mismatches = {
            key: (baseline_config.get(key), adaptive_config.get(key))
            for key in MATCHED_CONFIGURATION_KEYS
            if baseline_config.get(key) != adaptive_config.get(key)
        }
        if mismatches:
            summaries["matched_pairs"].append({
                "baseline_run": baseline_id,
                "adaptive_run": adaptive_id,
                "matched": False,
                "configuration_mismatches": mismatches,
            })
            continue

        baseline_panel = select_nonredundant_panel(
            baseline_records,
            run_id=baseline_id,
            arm="baseline",
            panel_size=panel_size,
        )
        adaptive_panel = select_nonredundant_panel(
            adaptive_records,
            run_id=adaptive_id,
            arm="adaptive",
            panel_size=panel_size,
        )
        panels.extend(baseline_panel)
        panels.extend(adaptive_panel)

        baseline_fitness = [row["original_fitness_score"] for row in baseline_panel]
        adaptive_fitness = [row["original_fitness_score"] for row in adaptive_panel]
        summaries["matched_pairs"].append({
            "baseline_run": baseline_id,
            "adaptive_run": adaptive_id,
            "matched": True,
            "random_seed": baseline_config["random_seed"],
            "population_size": baseline_config["population_size"],
            "generations": baseline_config["generations"],
            "mutation_rate": baseline_config["mutation_rate"],
            "source_candidate_count": {
                "baseline": len(baseline_records),
                "adaptive": len(adaptive_records),
            },
            "source_unique_sequence_count": {
                "baseline": len({r["sequence"] for r in baseline_records}),
                "adaptive": len({r["sequence"] for r in adaptive_records}),
            },
            "panel_unique_sequence_count": {
                "baseline": len({r["sequence"] for r in baseline_panel}),
                "adaptive": len({r["sequence"] for r in adaptive_panel}),
            },
            "fitness_median": {
                "baseline": statistics.median(baseline_fitness),
                "adaptive": statistics.median(adaptive_fitness),
                "adaptive_minus_baseline": (
                    statistics.median(adaptive_fitness)
                    - statistics.median(baseline_fitness)
                ),
            },
            "fitness_mean": {
                "baseline": statistics.mean(baseline_fitness),
                "adaptive": statistics.mean(adaptive_fitness),
                "adaptive_minus_baseline": (
                    statistics.mean(adaptive_fitness)
                    - statistics.mean(baseline_fitness)
                ),
            },
            "independent_run_pair_count": 1,
            "inferential_uncertainty_status": (
                "descriptive_only_one_matched_run_pair"
            ),
        })
    return panels, summaries


def write_outputs(panels: list[dict], summaries: dict, output_dir: Path) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    fields = list(panels[0]) if panels else []
    with (output_dir / "selected_candidates.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(panels)
    (output_dir / "panel_summary.json").write_text(
        json.dumps(summaries, indent=2) + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "prototype2" / "validation_checkpoint_20261016",
    )
    parser.add_argument("--panel-size", type=int, default=PANEL_SIZE)
    args = parser.parse_args()
    panels, summaries = build_panels(args.panel_size)
    write_outputs(panels, summaries, args.output_dir)
    print(f"Selected candidates: {len(panels)}")
    print(f"Wrote {args.output_dir / 'selected_candidates.csv'}")
    print(f"Wrote {args.output_dir / 'panel_summary.json'}")
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
