"""P4 target-aware candidate evaluation and optional Pareto selection.

Legacy CandidateRanker scores and historical fitness weights are untouched.
Only explicitly supplied, sequence-verified 9SYA receptor-peptide complexes
are eligible for target-aware evaluation.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
from typing import Any


TARGET_VERSION = "HERVK_9SYA_target_geometry_v2"


def _proximity(d203: float, d205: float) -> float:
    values = [
        math.exp(-max(float(d) - 4.0, 0.0) / 2.0)
        for d in (d203, d205)
    ]
    return sum(values) / 2.0


def register_target_complex(
    candidate: Any,
    complex_path: str | Path,
    *,
    target_structure_id: str = "9SYA",
    receptor_chain: str = "A",
    peptide_chain: str = "B",
    provenance: str = "HADDOCK3_predicted_complex",
) -> None:
    """Attach a verified complex input path to a Candidate's metadata."""
    candidate.metadata["target_complex_path"] = str(Path(complex_path).resolve())
    candidate.metadata["target_structure_id"] = target_structure_id
    candidate.metadata["target_receptor_chain"] = receptor_chain
    candidate.metadata["target_peptide_chain"] = peptide_chain
    candidate.metadata["target_complex_provenance"] = provenance


def _score_one_complex(candidate: Any, entry: dict[str, Any]) -> dict[str, Any]:
    # Lazy import keeps the original pipeline startup path unchanged.
    from src.target_aware_scoring import (
        HotspotSpec,
        StructuralScoreConfig,
        score_peptide_receptor_complex,
    )

    path = Path(entry["path"]).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)

    target_id = str(entry.get("target_structure_id", ""))
    if target_id != "9SYA":
        raise ValueError(f"Expected target 9SYA, received {target_id!r}")

    receptor_chain = str(entry.get("receptor_chain", "A"))
    peptide_chain = str(entry.get("peptide_chain", "B"))

    cfg = StructuralScoreConfig(
        hotspot_residues=(
            HotspotSpec(receptor_chain, 203, "ARG"),
            HotspotSpec(receptor_chain, 205, "ARG"),
        )
    )

    result = score_peptide_receptor_complex(
        path,
        config=cfg,
        receptor_chain=receptor_chain,
        peptide_chain=peptide_chain,
        target_structure_id=target_id,
        structure_provenance=str(entry.get("provenance", "external_complex")),
    )

    if result.get("peptide_sequence", "").upper() != candidate.sequence.upper():
        raise ValueError(
            "Candidate/complex sequence mismatch: "
            f"{candidate.sequence} != {result.get('peptide_sequence')}"
        )

    if result.get("score_valid") is not True:
        raise ValueError(
            "Target complex failed structural validation: "
            + ";".join(result.get("eligibility_reasons", []))
        )

    metrics = result.get("structural_metrics", {})
    hotspot_rows = metrics.get("hotspot_metrics", [])

    distances = {}
    for row in hotspot_rows:
        label = row.get("hotspot", "")
        if ":203:" in label:
            distances["arg203_min_A"] = float(row["minimum_distance_A"])
        elif ":205:" in label:
            distances["arg205_min_A"] = float(row["minimum_distance_A"])

    if set(distances) != {"arg203_min_A", "arg205_min_A"}:
        raise ValueError("Both Arg203 and Arg205 hotspot metrics are required")

    return {
        "path": str(path),
        "target_structure_id": target_id,
        "score_version": result.get("score_version"),
        "complex_sha256": result.get("structure_provenance", {}).get("sha256"),
        "target_geometry_score": float(
            result["target_aware_structural_score"]
        ),
        **distances,
        "minimum_interchain_distance_A": metrics.get(
            "minimum_interchain_heavy_atom_distance_A"
        ),
        "severe_clash_count": metrics.get(
            "severe_clash_heavy_atom_pairs"
        ),
        "interface_residue_pairs_5A": metrics.get(
            "distinct_interchain_residue_pairs", {}
        ).get("5A"),
        "sulfate_analysis": "not_included_in_v1_1_target_score",
    }


def _pareto_fronts(objectives: list[tuple[float, ...]]) -> list[int]:
    remaining = set(range(len(objectives)))
    ranks = [-1] * len(objectives)
    rank = 0

    while remaining:
        front = []
        for i in sorted(remaining):
            dominated = False
            for j in remaining:
                if i == j:
                    continue
                a, b = objectives[j], objectives[i]
                if all(x >= y for x, y in zip(a, b)) and any(
                    x > y for x, y in zip(a, b)
                ):
                    dominated = True
                    break
            if not dominated:
                front.append(i)

        if not front:
            raise RuntimeError("Pareto front construction failed")

        for i in front:
            ranks[i] = rank
            remaining.remove(i)
        rank += 1

    return ranks


def apply_p4_target_objective(
    candidates: list[Any],
    *,
    generation: int,
    output_directory: str | Path,
    pareto_enabled: bool | None = None,
) -> tuple[list[Any], dict[str, Any]]:
    """Evaluate explicitly registered complexes and optionally reorder.

    Pareto reordering is fail-closed: every candidate in the generation must
    have valid 9SYA evidence. Legacy overall_score is never overwritten.
    """
    if pareto_enabled is None:
        pareto_enabled = (
            os.environ.get("SILICON_VIROVORE_P4_TARGET_PARETO", "0") == "1"
        )

    records = []
    eligible = []

    for candidate in candidates:
        meta = candidate.metadata
        raw_entries = meta.get("target_complexes")

        if raw_entries is None:
            raw_path = meta.get("target_complex_path")
            raw_entries = [{"path": raw_path}] if raw_path else []

        if isinstance(raw_entries, dict):
            raw_entries = [raw_entries]

        entry_results = []
        errors = []

        for raw in raw_entries:
            try:
                if isinstance(raw, str):
                    raw = {"path": raw}

                entry = {
                    "path": raw["path"],
                    "target_structure_id": raw.get(
                        "target_structure_id",
                        meta.get("target_structure_id", ""),
                    ),
                    "receptor_chain": raw.get(
                        "receptor_chain",
                        meta.get("target_receptor_chain", "A"),
                    ),
                    "peptide_chain": raw.get(
                        "peptide_chain",
                        meta.get("target_peptide_chain", "B"),
                    ),
                    "provenance": raw.get(
                        "provenance",
                        meta.get(
                            "target_complex_provenance",
                            "external_complex",
                        ),
                    ),
                }
                entry_results.append(_score_one_complex(candidate, entry))
            except Exception as exc:
                errors.append(f"{type(exc).__name__}: {exc}")

        if entry_results and not errors:
            geometry = sum(
                x["target_geometry_score"] for x in entry_results
            ) / len(entry_results)
            d203 = sum(
                x["arg203_min_A"] for x in entry_results
            ) / len(entry_results)
            d205 = sum(
                x["arg205_min_A"] for x in entry_results
            ) / len(entry_results)
            proximity = _proximity(d203, d205)

            meta["p4_target_evidence"] = {
                "version": TARGET_VERSION,
                "status": "success",
                "replicate_count": len(entry_results),
                "target_geometry_score_mean": geometry,
                "arg203_min_A_mean": d203,
                "arg205_min_A_mean": d205,
                "hotspot_proximity_score": proximity,
                "complex_results": entry_results,
                "reliability_calibrated": False,
            }
            eligible.append(candidate)
        else:
            status = "missing_complex" if not raw_entries else "invalid_complex"
            meta["p4_target_evidence"] = {
                "version": TARGET_VERSION,
                "status": status,
                "errors": errors,
                "reliability_calibrated": False,
            }

        records.append({
            "sequence": candidate.sequence,
            "candidate_id": meta.get("candidate_id"),
            **meta["p4_target_evidence"],
        })

    coverage = len(eligible) / len(candidates) if candidates else 0.0
    reordered = False
    reason = "target evidence recorded; legacy ranking retained"

    if pareto_enabled:
        if candidates and len(eligible) == len(candidates):
            objectives = []
            for candidate in candidates:
                target = candidate.metadata["p4_target_evidence"]
                fitness = max(0.0, min(1.0, float(candidate.c_score)))
                objectives.append((
                    fitness,
                    float(target["target_geometry_score_mean"]),
                    float(target["hotspot_proximity_score"]),
                ))

            fronts = _pareto_fronts(objectives)
            original_order = {id(c): i for i, c in enumerate(candidates)}

            for candidate, front in zip(candidates, fronts):
                candidate.metadata["p4_target_pareto_front"] = front

            candidates = sorted(
                candidates,
                key=lambda c: (
                    c.metadata["p4_target_pareto_front"],
                    -float(c.overall_score),
                    original_order[id(c)],
                ),
            )
            for rank, candidate in enumerate(candidates, start=1):
                candidate.rank = rank

            reordered = True
            reason = "all candidates have valid target evidence; Pareto selection applied"
        else:
            reason = (
                "Pareto selection withheld: complete valid target evidence "
                "is required for every candidate"
            )

    report = {
        "generation": generation,
        "version": TARGET_VERSION,
        "candidate_count": len(candidates),
        "target_evidence_count": len(eligible),
        "target_evidence_coverage": coverage,
        "pareto_requested": bool(pareto_enabled),
        "pareto_reordered_candidates": reordered,
        "reason": reason,
        "legacy_fitness_scores_modified": False,
        "reliability_feedback_active": False,
        "candidates": records,
    }

    out = Path(output_directory) / "p4_target_evidence"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"generation_{generation:04d}.json"
    path.write_text(json.dumps(report, indent=2, default=str) + "\n")
    report["report_path"] = str(path)

    return candidates, report


def self_test() -> None:
    from tempfile import TemporaryDirectory
    from types import SimpleNamespace

    candidate = SimpleNamespace(
        sequence="MKLAVDALLRTFAESSDKLRR",
        c_score=0.8,
        overall_score=0.7,
        rank=1,
        metadata={},
    )

    with TemporaryDirectory() as tmp:
        ranked, report = apply_p4_target_objective(
            [candidate],
            generation=1,
            output_directory=tmp,
            pareto_enabled=True,
        )
        assert ranked[0] is candidate
        assert report["target_evidence_coverage"] == 0.0
        assert report["pareto_reordered_candidates"] is False
        assert candidate.metadata["p4_target_evidence"]["status"] == "missing_complex"

    print("P4 target integration self-test passed.")


if __name__ == "__main__":
    self_test()
