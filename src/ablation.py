"""
ablation.py

Scientific ablation study for Silicon Virovore.

Evaluates the contribution of each computational module by
systematically removing one component at a time and measuring
the impact on benchmark performance.
"""

from __future__ import annotations

from dataclasses import dataclass
from copy import deepcopy
from pathlib import Path
import csv

from src.models import Candidate
from src.benchmark import BenchmarkRunner, BenchmarkSample
from src.ranking import CandidateRanker


# ============================================================
# Ablation Result
# ============================================================

@dataclass(slots=True)
class AblationResult:

    component: str

    top1: float

    top5: float

    mae: float

    score_drop: float


# ============================================================
# Ablation Engine
# ============================================================

class AblationStudy:

    def __init__(self):

        self.benchmark = BenchmarkRunner()

    # --------------------------------------------------------

    def _disable_component(
        self,
        candidates: list[Candidate],
        component: str,
    ):

        modified = deepcopy(candidates)

        for c in modified:

            if component == "fitness":
                c.c_score = 0

            elif component == "docking":
                c.strongest_anchor_delta_g = 0

            elif component == "structure":
                c.structure_confidence = 0

            elif component == "helix":
                c.helix_propensity = 0

            elif component == "hydrophobicity":
                c.hydrophobic_moment = 0

            elif component == "solvation":
                c.solvation_energy = 0

            elif component == "md":
                c.md_stability_score = 0

            elif component == "safety":
                c.toxicity_score = 0

        return modified

    # --------------------------------------------------------

    def run(
        self,
        candidates: list[Candidate],
        benchmark: list[BenchmarkSample],
    ) -> list[AblationResult]:

        ranker = CandidateRanker()

        baseline_candidates = ranker.rank(
            deepcopy(candidates)
        )

        baseline_results = self.benchmark.evaluate(
            baseline_candidates,
            benchmark,
        )

        baseline = self.benchmark.summarize(
            baseline_results
        )

        baseline_top1 = baseline["top1"]

        components = [

            "fitness",

            "docking",

            "structure",

            "helix",

            "hydrophobicity",

            "solvation",

            "md",

            "safety",
        ]

        results = []

        for component in components:

            modified = self._disable_component(
                candidates,
                component,
            )

            ranked = ranker.rank(modified)

            benchmark_results = self.benchmark.evaluate(
                ranked,
                benchmark,
            )

            summary = self.benchmark.summarize(
                benchmark_results
            )

            results.append(

                AblationResult(

                    component=component,

                    top1=summary["top1"],

                    top5=summary["top5"],

                    mae=summary["mae"],

                    score_drop=baseline_top1 -
                    summary["top1"],
                )
            )

        return results

    # --------------------------------------------------------

    def save_csv(
        self,
        results: list[AblationResult],
        output_file,
    ):

        output_file = Path(output_file)

        with open(
            output_file,
            "w",
            newline="",
        ) as f:

            writer = csv.writer(f)

            writer.writerow([

                "Component",

                "Top1",

                "Top5",

                "MAE",

                "Top1 Drop",

            ])

            for r in results:

                writer.writerow([

                    r.component,

                    round(r.top1,4),

                    round(r.top5,4),

                    round(r.mae,4),

                    round(r.score_drop,4),

                ])

    # --------------------------------------------------------

    def print_summary(
        self,
        results,
    ):

        print()

        print("="*60)

        print("ABLATION STUDY")

        print("="*60)

        print()

        for r in sorted(
            results,
            key=lambda x: x.score_drop,
            reverse=True,
        ):

            print(

                f"{r.component:18}"

                f"Top1={r.top1:.3f}   "

                f"Drop={r.score_drop:.3f}"

            )

        print()

        print("="*60)