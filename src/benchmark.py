"""
benchmark.py

Scientific benchmarking framework for Silicon Virovore.

Runs known peptide datasets through the pipeline and evaluates
ranking performance, reproducibility, and predictive quality.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import statistics

from src.models import Candidate
from src.ranking import CandidateRanker


# ==========================================================
# Benchmark Sample
# ==========================================================

@dataclass(slots=True)
class BenchmarkSample:

    sequence: str

    expected_activity: float

    label: str = ""


# ==========================================================
# Benchmark Result
# ==========================================================

@dataclass(slots=True)
class BenchmarkResult:

    sequence: str

    expected_activity: float

    predicted_score: float

    predicted_rank: int


# ==========================================================
# Benchmark Runner
# ==========================================================

class BenchmarkRunner:

    def __init__(self):

        self.ranker = CandidateRanker()

    # -----------------------------------------------------

    def evaluate(
        self,
        candidates: list[Candidate],
        benchmark: list[BenchmarkSample],
    ) -> list[BenchmarkResult]:

        ranked = self.ranker.rank(candidates)

        lookup = {
            c.sequence: c
            for c in ranked
        }

        results = []

        for sample in benchmark:

            if sample.sequence not in lookup:
                continue

            candidate = lookup[sample.sequence]

            results.append(

                BenchmarkResult(

                    sequence=sample.sequence,

                    expected_activity=sample.expected_activity,

                    predicted_score=candidate.overall_score,

                    predicted_rank=candidate.rank,
                )
            )

        return results

    # -----------------------------------------------------

    def mean_absolute_error(
        self,
        results: list[BenchmarkResult],
    ):

        if not results:
            return 0.0

        errors = [

            abs(
                r.expected_activity -
                r.predicted_score
            )

            for r in results
        ]

        return statistics.mean(errors)

    # -----------------------------------------------------

    def top1_accuracy(
        self,
        results: list[BenchmarkResult],
    ):

        if not results:
            return 0.0

        expected_best = max(
            results,
            key=lambda r: r.expected_activity,
        )

        predicted_best = min(
            results,
            key=lambda r: r.predicted_rank,
        )

        return float(
            expected_best.sequence ==
            predicted_best.sequence
        )

    # -----------------------------------------------------

    def top5_accuracy(
        self,
        results: list[BenchmarkResult],
    ):

        expected = sorted(

            results,

            key=lambda r: r.expected_activity,

            reverse=True,

        )[:5]

        predicted = sorted(

            results,

            key=lambda r: r.predicted_rank,

        )[:5]

        expected = {

            x.sequence

            for x in expected
        }

        predicted = {

            x.sequence

            for x in predicted
        }

        return len(
            expected &
            predicted
        ) / len(expected)

    # -----------------------------------------------------

    def save_csv(
        self,
        results,
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

                "Sequence",

                "Expected Activity",

                "Predicted Score",

                "Predicted Rank",

            ])

            for r in results:

                writer.writerow([

                    r.sequence,

                    r.expected_activity,

                    r.predicted_score,

                    r.predicted_rank,

                ])

    # -----------------------------------------------------

    def summarize(
        self,
        results,
    ):

        return {

            "samples":

                len(results),

            "mae":

                self.mean_absolute_error(
                    results
                ),

            "top1":

                self.top1_accuracy(
                    results
                ),

            "top5":

                self.top5_accuracy(
                    results
                ),

        }