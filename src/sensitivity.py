"""
sensitivity.py

Performs parameter sensitivity analysis for Silicon Virovore.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
from statistics import mean

from src.models import Candidate
from src.population_runner import generate_candidates
from src.ranking import CandidateRanker


# ==========================================================
# Sensitivity Result
# ==========================================================

@dataclass(slots=True)
class SensitivityResult:

    parameter: str

    value: float

    mean_score: float

    best_score: float

    std_score: float

    top_sequence: str


# ==========================================================
# Sensitivity Engine
# ==========================================================

class SensitivityStudy:

    def __init__(self):

        self.ranker = CandidateRanker()

    # ------------------------------------------------------

    def evaluate_population(
        self,
        candidates: list[Candidate],
    ):

        ranked = self.ranker.rank(candidates)

        scores = [
            c.overall_score
            for c in ranked
        ]

        return SensitivityResult(

            parameter="",

            value=0,

            mean_score=mean(scores),

            best_score=max(scores),

            std_score=(
                0
                if len(scores) < 2
                else (
                    sum(
                        (x - mean(scores)) ** 2
                        for x in scores
                    )
                    / (len(scores) - 1)
                ) ** 0.5
            ),

            top_sequence=ranked[0].sequence,
        )

    # ------------------------------------------------------

    def mutation_rate_study(

        self,

        seed_sequence,

        mutation_rates,

        population_size=100,

    ):

        results = []

        for rate in mutation_rates:

            candidates = generate_candidates(

                seed_sequence,

                population_size=population_size,

                mutation_rate=rate,

            )

            result = self.evaluate_population(
                candidates
            )

            result.parameter = "mutation_rate"

            result.value = rate

            results.append(result)

        return results

    # ------------------------------------------------------

    def population_size_study(

        self,

        seed_sequence,

        population_sizes,

        mutation_rate=0.05,

    ):

        results = []

        for size in population_sizes:

            candidates = generate_candidates(

                seed_sequence,

                population_size=size,

                mutation_rate=mutation_rate,

            )

            result = self.evaluate_population(
                candidates
            )

            result.parameter = "population_size"

            result.value = size

            results.append(result)

        return results

    # ------------------------------------------------------

    def save_csv(

        self,

        results,

        filename,

    ):

        filename = Path(filename)

        with open(
            filename,
            "w",
            newline="",
        ) as f:

            writer = csv.writer(f)

            writer.writerow([

                "Parameter",

                "Value",

                "Mean Score",

                "Best Score",

                "Std",

                "Top Sequence",

            ])

            for r in results:

                writer.writerow([

                    r.parameter,

                    r.value,

                    round(r.mean_score,4),

                    round(r.best_score,4),

                    round(r.std_score,4),

                    r.top_sequence,

                ])