"""
bootstrap.py

Bootstrap and repeated-run stability analysis for Silicon Virovore.

This module quantifies how consistently candidate peptides
rank highly across repeated optimization runs.

Future upgrades:
    - Confidence intervals
    - Bayesian uncertainty estimation
    - Jackknife resampling
    - Cross-validation across datasets
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import random
import statistics


# ------------------------------------------------------------

@dataclass(slots=True)
class BootstrapResult:

    sequence: str

    wins: int

    frequency: float

    mean_score: float

    std_score: float


# ------------------------------------------------------------

class BootstrapAnalyzer:

    def __init__(self, iterations=100):

        self.iterations = iterations

    # --------------------------------------------------------

    def resample(self, candidates):

        """
        Sample with replacement.
        """

        return random.choices(

            candidates,

            k=len(candidates),

        )

    # --------------------------------------------------------

    def evaluate(self, candidates):

        """
        Returns the candidate with the highest overall score.
        """

        return max(

            candidates,

            key=lambda c: c.overall_score,

        )

    # --------------------------------------------------------

    def run(self, candidates):

        winners = []

        score_history = {}

        for _ in range(self.iterations):

            sample = self.resample(candidates)

            winner = self.evaluate(sample)

            winners.append(winner.sequence)

            score_history.setdefault(

                winner.sequence,

                []

            ).append(

                winner.overall_score

            )

        counts = Counter(winners)

        results = []

        for seq, wins in counts.items():

            scores = score_history[seq]

            results.append(

                BootstrapResult(

                    sequence=seq,

                    wins=wins,

                    frequency=round(

                        wins / self.iterations,

                        4,

                    ),

                    mean_score=round(

                        statistics.mean(scores),

                        4,

                    ),

                    std_score=round(

                        statistics.pstdev(scores),

                        4,

                    ),

                )

            )

        results.sort(

            key=lambda x: x.frequency,

            reverse=True,

        )

        return results

    # --------------------------------------------------------

    def summary(self, candidates):

        results = self.run(candidates)

        return {

            "iterations":

                self.iterations,

            "num_unique_winners":

                len(results),

            "results":

                results,

        }