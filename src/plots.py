"""
plots.py

Publication-quality plotting utilities for Silicon Virovore.
"""

from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt

from src.models import Candidate


class PlotGenerator:

    def __init__(self, output_dir: Path):

        self.output_dir = Path(output_dir)

        self.output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

    # -----------------------------------------------------

    def fitness_histogram(
        self,
        candidates: list[Candidate],
    ):

        values = [
            c.c_score
            for c in candidates
        ]

        plt.figure(figsize=(8,5))

        plt.hist(
            values,
            bins=20,
        )

        plt.xlabel("Native Fitness")

        plt.ylabel("Count")

        plt.title("Distribution of Native Fitness Scores")

        plt.tight_layout()

        output = self.output_dir / "fitness_histogram.png"

        plt.savefig(output, dpi=300)

        plt.close()

        return output

    # -----------------------------------------------------

    def docking_histogram(
        self,
        candidates: list[Candidate],
    ):

        values = [
            c.strongest_anchor_delta_g
            for c in candidates
            if c.strongest_anchor_delta_g is not None
        ]

        plt.figure(figsize=(8,5))

        plt.hist(
            values,
            bins=20,
        )

        plt.xlabel("Strongest Anchor ΔG")

        plt.ylabel("Count")

        plt.title("Docking Score Distribution")

        plt.tight_layout()

        output = self.output_dir / "docking_histogram.png"

        plt.savefig(output, dpi=300)

        plt.close()

        return output

    # -----------------------------------------------------

    def score_scatter(
        self,
        candidates: list[Candidate],
    ):

        x = [
            c.c_score
            for c in candidates
        ]

        y = [
            c.strongest_anchor_delta_g
            for c in candidates
            if c.strongest_anchor_delta_g is not None
        ]

        valid = [
            c
            for c in candidates
            if c.strongest_anchor_delta_g is not None
        ]

        x = [c.c_score for c in valid]
        y = [c.strongest_anchor_delta_g for c in valid]

        plt.figure(figsize=(6,6))

        plt.scatter(
            x,
            y,
        )

        plt.xlabel("Native Fitness")

        plt.ylabel("Docking ΔG")

        plt.title("Fitness vs Docking")

        plt.tight_layout()

        output = self.output_dir / "fitness_vs_docking.png"

        plt.savefig(output, dpi=300)

        plt.close()

        return output

    # -----------------------------------------------------

    def convergence_plot(
        self,
        best_scores: list[float],
    ):

        plt.figure(figsize=(8,5))

        plt.plot(best_scores)

        plt.xlabel("Generation")

        plt.ylabel("Best Score")

        plt.title("Evolutionary Convergence")

        plt.tight_layout()

        output = self.output_dir / "convergence.png"

        plt.savefig(output, dpi=300)

        plt.close()

        return output

    # -----------------------------------------------------

    def generate_all(
        self,
        candidates: list[Candidate],
        best_scores: list[float],
    ):

        return {

            "fitness":

                self.fitness_histogram(
                    candidates
                ),

            "docking":

                self.docking_histogram(
                    candidates
                ),

            "scatter":

                self.score_scatter(
                    candidates
                ),

            "convergence":

                self.convergence_plot(
                    best_scores
                ),
        }