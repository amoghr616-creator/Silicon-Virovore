"""
correlation.py

Statistical correlation analysis for Silicon Virovore.

This module evaluates agreement between
multiple scientific scoring metrics.

Future versions may include:
    - Partial correlations
    - Mutual information
    - PCA
    - Feature importance
"""

from __future__ import annotations

from itertools import combinations
from math import sqrt


class CorrelationAnalyzer:

    # ---------------------------------------------------------

    @staticmethod
    def pearson(x, y):
        """
        Pearson correlation coefficient.
        """

        if len(x) != len(y):
            raise ValueError("Vectors must have equal length.")

        n = len(x)

        if n < 2:
            return 0.0

        mean_x = sum(x) / n
        mean_y = sum(y) / n

        numerator = sum(
            (a - mean_x) * (b - mean_y)
            for a, b in zip(x, y)
        )

        denominator_x = sqrt(
            sum((a - mean_x) ** 2 for a in x)
        )

        denominator_y = sqrt(
            sum((b - mean_y) ** 2 for b in y)
        )

        if denominator_x == 0 or denominator_y == 0:
            return 0.0

        return numerator / (denominator_x * denominator_y)

    # ---------------------------------------------------------

    def metric_vectors(self, ranked_candidates):

        return {

            "overall_score": [
                rc.overall_score
                for rc in ranked_candidates
            ],

            "confidence": [
                rc.confidence
                for rc in ranked_candidates
            ],

            "c_score": [
                rc.candidate.c_score
                for rc in ranked_candidates
            ],

            "docking": [
                rc.candidate.strongest_anchor_delta_g
                for rc in ranked_candidates
            ],

            "helix": [
                rc.candidate.helix_propensity
                for rc in ranked_candidates
            ],

            "hydrophobicity": [
                rc.candidate.hydrophobic_moment
                for rc in ranked_candidates
            ],

            "solvation": [
                rc.candidate.solvation_energy
                for rc in ranked_candidates
            ],
        }

    # ---------------------------------------------------------

    def correlation_matrix(self, ranked_candidates):

        metrics = self.metric_vectors(
            ranked_candidates
        )

        matrix = {}

        for a, b in combinations(metrics.keys(), 2):

            r = self.pearson(
                metrics[a],
                metrics[b],
            )

            matrix[f"{a} vs {b}"] = round(r, 4)

        return matrix

    # ---------------------------------------------------------

    def strongest_correlations(self, ranked_candidates, top_n=5):

        corr = self.correlation_matrix(
            ranked_candidates
        )

        ranked = sorted(

            corr.items(),

            key=lambda x: abs(x[1]),

            reverse=True,

        )

        return ranked[:top_n]

    # ---------------------------------------------------------

    def summary(self, ranked_candidates):

        return {

            "correlations":

                self.correlation_matrix(
                    ranked_candidates
                ),

            "strongest":

                self.strongest_correlations(
                    ranked_candidates
                ),
        }