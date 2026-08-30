"""
diversity.py

Sequence diversity analysis for Silicon Virovore.

Provides quantitative measurements describing
how broadly the genetic algorithm explored
peptide sequence space.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from math import log2


class DiversityAnalyzer:

    # ---------------------------------------------------------

    @staticmethod
    def sequence_identity(seq1: str, seq2: str) -> float:
        """
        Pairwise sequence identity.
        """

        matches = sum(
            a == b
            for a, b in zip(seq1, seq2)
        )

        return matches / len(seq1)

    # ---------------------------------------------------------

    def mean_pairwise_identity(self, sequences):

        if len(sequences) < 2:
            return 1.0

        identities = [

            self.sequence_identity(a, b)

            for a, b in combinations(sequences, 2)
        ]

        return sum(identities) / len(identities)

    # ---------------------------------------------------------

    def mutation_frequency(
        self,
        sequences,
        reference,
    ):

        freq = []

        for i in range(len(reference)):

            mutations = sum(

                seq[i] != reference[i]

                for seq in sequences
            )

            freq.append(
                mutations / len(sequences)
            )

        return freq

    # ---------------------------------------------------------

    def amino_acid_frequency(
        self,
        sequences,
    ):

        result = []

        L = len(sequences[0])

        for i in range(L):

            counter = Counter(

                seq[i]

                for seq in sequences
            )

            result.append(counter)

        return result

    # ---------------------------------------------------------

    def shannon_entropy(
        self,
        sequences,
    ):

        entropy = []

        L = len(sequences[0])

        for i in range(L):

            counter = Counter(

                seq[i]

                for seq in sequences
            )

            total = sum(counter.values())

            H = 0

            for count in counter.values():

                p = count / total

                H -= p * log2(p)

            entropy.append(H)

        return entropy

    # ---------------------------------------------------------

    def summary(
        self,
        sequences,
        reference,
    ):

        return {

            "num_sequences":

                len(sequences),

            "unique_sequences":

                len(set(sequences)),

            "mean_pairwise_identity":

                round(
                    self.mean_pairwise_identity(
                        sequences
                    ),
                    4,
                ),

            "mutation_frequency":

                self.mutation_frequency(
                    sequences,
                    reference,
                ),

            "entropy":

                self.shannon_entropy(
                    sequences,
                ),
        }