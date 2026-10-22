"""Pairwise mutation search and Shannon-entropy gating.

E_ij is a finite-difference interaction on Candidate.overall_score,
not proof of a biochemical interaction. In the current default
pipeline that score is sequence-level fitness, not validated
target-binding evidence.
"""
from __future__ import annotations

from collections import Counter
from itertools import combinations
from math import exp, isfinite, log2
import random
from typing import Any, Iterable

AA_POOL = "ACDEFGHIKLMNPQRSTVWY"

PAIRWISE_MODES = frozenset({
    "uniform_pair_matched",
    "entropy_pair_matched",
    "epistasis_pair_matched",
    "epistasis_entropy",
})


class EpistasisEntropyController:
    def __init__(
        self,
        *,
        entropy_floor: float = 0.5,
        entropy_boost: float = 2.0,
        pair_exploration: float = 0.25,
        pair_strength: float = 1.25,
        pair_prior_count: float = 3.0,
    ):
        if not isfinite(entropy_floor) or entropy_floor <= 0:
            raise ValueError(
                "entropy_floor must be finite and positive"
            )
        if not isfinite(entropy_boost) or entropy_boost < 0:
            raise ValueError(
                "entropy_boost must be finite and nonnegative"
            )
        if (
            not 0 <= pair_exploration <= 1
            or pair_strength < 0
            or pair_prior_count < 0
        ):
            raise ValueError(
                "invalid pairwise exploration parameters"
            )

        self.entropy_floor = entropy_floor
        self.entropy_boost = entropy_boost
        self.pair_exploration = pair_exploration
        self.pair_strength = pair_strength
        self.pair_prior_count = pair_prior_count

        self.entropy_by_position: list[float] = []
        self.pair_mean: list[list[float]] = []
        self.pair_count: list[list[int]] = []
        self.last_diagnostics: dict[str, Any] = {}

        # Telemetry for the latest population-generation call.
        self.last_selected_pairs: list[tuple[int, int]] = []
        self.last_pair_weights: list[
            tuple[tuple[int, int], float]
        ] = []
    def reset(self) -> None:
        """Reset learned state and generation-level telemetry."""
        self.entropy_by_position.clear()
        self.pair_mean.clear()
        self.pair_count.clear()
        self.last_diagnostics.clear()
        self.last_selected_pairs.clear()
        self.last_pair_weights.clear()

    def _ensure_length(self, n: int) -> None:
        if n < 1:
            raise ValueError("sequence length must be positive")
        if self.pair_mean and len(self.pair_mean) != n:
            raise ValueError("sequence length changed during this run")
        if not self.pair_mean:
            self.pair_mean = [[0.0] * n for _ in range(n)]
            self.pair_count = [[0] * n for _ in range(n)]

    @staticmethod
    def shannon_entropy(sequences: Iterable[str]) -> list[float]:
        """Per-position Shannon entropy in bits."""
        seqs = [s.upper() for s in sequences]
        if not seqs:
            return []
        n, length = len(seqs), len(seqs[0])
        if any(len(s) != length for s in seqs):
            raise ValueError("all sequences must have equal length")

        out = []
        for i in range(length):
            h = 0.0
            for count in Counter(s[i] for s in seqs).values():
                p = count / n
                h -= p * log2(p)
            out.append(h)
        return out

    def position_weights(
        self,
        n: int,
        *,
        use_entropy: bool = True,
        protected_positions: Iterable[int] = (),
    ) -> list[float]:
        """Boost low-entropy positions while retaining nonzero sampling mass."""
        self._ensure_length(n)
        protected = set(protected_positions)
        if any(i < 0 or i >= n for i in protected):
            raise ValueError("protected position out of range")

        raw = []
        for i in range(n):
            if (
                not use_entropy
                or i in protected
                or len(self.entropy_by_position) != n
            ):
                raw.append(1.0)
            else:
                deficit = max(
                    0.0,
                    self.entropy_floor - self.entropy_by_position[i],
                ) / self.entropy_floor
                raw.append(1.0 + self.entropy_boost * deficit)

        total = sum(raw)
        return [x / total for x in raw]

    @staticmethod
    def _draw(weights: list[float], rng: random.Random) -> int:
        total = sum(weights)
        if not weights or not isfinite(total) or total <= 0:
            raise ValueError("weights must have positive finite mass")

        draw, acc = rng.random() * total, 0.0
        for i, weight in enumerate(weights):
            acc += weight
            if draw < acc:
                return i
        return len(weights) - 1

    def _pair_weights(
        self,
        n: int,
        positions: list[float],
        use_epistasis: bool,
    ):
        pairs = list(combinations(range(n), 2))
        observed = [
            abs(self.pair_mean[i][j])
            for i, j in pairs
            if self.pair_count[i][j]
        ]
        scale = (
            max(sum(observed) / len(observed), 0.01)
            if observed
            else 0.02
        )

        base, guided = [], []
        for i, j in pairs:
            w = positions[i] * positions[j]
            base.append(w)

            multiplier = 1.0
            count = self.pair_count[i][j]
            if use_epistasis and count:
                confidence = count / (
                    count + self.pair_prior_count + 1e-12
                )
                z = max(
                    -2.0,
                    min(2.0, self.pair_mean[i][j] / scale),
                )
                multiplier = exp(
                    self.pair_strength * confidence * z
                )
            guided.append(w * multiplier)

        if use_epistasis:
            # Preserve exploration of pairs with weak or no evidence.
            weights = [
                (1 - self.pair_exploration) * g
                + self.pair_exploration * b
                for g, b in zip(guided, base)
            ]
        else:
            weights = base

        return list(zip(pairs, weights))

    @staticmethod
    def _mutate(parent: str, changes: dict[int, str]) -> str:
        child = list(parent)
        for i, aa in changes.items():
            child[i] = aa
        return "".join(child)

    def generate_population(
        self,
        parent: str,
        *,
        pop_size: int,
        seed: int | None,
        mode: str,
        protected_positions: Iterable[int] = (),
    ) -> list[str]:
        """Generate a factorial population and record pair-selection telemetry."""
        parent = parent.upper()

        if not parent or any(aa not in AA_POOL for aa in parent):
            raise ValueError(
                "parent must contain canonical amino acids"
            )
        if pop_size < 1 or mode not in PAIRWISE_MODES:
            raise ValueError(
                "invalid population size or pairwise mode"
            )

        n = len(parent)
        self._ensure_length(n)
        rng = random.Random(seed)

        # Clear only the most recent call's telemetry.
        # Do NOT clear pair_mean or pair_count here.
        self.last_selected_pairs = []
        self.last_pair_weights = []

        use_entropy = mode in {
            "entropy_pair_matched",
            "epistasis_entropy",
        }
        use_epistasis = mode in {
            "epistasis_pair_matched",
            "epistasis_entropy",
        }

        pos_weights = self.position_weights(
            n,
            use_entropy=use_entropy,
            protected_positions=protected_positions,
        )

        weighted_pairs = self._pair_weights(
            n,
            pos_weights,
            use_epistasis,
        )

        # Record the actual pair weights used in this generation.
        self.last_pair_weights = sorted(
            [
                (pair, float(weight))
                for pair, weight in weighted_pairs
            ],
            key=lambda item: item[1],
            reverse=True,
        )[:5]

        result = [parent]
        seen = {parent}
        available = list(weighted_pairs)

        # Each selected pair produces two single mutants and one
        # corresponding double mutant.
        pair_budget = min(
            (pop_size - 1) // 3,
            len(available),
        )

        for _ in range(pair_budget):
            k = self._draw(
                [weight for _, weight in available],
                rng,
            )
            (i, j), selected_weight = available.pop(k)

            self.last_selected_pairs.append((i, j))

            aa_i = [
                aa for aa in AA_POOL
                if aa != parent[i]
            ]
            aa_j = [
                aa for aa in AA_POOL
                if aa != parent[j]
            ]

            trio = None

            for _attempt in range(100):
                a = rng.choice(aa_i)
                b = rng.choice(aa_j)

                test = [
                    self._mutate(parent, {i: a}),
                    self._mutate(parent, {j: b}),
                    self._mutate(parent, {i: a, j: b}),
                ]

                if (
                    len(set(test)) == 3
                    and not any(sequence in seen for sequence in test)
                ):
                    trio = test
                    break

            if trio is None:
                raise RuntimeError(
                    "could not create unique factorial pair variants"
                )

            result.extend(trio)
            seen.update(trio)

        # Fill remaining slots with single mutants. These do not
        # count as additional pair-selection observations.
        while len(result) < pop_size:
            i = self._draw(pos_weights, rng)
            aa = rng.choice(
                [residue for residue in AA_POOL if residue != parent[i]]
            )
            sequence = self._mutate(parent, {i: aa})

            if sequence not in seen:
                result.append(sequence)
                seen.add(sequence)

        return result[:pop_size]
    def observe_generation(
        self,
        candidates: list[Any],
        *,
        parent_sequence: str,
    ) -> dict[str, Any]:
        """Update entropy and pairwise finite-difference estimates."""
        parent = parent_sequence.upper()
        n = len(parent)
        self._ensure_length(n)

        valid = [
            c for c in candidates
            if isinstance(getattr(c, "sequence", None), str)
            and len(c.sequence) == n
        ]
        if not valid:
            raise ValueError("generation has no valid sequences")

        self.entropy_by_position = self.shannon_entropy(
            c.sequence for c in valid
        )

        scores: dict[str, float] = {}

        for candidate in valid:
            try:
                raw_score = getattr(candidate, "c_score", None)
                if raw_score is None:
                    continue

                score = float(raw_score)

                if isfinite(score):
                    scores[candidate.sequence.upper()] = score

            except (TypeError, ValueError, AttributeError):
                continue

        base = scores.get(parent)
        updates = 0

        if base is not None:
            singles, doubles = {}, []
            for seq, score in scores.items():
                changed = [
                    i for i, (a, b) in enumerate(zip(parent, seq))
                    if a != b
                ]
                if len(changed) == 1:
                    i = changed[0]
                    singles[(i, seq[i])] = score
                elif len(changed) == 2:
                    doubles.append((changed[0], changed[1], seq, score))

            for i, j, seq, double_score in doubles:
                left = singles.get((i, seq[i]))
                right = singles.get((j, seq[j]))
                if left is None or right is None:
                    continue

                # E_ij = f(double) - f(single_i) - f(single_j) + f(parent)
                effect = double_score - left - right + base
                count = self.pair_count[i][j]
                mean = self.pair_mean[i][j] + (
                    effect - self.pair_mean[i][j]
                ) / (count + 1)

                self.pair_count[i][j] = self.pair_count[j][i] = count + 1
                self.pair_mean[i][j] = self.pair_mean[j][i] = mean
                updates += 1

        self.last_diagnostics = {
            "mean_entropy_bits": sum(self.entropy_by_position) / n,
            "low_entropy_positions": sum(
                h < self.entropy_floor for h in self.entropy_by_position
            ),
            "epistasis_updates": updates,
            "observed_pair_count": sum(
                1 for i in range(n) for j in range(i + 1, n)
                if self.pair_count[i][j] > 0
            ),
            "parent_score_available": base is not None,
            "interaction_matrix_kind": "sequence_objective_finite_difference",
        }
        return dict(self.last_diagnostics)