"""Pairwise mutation search and Shannon-entropy gating.

E_ij is a finite-difference interaction on Candidate.overall_score,
not proof of a biochemical interaction. In the current default
pipeline that score is sequence-level fitness, not validated
target-binding evidence.
"""
from __future__ import annotations

from collections import Counter
from itertools import combinations
from math import exp, isfinite, log2, sqrt
import random
from typing import Any, Iterable

AA_POOL = "ACDEFGHIKLMNPQRSTVWY"

PAIRWISE_MODES = frozenset({
    "uniform_pair_matched",
    "entropy_pair_matched",
    "epistasis_pair_matched",
    "epistasis_entropy",
    "hotspot_epistasis_entropy",
    "hotspot_epistasis_entropy_exploratory",
    "hotspot_entropy_substitution_epistasis_exploratory",
    "hotspot_entropy_substitution_epistasis_exploratory_v2",
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

        # Learned per-substitution effects and exact substitution-pair
        # interaction estimates. These are empirical sequence-score
        # estimates, not biochemical interaction energies.
        self.substitution_mean: dict[tuple[int, str, str], float] = {}
        self.substitution_count: dict[tuple[int, str, str], int] = {}
        self.substitution_pair_mean: dict[
            tuple[int, str, str, int, str, str], float
        ] = {}
        self.substitution_pair_count: dict[
            tuple[int, str, str, int, str, str], int
        ] = {}

        self.substitution_prior_count = 3.0
        self.substitution_temperature = 0.01
        self.substitution_random_probability = 0.25
        self.substitution_ucb_bonus = 0.002
        self.last_diagnostics: dict[str, Any] = {}

        # Telemetry for the latest population-generation call.
        self.last_selected_pairs: list[tuple[int, int]] = []
        self.last_pair_weights: list[
            tuple[tuple[int, int], float]
        ] = []
        self.last_selected_substitutions: list[
            tuple[int, str, str, int, str, str]
        ] = []
    def reset(self) -> None:
        """Reset learned state and generation-level telemetry."""
        self.entropy_by_position.clear()
        self.pair_mean.clear()
        self.pair_count.clear()
        self.substitution_mean.clear()
        self.substitution_count.clear()
        self.substitution_pair_mean.clear()
        self.substitution_pair_count.clear()
        self.last_diagnostics.clear()
        self.last_selected_pairs.clear()
        self.last_pair_weights.clear()
        self.last_selected_substitutions.clear()

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
        hotspot_positions: Iterable[int] = (),
        hotspot_boost: float = 2.0,
    ) -> list[float]:
        """Boost low-entropy positions while retaining nonzero sampling mass."""
        self._ensure_length(n)
        protected = set(protected_positions)
        if any(i < 0 or i >= n for i in protected):
            raise ValueError("protected position out of range")

        hotspots = set(hotspot_positions)
        if any(i < 0 or i >= n for i in hotspots):
            raise ValueError("hotspot position out of range")
        if not isfinite(hotspot_boost) or hotspot_boost <= 0:
            raise ValueError("hotspot_boost must be finite and positive")

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

        # A validated hotspot increases sampling weight, not certainty:
        # every position retains nonzero probability.
        for i in hotspots:
            raw[i] *= hotspot_boost

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

    def _update_running_mean(
        self,
        means: dict,
        counts: dict,
        key: tuple,
        value: float,
    ) -> None:
        """Online mean update for one substitution effect."""
        count = counts.get(key, 0)
        old_mean = means.get(key, 0.0)
        means[key] = old_mean + (value - old_mean) / (count + 1)
        counts[key] = count + 1

    def _substitution_pair_scores(
        self,
        parent: str,
        i: int,
        j: int,
    ) -> list[tuple[tuple[str, str], float]]:
        """Predict pair score changes from learned substitution evidence."""
        scored = []
        prior = self.substitution_prior_count

        for aa_i in AA_POOL:
            if aa_i == parent[i]:
                continue
            key_i = (i, parent[i], aa_i)
            count_i = self.substitution_count.get(key_i, 0)
            mean_i = self.substitution_mean.get(key_i, 0.0)
            estimate_i = mean_i * count_i / (count_i + prior)

            for aa_j in AA_POOL:
                if aa_j == parent[j]:
                    continue

                key_j = (j, parent[j], aa_j)
                count_j = self.substitution_count.get(key_j, 0)
                mean_j = self.substitution_mean.get(key_j, 0.0)
                estimate_j = mean_j * count_j / (count_j + prior)

                pair_key = (
                    i, parent[i], aa_i,
                    j, parent[j], aa_j,
                )
                pair_n = self.substitution_pair_count.get(pair_key, 0)
                pair_mean = self.substitution_pair_mean.get(pair_key, 0.0)
                pair_estimate = pair_mean * pair_n / (pair_n + prior)

                # A small, bounded optimism bonus keeps less-tested choices
                # in contention. It is not evidence that they are beneficial.
                support = count_i + count_j + 2 * pair_n
                bonus = self.substitution_ucb_bonus / sqrt(1 + support)
                predicted_delta = (
                    estimate_i + estimate_j + pair_estimate + bonus
                )
                scored.append(((aa_i, aa_j), predicted_delta))

        return scored

    def _single_substitution_scores(
        self,
        parent: str,
        position: int,
    ) -> list[tuple[str, float]]:
        """Return shrunken single-substitution gains plus a small UCB bonus."""
        prior = self.substitution_prior_count
        scored = []
        for aa in AA_POOL:
            if aa == parent[position]:
                continue
            key = (position, parent[position], aa)
            count = self.substitution_count.get(key, 0)
            mean = self.substitution_mean.get(key, 0.0)
            estimate = mean * count / (count + prior)
            bonus = self.substitution_ucb_bonus / sqrt(1 + count)
            scored.append((aa, estimate + bonus))
        return scored

    def _choose_single_substitution(
        self,
        parent: str,
        position: int,
        rng: random.Random,
    ) -> str:
        """Select a learned single mutation while retaining random exploration."""
        options = [aa for aa in AA_POOL if aa != parent[position]]
        if rng.random() < self.substitution_random_probability:
            return rng.choice(options)
        scored = self._single_substitution_scores(parent, position)
        best = max(score for _, score in scored)
        temperature = max(self.substitution_temperature, 1e-9)
        weights = [
            exp(max(-50.0, min(0.0, (score - best) / temperature)))
            for _, score in scored
        ]
        selected = self._draw(weights, rng)
        return scored[selected][0]

    def _choose_substitution_pair(
        self,
        parent: str,
        i: int,
        j: int,
        rng: random.Random,
    ) -> tuple[str, str]:
        """Use learned substitutions most of the time; retain random exploration."""
        options_i = [aa for aa in AA_POOL if aa != parent[i]]
        options_j = [aa for aa in AA_POOL if aa != parent[j]]

        if rng.random() < self.substitution_random_probability:
            return rng.choice(options_i), rng.choice(options_j)

        scored = self._substitution_pair_scores(parent, i, j)
        best = max(score for _, score in scored)
        temperature = max(self.substitution_temperature, 1e-9)

        weights = [
            exp(max(-50.0, min(0.0, (score - best) / temperature)))
            for _, score in scored
        ]
        selected_index = self._draw(weights, rng)
        return scored[selected_index][0]

    def generate_population(
        self,
        parent: str,
        *,
        pop_size: int,
        seed: int | None,
        mode: str,
        protected_positions: Iterable[int] = (),
        hotspot_positions: Iterable[int] = (),
        hotspot_boost: float = 2.0,
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
        self.last_selected_substitutions = []

        use_entropy = mode in {
            "entropy_pair_matched",
            "epistasis_entropy",
            "hotspot_epistasis_entropy",
            "hotspot_epistasis_entropy_exploratory",
            "hotspot_entropy_substitution_epistasis_exploratory",
            "hotspot_entropy_substitution_epistasis_exploratory_v2",
        }
        use_epistasis = mode in {
            "epistasis_pair_matched",
            "epistasis_entropy",
            "hotspot_epistasis_entropy",
            "hotspot_epistasis_entropy_exploratory",
            "hotspot_entropy_substitution_epistasis_exploratory",
            "hotspot_entropy_substitution_epistasis_exploratory_v2",
        }

        active_hotspots = (
            set(hotspot_positions)
            if mode in {
                "hotspot_epistasis_entropy",
                "hotspot_epistasis_entropy_exploratory",
                "hotspot_entropy_substitution_epistasis_exploratory",
                "hotspot_entropy_substitution_epistasis_exploratory_v2",
            }
            else set()
        )
        pos_weights = self.position_weights(
            n,
            use_entropy=use_entropy,
            protected_positions=protected_positions,
            hotspot_positions=active_hotspots,
            hotspot_boost=hotspot_boost,
        )

        if mode in {"hotspot_entropy_substitution_epistasis_exploratory", "hotspot_entropy_substitution_epistasis_exploratory_v2"}:
            # Do not let entropy alone decide where to search. Mix its
            # exploration distribution with a softmax over learned single-AA
            # gains, while retaining the explicit hotspot boost in both paths.
            position_scores = [
                max(score for _, score in self._single_substitution_scores(parent, i))
                for i in range(n)
            ]
            temperature = max(2.0 * self.substitution_temperature, 0.02)
            best_position_score = max(position_scores)
            learned_weights = [
                exp(max(-50.0, min(0.0, (score - best_position_score) / temperature)))
                for score in position_scores
            ]
            for i in active_hotspots:
                learned_weights[i] *= hotspot_boost
            learned_total = sum(learned_weights)
            learned_weights = [weight / learned_total for weight in learned_weights]
            pos_weights = [
                0.25 * exploration + 0.75 * learned
                for exploration, learned in zip(pos_weights, learned_weights)
            ]

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
        # Keep double-mutant interaction probes to roughly 20% of offspring
        # slots; the remaining budget is filled with additional singles.
        # This preserves matched factorial controls without spending one-third
        # of all evaluations on double mutants.
        requested_pair_budget = (
            2 if mode == "hotspot_entropy_substitution_epistasis_exploratory_v2" else max(1, (pop_size - 1) // 10)
        )
        pair_budget = min(
            (pop_size - 1) // 3,
            requested_pair_budget,
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
                if mode in {"hotspot_entropy_substitution_epistasis_exploratory", "hotspot_entropy_substitution_epistasis_exploratory_v2"}:
                    a, b = self._choose_substitution_pair(
                        parent, i, j, rng
                    )
                    self.last_selected_substitutions.append(
                        (i, parent[i], a, j, parent[j], b)
                    )
                else:
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
            if mode in {"hotspot_entropy_substitution_epistasis_exploratory", "hotspot_entropy_substitution_epistasis_exploratory_v2"}:
                aa = self._choose_single_substitution(parent, i, rng)
            else:
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

        substitution_updates = 0
        substitution_pair_updates = 0

        # Learn actual single-substitution score changes relative to this
        # generation's parent; pair interactions require matching singles
        # from the same parent background and generation.
        if base is not None:
            single_scores = {}
            double_records = []

            for seq, score in scores.items():
                changed = [
                    pos for pos, (old, new) in enumerate(zip(parent, seq))
                    if old != new
                ]
                if len(changed) == 1:
                    i = changed[0]
                    key = (i, parent[i], seq[i])
                    self._update_running_mean(
                        self.substitution_mean,
                        self.substitution_count,
                        key,
                        score - base,
                    )
                    single_scores[key] = score
                    substitution_updates += 1
                elif len(changed) == 2:
                    double_records.append(
                        (changed[0], changed[1], seq, score)
                    )

            for i, j, seq, double_score in double_records:
                key_i = (i, parent[i], seq[i])
                key_j = (j, parent[j], seq[j])
                score_i = single_scores.get(key_i)
                score_j = single_scores.get(key_j)

                # Do not invent a pair effect when either matched single is
                # missing from this generation.
                if score_i is None or score_j is None:
                    continue

                pair_key = (
                    i, parent[i], seq[i],
                    j, parent[j], seq[j],
                )
                interaction = double_score - score_i - score_j + base
                self._update_running_mean(
                    self.substitution_pair_mean,
                    self.substitution_pair_count,
                    pair_key,
                    interaction,
                )
                substitution_pair_updates += 1

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
            "substitution_updates": substitution_updates,
            "substitution_pair_updates": substitution_pair_updates,
            "learned_substitution_count": len(self.substitution_mean),
            "learned_substitution_pair_count": len(self.substitution_pair_mean),
        }
        return dict(self.last_diagnostics)