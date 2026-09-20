"""
ranking.py

Multi-objective ranking engine for Silicon Virovore.

Combines every computational stage into one overall score.
"""

from __future__ import annotations
import logging

from src.models import Candidate, candidate_has_valid_structure

logger = logging.getLogger(__name__)

class CandidateRanker:
    """Rank candidates with an explicit missing-evidence policy.

    Missing weighted metrics are excluded from the candidate's denominator;
    they are never converted into a favorable numeric measurement. The
    candidate remains eligible for ranking on the metrics that were actually
    computed, while ``ranking_missing_evidence`` records what was absent.
    """

    def __init__(
        self,
        fitness_weight=0.30,
        docking_weight=0.30,
        structure_weight=0.10,
        helix_weight=0.08,
        hydro_weight=0.0,
        solvation_weight=0.06,
        md_weight=0.05,
        safety_weight=0.0,
    ):

        self.weights = {
            "fitness": fitness_weight,
            "docking": docking_weight,
            "structure": structure_weight,
            "helix": helix_weight,
            "hydrophobicity": hydro_weight,
            "solvation": solvation_weight,
            "md": md_weight,
            "safety": safety_weight,
        }

    def normalize(
        self,
        values,
        reverse=False,
    ):
        """
        Min-max normalize available values to [0, 1].

        Missing values remain missing and are represented by None.
        A metric with no variation among available observations
        receives 0.5 rather than falsely becoming perfect.
        """

        if not values:
            return []

        cleaned = [
            None if value is None else float(value)
            for value in values
        ]

        available = [
            value for value in cleaned
            if value is not None
        ]

        if not available:
            return [None] * len(cleaned)

        minimum = min(available)
        maximum = max(available)

        if abs(maximum - minimum) < 1e-9:
            normalized = [
                0.5 if value is not None else None
                for value in cleaned
            ]
        else:
            normalized = [
                None
                if value is None
                else (value - minimum) / (maximum - minimum)
                for value in cleaned
            ]

        if reverse:
            normalized = [
                None
                if value is None
                else 1.0 - value
                for value in normalized
            ]

        return normalized

    def score_raw_objective_records(
        self,
        raw_records: list[dict],
    ) -> list[float]:
        """Score raw objective records on one shared normalization scope.

        This is used by ARISE when it pools candidates from multiple
        generations.  The objective weights are the same as ``rank``; only
        the normalization scope differs.
        """

        if not raw_records:
            return []

        names = tuple(self.weights)
        raw_values = {
            name: [
                record.get(name)
                if isinstance(record, dict)
                else None
                for record in raw_records
            ]
            for name in names
        }
        # Match rank(): these native metrics are defined for every candidate
        # and use zero as their defensive fallback.
        for name in ("fitness", "helix", "hydrophobicity", "solvation"):
            raw_values[name] = [
                0.0 if value is None else value
                for value in raw_values[name]
            ]

        normalized = {
            "fitness": self.normalize(raw_values["fitness"]),
            "docking": self.normalize(raw_values["docking"], reverse=True),
            "structure": self.normalize(raw_values["structure"]),
            "helix": self.normalize(raw_values["helix"]),
            "hydrophobicity": self.normalize(raw_values["hydrophobicity"]),
            "solvation": self.normalize(
                raw_values["solvation"],
                reverse=True,
            ),
            "md": self.normalize(raw_values["md"]),
            "safety": self.normalize(raw_values["safety"], reverse=True),
        }

        scores = []
        for index in range(len(raw_records)):
            weighted_total = 0.0
            available_weight = 0.0
            for name, weight in self.weights.items():
                value = normalized[name][index]
                if value is None or weight <= 0.0:
                    continue
                weighted_total += value * weight
                available_weight += weight
            scores.append(
                round(weighted_total / available_weight, 6)
                if available_weight
                else 0.0
            )
        return scores
    # --------------------------------------------------------

    @staticmethod
    def _safe(
        value,
        default=0.0,
    ):

        if value is None:
            return default

        return value

    # --------------------------------------------------------

    @staticmethod
    def _metadata(
        candidate: Candidate,
        key: str,
        default=0.0,
    ):

        return candidate.metadata.get(
            key,
            default,
        )

    # --------------------------------------------------------

    @staticmethod
    def _confidence(
        candidate: Candidate,
    ) -> float:

        confidence = 0.0

        if (
            candidate_has_valid_structure(candidate)
            and candidate.structure_confidence is not None
        ):

            confidence += (
                candidate.structure_confidence / 100.0
            ) * 0.40

        if candidate.vina_delta_g is not None:

            confidence += 0.20

        if candidate.md_stability_score is not None:

            confidence += 0.20

        if candidate.toxicity_score is not None:

            confidence += 0.20

        return round(
            min(confidence, 1.0),
            3,
        )

    def rank(
        self,
        candidates: list[Candidate],
        *,
        score_field: str = "overall_score",
        metadata_prefix: str = "",
        normalization_scope: str = "current_candidate_set",
    ) -> list[Candidate]:

        if not candidates:
            return []

        # ----------------------------------------------------
        # Extract raw metrics
        # ----------------------------------------------------

        fitness = [
            self._safe(c.c_score)
            for c in candidates
        ]

        docking = [
            c.metadata.get("consensus_docking")
            if c.metadata.get("docking_status") in {
                "complete",
                "incomplete",
            }
            else None
            for c in candidates
        ]

        docking_consistency = [
            self._metadata(
                c,
                "docking_std",
                0.0,
            )
            for c in candidates
        ]

        structure = [
            c.structure_confidence
            if candidate_has_valid_structure(c)
            else None
            for c in candidates
        ]

        helix = [
            self._safe(
                c.helix_propensity
            )
            for c in candidates
        ]

        hydro = [
            self._safe(
                c.hydrophobic_moment
            )
            for c in candidates
        ]

        solvation = [
            self._safe(
                c.solvation_energy
            )
            for c in candidates
        ]

        md = [
        
            c.md_stability_score
            for c in candidates
        ]

        safety = [
          
            c.toxicity_score
            for c in candidates
        ]

        # ----------------------------------------------------
        # Normalize metrics
        # ----------------------------------------------------

        fitness_norm = self.normalize(
            fitness
        )

        docking_norm = self.normalize(
            docking,
            reverse=True,
        )

        structure_norm = self.normalize(
            structure
        )

        helix_norm = self.normalize(
            helix
        )

        hydro_norm = self.normalize(
            hydro
        )

        solvation_norm = self.normalize(
            solvation,
            reverse=True,
        )

        md_norm = self.normalize(
            md
        )

        safety_norm = self.normalize(
            safety,
            reverse=True,
        )

        # Lower docking standard deviation is better.
        consistency_norm = self.normalize(
            docking_consistency,
            reverse=True,
        )

        # ----------------------------------------------------
        # Candidate scoring
        # ----------------------------------------------------

        raw_objectives = {
            "fitness": fitness,
            "docking": docking,
            "structure": structure,
            "helix": helix,
            "hydrophobicity": hydro,
            "solvation": solvation,
            "md": md,
            "safety": safety,
        }
        normalized_objectives = {
            "fitness": fitness_norm,
            "docking": docking_norm,
            "structure": structure_norm,
            "helix": helix_norm,
            "hydrophobicity": hydro_norm,
            "solvation": solvation_norm,
            "md": md_norm,
            "safety": safety_norm,
        }
        computed_scores: dict[int, float] = {}

        for i, candidate in enumerate(candidates):
            breakdown = {}

            def weighted_component(normalized_value, weight):
                if normalized_value is None:
                    return None
                return normalized_value * weight

            breakdown["fitness"] = weighted_component(fitness_norm[i], self.weights["fitness"])
            breakdown["docking"] = weighted_component(docking_norm[i], self.weights["docking"])
            breakdown["structure"] = weighted_component(structure_norm[i], self.weights["structure"])
            breakdown["helix"] = weighted_component(helix_norm[i], self.weights["helix"])
            breakdown["hydrophobicity"] = weighted_component(hydro_norm[i], self.weights["hydrophobicity"])
            breakdown["solvation"] = weighted_component(solvation_norm[i], self.weights["solvation"])
            breakdown["md"] = weighted_component(md_norm[i], self.weights["md"])
            breakdown["safety"] = weighted_component(safety_norm[i], self.weights["safety"])

            available_weight = 0.0
            weighted_total = 0.0
            for key, value in breakdown.items():
                if value is None or self.weights[key] <= 0.0:
                    continue
                weighted_total += value
                available_weight += self.weights[key]

            score = round(weighted_total / available_weight, 6) if available_weight else 0.0
            computed_scores[id(candidate)] = score
            candidate.confidence = self._confidence(candidate)

            prefix = f"{metadata_prefix}_" if metadata_prefix else ""
            candidate.metadata[f"{prefix}raw_objectives"] = {
                key: values[i] for key, values in raw_objectives.items()
            }
            candidate.metadata[f"{prefix}normalized_objectives"] = {
                key: values[i] for key, values in normalized_objectives.items()
            }
            candidate.metadata[f"{prefix}ranking_total"] = round(
                sum(value for value in breakdown.values() if value is not None),
                6,
            )
            candidate.metadata[f"{prefix}ranking_components"] = {
                key: round(value, 6)
                for key, value in breakdown.items()
                if value is not None
            }
            candidate.metadata[f"{prefix}ranking_missing_evidence"] = [
                key
                for key, value in breakdown.items()
                if value is None and self.weights[key] > 0.0
            ]
            candidate.metadata[f"{prefix}available_weight"] = available_weight
            candidate.metadata[f"{prefix}normalization_scope"] = normalization_scope
            candidate.metadata[f"{prefix}score"] = score

            setattr(candidate, score_field, score)
            if metadata_prefix or score_field != "overall_score":
                candidate.posthoc_recomputed_score = score
            else:
                candidate.ranking_breakdown = {
                    key: value
                    for key, value in breakdown.items()
                    if value is not None
                }

            logger.debug(
                "RANK %s | scope=%s | score=%.4f | breakdown=%s",
                candidate.sequence,
                metadata_prefix or "generation",
                score,
                candidate.metadata[f"{prefix}ranking_components"],
            )

        candidates.sort(
            key=lambda c: (computed_scores[id(c)], c.confidence),
            reverse=True,
        )

        for rank, candidate in enumerate(candidates, start=1):
            if metadata_prefix:
                candidate.metadata[f"{metadata_prefix}_rank"] = rank
            else:
                candidate.rank = rank

        return candidates
if __name__ == "__main__":

    candidates = []

    for i in range(5):

        c = Candidate(
            sequence=f"SEQ{i}",
            c_score=i + 5,
        )

        c.strongest_anchor_delta_g = -6 - i

        c.helix_propensity = i * 0.3

        c.hydrophobic_moment = i * 0.4

        c.solvation_energy = -10 + i

        candidates.append(c)

    ranked = CandidateRanker().rank(candidates)

    for c in ranked:

        print(
            c.rank,
            c.sequence,
            round(c.overall_score, 3),
        )
