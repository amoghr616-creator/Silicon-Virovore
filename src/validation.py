"""Consistency checks for computational pipeline evidence."""

from __future__ import annotations

from dataclasses import dataclass
import datetime
import math
from pathlib import Path
import platform
import sys

from src.config import AMINO_ACIDS, PEPTIDE_LENGTH
from src.models import candidate_has_valid_structure


@dataclass(slots=True)
class ValidationResult:
    passed: bool
    errors: list[str]
    warnings: list[str]
    metadata: dict


class Validator:
    def __init__(self):
        self.errors: list[str] = []
        self.warnings: list[str] = []

    def validate_sequences(self, candidates) -> None:
        seen = set()
        for candidate in candidates:
            sequence = candidate.sequence
            if not sequence:
                self.errors.append("Empty sequence detected.")
            if len(sequence) != PEPTIDE_LENGTH:
                self.errors.append(f"Invalid sequence length: {sequence}")
            if any(residue not in AMINO_ACIDS for residue in sequence):
                self.errors.append(f"Invalid amino acid: {sequence}")
            if sequence in seen:
                self.warnings.append(f"Duplicate sequence: {sequence}")
            seen.add(sequence)

    def validate_scores(self, candidates) -> None:
        for candidate in candidates:
            if not math.isfinite(candidate.overall_score):
                self.errors.append(f"Non-finite overall score: {candidate.sequence}")

    def validate_docking(self, candidates) -> None:
        for candidate in candidates:
            strongest = candidate.strongest_anchor_delta_g
            if strongest is not None and not math.isfinite(strongest):
                self.errors.append(f"Invalid docking score: {candidate.sequence}")
            if (
                candidate.mean_delta_g is not None
                and not math.isfinite(candidate.mean_delta_g)
            ):
                self.errors.append(f"Invalid mean docking score: {candidate.sequence}")

            validated = bool(candidate.metadata.get("tier2_validated", False))
            if candidate.passed_tier_2 != validated:
                self.errors.append(f"Tier-2 state mismatch: {candidate.sequence}")
            vina_status = candidate.metadata.get("vina_status")
            if validated and vina_status != "success":
                self.errors.append(
                    f"Tier-2 validated without successful Vina status: {candidate.sequence}"
                )
            if vina_status == "success" and candidate.vina_delta_g is None:
                self.errors.append(
                    f"Vina marked successful without a score: {candidate.sequence}"
                )
            if vina_status != "success" and candidate.vina_delta_g is not None:
                self.errors.append(
                    f"Vina score exists without successful status: {candidate.sequence}"
                )
            if validated and not candidate_has_valid_structure(candidate):
                self.errors.append(
                    f"Tier-2 passed without valid structure: {candidate.sequence}"
                )

            structure_available = candidate.metadata.get("structure_available")
            if structure_available is True and not candidate_has_valid_structure(candidate):
                self.errors.append(
                    f"Structure marked available but is invalid: {candidate.sequence}"
                )
            if structure_available is False and candidate.structure_confidence is not None:
                self.errors.append(
                    f"pLDDT exists while structure is unavailable: {candidate.sequence}"
                )

    def validate_structures(self, structure_directory) -> None:
        directory = Path(structure_directory)
        if not directory.exists():
            self.warnings.append("Structure directory missing.")
            return
        if not list(directory.glob("*.pdb")):
            self.warnings.append("No predicted structures found.")

    def metadata(self) -> dict:
        return {
            "timestamp": datetime.datetime.now().isoformat(),
            "python": sys.version,
            "platform": platform.platform(),
        }

    def validate(self, candidates, structure_directory) -> ValidationResult:
        self.errors = []
        self.warnings = []
        self.validate_sequences(candidates)
        self.validate_scores(candidates)
        self.validate_docking(candidates)
        self.validate_structures(structure_directory)
        return ValidationResult(
            passed=not self.errors,
            errors=self.errors,
            warnings=self.warnings,
            metadata=self.metadata(),
        )
