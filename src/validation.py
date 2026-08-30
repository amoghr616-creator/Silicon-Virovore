"""
validation.py

Scientific validation framework for Silicon Virovore.

Ensures every computational stage produced
scientifically valid outputs before reports
or downstream analyses are generated.

Future upgrades:
    • Experimental benchmark validation
    • AlphaFold confidence validation
    • Docking reproducibility statistics
    • MD convergence validation
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import platform
import sys
import datetime


# ------------------------------------------------------------

@dataclass(slots=True)
class ValidationResult:

    passed: bool

    errors: list[str]

    warnings: list[str]

    metadata: dict


# ------------------------------------------------------------

class Validator:

    def __init__(self):

        self.errors = []

        self.warnings = []

    # ---------------------------------------------------------

    def validate_sequences(self, candidates):

        sequences = set()

        for c in candidates:

            seq = c.sequence

            if not seq:

                self.errors.append(
                    "Empty sequence detected."
                )

            if len(seq) != 21:

                self.errors.append(
                    f"Invalid sequence length: {seq}"
                )

            if seq in sequences:

                self.warnings.append(
                    f"Duplicate sequence: {seq}"
                )

            sequences.add(seq)

    # ---------------------------------------------------------

    def validate_scores(self, candidates):

        for c in candidates:

            if c.overall_score != c.overall_score:

                self.errors.append(
                    f"NaN score: {c.sequence}"
                )

            if abs(c.overall_score) > 1e6:

                self.errors.append(
                    f"Impossible score: {c.sequence}"
                )

    # ---------------------------------------------------------

    def validate_docking(self, candidates):

        for c in candidates:

            if c.strongest_anchor_delta_g > 100:

                self.errors.append(
                    f"Invalid docking score: {c.sequence}"
                )

    # ---------------------------------------------------------

    def validate_structures(self, structure_directory):

        structure_directory = Path(structure_directory)

        if not structure_directory.exists():

            self.errors.append(
                "Structure directory missing."
            )

            return

        pdbs = list(
            structure_directory.glob("*.pdb")
        )

        if len(pdbs) == 0:

            self.errors.append(
                "No predicted structures found."
            )

    # ---------------------------------------------------------

    def metadata(self):

        return {

            "timestamp":

                datetime.datetime.now().isoformat(),

            "python":

                sys.version,

            "platform":

                platform.platform(),
        }

    # ---------------------------------------------------------

    def validate(

        self,

        candidates,

        structure_directory,

    ):

        self.validate_sequences(candidates)

        self.validate_scores(candidates)

        self.validate_docking(candidates)

        self.validate_structures(
            structure_directory
        )

        return ValidationResult(

            passed=len(self.errors) == 0,

            errors=self.errors,

            warnings=self.warnings,

            metadata=self.metadata(),
        )