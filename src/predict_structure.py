"""
predict_structure.py

Predict peptide structures using ESMFold and attach structural
metadata to Candidate objects.
"""

from __future__ import annotations

import logging
from pathlib import Path

import torch

from src import candidate
from src.config import STRUCTURE_DIR
from src.models import Candidate

logger = logging.getLogger(__name__)

STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)


class StructurePredictor:
    """
    Predicts peptide structures using ESMFold.

    Features
    --------
    • Loads the model only once.
    • Caches predicted PDB files.
    • Extracts mean pLDDT confidence.
    • Attaches structural metadata to Candidate.
    """

    _model = None
    _device = None

    def __init__(self):

        if StructurePredictor._model is None:

            logger.info("Loading ESMFold model...")

            import esm

            device = (
                "cuda"
                if torch.cuda.is_available()
                else "cpu"
            )

            model = esm.pretrained.esmfold_v1()

            model = model.eval().to(device)

            StructurePredictor._model = model
            StructurePredictor._device = device

            logger.info(
                "ESMFold loaded on %s",
                device,
            )

        self.model = StructurePredictor._model
        self.device = StructurePredictor._device

    @staticmethod
    def _extract_mean_plddt(
        pdb_string: str,
    ) -> float | None:
        """
        Extract mean pLDDT from the B-factor column
        of the generated PDB.
        """

        values = []

        for line in pdb_string.splitlines():

            if line.startswith("ATOM"):

                try:
                    values.append(
                        float(line[60:66])
                    )
                except ValueError:
                    continue

        if not values:
            return None

        return sum(values) / len(values)

    def predict(
        self,
        candidate: Candidate,
    ) -> Candidate:

        pdb_path = STRUCTURE_DIR / f"{candidate.sequence}.pdb"

        #######################################################
        # Cached structure
        #######################################################

        if pdb_path.exists():

            candidate.structure_path = pdb_path
            candidate.metadata["structure_available"] = True
            candidate.metadata["structure_backend"] = "ESMFold"
            candidate.metadata["structure_status"] = "cached"

            try:
                pdb_string = pdb_path.read_text()

                candidate.structure_confidence = (
                    self._extract_mean_plddt(
                        pdb_string
                    )
                )

            except Exception:

                pass

            return candidate

        #######################################################
        # Run ESMFold
        #######################################################

        try:

            with torch.inference_mode():

                pdb_string = self.model.infer_pdb(
                    candidate.sequence
                )

            pdb_path.write_text(pdb_string)

            candidate.structure_path = pdb_path

            candidate.metadata["structure_available"] = True
            candidate.metadata["structure_backend"] = "ESMFold"
            candidate.metadata["structure_status"] = "success"

            candidate.structure_confidence = (
                self._extract_mean_plddt(
                    pdb_string
                )
            )

            logger.info(
                "Predicted structure for %s (mean pLDDT %.2f)",
                candidate.sequence,
                candidate.structure_confidence
                if candidate.structure_confidence
                else -1.0,
            )

        except Exception as exc:

            logger.warning(
                "Structure prediction failed for %s: %s",
                candidate.sequence,
                exc,
            )

            candidate.structure_path = None
            candidate.structure_confidence = None

            candidate.metadata["structure_available"] = False
            candidate.metadata["structure_status"] = "failed"

        return candidate


def predict_population_structures(
    candidates: list[Candidate],
) -> list[Candidate]:
    """
    Predict structures for a population when the ESMFold
    environment is available.

    If ESMFold/OpenFold is unavailable, explicitly mark the
    structural stage as unavailable and return candidates
    without fabricated structural metrics.
    """

    try:
        predictor = StructurePredictor()

    except ImportError as exc:

        logger.warning(
            "ESMFold unavailable: required dependency missing (%s). "
            "Structural prediction disabled for this run.",
            exc,
        )

        for candidate in candidates:
            candidate.structure_path = None
            candidate.structure_confidence = None
            candidate.metadata["structure_available"] = False
            candidate.metadata["structure_backend"] = "ESMFold"
            candidate.metadata["structure_status"] = "unavailable"

        return candidates

    except Exception as exc:

        logger.warning(
            "ESMFold initialization failed: %s. "
            "Structural prediction disabled for this run.",
            exc,
        )

        for candidate in candidates:
            candidate.structure_path = None
            candidate.structure_confidence = None
            candidate.metadata["structure_available"] = False
            candidate.metadata["structure_backend"] = "ESMFold"
            candidate.metadata["structure_status"] = "failed"

        return candidates

    predicted = []

    for candidate in candidates:

        candidate = predictor.predict(candidate)

        candidate.metadata["structure_backend"] = "ESMFold"

        if candidate.structure_path is not None:
            candidate.metadata["structure_available"] = True
            candidate.metadata["structure_status"] = "success"
        else:
            candidate.metadata["structure_available"] = False
            candidate.metadata["structure_status"] = "failed"

        predicted.append(candidate)

    return predicted