"""
predict_structure.py

Predict peptide structures with ESMFold and attach the
generated PDB path to Candidate objects.
"""

from __future__ import annotations

import logging
from pathlib import Path

import torch

from src.config import STRUCTURE_DIR
from src.models import Candidate

logger = logging.getLogger(__name__)

STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)


class StructurePredictor:
    """
    Uses ESMFold to predict peptide structures.

    Structures are cached on disk so they are only
    generated once.
    """

    def __init__(self):

        logger.info("Loading ESMFold model...")

        self.device = (
            "cuda"
            if torch.cuda.is_available()
            else "cpu"
        )

        import esm

        self.model = esm.pretrained.esmfold_v1()

        self.model = self.model.eval().to(self.device)

        logger.info(
            "ESMFold loaded on %s",
            self.device,
        )

    def predict(
        self,
        candidate: Candidate,
    ) -> Candidate:

        pdb_path = STRUCTURE_DIR / f"{candidate.sequence}.pdb"

        # Already predicted
        if pdb_path.exists():
            candidate.structure_path = pdb_path
            return candidate

        with torch.no_grad():

            pdb_string = self.model.infer_pdb(
                candidate.sequence
            )

        pdb_path.write_text(pdb_string)

        candidate.structure_path = pdb_path

        logger.info(
            "Predicted structure for %s",
            candidate.sequence,
        )

        return candidate


def predict_population_structures(
    candidates: list[Candidate],
) -> list[Candidate]:

    predictor = StructurePredictor()

    results = []

    for candidate in candidates:
        results.append(
            predictor.predict(candidate)
        )

    return results