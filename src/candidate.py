"""
candidate.py

Unified data model for Silicon Virovore.

Every pipeline stage updates the same Candidate object.
"""

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Candidate:

    # Identity
    sequence: str

    # Generated peptide structure
    structure_path: Path | None = None

    # Native C engine
    c_score: float = 0.0

    hydrophobic_moment: float = 0.0
    solvation_energy: float = 0.0
    helix_propensity: float = 0.0

    # Docking
    mean_delta_g: float = 0.0
    strongest_anchor_delta_g: float = 0.0
    vina_delta_g: float | None = None

    # Safety
    toxicity_score: float | None = None
    aggregation_score: float | None = None
    solubility_score: float | None = None

    # Molecular Dynamics
    rmsd: float | None = None
    rmsf: float | None = None

    # Final ranking
    overall_score: float = 0.0
    confidence: float = 0.0

    notes: list[str] = field(default_factory=list)

    def add_note(self, text: str):
        self.notes.append(text)