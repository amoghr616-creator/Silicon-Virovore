"""
models.py

Shared data models for Silicon Virovore.

Every stage of the pipeline exchanges these dataclasses,
ensuring a consistent interface across receptor preparation,
population generation, docking, ranking, analysis, validation,
and reporting.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


# ============================================================
# Candidate
# ============================================================

@dataclass(slots=True)
class Candidate:
    """
    Represents one peptide throughout the entire pipeline.
    """

    # ============================================================
    # Identity
    # ============================================================

    sequence: str

    # ============================================================
    # Native Backend
    # ============================================================

    c_score: float = 0.0

    # ============================================================
    # Pipeline Notes
    # ============================================================

    notes: list[str] = field(default_factory=list)

    # ============================================================
    # Docking
    # ============================================================

    mean_delta_g: float | None = None
    strongest_anchor_delta_g: float | None = None
    passed_tier_2: bool = False
    vina_delta_g: float | None = None

    docking_scores: list[float] = field(default_factory=list)
    fragment_scores: dict[str, float] = field(default_factory=dict)
    best_fragment: str | None = None

    # ============================================================
    # Native Biophysical Metrics
    # ============================================================

    hydrophobic_moment: float = 0.0
    helix_propensity: float = 0.0
    solvation_energy: float = 0.0

    net_charge: float | None = None
    hydropathy_index: float | None = None

    # ============================================================
    # Molecular Dynamics
    # ============================================================

    rmsd: float | None = None
    rmsf: float | None = None
    radius_of_gyration: float | None = None
    md_stability_score: float | None = None

    # ============================================================
    # Safety Prediction
    # ============================================================

    toxicity_score: float | None = None
    hemolysis_score: float | None = None
    immunogenicity_score: float | None = None

    # ============================================================
    # Ranking
    # ============================================================

    confidence: float = 1.0
    overall_score: float = 0.0
    posthoc_recomputed_score: float | None = None
    rank: int = 0

    ranking_breakdown: dict[str, float] = field(default_factory=dict)

    # ============================================================
    # Structure Prediction
    # ============================================================

    structure_path: Path | None = None
    structure_confidence: float | None = None
    predicted_tm_score: float | None = None

    fragments: list[str] = field(default_factory=list)

    # ============================================================
    # Provenance
    # ============================================================

    created_at: str = ""
    updated_at: str = ""

    metadata: dict[str, Any] = field(default_factory=dict)

    # ============================================================
    # Methods
    # ============================================================

    def add_note(self, text: str) -> None:
        """Append a pipeline note."""
        self.notes.append(text)


def candidate_has_valid_structure(candidate: Candidate) -> bool:
    """Return whether a candidate has verified structural evidence."""

    path = candidate.structure_path

    return bool(
        candidate.metadata.get("structure_available") is True
        and path is not None
        and Path(path).is_file()
        and candidate.structure_confidence is not None
    )


def mark_structure_unavailable(
    candidate: Candidate,
    status: str = "unavailable",
    error: str | None = None,
) -> None:
    """Clear structural evidence without fabricating replacement values."""

    candidate.structure_path = None
    candidate.structure_confidence = None
    candidate.metadata["structure_available"] = False
    candidate.metadata["structure_status"] = status

    if error:
        candidate.metadata["structure_error"] = error
    else:
        candidate.metadata.pop("structure_error", None)
# ============================================================
# Ranked Candidate
# ============================================================

@dataclass(slots=True)
class RankedCandidate:
    """
    Candidate after ranking.
    """

    candidate: Candidate

    overall_score: float

    confidence: float

    rank: int = 0


# ============================================================
# Diversity Result
# ============================================================

@dataclass(slots=True)
class DiversityResult:

    entropy: list[float]

    mean_entropy: float

    pairwise_identity: float

    unique_sequences: int


# ============================================================
# Correlation Result
# ============================================================

@dataclass(slots=True)
class CorrelationResult:

    metric_a: str

    metric_b: str

    coefficient: float


# ============================================================
# Pareto Result
# ============================================================

@dataclass(slots=True)
class ParetoResult:

    sequence: str

    objectives: tuple[float, ...]

    dominates: int

    dominated_by: int

    front: int


# ============================================================
# Bootstrap Result
# ============================================================

@dataclass(slots=True)
class BootstrapResult:

    sequence: str

    wins: int

    frequency: float

    mean_score: float

    std_score: float


# ============================================================
# Validation Result
# ============================================================

@dataclass(slots=True)
class ValidationResult:

    passed: bool

    errors: list[str]

    warnings: list[str]

    metadata: dict[str, Any]


# ============================================================
# Pipeline Report
# ============================================================

@dataclass(slots=True)
class PipelineReport:

    top_candidates: list[RankedCandidate]

    diversity: DiversityResult | None = None

    correlations: list[CorrelationResult] = field(default_factory=list)

    pareto_front: list[ParetoResult] = field(default_factory=list)

    bootstrap: list[BootstrapResult] = field(default_factory=list)

    validation: ValidationResult | None = None

    runtime_seconds: float = 0.0

    timestamp: str = ""

    audit: dict[str, Any] = field(default_factory=dict)


# ============================================================
# Pipeline Configuration
# ============================================================

@dataclass(slots=True)
class PipelineConfig:

    population_size: int = 100

    generations: int = 20

    mutation_rate: float = 0.05

    elite_count: int = 10

    tier2_threshold: float = -7.0

    bootstrap_iterations: int = 100

    random_seed: int | None = None

    receptor_path: Path | None = None

    output_directory: Path | None = None

    def __post_init__(self) -> None:
        from src.config import validate_configuration

        validate_configuration(self)

@dataclass(slots=True)
class PositionStatistics:

    observations: int = 0
    cumulative_score: float = 0.0
    importance: float = 0.0

    mutations: int = 0
    successful_mutations: int = 0
