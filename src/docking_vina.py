"""
docking_vina.py

Docking stage for Silicon Virovore.

Evaluates peptide candidates using a lightweight surrogate model
and optionally escalates promising candidates to AutoDock Vina.

This module is intentionally backend-agnostic so future docking
engines (DiffDock, Boltz-2, Rosetta, GNINA, etc.) can be swapped
without changing the remainder of the pipeline.
"""

import logging
import math
import os
import re
import shutil
import statistics
import subprocess
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from src.config import (
    DOCKING_DIR,
    PEPTIDE_FRAGMENT_SIZE,
    RECEPTOR_PDBQT,
    TIER2_DOCKING_THRESHOLD,
    VINA_BOX_PADDING,
    VINA_CPU,
    VINA_EXECUTABLE,
    VINA_EXHAUSTIVENESS,
    VINA_MAX_EVALS,
    VINA_NUM_MODES,
)
from src.models import Candidate, candidate_has_valid_structure

logger = logging.getLogger(__name__)

TIER2_THRESHOLD = TIER2_DOCKING_THRESHOLD


# ============================================================
# Docking Backends
# ============================================================

class DockingBackend:
    """
    Abstract docking backend.
    """

    def score_fragment(self, fragment: str) -> float:
        raise NotImplementedError


class MLSurrogateBackend(DockingBackend):
    """
    Lightweight heuristic docking surrogate.

    This backend is deterministic, local, and used only for Tier-1
    computational screening. It is not an experimental affinity.
    """

    HYDROPHOBIC = set("AILMFWYV")
    CHARGED = set("RHKDE")

    def __init__(self, cache_enabled: bool = True):
        self.cache_enabled = cache_enabled
        self._cache: dict[str, float] = {}
        self.cache_hits = 0
        self.cache_misses = 0

    def score_fragment(self, fragment: str) -> float:
        if not fragment:
            raise ValueError("Cannot score an empty fragment.")
        if self.cache_enabled and fragment in self._cache:
            self.cache_hits += 1
            return self._cache[fragment]

        self.cache_misses += 1

        score = -5.0

        for aa in fragment:
            if aa in self.HYDROPHOBIC:
                score -= 0.30
            elif aa in self.CHARGED:
                score += 0.10

        result = round(score, 2)
        if self.cache_enabled:
            self._cache[fragment] = result
        return result


class AutoDockBackend(DockingBackend):
    """
    Optional structure-validation backend.

    This class is intentionally disabled unless a valid executable,
    receptor, and candidate structure are available.
    """

    def __init__(self, docking_directory: str | Path = DOCKING_DIR):

        self.vina = Path(
            os.environ.get(
                "VINA_EXECUTABLE",
                str(VINA_EXECUTABLE),
            )
        )

        self.obabel = (
            shutil.which("obabel")
            or shutil.which("babel")
        )

        self.receptor = RECEPTOR_PDBQT

        if not self.vina.exists():
            raise RuntimeError(
                f"Vina executable not found: {self.vina}"
            )

        if not os.access(self.vina, os.X_OK):
            raise RuntimeError(
                f"Vina executable is not executable: {self.vina}"
            )

        if self.obabel is None:
            raise RuntimeError(
                "Open Babel executable not found."
            )

        if not self.receptor.exists():
            raise RuntimeError(
                f"Missing receptor: {self.receptor}"
            )

        self.docking_directory = Path(docking_directory)
        self.docking_directory.mkdir(parents=True, exist_ok=True)

        self.receptor_clean = (
            self.docking_directory / "receptor_clean.pdbqt"
        )

        self._sanitize_pdbqt(
            self.receptor,
            self.receptor_clean,
            ligand=False,
        )

        self.center, self.size = self._receptor_box(
            self.receptor,
            padding=VINA_BOX_PADDING,
        )

    @staticmethod
    def _sanitize_pdbqt(
        source: Path,
        destination: Path,
        ligand: bool,
    ) -> None:
        """Keep only records accepted by the bundled Vina parser."""

        allowed_prefixes = (
            "ATOM  ",
            "HETATM",
            "TER",
            "END",
        )

        if ligand:
            allowed_prefixes += (
                "ROOT",
                "ENDROOT",
                "BRANCH",
                "ENDBRANCH",
                "TORSDOF",
            )

        lines = [
            line
            for line in source.read_text().splitlines()
            if line.startswith(allowed_prefixes)
        ]

        if not any(
            line.startswith(("ATOM  ", "HETATM"))
            for line in lines
        ):
            raise RuntimeError(
                f"No atom records found in PDBQT input: {source}"
            )

        destination.write_text("\n".join(lines) + "\n")

    @staticmethod
    def _receptor_box(
        receptor: Path,
        padding: float,
    ) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        """Calculate a reproducible whole-receptor Vina search box."""

        coordinates = []

        for line in receptor.read_text().splitlines():
            if not line.startswith(("ATOM  ", "HETATM")):
                continue

            try:
                coordinates.append(
                    tuple(
                        float(line[start:start + 8])
                        for start in (30, 38, 46)
                    )
                )
            except ValueError:
                continue

        if not coordinates:
            raise RuntimeError(
                f"No receptor coordinates found: {receptor}"
            )

        minimum = tuple(
            min(point[index] for point in coordinates)
            for index in range(3)
        )
        maximum = tuple(
            max(point[index] for point in coordinates)
            for index in range(3)
        )

        center = tuple(
            (minimum[index] + maximum[index]) / 2
            for index in range(3)
        )
        size = tuple(
            maximum[index] - minimum[index] + 2 * padding
            for index in range(3)
        )

        return center, size

    def score_fragment(
        self,
        fragment: str,
        structure_path: str | Path | None = None,
    ) -> float:
        """
        Dock a predicted peptide structure with AutoDock Vina.

        The fragment is used for Tier-1 eligibility, while the complete
        predicted candidate structure is docked for Tier-2 validation.
        """

        if structure_path is None:
            raise RuntimeError(
                "Structure path required for validation."
            )

        structure_path = Path(structure_path)

        if not structure_path.exists():
            raise RuntimeError(
                f"Missing candidate structure: {structure_path}"
            )

        with tempfile.TemporaryDirectory(
            prefix="virovore-vina-",
        ) as temporary_directory:
            temporary_directory = Path(temporary_directory)
            ligand_raw = temporary_directory / "ligand.pdbqt"
            ligand_clean = temporary_directory / "ligand_clean.pdbqt"
            pose_path = self.docking_directory / f"{structure_path.stem}.pdbqt"

            conversion = subprocess.run(
                [
                    self.obabel,
                    str(structure_path),
                    "-O",
                    str(ligand_raw),
                    "--partialcharge",
                    "gasteiger",
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            if conversion.returncode != 0:
                raise RuntimeError(
                    "Open Babel failed: "
                    f"{conversion.stderr.strip()[-500:]}"
                )

            self._sanitize_pdbqt(
                ligand_raw,
                ligand_clean,
                ligand=True,
            )

            command = [
                str(self.vina),
                "--receptor",
                str(self.receptor_clean),
                "--ligand",
                str(ligand_clean),
                "--center_x",
                f"{self.center[0]:.3f}",
                "--center_y",
                f"{self.center[1]:.3f}",
                "--center_z",
                f"{self.center[2]:.3f}",
                "--size_x",
                f"{self.size[0]:.3f}",
                "--size_y",
                f"{self.size[1]:.3f}",
                "--size_z",
                f"{self.size[2]:.3f}",
                "--exhaustiveness",
                str(VINA_EXHAUSTIVENESS),
                "--max_evals",
                str(VINA_MAX_EVALS),
                "--num_modes",
                str(VINA_NUM_MODES),
                "--cpu",
                str(VINA_CPU),
                "--out",
                str(pose_path),
            ]

            docking = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False,
            )

            if docking.returncode != 0:
                raise RuntimeError(
                    "AutoDock Vina failed: "
                    f"{docking.stderr.strip()[-500:]}"
                )

            match = re.search(
                r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+",
                docking.stdout,
                re.MULTILINE,
            )

            if match is None:
                raise RuntimeError(
                    "AutoDock Vina returned no affinity score."
                )

            return float(match.group(1))
# ============================================================
# Docking Engine
# ============================================================

class PeptideDockingScorer:
    """
    Performs fragment-based docking evaluation.

    Computes:

    • Best fragment ΔG
    • Mean ΔG
    • Docking consistency
    • Consensus docking score
    • Optional Tier-2 AutoDock validation
    """

    def __init__(
        self,
        backend: DockingBackend | None = None,
        docking_directory: str | Path = DOCKING_DIR,
        workers: int = 1,
        cache_enabled: bool = True,
        tier2_threshold: float = TIER2_THRESHOLD,
    ):

        self.backend = backend or MLSurrogateBackend(
            cache_enabled=cache_enabled
        )
        self.workers = workers
        self.cache_enabled = cache_enabled
        self.tier2_threshold = tier2_threshold
        self._vina_cache: dict[tuple, float] = {}
        self._cache_lock = threading.Lock()
        self.vina_cache_hits = 0
        self.vina_cache_misses = 0

        try:

            self.vina_backend = AutoDockBackend(docking_directory)

        except Exception as exc:

            logger.warning(
                "AutoDock disabled: %s",
                exc,
            )

            self.vina_backend = None

    def cache_stats(self) -> dict[str, int]:
        """Return deterministic surrogate and Tier-2 cache counters."""

        return {
            "surrogate_hits": getattr(self.backend, "cache_hits", 0),
            "surrogate_misses": getattr(self.backend, "cache_misses", 0),
            "vina_hits": self.vina_cache_hits,
            "vina_misses": self.vina_cache_misses,
        }

    def evaluate_candidates(self, candidates: list[Candidate]) -> list[Candidate]:
        """Evaluate candidates serially or with bounded candidate-level workers."""

        if self.workers == 1:
            return [self.evaluate_candidate(candidate) for candidate in candidates]

        with ThreadPoolExecutor(max_workers=self.workers) as executor:
            return list(executor.map(self.evaluate_candidate, candidates))

    def _score_vina(
        self,
        fragment: str,
        structure_path: str | Path,
    ) -> float:
        if not self.cache_enabled or self.vina_backend is None:
            return self.vina_backend.score_fragment(fragment, structure_path)

        structure_path = Path(structure_path)
        stat = structure_path.stat()
        key = (
            fragment,
            str(structure_path.resolve()),
            stat.st_mtime_ns,
            VINA_EXHAUSTIVENESS,
            VINA_NUM_MODES,
            VINA_MAX_EVALS,
            VINA_CPU,
            VINA_BOX_PADDING,
        )
        with self._cache_lock:
            cached = self._vina_cache.get(key)
        if cached is not None:
            self.vina_cache_hits += 1
            return cached

        self.vina_cache_misses += 1
        score = float(self.vina_backend.score_fragment(fragment, structure_path))
        with self._cache_lock:
            self._vina_cache[key] = score
        return score

    # --------------------------------------------------------

    @staticmethod
    def _fragments(sequence: str) -> list[str]:
        """
        Generate overlapping peptide fragments.
        """

        if len(sequence) <= PEPTIDE_FRAGMENT_SIZE:
            return [sequence]

        return [
            sequence[i:i + PEPTIDE_FRAGMENT_SIZE]
            for i in range(
                len(sequence) - PEPTIDE_FRAGMENT_SIZE + 1
            )
        ]

    # --------------------------------------------------------

    @staticmethod
    def _consensus_score(
        best_delta_g: float,
        mean_delta_g: float,
    ) -> float:
        """
        Composite surrogate docking score.

        Lower (more negative) remains better, matching the
        sign convention of the underlying docking scores.
        """

        return round(
            (
                0.60 * best_delta_g
                + 0.40 * mean_delta_g
            ),
            3,
        )

    # --------------------------------------------------------

    def evaluate_candidate(
        self,
        candidate: Candidate,
    ) -> Candidate:

        logger.info(
            "Docking %s",
            candidate.sequence,
        )

        fragments = self._fragments(candidate.sequence)
        scores: list[float] = []
        valid_fragments: list[str] = []
        failed_fragments: dict[str, str] = {}

        candidate.fragment_scores.clear()
        candidate.docking_scores = []
        candidate.fragments = fragments
        candidate.best_fragment = None
        candidate.mean_delta_g = None
        candidate.strongest_anchor_delta_g = None
        candidate.vina_delta_g = None
        candidate.passed_tier_2 = False

        # Reset all evidence flags on every evaluation.
        candidate.metadata.update({
            "docking_backend": (
                "MLSurrogate"
                if isinstance(self.backend, MLSurrogateBackend)
                else type(self.backend).__name__
            ),
            "docking_is_surrogate": isinstance(
                self.backend,
                MLSurrogateBackend,
            ),
            "docking_validated": False,
            "docking_status": "pending",
            "tier2_eligible": False,
            "tier2_attempted": False,
            "tier2_validated": False,
            "tier2_status": "not_eligible",
            "vina_status": "not_attempted",
            "vina_score_available": False,
        })
        for key in (
            "docking_error",
            "docking_failed_fragments",
            "tier2_error",
            "tier2_backend",
        ):
            candidate.metadata.pop(key, None)

        for fragment in fragments:
            try:
                delta_g = float(self.backend.score_fragment(fragment))
                if not math.isfinite(delta_g):
                    raise ValueError("backend returned a non-finite score")
            except Exception as exc:
                failed_fragments[fragment] = str(exc)
                continue

            scores.append(delta_g)
            valid_fragments.append(fragment)
            candidate.fragment_scores[fragment] = delta_g

        candidate.docking_scores = scores

        if failed_fragments:
            candidate.metadata["docking_failed_fragments"] = failed_fragments

        if not scores:
            candidate.metadata["docking_status"] = "unavailable"
            candidate.metadata["docking_error"] = (
                "No fragment produced a valid docking score."
            )
            candidate.metadata["ranking_observation_valid"] = False
            logger.warning(
                "Docking unavailable | sequence=%s | failed_fragments=%d",
                candidate.sequence,
                len(failed_fragments),
            )
            return candidate

        candidate.mean_delta_g = round(sum(scores) / len(scores), 2)
        candidate.strongest_anchor_delta_g = min(scores)
        candidate.best_fragment = valid_fragments[scores.index(
            candidate.strongest_anchor_delta_g
        )]
        candidate.metadata["docking_std"] = (
            statistics.stdev(scores) if len(scores) > 1 else 0.0
        )
        consensus_score = self._consensus_score(
            candidate.strongest_anchor_delta_g,
            candidate.mean_delta_g,
        )
        candidate.metadata["consensus_docking"] = consensus_score
        candidate.metadata["docking_status"] = (
            "incomplete" if failed_fragments else "complete"
        )
        candidate.metadata["ranking_observation_valid"] = not failed_fragments

        # Tier 2 eligibility is a surrogate gate; it is not validation.
        tier2_eligible = (
            candidate.strongest_anchor_delta_g <= self.tier2_threshold
        )
        candidate.metadata["tier2_eligible"] = tier2_eligible

        if not tier2_eligible:
            candidate.metadata["tier2_status"] = "not_eligible"
        elif not candidate_has_valid_structure(candidate):
            candidate.metadata["tier2_status"] = "eligible_no_structure"
            candidate.add_note(
                "Tier-2 eligibility reached, but validation was skipped "
                "because no verified structure was available."
            )
        elif self.vina_backend is None:
            candidate.metadata["tier2_status"] = "unavailable"
            candidate.metadata["vina_status"] = "unavailable"
            candidate.add_note("Tier-2 validation unavailable.")
        else:
            candidate.metadata["tier2_attempted"] = True
            candidate.metadata["vina_status"] = "attempted"
            logger.info("Tier-2 validation triggered.")
            try:
                vina_score = self._score_vina(
                    candidate.best_fragment,
                    candidate.structure_path,
                )
                if not math.isfinite(vina_score):
                    raise ValueError("Vina returned a non-finite score")

                candidate.vina_delta_g = vina_score
                candidate.passed_tier_2 = True
                candidate.metadata["tier2_validated"] = True
                candidate.metadata["docking_validated"] = True
                candidate.metadata["tier2_status"] = "validated"
                candidate.metadata["tier2_backend"] = "AutoDock Vina"
                candidate.metadata["vina_status"] = "success"
                candidate.metadata["vina_score_available"] = True
                candidate.add_note(
                    "Tier-2 structure-based validation completed."
                )
            except Exception as exc:
                candidate.vina_delta_g = None
                candidate.metadata["tier2_status"] = "failed"
                candidate.metadata["tier2_error"] = str(exc)
                candidate.metadata["vina_status"] = "failed"
                candidate.add_note("Tier-2 validation failed.")
                logger.warning(
                    "Tier-2 validation failed for %s: %s",
                    candidate.sequence,
                    exc,
                )

        logger.info(
            (
                "Docking complete | "
                "Best %.2f | "
                "Mean %.2f | "
                "Consensus %.2f"
            ),
            candidate.strongest_anchor_delta_g,
            candidate.mean_delta_g,
            consensus_score,
        )

        return candidate

    # --------------------------------------------------------

    def evaluate(
        self,
        candidate: Candidate,
    ) -> Candidate:
        """
        Backwards-compatible alias.
        """

        return self.evaluate_candidate(
            candidate
        )
