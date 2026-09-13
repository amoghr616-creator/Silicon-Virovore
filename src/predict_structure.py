"""Remote ESMFold prediction with verified caching and per-candidate isolation."""

from __future__ import annotations

import logging
import os
from pathlib import Path
import socket
import time

from src.config import (
    ESMFOLD_API_URL,
    ESMFOLD_MAX_ATTEMPTS,
    ESMFOLD_MAX_RETRY_DELAY_SECONDS,
    ESMFOLD_RETRY_BACKOFF_SECONDS,
    ESMFOLD_TIMEOUT_SECONDS,
    STRUCTURE_DIR,
)
from src.models import (
    Candidate,
    candidate_has_valid_structure,
    mark_structure_unavailable,
)

logger = logging.getLogger(__name__)
STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)

RETRYABLE_STATUS_CODES = frozenset({429, 500, 502, 503, 504})


class StructureServiceUnavailable(RuntimeError):
    """Raised when the remote structure service cannot be reached."""


class PermanentStructurePredictionError(RuntimeError):
    """Raised for an input or response error that should not be retried."""


def _is_dns_failure(error: Exception) -> bool:
    """Identify host-resolution failures wrapped by requests/urllib3."""

    current: BaseException | None = error
    while current is not None:
        if isinstance(current, socket.gaierror):
            return True
        current = current.__cause__ or current.__context__

    message = str(error)
    return any(
        marker in message
        for marker in (
            "NameResolutionError",
            "Failed to resolve",
            "Could not resolve host",
        )
    )


class StructurePredictor:
    """Predict peptide structures through the remote ESMFold API.

    Each candidate is independent: a failed request is recorded as missing
    evidence and does not stop the remainder of the population.
    """

    def __init__(
        self,
        api_url: str | None = None,
        timeout: int = ESMFOLD_TIMEOUT_SECONDS,
    ):
        try:
            import requests
        except ImportError as exc:
            raise ImportError(
                "The remote ESMFold backend requires the 'requests' package."
            ) from exc

        self._requests = requests
        self.api_url = api_url or os.environ.get(
            "ESMFOLD_API_URL",
            ESMFOLD_API_URL,
        )
        self.timeout = timeout

    def _request_pdb(self, sequence: str) -> str:
        """Request one PDB with bounded retries and backoff."""

        headers = {
            "Content-Type": "text/plain",
            "User-Agent": "Silicon-Virovore/1.0",
        }
        last_error: Exception | None = None

        for attempt in range(1, ESMFOLD_MAX_ATTEMPTS + 1):
            try:
                response = self._requests.post(
                    self.api_url,
                    headers=headers,
                    data=sequence,
                    timeout=self.timeout,
                )

                if response.status_code >= 400:
                    if response.status_code not in RETRYABLE_STATUS_CODES:
                        raise PermanentStructurePredictionError(
                            f"ESMFold API returned HTTP {response.status_code}"
                        )

                    last_error = RuntimeError(
                        f"ESMFold API returned HTTP {response.status_code}"
                    )
                    if attempt < ESMFOLD_MAX_ATTEMPTS:
                        delay = min(
                            ESMFOLD_MAX_RETRY_DELAY_SECONDS,
                            ESMFOLD_RETRY_BACKOFF_SECONDS * 2 ** (attempt - 1),
                        )
                        logger.warning(
                            "ESMFold request %d/%d failed with HTTP %d; "
                            "retrying in %ss.",
                            attempt,
                            ESMFOLD_MAX_ATTEMPTS,
                            response.status_code,
                            delay,
                        )
                        time.sleep(delay)
                        continue

                    raise StructureServiceUnavailable(str(last_error))

                return response.text

            except PermanentStructurePredictionError:
                raise
            except self._requests.RequestException as exc:
                last_error = exc
                if attempt < ESMFOLD_MAX_ATTEMPTS:
                    delay = min(
                        ESMFOLD_MAX_RETRY_DELAY_SECONDS,
                        ESMFOLD_RETRY_BACKOFF_SECONDS * 2 ** (attempt - 1),
                    )
                    logger.warning(
                        "ESMFold request %d/%d failed: %s; retrying in %ss.",
                        attempt,
                        ESMFOLD_MAX_ATTEMPTS,
                        exc,
                        delay,
                    )
                    time.sleep(delay)
                    continue

                if _is_dns_failure(exc):
                    raise StructureServiceUnavailable(
                        f"Could not resolve ESMFold host: {self.api_url}"
                    ) from exc
                raise StructureServiceUnavailable(str(exc)) from exc

        raise StructureServiceUnavailable(
            "ESMFold request failed after all attempts."
        ) from last_error

    @staticmethod
    def _extract_mean_plddt(pdb_string: str) -> float | None:
        """Extract and normalize mean pLDDT from PDB B-factors."""

        values = []
        for line in pdb_string.splitlines():
            if line.startswith("ATOM"):
                try:
                    values.append(float(line[60:66]))
                except ValueError:
                    continue

        if not values:
            return None

        mean_plddt = sum(values) / len(values)
        if 0 < mean_plddt <= 1:
            mean_plddt *= 100
        return round(mean_plddt, 2)

    @staticmethod
    def _valid_pdb(pdb_string: str) -> bool:
        return any(
            line.startswith(("ATOM  ", "HETATM"))
            for line in pdb_string.splitlines()
        )

    @classmethod
    def _set_success(
        cls,
        candidate: Candidate,
        pdb_path: Path,
        pdb_string: str,
        status: str,
    ) -> None:
        if not cls._valid_pdb(pdb_string):
            raise ValueError("PDB contains no atom records.")

        confidence = cls._extract_mean_plddt(pdb_string)
        if confidence is None:
            raise ValueError("PDB contains no usable pLDDT values.")

        candidate.structure_path = pdb_path
        candidate.structure_confidence = confidence
        candidate.metadata["structure_available"] = True
        candidate.metadata["structure_backend"] = "ESMFold"
        candidate.metadata["structure_status"] = status
        candidate.metadata.pop("structure_error", None)

        if not candidate_has_valid_structure(candidate):
            raise ValueError("Structure state failed validation.")

    def predict(self, candidate: Candidate) -> Candidate:
        """Predict one candidate without affecting other candidates."""

        pdb_path = STRUCTURE_DIR / f"{candidate.sequence}.pdb"
        candidate.structure_path = None
        candidate.structure_confidence = None
        candidate.metadata["structure_available"] = False
        candidate.metadata["structure_backend"] = "ESMFold"
        candidate.metadata["structure_status"] = "pending"
        candidate.metadata.pop("structure_error", None)

        if pdb_path.is_file():
            try:
                self._set_success(
                    candidate,
                    pdb_path,
                    pdb_path.read_text(),
                    "cached",
                )
                logger.info(
                    "ESMFold cached | sequence=%s | pLDDT=%.2f",
                    candidate.sequence,
                    candidate.structure_confidence,
                )
                return candidate
            except Exception as exc:
                logger.warning(
                    "ESMFold cache invalid | sequence=%s | reason=%s",
                    candidate.sequence,
                    exc,
                )

        try:
            pdb_string = self._request_pdb(candidate.sequence)
            pdb_path.write_text(pdb_string)
            self._set_success(candidate, pdb_path, pdb_string, "success")
            logger.info(
                "ESMFold success | sequence=%s | pLDDT=%.2f",
                candidate.sequence,
                candidate.structure_confidence,
            )
        except StructureServiceUnavailable as exc:
            mark_structure_unavailable(
                candidate,
                status="unavailable",
                error=str(exc),
            )
            candidate.metadata["structure_backend"] = "ESMFold"
            logger.warning(
                "ESMFold failure | sequence=%s | attempt=%d/%d | reason=%s",
                candidate.sequence,
                ESMFOLD_MAX_ATTEMPTS,
                ESMFOLD_MAX_ATTEMPTS,
                exc,
            )
        except Exception as exc:
            mark_structure_unavailable(
                candidate,
                status="failed",
                error=str(exc),
            )
            candidate.metadata["structure_backend"] = "ESMFold"
            logger.warning(
                "ESMFold failure | sequence=%s | attempt=%d/%d | reason=%s",
                candidate.sequence,
                ESMFOLD_MAX_ATTEMPTS,
                ESMFOLD_MAX_ATTEMPTS,
                exc,
            )

        return candidate


def _mark_population_unavailable(
    candidates: list[Candidate],
    error: str,
) -> None:
    for candidate in candidates:
        mark_structure_unavailable(
            candidate,
            status="unavailable",
            error=error,
        )
        candidate.metadata["structure_backend"] = "ESMFold"


def predict_population_structures(
    candidates: list[Candidate],
) -> list[Candidate]:
    """Predict each candidate independently and report evidence counts."""

    try:
        predictor = StructurePredictor()
    except ImportError as exc:
        logger.warning("ESMFold unavailable: %s", exc)
        _mark_population_unavailable(candidates, str(exc))
        return candidates
    except Exception as exc:
        logger.warning("ESMFold initialization failed: %s", exc)
        _mark_population_unavailable(candidates, str(exc))
        return candidates

    counts = {"success": 0, "cached": 0, "failed": 0, "unavailable": 0}
    total = len(candidates)

    for index, candidate in enumerate(candidates, start=1):
        logger.info(
            "ESMFold progress: %d/%d | %s",
            index,
            total,
            candidate.sequence,
        )
        predictor.predict(candidate)
        status = candidate.metadata.get("structure_status", "failed")
        counts[status] = counts.get(status, 0) + 1

    logger.info(
        "Structure prediction status: successful=%d | unavailable=%d | "
        "cached=%d | failed=%d | total=%d",
        counts.get("success", 0),
        counts.get("unavailable", 0),
        counts.get("cached", 0),
        counts.get("failed", 0),
        total,
    )
    return candidates
