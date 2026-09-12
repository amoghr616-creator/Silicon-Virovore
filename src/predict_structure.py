"""
predict_structure.py

Predict peptide structures using ESMFold and attach structural
metadata to Candidate objects.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
import socket
import time

from src.config import (
    ESMFOLD_API_URL,
    ESMFOLD_MAX_ATTEMPTS,
    ESMFOLD_MAX_CONSECUTIVE_FAILURES,
    ESMFOLD_RETRY_BACKOFF_SECONDS,
    ESMFOLD_TIMEOUT_SECONDS,
    STRUCTURE_DIR,
)
from src.models import Candidate

logger = logging.getLogger(__name__)

STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)

RETRYABLE_STATUS_CODES = frozenset({
    429,
    500,
    502,
    503,
    504,
})


class StructureServiceUnavailable(RuntimeError):
    """Raised when the remote structure service cannot be reached."""


def _is_dns_failure(error: Exception) -> bool:
    """Identify host-resolution failures wrapped by requests/urllib3."""

    current: BaseException | None = error

    while current is not None:
        if isinstance(current, socket.gaierror):
            return True

        current = current.__cause__ or current.__context__

    message = str(error)
    return (
        "NameResolutionError" in message
        or "Failed to resolve" in message
        or "Could not resolve host" in message
    )


class StructurePredictor:
    """
    Predicts peptide structures using the remote ESMFold API.

    Features
    --------
    • Uses the remote ESMFold API.
    • Caches predicted PDB files.
    • Extracts mean pLDDT confidence.
    • Attaches structural metadata to Candidate.
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

    def _request_pdb(
        self,
        sequence: str,
    ) -> str:
        """Request one structure with bounded retry handling."""

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

                if response.status_code in RETRYABLE_STATUS_CODES:
                    last_error = RuntimeError(
                        f"ESMFold API returned HTTP {response.status_code}"
                    )

                    if attempt < ESMFOLD_MAX_ATTEMPTS:
                        delay = ESMFOLD_RETRY_BACKOFF_SECONDS * (
                            2 ** (attempt - 1)
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

                response.raise_for_status()
                return response.text

            except self._requests.RequestException as exc:
                last_error = exc

                if _is_dns_failure(exc):
                    raise StructureServiceUnavailable(
                        f"Could not resolve ESMFold host: {self.api_url}"
                    ) from exc

                if attempt < ESMFOLD_MAX_ATTEMPTS:
                    delay = ESMFOLD_RETRY_BACKOFF_SECONDS * (
                        2 ** (attempt - 1)
                    )
                    logger.warning(
                        "ESMFold request %d/%d failed: %s; "
                        "retrying in %ss.",
                        attempt,
                        ESMFOLD_MAX_ATTEMPTS,
                        exc,
                        delay,
                    )
                    time.sleep(delay)
                    continue

                raise

        raise RuntimeError(
            "ESMFold request failed after all attempts."
        ) from last_error

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

        mean_plddt = sum(values) / len(values)

        # ESM Atlas responses commonly store pLDDT in [0, 1], while
        # the rest of this project uses the conventional [0, 100] scale.
        if 0 < mean_plddt <= 1:
            mean_plddt *= 100

        return round(mean_plddt, 2)

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

            pdb_string = self._request_pdb(
                candidate.sequence,
            )

            if "ATOM" not in pdb_string:
                raise RuntimeError(
                    "ESMFold API returned no ATOM records."
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
                if candidate.structure_confidence is not None
                else -1.0,
            )

        except StructureServiceUnavailable as exc:

            logger.error(
                "ESMFold service unavailable for %s: %s",
                candidate.sequence,
                exc,
            )

            candidate.structure_path = None
            candidate.structure_confidence = None
            candidate.metadata["structure_available"] = False
            candidate.metadata["structure_status"] = "unavailable"

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
    Predict structures for a population through the remote ESMFold API.

    If the API client or remote service is unavailable, explicitly mark the
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
    consecutive_failures = 0
    total = len(candidates)

    for index, candidate in enumerate(candidates, start=1):

        logger.info(
            "ESMFold progress: %d/%d | %s",
            index,
            total,
            candidate.sequence,
        )

        candidate = predictor.predict(candidate)

        candidate.metadata["structure_backend"] = "ESMFold"

        if candidate.structure_path is not None:
            candidate.metadata["structure_available"] = True
            candidate.metadata["structure_status"] = "success"
            consecutive_failures = 0
        elif candidate.metadata.get("structure_status") == "unavailable":
            predicted.append(candidate)

            remaining = candidates[index:]
            logger.warning(
                "Stopping ESMFold requests because the service is "
                "unavailable; %d candidates marked unavailable.",
                len(remaining),
            )

            for remaining_candidate in remaining:
                remaining_candidate.structure_path = None
                remaining_candidate.structure_confidence = None
                remaining_candidate.metadata[
                    "structure_available"
                ] = False
                remaining_candidate.metadata[
                    "structure_backend"
                ] = "ESMFold"
                remaining_candidate.metadata[
                    "structure_status"
                ] = "unavailable"

            predicted.extend(remaining)
            break
        else:
            candidate.metadata["structure_available"] = False
            candidate.metadata["structure_status"] = "failed"
            consecutive_failures += 1

        predicted.append(candidate)

        if consecutive_failures >= ESMFOLD_MAX_CONSECUTIVE_FAILURES:
            remaining = candidates[index:]
            logger.warning(
                "Stopping ESMFold requests after %d consecutive failures; "
                "%d candidates marked unavailable.",
                consecutive_failures,
                len(remaining),
            )

            for remaining_candidate in remaining:
                remaining_candidate.structure_path = None
                remaining_candidate.structure_confidence = None
                remaining_candidate.metadata[
                    "structure_available"
                ] = False
                remaining_candidate.metadata[
                    "structure_backend"
                ] = "ESMFold"
                remaining_candidate.metadata[
                    "structure_status"
                ] = "unavailable"

            predicted.extend(remaining)
            break

    return predicted
