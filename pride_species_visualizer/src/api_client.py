"""PRIDE Archive REST API client for species-level dataset counts."""

from __future__ import annotations

import logging
import time
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from requests import Response, Session

logger = logging.getLogger(__name__)

BASE_URL = "https://www.ebi.ac.uk/pride/ws/archive/v3"
SEARCH_PROJECTS = f"{BASE_URL}/search/projects"
PROJECTS_COUNT = f"{BASE_URL}/projects/count"
ORGANISMS_COUNT = f"{BASE_URL}/findAllOrganismsCount"

DEFAULT_PAGE_SIZE = 100
DEFAULT_TIMEOUT = (10, 120)
USER_AGENT = "pride-species-visualizer/1.0.0"
MAX_RETRIES = 5
RETRY_BACKOFF_SECONDS = 2.0
EXPECTED_SPECIES_COUNT = 4604
SPECIES_COUNT_TOLERANCE = 500

UNKNOWN_SPECIES_LABEL = "Unknown"


class PrideApiError(Exception):
    """Raised when the PRIDE Archive API returns an unrecoverable error."""


def create_session() -> Session:
    """Create a configured HTTP session for PRIDE API requests.

    Returns:
        A ``requests.Session`` with timeouts and a descriptive User-Agent.
    """
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT, "Accept": "application/json"})
    return session


def _request_with_retry(
    session: Session,
    url: str,
    params: dict[str, Any] | None = None,
) -> Response:
    """Perform a GET request with exponential backoff on transient failures.

    Args:
        session: Active HTTP session.
        url: Request URL.
        params: Optional query parameters.

    Returns:
        Successful HTTP response.

    Raises:
        PrideApiError: If all retry attempts fail or the response is not OK.
    """
    last_error: Exception | None = None
    for attempt in range(MAX_RETRIES):
        try:
            response = session.get(url, params=params, timeout=DEFAULT_TIMEOUT)
            if response.status_code in (429, 500, 502, 503, 504):
                wait = RETRY_BACKOFF_SECONDS * (2**attempt)
                logger.warning(
                    "HTTP %s from %s; retrying in %.1fs (attempt %d/%d)",
                    response.status_code,
                    url,
                    wait,
                    attempt + 1,
                    MAX_RETRIES,
                )
                time.sleep(wait)
                continue
            if 400 <= response.status_code < 500:
                response.raise_for_status()
            response.raise_for_status()
            return response
        except requests.RequestException as exc:
            last_error = exc
            wait = RETRY_BACKOFF_SECONDS * (2**attempt)
            logger.warning(
                "Request failed for %s: %s; retrying in %.1fs (attempt %d/%d)",
                url,
                exc,
                wait,
                attempt + 1,
                MAX_RETRIES,
            )
            time.sleep(wait)

    message = f"Failed to fetch {url} after {MAX_RETRIES} attempts"
    if last_error is not None:
        raise PrideApiError(message) from last_error
    raise PrideApiError(message)


def get_total_project_count(session: Session) -> int:
    """Return the total number of projects in the PRIDE Archive.

    Args:
        session: Active HTTP session.

    Returns:
        Integer project count from ``/projects/count``.

    Raises:
        PrideApiError: If the count endpoint cannot be reached.
    """
    response = _request_with_retry(session, PROJECTS_COUNT)
    try:
        return int(response.text.strip())
    except ValueError as exc:
        raise PrideApiError(f"Invalid project count response: {response.text!r}") from exc


def get_expected_organism_count(session: Session) -> int | None:
    """Return PRIDE's reported distinct organism count, if available.

    Args:
        session: Active HTTP session.

    Returns:
        Distinct organism count, or ``None`` if the endpoint fails.
    """
    try:
        response = session.get(
            ORGANISMS_COUNT,
            timeout=DEFAULT_TIMEOUT,
            headers={"Accept": "text/plain, */*"},
        )
        response.raise_for_status()
        return int(response.text.strip())
    except (requests.RequestException, ValueError) as exc:
        logger.warning("Could not fetch organism count sanity check: %s", exc)
        return None


def fetch_projects_page(
    session: Session,
    page: int,
    page_size: int = DEFAULT_PAGE_SIZE,
) -> list[dict[str, Any]]:
    """Fetch one page of project metadata from the search endpoint.

    Args:
        session: Active HTTP session.
        page: Zero-based page index.
        page_size: Number of projects per page.

    Returns:
        List of project dictionaries, or an empty list when no results remain.

    Raises:
        PrideApiError: If the response is not a JSON list.
    """
    params = {"page": page, "pageSize": page_size}
    response = _request_with_retry(session, SEARCH_PROJECTS, params=params)
    data = response.json()
    if not isinstance(data, list):
        raise PrideApiError(
            f"Expected JSON list from search/projects, got {type(data).__name__}"
        )
    return data


def aggregate_organism_counts(
    projects: list[dict[str, Any]],
    *,
    bucket_unknown: bool = True,
) -> Counter[str]:
    """Count dataset occurrences per organism label from project metadata.

    Each organism listed on a project increments that species count by one.
    Multi-organism projects therefore contribute to every listed species.

    Args:
        projects: Project metadata dicts from the search API.
        bucket_unknown: If True, count projects with no organisms as ``Unknown``.

    Returns:
        Counter mapping species name to dataset count for this batch.
    """
    counts: Counter[str] = Counter()
    for project in projects:
        organisms = project.get("organisms") or []
        if not organisms:
            if bucket_unknown:
                counts[UNKNOWN_SPECIES_LABEL] += 1
            continue
        for organism in organisms:
            label = str(organism).strip()
            if label:
                counts[label] += 1
    return counts


def _counts_to_dataframe(counts: Counter[str]) -> pd.DataFrame:
    """Convert a species counter to a sorted DataFrame.

    Args:
        counts: Species name to dataset count mapping.

    Returns:
        DataFrame with columns ``species`` and ``dataset_count``, sorted descending.
    """
    if not counts:
        return pd.DataFrame(columns=["species", "dataset_count"])
    df = pd.DataFrame(
        [
            {"species": species, "dataset_count": int(count)}
            for species, count in counts.items()
        ]
    )
    return df.sort_values("dataset_count", ascending=False).reset_index(drop=True)


def fetch_all_species_counts(
    page_size: int = DEFAULT_PAGE_SIZE,
    *,
    show_progress: bool = True,
    bucket_unknown: bool = True,
) -> pd.DataFrame:
    """Paginate all PRIDE projects and aggregate dataset counts per species.

    Args:
        page_size: Projects requested per API page.
        show_progress: Log pagination progress at INFO level.
        bucket_unknown: Count projects without organisms as ``Unknown``.

    Returns:
        DataFrame with columns ``species`` and ``dataset_count``.

    Raises:
        PrideApiError: On unrecoverable API failures.
    """
    session = create_session()
    total_projects = get_total_project_count(session)
    expected_organisms = get_expected_organism_count(session)

    if show_progress:
        logger.info("Total PRIDE projects to scan: %d", total_projects)

    counts: Counter[str] = Counter()
    page = 0
    projects_processed = 0

    while True:
        batch = fetch_projects_page(session, page, page_size)
        if not batch:
            break

        counts.update(
            aggregate_organism_counts(batch, bucket_unknown=bucket_unknown)
        )
        projects_processed += len(batch)
        page += 1

        if show_progress and (page == 1 or page % 10 == 0):
            logger.info(
                "Fetched page %d (%d/%d projects, %d species so far)",
                page,
                min(projects_processed, total_projects),
                total_projects,
                len(counts),
            )

        if projects_processed >= total_projects:
            break

    if show_progress:
        logger.info(
            "Finished: %d projects processed, %d distinct species",
            projects_processed,
            len(counts),
        )

    df = _counts_to_dataframe(counts)

    if expected_organisms is not None:
        if abs(len(df) - expected_organisms) > SPECIES_COUNT_TOLERANCE:
            logger.warning(
                "Species count %d differs from PRIDE organism count %d by more than %d",
                len(df),
                expected_organisms,
                SPECIES_COUNT_TOLERANCE,
            )
    elif abs(len(df) - EXPECTED_SPECIES_COUNT) > SPECIES_COUNT_TOLERANCE:
        logger.warning(
            "Species count %d differs from expected ~%d",
            len(df),
            EXPECTED_SPECIES_COUNT,
        )

    return df


def save_counts_cache(df: pd.DataFrame, path: str | Path) -> None:
    """Write aggregated species counts to a CSV cache file.

    Args:
        df: DataFrame with ``species`` and ``dataset_count`` columns.
        path: Output CSV path (parent directories are created).
    """
    cache_path = Path(path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache_path, index=False)
    logger.info("Wrote species counts cache to %s", cache_path)


def load_counts_cache(path: str | Path) -> pd.DataFrame:
    """Load aggregated species counts from a CSV cache file.

    Args:
        path: Path to a CSV written by ``save_counts_cache``.

    Returns:
        DataFrame with ``species`` and ``dataset_count`` columns.

    Raises:
        FileNotFoundError: If the cache file does not exist.
        ValueError: If required columns are missing.
    """
    cache_path = Path(path)
    if not cache_path.is_file():
        raise FileNotFoundError(f"Cache file not found: {cache_path}")

    df = pd.read_csv(cache_path)
    required = {"species", "dataset_count"}
    if not required.issubset(df.columns):
        raise ValueError(
            f"Cache must contain columns {required}, got {list(df.columns)}"
        )
    df["dataset_count"] = df["dataset_count"].astype(int)
    return df.sort_values("dataset_count", ascending=False).reset_index(drop=True)
