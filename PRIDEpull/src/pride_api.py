"""PRIDE Archive REST API client (vendored patterns from pride_species_visualizer)."""

# Enable postponed evaluation of type annotations for cleaner forward refs
from __future__ import annotations

# Import logging so retry warnings and pagination progress are visible
import logging

# Import time so the retry loop can sleep between backoff attempts
import time

# Import Path for optional cache file locations passed as strings
from pathlib import Path

# Import Any because PRIDE JSON payloads use heterogeneous dict values
from typing import Any

# Import requests for HTTP session management and GET calls
import requests

# Import Response and Session types for function signatures
from requests import Response, Session

# Import config constants so no URLs or retry limits are hardcoded here
from config import config as cfg

# Create a module-level logger named after this file
logger = logging.getLogger(__name__)


class PrideApiError(Exception):
    """Raised when the PRIDE Archive API returns an unrecoverable error."""


def create_session() -> Session:
    """Create a configured HTTP session for PRIDE API requests."""

    # Instantiate a persistent requests Session object for connection reuse
    session = requests.Session()

    # Attach User-Agent and Accept headers expected by the PRIDE Archive API
    session.headers.update(
        {
            "User-Agent": cfg.PRIDE_USER_AGENT,
            "Accept": "application/json",
        }
    )

    # Return the configured session to the caller
    return session


def _request_with_retry(
    session: Session,
    url: str,
    params: dict[str, Any] | None = None,
) -> Response:
    """Perform a GET request with exponential backoff on transient failures."""

    # Track the last exception so we can chain it in the final PrideApiError
    last_error: Exception | None = None

    # Attempt the request up to PRIDE_MAX_RETRIES times
    for attempt in range(cfg.PRIDE_MAX_RETRIES):
        try:
            # Issue a GET request with the configured connect/read timeout tuple
            response = session.get(url, params=params, timeout=cfg.PRIDE_TIMEOUT)

            # Retry on rate-limit and server-side error status codes
            if response.status_code in (429, 500, 502, 503, 504):
                # Compute exponential backoff wait time for this attempt index
                wait = cfg.PRIDE_RETRY_BACKOFF_SECONDS * (2**attempt)

                # Log the HTTP status and planned wait before retrying
                logger.warning(
                    "HTTP %s from %s; retrying in %.1fs (attempt %d/%d)",
                    response.status_code,
                    url,
                    wait,
                    attempt + 1,
                    cfg.PRIDE_MAX_RETRIES,
                )

                # Pause execution before the next retry attempt
                time.sleep(wait)

                # Skip to the next loop iteration without raising
                continue

            # Fail fast on other client errors (4xx except 429 handled above)
            if 400 <= response.status_code < 500:
                response.raise_for_status()

            # Raise for any remaining non-success status codes
            response.raise_for_status()

            # Return the successful response object to the caller
            return response

        except requests.RequestException as exc:
            # Store the exception for the final error message if all retries fail
            last_error = exc

            # Compute backoff wait after a connection-level failure
            wait = cfg.PRIDE_RETRY_BACKOFF_SECONDS * (2**attempt)

            # Log the connection error and planned wait before retrying
            logger.warning(
                "Request failed for %s: %s; retrying in %.1fs (attempt %d/%d)",
                url,
                exc,
                wait,
                attempt + 1,
                cfg.PRIDE_MAX_RETRIES,
            )

            # Pause execution before the next retry attempt
            time.sleep(wait)

    # Build a human-readable message after all retry attempts are exhausted
    message = f"Failed to fetch {url} after {cfg.PRIDE_MAX_RETRIES} attempts"

    # Chain the last caught exception when one exists for easier debugging
    if last_error is not None:
        raise PrideApiError(message) from last_error

    # Raise without a chained cause when no RequestException was captured
    raise PrideApiError(message)


def search_projects_page(
    session: Session,
    page: int,
    page_size: int | None = None,
    *,
    keyword: str | None = None,
    organism_facet: str | None = None,
) -> list[dict[str, Any]]:
    """Fetch one page of project metadata from /search/projects."""

    # Default page size from config when the caller does not override it
    if page_size is None:
        page_size = cfg.PRIDE_PAGE_SIZE

    # Start building query parameters with zero-based page index and page size
    params: dict[str, Any] = {"page": page, "pageSize": page_size}

    # Add free-text keyword to params when provided by the caller
    if keyword:
        params["keyword"] = keyword

    # Add organisms_facet filter when a species facet string is provided
    if organism_facet:
        params["filter"] = f"organisms_facet=={organism_facet}"

    # Compose the full search/projects URL from the configured base URL
    url = f"{cfg.PRIDE_BASE_URL}/search/projects"

    # Execute the GET with retry logic and return the parsed JSON body
    response = _request_with_retry(session, url, params=params)

    # Parse the response body as JSON (expected to be a list of project dicts)
    data = response.json()

    # Validate that the API returned a JSON list as documented
    if not isinstance(data, list):
        raise PrideApiError(
            f"Expected JSON list from search/projects, got {type(data).__name__}"
        )

    # Return the list of project metadata dictionaries for this page
    return data


def fetch_project_files_page(
    session: Session,
    accession: str,
    page: int,
    page_size: int = 100,
) -> list[dict[str, Any]]:
    """Fetch one page of file metadata for a PRIDE project accession."""

    # Build query parameters with zero-based page index and page size
    params: dict[str, Any] = {"page": page, "pageSize": page_size}

    # Compose the per-project files endpoint URL using the PXD accession
    url = f"{cfg.PRIDE_BASE_URL}/projects/{accession}/files"

    # Execute the GET with retry logic
    response = _request_with_retry(session, url, params=params)

    # Parse the response body as JSON
    data = response.json()

    # Validate that the API returned a JSON list
    if not isinstance(data, list):
        raise PrideApiError(
            f"Expected JSON list from projects/{accession}/files, "
            f"got {type(data).__name__}"
        )

    # Return the list of file metadata dictionaries for this page
    return data


def fetch_all_project_files(
    session: Session,
    accession: str,
    page_size: int = 100,
) -> list[dict[str, Any]]:
    """Paginate all file records for a single PRIDE project accession."""

    # Accumulate file dicts from every page into this list
    all_files: list[dict[str, Any]] = []

    # Start pagination at page zero (PRIDE API uses zero-based page index)
    page = 0

    # Loop until an empty page is returned
    while True:
        # Fetch one page of file metadata for this project accession
        batch = fetch_project_files_page(session, accession, page, page_size)

        # Stop pagination when the API returns no more files
        if not batch:
            break

        # Append this page's file records to the accumulated list
        all_files.extend(batch)

        # Advance to the next page index
        page += 1

    # Return the complete file inventory for the project
    return all_files


def extract_download_url(file_record: dict[str, Any]) -> str | None:
    """Pick the best public download URL from a PRIDE file record."""

    # Read the list of publicFileLocations dicts from the file record
    locations = file_record.get("publicFileLocations") or []

    # Iterate preferred protocol names in config order (FTP first)
    for protocol_name in cfg.DOWNLOAD_PROTOCOL_PREFERENCE:
        # Scan each location entry for a matching protocol name
        for location in locations:
            # Compare the location's name field to the preferred protocol
            if location.get("name") == protocol_name:
                # Return the URL/value string for the first matching protocol
                url = location.get("value")
                if url:
                    return str(url)

    # Return None when no preferred download URL could be found
    return None
