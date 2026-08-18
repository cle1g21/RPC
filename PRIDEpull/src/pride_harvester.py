"""PRIDE immunopeptidomics project harvester and multi-format file inventory."""

# Enable postponed evaluation of type annotations
from __future__ import annotations

# Import json so harvest results can be cached to disk
import json

# Import logging for harvest progress messages
import logging

# Import re for case-insensitive keyword boundary checks when needed
import re

# Import Path for cache file read/write
from pathlib import Path

# Import Any for PRIDE JSON dict typing
from typing import Any

# Import config constants for keywords, species, and cache paths
from config import config as cfg

# Import PRIDE API helpers for search, file listing, and download URL extraction
from src.pride_api import (
    PrideApiError,
    create_session,
    extract_download_url,
    fetch_all_project_files,
    search_projects_page,
)

# Create a module-level logger
logger = logging.getLogger(__name__)


def _build_searchable_text(project: dict[str, Any]) -> str:
    """Concatenate project metadata fields into one lowercase searchable string."""

    # Collect text fragments from multiple metadata fields into a list
    parts: list[str] = []

    # Add the project title when present
    if project.get("title"):
        parts.append(str(project["title"]))

    # Add the project description / abstract when present
    if project.get("projectDescription"):
        parts.append(str(project["projectDescription"]))

    # Add each keyword tag when the keywords list is non-empty
    for kw in project.get("keywords") or []:
        parts.append(str(kw))

    # Add each project tag when the projectTags list is non-empty
    for tag in project.get("projectTags") or []:
        parts.append(str(tag))

    # Add each experiment type label when present
    for exp_type in project.get("experimentTypes") or []:
        parts.append(str(exp_type))

    # Join all parts with spaces and normalize to lowercase for substring matching
    return " ".join(parts).lower()


def match_immuno_keywords(project: dict[str, Any]) -> list[str]:
    """Return IMMUNO_KEYWORDS terms found in project metadata (case-insensitive)."""

    # Build the combined searchable text blob for this project
    searchable = _build_searchable_text(project)

    # Collect keyword terms that appear as substrings in the searchable text
    matched: list[str] = []

    # Check each configured immunopeptidomics keyword against the searchable text
    for keyword in cfg.IMMUNO_KEYWORDS:
        # Normalize the keyword to lowercase for consistent matching
        kw_lower = keyword.lower()

        # Record the keyword when it appears anywhere in the metadata text
        if kw_lower in searchable:
            matched.append(keyword)

    # Return the list of matched keyword strings (may be empty)
    return matched


def is_human_project(project: dict[str, Any]) -> bool:
    """Return True when the project organisms list includes Homo sapiens."""

    # Read the organisms list from project metadata (may be empty)
    organisms = project.get("organisms") or []

    # Check each organism label for the configured human substring
    for organism in organisms:
        # Accept the project when any organism string contains the human label
        if cfg.HUMAN_ORGANISM_LABEL.lower() in str(organism).lower():
            return True

    # Reject the project when no human organism label was found
    return False


def _classify_file_extension(file_name: str) -> str | None:
    """Map a PRIDE fileName to native_mgf, thermo_raw, bruker_d, or None."""

    # Normalize the file name to lowercase for extension comparison
    lower_name = file_name.lower()

    # Skip mzid-prefixed MGF duplicates (e.g. P1028.mzid_P1028.MGF)
    if ".mzid_" in lower_name and lower_name.endswith(".mgf"):
        return None

    # Classify native MGF peak list files
    if lower_name.endswith(cfg.MGF_EXTENSIONS):
        return "native_mgf"

    # Classify Thermo .raw vendor files
    if lower_name.endswith(cfg.RAW_EXTENSIONS):
        return "thermo_raw"

    # Classify Bruker .d folders and common archive wrappers
    for ext in cfg.BRUKER_D_EXTENSIONS:
        if lower_name.endswith(ext):
            return "bruker_d"

    # Return None for file types outside mgf/raw/d scope
    return None


def inventory_project_files(
    session: Any,
    accession: str,
) -> dict[str, list[dict[str, Any]]]:
    """Query PRIDE and classify .mgf, .raw, and .d files for one project."""

    # Initialize empty buckets for each supported file category
    buckets: dict[str, list[dict[str, Any]]] = {
        "native_mgf": [],
        "thermo_raw": [],
        "bruker_d": [],
    }

    # Fetch the complete paginated file list for this accession
    raw_files = fetch_all_project_files(session, accession)

    # Iterate every file record returned by the PRIDE files API
    for file_record in raw_files:
        # Read the original PRIDE file name from the record
        file_name = file_record.get("fileName") or ""

        # Skip records with no file name
        if not file_name:
            continue

        # Determine the bucket key from the file extension
        category = _classify_file_extension(file_name)

        # Skip files that do not match mgf, raw, or d extensions
        if category is None:
            continue

        # Resolve the best FTP/HTTPS download URL for this file
        download_url = extract_download_url(file_record)

        # Skip files that have no public download location
        if not download_url:
            continue

        # Build a normalized inventory entry for downstream download/routing
        entry = {
            "fileName": file_name,
            "downloadUrl": download_url,
            "fileSizeBytes": file_record.get("fileSizeBytes") or 0,
            "fileCategory": (file_record.get("fileCategory") or {}).get("value", ""),
        }

        # Append the entry to the appropriate category bucket
        buckets[category].append(entry)

    # Return the classified file buckets for this project
    return buckets


def select_project_files(
    buckets: dict[str, list[dict[str, Any]]],
    max_files: int | None = None,
) -> tuple[str, list[dict[str, Any]]]:
    """
    Select up to max_files candidates, preferring largest native MGF first.

    Returns:
        (routing, selected_files) where routing is native_mgf_direct,
        vendor_raw_msconvert, or skip.
    """

    # Default to config MAX_FILES_PER_PROJECT when max_files is not specified
    if max_files is None:
        max_files = cfg.MAX_FILES_PER_PROJECT

    # Prefer native MGF when any .mgf files are available on PRIDE
    if buckets.get("native_mgf"):
        # Sort native MGF entries by fileSizeBytes descending (largest first)
        sorted_mgf = sorted(
            buckets["native_mgf"],
            key=lambda f: f.get("fileSizeBytes") or 0,
            reverse=True,
        )

        # Return the top N native MGF files with direct routing label
        return "native_mgf_direct", sorted_mgf[:max_files]

    # Combine Thermo RAW and Bruker .d into one vendor raw candidate list
    vendor_files = (buckets.get("thermo_raw") or []) + (buckets.get("bruker_d") or [])

    # Fall back to vendor conversion when no native MGF exists
    if vendor_files:
        # Sort vendor files by fileSizeBytes descending (largest first)
        sorted_vendor = sorted(
            vendor_files,
            key=lambda f: f.get("fileSizeBytes") or 0,
            reverse=True,
        )

        # Return the top N vendor files with msconvert routing label
        return "vendor_raw_msconvert", sorted_vendor[:max_files]

    # Return skip when no supported file types were found
    return "skip", []


def harvest_immunopeptidomics_projects(
    *,
    refresh: bool = False,
    max_projects: int | None = None,
) -> list[dict[str, Any]]:
    """
    Search PRIDE for human immunopeptidomics projects and inventory their files.

    Returns a list of manifest dicts ready for pipeline_runner consumption.
    """

    # Resolve the harvest cache path from config
    cache_path = Path(cfg.HARVEST_CACHE_PATH)

    # Load cached manifest when refresh is False and the cache file exists
    if not refresh and cache_path.is_file():
        logger.info("Loading harvest cache from %s", cache_path)

        # Read and parse the JSON cache file
        with cache_path.open(encoding="utf-8") as handle:
            cached = json.load(handle)

        # Return cached manifest when it is a non-empty list
        if isinstance(cached, list) and cached:
            return cached

    # Create a persistent HTTP session for all PRIDE API calls in this harvest
    session = create_session()

    # Accumulate unique projects keyed by PXD accession
    projects_by_accession: dict[str, dict[str, Any]] = {}

    # Iterate each server-side keyword to cast a wide net with fewer missed studies
    for server_keyword in cfg.PRIDE_SERVER_KEYWORDS:
        logger.info("Searching PRIDE with keyword=%r", server_keyword)

        # Start pagination at page zero for this keyword query
        page = 0

        # Paginate until an empty page is returned
        while True:
            # Fetch one page of human-filtered search results for this keyword
            batch = search_projects_page(
                session,
                page,
                keyword=server_keyword,
                organism_facet=cfg.HUMAN_ORGANISM_FACET,
            )

            # Stop pagination when no more projects are returned
            if not batch:
                break

            # Process each project dict in this page
            for project in batch:
                # Read the PXD accession identifier
                accession = project.get("accession")

                # Skip records without an accession field
                if not accession:
                    continue

                # Skip non-human projects even if the server filter misfired
                if not is_human_project(project):
                    continue

                # Apply client-side immunopeptidomics keyword matrix filter
                matched_keywords = match_immuno_keywords(project)

                # Skip projects that match no immunopeptidomics keywords
                if not matched_keywords:
                    continue

                # Store or update the project in the deduplicated accession map
                projects_by_accession[accession] = project

            # Advance to the next page index
            page += 1

    logger.info(
        "Found %d unique human immunopeptidomics candidate projects",
        len(projects_by_accession),
    )

    # Build the final manifest list with file inventory per project
    manifest: list[dict[str, Any]] = []

    # Iterate each deduplicated project in arbitrary dict order
    for accession, project in projects_by_accession.items():
        # Apply optional max_projects cap during manifest construction
        if max_projects is not None and len(manifest) >= max_projects:
            break

        logger.info("Inventorying files for %s", accession)

        try:
            # Query and classify all mgf/raw/d files for this accession
            buckets = inventory_project_files(session, accession)

            # Select up to MAX_FILES_PER_PROJECT files with routing decision
            routing, selected_files = select_project_files(buckets)

            # Skip projects with no processable files
            if routing == "skip" or not selected_files:
                logger.info("Skipping %s: no .mgf, .raw, or .d files found", accession)
                continue

            # Append a manifest entry for the pipeline runner
            manifest.append(
                {
                    "accession": accession,
                    "title": project.get("title") or "",
                    "matched_keywords": match_immuno_keywords(project),
                    "routing": routing,
                    "selected_files": selected_files,
                    "status": "pending",
                }
            )

        except PrideApiError as exc:
            # Log file inventory failures without aborting the entire harvest
            logger.warning("Failed to inventory %s: %s", accession, exc)

    # Ensure the cache directory exists before writing
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    # Write the manifest to the harvest cache JSON file
    with cache_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    logger.info("Wrote harvest cache with %d projects to %s", len(manifest), cache_path)

    # Return the completed manifest list
    return manifest


def load_harvest_cache() -> list[dict[str, Any]]:
    """Load a previously saved harvest manifest from disk."""

    # Resolve the cache path from config
    cache_path = Path(cfg.HARVEST_CACHE_PATH)

    # Raise when the cache file does not exist
    if not cache_path.is_file():
        raise FileNotFoundError(f"Harvest cache not found: {cache_path}")

    # Read and parse the JSON cache
    with cache_path.open(encoding="utf-8") as handle:
        data = json.load(handle)

    # Validate that the cache contains a list
    if not isinstance(data, list):
        raise ValueError(f"Expected list in harvest cache, got {type(data).__name__}")

    # Return the manifest list
    return data
