"""End-to-end orchestration: download, convert, predict, cleanup, ledger."""

# Enable postponed evaluation of type annotations
from __future__ import annotations

# Import logging for pipeline progress
import logging

# Import sys so the repo root can be added to the module search path
import sys

# Import Path for filesystem operations
from pathlib import Path

# Import Any for manifest and summary dict typing
from typing import Any

# Resolve PRIDEpull repo root (parent of src/)
_REPO_ROOT = Path(__file__).resolve().parent.parent

# Insert repo root at front of sys.path for config and src imports
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

# Import centralized configuration
from config import config as cfg

# Import audit ledger helpers
from src.audit_logger import (
    append_ledger_entry,
    cleanup_data_files,
    compute_prediction_stats,
)

# Import msconvert wrapper and MGF subset slicer
from src.file_converter import convert_raw_to_mgf, slice_mgf

# Import chunked PRIDE file downloader
from src.file_downloader import download_pride_file, project_download_path

# Import InstaNovo prediction launcher
from src.instanovo_predictor import run_prediction

# Import harvester for manifest loading
from src.pride_harvester import harvest_immunopeptidomics_projects, load_harvest_cache

# Create module-level logger
logger = logging.getLogger(__name__)


def _file_size(path: Path | None) -> int | None:
    """Return file size in bytes or None when path is missing."""

    # Return None when no path was provided
    if path is None:
        return None

    # Return None when the path does not exist
    if not path.is_file():
        return None

    # Return stat size for existing files
    return path.stat().st_size


def _route_label(routing: str) -> str:
    """Map internal routing key to human-readable ledger PROCESSING_ROUTE string."""

    # Map native MGF direct routing
    if routing == "native_mgf_direct":
        return "Native MGF streamed direct"

    # Map vendor raw msconvert routing
    if routing == "vendor_raw_msconvert":
        return "Vendor RAW converted via msconvert"

    # Default label for unknown routing keys
    return routing


def process_single_file(
    accession: str,
    file_entry: dict[str, Any],
    routing: str,
    *,
    skip_download: bool = False,
    delete_on_success: bool | None = None,
    run_subset_only: bool | None = None,
) -> dict[str, Any]:
    """
    Process one PRIDE file through download → convert → subset → predict → cleanup.

    Returns a summary dict suitable for audit_logger.append_ledger_entry.
    """

    # Read original PRIDE file name from inventory entry
    original_filename = file_entry["fileName"]

    # Resolve delete_on_success from config when not overridden
    if delete_on_success is None:
        delete_on_success = cfg.DELETE_ON_SUCCESS

    # Resolve run_subset_only from config when not overridden
    if run_subset_only is None:
        run_subset_only = cfg.RUN_SUBSET_ONLY

    # Build standard landing path under data/{accession}/
    dest_path = project_download_path(accession, original_filename)

    # Track paths to delete after successful prediction
    cleanup_paths: list[Path] = []

    # Track intermediate paths for ledger size reporting
    downloaded_path: Path | None = None
    converted_mgf_path: Path | None = None
    subset_mgf_path: Path | None = None
    active_mgf_path: Path | None = None
    output_csv: Path | None = None

    # Initialize ledger entry with mandatory identifiers first
    ledger_entry: dict[str, Any] = {
        "pride_accession": accession,
        "original_filename": original_filename,
        "processing_route": _route_label(routing),
        "run_mode": f"subset ({cfg.SUBSET_SPECTRUM_COUNT} spectra)"
        if run_subset_only
        else "full",
        "download_size_bytes": file_entry.get("fileSizeBytes"),
        "status": "FAILED",
    }

    try:
        # Step 1: Download (unless skip_download and file already in data/)
        if skip_download and dest_path.is_file():
            logger.info("Using existing file: %s", dest_path)
            downloaded_path = dest_path
        else:
            logger.info("Downloading %s → %s", original_filename, dest_path)
            downloaded_path = download_pride_file(
                file_entry["downloadUrl"],
                dest_path,
                expected_size=file_entry.get("fileSizeBytes"),
            )

        # Register downloaded file for post-success cleanup
        cleanup_paths.append(downloaded_path)

        # Step 2: Route to native MGF or msconvert conversion
        if routing == "native_mgf_direct":
            # Use downloaded MGF directly without conversion
            active_mgf_path = downloaded_path
        else:
            # Convert vendor raw/.d to MGF in the same accession data folder
            logger.info("Converting %s via msconvert", downloaded_path)
            converted_mgf_path = convert_raw_to_mgf(
                downloaded_path,
                output_dir=downloaded_path.parent,
            )
            active_mgf_path = converted_mgf_path
            cleanup_paths.append(converted_mgf_path)

        # Record converted MGF size in ledger when conversion occurred
        ledger_entry["converted_mgf_size_bytes"] = _file_size(converted_mgf_path)

        # Step 3: Optional subset slice when RUN_SUBSET_ONLY is enabled
        if run_subset_only and active_mgf_path is not None:
            stem = active_mgf_path.stem
            subset_name = f"{stem}_subset{cfg.SUBSET_SPECTRUM_COUNT}.mgf"
            subset_mgf_path = active_mgf_path.parent / subset_name
            logger.info("Slicing subset MGF: %s", subset_mgf_path)
            active_mgf_path = slice_mgf(
                active_mgf_path,
                subset_mgf_path,
                max_spectra=cfg.SUBSET_SPECTRUM_COUNT,
            )
            cleanup_paths.append(subset_mgf_path)

        # Record subset MGF size when a subset was created
        ledger_entry["subset_mgf_size_bytes"] = _file_size(subset_mgf_path)

        # Step 4: Run InstaNovo prediction on the active MGF path
        if active_mgf_path is None:
            raise RuntimeError("No MGF path available for InstaNovo")

        logger.info("Running InstaNovo on %s", active_mgf_path)
        output_csv = run_prediction(
            accession,
            active_mgf_path,
            run_subset_only=run_subset_only,
        )

        # Step 5: Compute prediction statistics from output CSV
        stats = compute_prediction_stats(output_csv)

        # Merge stats into ledger entry
        ledger_entry.update(stats)
        ledger_entry["output_csv"] = str(output_csv)
        ledger_entry["instanovo_batch_size"] = cfg.INSTANOVO_BATCH_SIZE
        ledger_entry["instanovo_num_workers"] = cfg.INSTANOVO_NUM_WORKERS
        ledger_entry["status"] = "COMPLETED"

        # Step 6: Delete-on-success cleanup of data/ intermediates
        if delete_on_success:
            confirmations = cleanup_data_files(cleanup_paths)
            ledger_entry["cleanup_confirmations"] = confirmations
            ledger_entry["cleanup_status"] = "SUCCESS"
        else:
            ledger_entry["cleanup_confirmations"] = ["(cleanup disabled)"]
            ledger_entry["cleanup_status"] = "SKIPPED"

        # Append completed ledger block
        append_ledger_entry(ledger_entry)

        # Return summary for caller
        return ledger_entry

    except Exception as exc:
        # Record failure in ledger without deleting data/ files
        ledger_entry["error_message"] = str(exc)
        ledger_entry["cleanup_confirmations"] = ["(skipped on failure)"]
        ledger_entry["cleanup_status"] = "SKIPPED"
        append_ledger_entry(ledger_entry)
        logger.exception("Failed processing %s / %s", accession, original_filename)
        raise


def run_pipeline(
    *,
    manifest: list[dict[str, Any]] | None = None,
    refresh_harvest: bool = False,
    accession_filter: str | None = None,
    max_projects: int | None = None,
    skip_download: bool = False,
    delete_on_success: bool | None = None,
    run_subset_only: bool | None = None,
) -> dict[str, Any]:
    """
    Run the full PRIDEpull pipeline over a harvest manifest.

    Returns a summary dict with completed and failed counts.
    """

    # Ensure data landing directory exists
    cfg.DATA_LANDING_DIR.mkdir(parents=True, exist_ok=True)

    # Load or harvest manifest when not provided by caller
    if manifest is None:
        if refresh_harvest:
            manifest = harvest_immunopeptidomics_projects(
                refresh=True,
                max_projects=max_projects or cfg.MAX_PROJECTS,
            )
        else:
            try:
                manifest = load_harvest_cache()
            except FileNotFoundError:
                manifest = harvest_immunopeptidomics_projects(
                    max_projects=max_projects or cfg.MAX_PROJECTS,
                )

    # Filter to a single accession when requested via CLI
    if accession_filter:
        manifest = [m for m in manifest if m.get("accession") == accession_filter]

    # Apply max_projects cap when specified
    if max_projects is not None:
        manifest = manifest[:max_projects]

    # Initialize run summary counters
    summary: dict[str, Any] = {
        "total_projects": len(manifest),
        "completed": 0,
        "failed": 0,
        "skipped": 0,
        "errors": [],
    }

    # Iterate each project manifest entry
    for project_entry in manifest:
        accession = project_entry.get("accession", "UNKNOWN")
        routing = project_entry.get("routing", "skip")
        selected_files = project_entry.get("selected_files") or []

        # Skip entries with no selected files
        if not selected_files or routing == "skip":
            summary["skipped"] += 1
            continue

        # Process each selected file for this project (up to MAX_FILES_PER_PROJECT)
        for file_entry in selected_files:
            try:
                process_single_file(
                    accession,
                    file_entry,
                    routing,
                    skip_download=skip_download,
                    delete_on_success=delete_on_success,
                    run_subset_only=run_subset_only,
                )
                summary["completed"] += 1
            except Exception as exc:
                summary["failed"] += 1
                summary["errors"].append(
                    {"accession": accession, "file": file_entry.get("fileName"), "error": str(exc)}
                )

    # Log final summary
    logger.info(
        "Pipeline finished: %d completed, %d failed, %d skipped",
        summary["completed"],
        summary["failed"],
        summary["skipped"],
    )

    return summary
