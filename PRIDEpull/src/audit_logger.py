"""Append processing history blocks and compute InstaNovo prediction statistics."""

# Enable postponed evaluation of type annotations
from __future__ import annotations

# Import logging for ledger write confirmation
import logging

# Import os for file deletion during cleanup
import os

# Import shutil for recursive directory removal (.d folders)
import shutil

# Import datetime for UTC timestamps in ledger blocks
from datetime import datetime, timezone

# Import Path for ledger and CSV paths
from pathlib import Path

# Import Any for flexible run metadata dicts
from typing import Any

# Import pandas for prediction CSV statistics
import pandas as pd

# Import config for ledger path and column names
from config import config as cfg

# Create a module-level logger
logger = logging.getLogger(__name__)

# Separator line used between ledger blocks
_BLOCK_SEPARATOR = "=" * 80

# Sub-separator within a ledger block
_SUB_SEPARATOR = "-" * 80


def _format_timestamp() -> str:
    """Return current UTC timestamp in ISO 8601 format."""

    # Get timezone-aware UTC datetime and format as ISO string
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def compute_prediction_stats(csv_path: Path | str) -> dict[str, Any]:
    """
    Parse InstaNovo output CSV and return spectrum/peptide/confidence metrics.
    """

    # Normalize csv_path to Path
    path = Path(csv_path)

    # Initialize default stats when file is missing or empty
    stats: dict[str, Any] = {
        "spectra_processed": 0,
        "peptides_sequenced": 0,
        "mean_prediction_log_prob": None,
        "median_prediction_log_prob": None,
        "output_csv_size_bytes": 0,
    }

    # Return defaults when the CSV does not exist
    if not path.is_file():
        return stats

    # Record output file size in bytes
    stats["output_csv_size_bytes"] = path.stat().st_size

    # Read the predictions CSV into a DataFrame
    df = pd.read_csv(path)

    # Count total rows as spectra processed
    stats["spectra_processed"] = len(df)

    # Count rows with non-empty predictions column as peptides sequenced
    if "predictions" in df.columns:
        non_empty = df["predictions"].astype(str).str.strip()
        stats["peptides_sequenced"] = int((non_empty != "").sum())

    # Compute mean/median log probability from InstaNovo confidence column
    prob_col = None
    for candidate in (
        "instanovo_prediction_log_probability",
        "log_probs",
    ):
        if candidate in df.columns:
            prob_col = candidate
            break

    # Calculate statistics when a probability column exists and has numeric values
    if prob_col is not None:
        numeric = pd.to_numeric(df[prob_col], errors="coerce").dropna()
        if len(numeric) > 0:
            stats["mean_prediction_log_prob"] = float(numeric.mean())
            stats["median_prediction_log_prob"] = float(numeric.median())

    # Return the computed statistics dictionary
    return stats


def cleanup_data_files(file_paths: list[Path | str]) -> list[str]:
    """
    Delete intermediate files in data/ using os.remove or shutil.rmtree.

    Returns a list of human-readable deletion confirmation strings.
    """

    # Accumulate confirmation messages for each deleted path
    confirmations: list[str] = []

    # Iterate each path requested for deletion
    for file_path in file_paths:
        # Normalize to Path
        path = Path(file_path)

        # Skip paths that do not exist on disk
        if not path.exists():
            confirmations.append(f"{path} [NOT FOUND - skipped]")
            continue

        try:
            # Remove directories recursively (Bruker .d folders)
            if path.is_dir():
                shutil.rmtree(path)
                confirmations.append(f"{path} [DELETED dir]")
            else:
                # Remove a single file with os.remove
                os.remove(path)
                confirmations.append(f"{path} [DELETED]")

            # Log each successful deletion
            logger.info("Deleted: %s", path)

        except OSError as exc:
            # Record failure without aborting cleanup of remaining paths
            confirmations.append(f"{path} [DELETE FAILED: {exc}]")
            logger.warning("Failed to delete %s: %s", path, exc)

    # Return the list of confirmation strings for the ledger
    return confirmations


def append_ledger_entry(entry: dict[str, Any]) -> None:
    """
    Append a structured processing block to processing_history_ledger.txt.

    Required entry keys:
        pride_accession, original_filename, processing_route, status
    Optional keys populate sizes, stats, cleanup, run_mode, etc.
    """

    # Resolve ledger path from config
    ledger_path = Path(cfg.PROCESSING_LEDGER_PATH)

    # Ensure parent directory exists (ledger lives at repo root)
    ledger_path.parent.mkdir(parents=True, exist_ok=True)

    # Read mandatory identifiers from entry dict
    accession = entry.get("pride_accession", "UNKNOWN")
    original_filename = entry.get("original_filename", "UNKNOWN")

    # Build ledger block lines in the specified format
    lines: list[str] = [
        "",
        _BLOCK_SEPARATOR,
        f"PRIDE_ACCESSION: {accession}",
        f"ORIGINAL_PRIDE_FILENAME: {original_filename}",
        f"TIMESTAMP: {entry.get('timestamp') or _format_timestamp()}",
        f"PROCESSING_ROUTE: {entry.get('processing_route', 'N/A')}",
        f"RUN_MODE: {entry.get('run_mode', 'N/A')}",
        _SUB_SEPARATOR,
        f"DOWNLOAD_SIZE_BYTES: {entry.get('download_size_bytes', 'N/A')}",
        f"CONVERTED_MGF_SIZE_BYTES: {entry.get('converted_mgf_size_bytes', 'N/A')}",
        f"SUBSET_MGF_SIZE_BYTES: {entry.get('subset_mgf_size_bytes', 'N/A')}",
        f"OUTPUT_CSV: {entry.get('output_csv', 'N/A')}",
        f"OUTPUT_CSV_SIZE_BYTES: {entry.get('output_csv_size_bytes', 'N/A')}",
        _SUB_SEPARATOR,
        f"SPECTRA_PROCESSED: {entry.get('spectra_processed', 'N/A')}",
        f"PEPTIDES_SEQUENCED: {entry.get('peptides_sequenced', 'N/A')}",
        f"MEAN_PREDICTION_LOG_PROB: {entry.get('mean_prediction_log_prob', 'N/A')}",
        f"MEDIAN_PREDICTION_LOG_PROB: {entry.get('median_prediction_log_prob', 'N/A')}",
        f"INSTANOVO_BATCH_SIZE: {entry.get('instanovo_batch_size', cfg.INSTANOVO_BATCH_SIZE)}",
        f"INSTANOVO_NUM_WORKERS: {entry.get('instanovo_num_workers', cfg.INSTANOVO_NUM_WORKERS)}",
        _SUB_SEPARATOR,
        "CLEANUP_DELETED_FILES:",
    ]

    # Append each cleanup confirmation line indented
    cleanup_lines = entry.get("cleanup_confirmations") or []
    if cleanup_lines:
        for conf in cleanup_lines:
            lines.append(f"  - {conf}")
    else:
        lines.append("  - (none)")

    # Append cleanup status and overall run status
    lines.extend(
        [
            f"CLEANUP_STATUS: {entry.get('cleanup_status', 'N/A')}",
            f"STATUS: {entry.get('status', 'UNKNOWN')}",
        ]
    )

    # Append error message line when the run failed
    if entry.get("error_message"):
        lines.append(f"ERROR: {entry['error_message']}")

    # Close the block with the separator
    lines.append(_BLOCK_SEPARATOR)

    # Join lines with newlines for append write
    block_text = "\n".join(lines) + "\n"

    # Append the block to the ledger file (create if missing)
    with ledger_path.open("a", encoding="utf-8") as handle:
        handle.write(block_text)

    # Log confirmation of ledger write
    logger.info("Appended ledger entry for %s / %s", accession, original_filename)
