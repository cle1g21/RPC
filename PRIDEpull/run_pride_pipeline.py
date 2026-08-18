#!/usr/bin/env python3
"""Master entry point for the PRIDEpull immunopeptidomics pipeline."""

# Enable postponed evaluation of type annotations
from __future__ import annotations

# Import argparse for command-line flag parsing
import argparse

# Import logging for console progress output
import logging

# Import sys so the repo root can be prepended to the module search path
import sys

# Import Path to resolve the PRIDEpull repository root directory
from pathlib import Path

# Resolve the directory containing this script (PRIDEpull repo root)
_ROOT = Path(__file__).resolve().parent

# Insert repo root at the front of sys.path for config and src imports
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

# Import centralized configuration defaults
from config import config as cfg

# Import harvester for harvest-only mode
from src.pride_harvester import harvest_immunopeptidomics_projects

# Import pipeline orchestrator
from src.pipeline_runner import run_pipeline


def parse_args() -> argparse.Namespace:
    """Define and parse command-line arguments for run_pride_pipeline.py."""

    # Create the argument parser with a short pipeline description
    parser = argparse.ArgumentParser(
        description="PRIDEpull: harvest human immunopeptidomics data and run InstaNovo.",
    )

    # Flag to run only the PRIDE harvester and write cache, then exit
    parser.add_argument(
        "--harvest-only",
        action="store_true",
        help="Query PRIDE and write harvest cache without processing files.",
    )

    # Flag to ignore existing harvest cache and re-query PRIDE
    parser.add_argument(
        "--refresh-harvest",
        action="store_true",
        help="Ignore harvest cache and re-query the PRIDE API.",
    )

    # Optional single PXD accession to process
    parser.add_argument(
        "--accession",
        type=str,
        default=None,
        help="Process only this PRIDE project accession (e.g. PXD077095).",
    )

    # Optional cap on number of projects to process in one run
    parser.add_argument(
        "--max-projects",
        type=int,
        default=None,
        help="Maximum number of projects to process (default: config MAX_PROJECTS).",
    )

    # Flag to disable subset mode and run full InstaNovo prediction
    parser.add_argument(
        "--full-run",
        action="store_true",
        help="Set RUN_SUBSET_ONLY=False for full-scale InstaNovo runs.",
    )

    # Flag to skip download and use files already present in data/
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Skip PRIDE download; use files already in PRIDEpull/data/.",
    )

    # Flag to disable delete-on-success cleanup (debugging)
    parser.add_argument(
        "--no-cleanup",
        action="store_true",
        help="Do not delete intermediate files in data/ after success.",
    )

    # Verbose logging flag
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable DEBUG level logging.",
    )

    # Parse argv and return the resulting namespace
    return parser.parse_args()


def configure_logging(verbose: bool) -> None:
    """Configure root logger level and format."""

    # Choose DEBUG when verbose, otherwise INFO
    level = logging.DEBUG if verbose else logging.INFO

    # Apply basicConfig with timestamp, level, and message format
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


def main() -> int:
    """Run harvest-only or full pipeline based on CLI flags."""

    # Parse command-line arguments
    args = parse_args()

    # Configure logging from verbose flag
    configure_logging(args.verbose)

    # Resolve run_subset_only: False when --full-run, else config default
    run_subset_only = False if args.full_run else cfg.RUN_SUBSET_ONLY

    # Resolve delete_on_success: False when --no-cleanup, else config default
    delete_on_success = False if args.no_cleanup else cfg.DELETE_ON_SUCCESS

    # Harvest-only mode: query PRIDE, write cache, exit
    if args.harvest_only:
        manifest = harvest_immunopeptidomics_projects(
            refresh=args.refresh_harvest,
            max_projects=args.max_projects,
        )
        print(f"Harvest complete: {len(manifest)} projects written to cache.")
        return 0

    # Full pipeline mode: harvest (if needed) + download + convert + predict
    summary = run_pipeline(
        refresh_harvest=args.refresh_harvest,
        accession_filter=args.accession,
        max_projects=args.max_projects,
        skip_download=args.skip_download,
        delete_on_success=delete_on_success,
        run_subset_only=run_subset_only,
    )

    # Print summary to stdout for operator visibility
    print(
        f"Pipeline complete: {summary['completed']} succeeded, "
        f"{summary['failed']} failed, {summary['skipped']} skipped."
    )

    # Return non-zero exit code when any file processing failed
    return 1 if summary["failed"] > 0 else 0


# Execute main() only when this file is run as a script
if __name__ == "__main__":
    sys.exit(main())
