#!/usr/bin/env python3
"""Master entry point for UniProt proteome download and KDE visualization."""

# Enable postponed evaluation of annotations used in the CLI helpers
from __future__ import annotations

# Import argparse to expose an optional force-refetch flag
import argparse
# Import logging to configure console progress output for the full pipeline
import logging
# Import sys so the repository root can be prepended to the import path
import sys
# Import Path for repository-root resolution used by import bootstrapping
from pathlib import Path

# Resolve the absolute path of this entry-point script
_SCRIPT_PATH = Path(__file__).resolve()
# Resolve the repository root as the parent directory of this script
_REPO_ROOT = _SCRIPT_PATH.parent
# Prepend the repository root to sys.path so config/ and src/ imports resolve
sys.path.insert(0, str(_REPO_ROOT))

# Import centralized directories and output filenames from config
from config.config import (
    DATA_DIR,
    ENABLE_LOG_PANEL,
    FORCE_REFETCH,
    OUTPUT_DIR,
    PROTEOMES,
    SUMMARY_CSV_NAME,
)
# Import the multi-proteome ingestion pipeline
from src.fetch_data import load_all_proteomes
# Import the single-panel and dual-panel Seaborn KDE plotting routines
from src.plot_distributions import plot_kde_density, plot_kde_linear_vs_log


# Configure root logging to emit INFO messages to the console
logging.basicConfig(
    # Set the minimum severity level shown during a pipeline run
    level=logging.INFO,
    # Format log lines with a timestamp, logger name, level, and message
    format="%(asctime)s | %(name)s | %(levelname)s | %(message)s",
)
# Create a module-level logger for entry-point status messages
logger = logging.getLogger("__main__")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the visualization pipeline."""
    # Create an argument parser describing the master entry-point script
    parser = argparse.ArgumentParser(
        # Provide a short description of the pipeline purpose
        description="Download UniProt proteomes and plot ORF-length KDE densities.",
    )
    # Add an optional flag that forces fresh UniProt downloads despite cache
    parser.add_argument(
        # Expose the long-form force-refetch switch
        "--force-refetch",
        # Store True when the flag is present on the command line
        action="store_true",
        # Default to the centralized FORCE_REFETCH configuration value
        default=FORCE_REFETCH,
        # Document that existing FASTA caches will be ignored when set
        help="Ignore local FASTA caches and re-download all proteomes from UniProt.",
    )
    # Parse argv and return the populated namespace
    return parser.parse_args()


def write_summary_csv(df, output_dir: Path) -> Path:
    """Write a per-species ORF-length summary CSV for auditability."""
    # Aggregate descriptive statistics of coding bp lengths for every species
    summary = (
        # Group on the species display name
        df.groupby("species", observed=False)["bp_length"]
        # Compute count plus common distribution summaries for auditing
        .agg(count="size", mean_bp="mean", median_bp="median", min_bp="min", max_bp="max")
        # Convert the grouped Series index into ordinary columns
        .reset_index()
    )
    # Build the destination path for the audit summary CSV
    summary_path = output_dir / SUMMARY_CSV_NAME
    # Persist the summary table without writing the pandas index
    summary.to_csv(summary_path, index=False)
    # Log the summary CSV destination for the operator
    logger.info("Wrote length summary to %s", summary_path)
    # Return the written summary path for the final artifact listing
    return summary_path


def main() -> int:
    """Run download, summary export, and KDE publication plotting end to end."""
    # Parse optional CLI flags such as --force-refetch
    args = parse_args()
    # Ensure the local FASTA cache directory exists before downloads begin
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    # Ensure the publication figure output directory exists before plotting
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    # Download, parse, and convert AA lengths to coding bp for all proteomes
    df = load_all_proteomes(PROTEOMES, force=args.force_refetch)
    # Write an audit CSV summarizing ORF-length statistics by species
    summary_path = write_summary_csv(df, OUTPUT_DIR)
    # Render the single-panel linear multi-species KDE and export PNG plus SVG
    kde_paths = plot_kde_density(df)
    # Initialize the artifact list with the summary CSV and single-panel figures
    artifacts = [summary_path, *kde_paths]
    # Generate the dual-panel linear-versus-log figure when enabled in config
    if ENABLE_LOG_PANEL:
        # Render Panels A/B comparative KDE and export PNG plus SVG
        dual_paths = plot_kde_linear_vs_log(df)
        # Append the dual-panel artifact paths to the console summary list
        artifacts.extend(dual_paths)
    # Log a header before the final artifact listing
    logger.info("Pipeline complete. Written artifacts:")
    # Print each artifact path so the operator can locate outputs quickly
    for artifact in artifacts:
        # Emit one log line per written file path
        logger.info("  %s", artifact)
    # Return a successful process exit code
    return 0


# Execute the master pipeline only when this module is run as a script
if __name__ == "__main__":
    # Propagate the main() return code as the process exit status
    raise SystemExit(main())
