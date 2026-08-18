#!/usr/bin/env python3
"""CLI entry point for PRIDE species dataset counting and visualization."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from src.api_client import (
    fetch_all_species_counts,
    load_counts_cache,
    save_counts_cache,
)
from src.visualizer import plot_species_distribution

DEFAULT_OUTPUT = "output/pride_species_counts.png"
DEFAULT_CACHE = "output/species_counts.csv"
DEFAULT_TOP = 30
DEFAULT_PAGE_SIZE = 100


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the PRIDE species visualizer.

    Args:
        argv: Optional argument list (defaults to ``sys.argv[1:]``).

    Returns:
        Parsed namespace with pipeline configuration.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Fetch PRIDE Archive dataset counts per species and generate "
            "a publication-ready vertical log-scale bar chart."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(DEFAULT_OUTPUT),
        help=f"Output figure path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=DEFAULT_TOP,
        help=f"Number of top species to plot (default: {DEFAULT_TOP})",
    )
    parser.add_argument(
        "--plot-all",
        action="store_true",
        help="Plot every species (produces a very wide figure)",
    )
    parser.add_argument(
        "--page-size",
        type=int,
        default=DEFAULT_PAGE_SIZE,
        help=f"API page size for project pagination (default: {DEFAULT_PAGE_SIZE})",
    )
    parser.add_argument(
        "--cache",
        type=Path,
        default=Path(DEFAULT_CACHE),
        help=f"CSV cache path for aggregated counts (default: {DEFAULT_CACHE})",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Ignore existing cache and re-fetch from the PRIDE API",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable pagination progress logging",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable debug logging",
    )
    return parser.parse_args(argv)


def run_pipeline(args: argparse.Namespace) -> int:
    """Execute fetch, cache, and plot steps.

    Args:
        args: Parsed CLI arguments from ``parse_args``.

    Returns:
        Process exit code (0 on success).
    """
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    log = logging.getLogger(__name__)

    cache_path: Path = args.cache
    use_cache = cache_path.is_file() and not args.refresh

    if use_cache:
        log.info("Loading species counts from cache: %s", cache_path)
        counts_df = load_counts_cache(cache_path)
    else:
        log.info("Fetching species counts from PRIDE API (page_size=%d)", args.page_size)
        counts_df = fetch_all_species_counts(
            page_size=args.page_size,
            show_progress=not args.no_progress,
        )
        save_counts_cache(counts_df, cache_path)

    total_datasets = int(counts_df["dataset_count"].sum())
    n_species = len(counts_df)

    log.info(
        "Data ready: %d species, %d total organism-project assignments",
        n_species,
        total_datasets,
    )

    figure_path = plot_species_distribution(
        counts_df,
        args.output,
        top=args.top,
        plot_all=args.plot_all,
    )

    log.info("Figure saved to %s", figure_path.resolve())
    print(f"Species in dataset: {n_species}")
    print(f"Organism-project assignments: {total_datasets:,}")
    print(f"Cache: {cache_path.resolve()}")
    print(f"Figure: {figure_path.resolve()}")
    return 0


def main() -> None:
    """Run the CLI and exit with the pipeline status code."""
    sys.exit(run_pipeline(parse_args()))


if __name__ == "__main__":
    main()
