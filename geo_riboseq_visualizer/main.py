#!/usr/bin/env python3
"""CLI entry point for GEO Ribo-seq species aggregation and visualization."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from src.geo_client import build_riboseq_query, create_entrez_config, fetch_species_counts
from src.visualizer import plot_species_distribution

DEFAULT_OUTPUT = Path("output/geo_riboseq_species_counts.png")
DEFAULT_CACHE = Path("output/species_counts.csv")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments.

    Args:
        argv: Optional argv list. If None, argparse uses sys.argv.

    Returns:
        Parsed arguments.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Query NCBI GEO via Entrez for Ribo-seq studies, aggregate dataset counts "
            "per species, and generate a publication-ready vertical log-scale plot."
        )
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output figure path.")
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE, help="CSV cache path.")
    parser.add_argument("--refresh", action="store_true", help="Re-query Entrez even if cache exists.")
    parser.add_argument("--top", type=int, default=30, help="Top species to show after 'Other'.")
    parser.add_argument("--batch-size", type=int, default=200, help="Entrez esummary batch size.")
    parser.add_argument("--sleep-seconds", type=float, default=0.34, help="Delay between Entrez calls.")
    parser.add_argument("--api-key", type=str, default=None, help="NCBI API key (optional).")
    parser.add_argument("--query", type=str, default=None, help="Optional override Entrez query term.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging.")
    return parser.parse_args(argv)


def load_or_fetch_counts(args: argparse.Namespace) -> pd.DataFrame:
    """Load cached counts or fetch fresh counts from GEO.

    Args:
        args: CLI args.

    Returns:
        DataFrame with columns species and dataset_count.
    """

    if args.cache.is_file() and not args.refresh:
        logging.getLogger(__name__).info("Loading cached counts from %s", args.cache)
        df = pd.read_csv(args.cache)
        df["dataset_count"] = df["dataset_count"].astype(int)
        return df

    create_entrez_config(api_key=args.api_key)
    term = args.query or build_riboseq_query()
    logging.getLogger(__name__).info("Entrez term: %s", term)
    df = fetch_species_counts(term, batch_size=args.batch_size, sleep_seconds=args.sleep_seconds)
    args.cache.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.cache, index=False)
    logging.getLogger(__name__).info("Wrote cache to %s", args.cache)
    return df


def main(argv: list[str] | None = None) -> int:
    """Run the CLI pipeline.

    Returns:
        Exit code (0 on success).
    """

    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    df = load_or_fetch_counts(args)
    out = plot_species_distribution(df, args.output, top=args.top)
    logging.getLogger(__name__).info("Saved figure to %s", out.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

