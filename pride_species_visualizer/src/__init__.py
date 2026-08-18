"""PRIDE Archive species dataset counter and visualizer."""

__version__ = "1.0.0"

from src.api_client import (
    PrideApiError,
    fetch_all_species_counts,
    load_counts_cache,
    save_counts_cache,
)
from src.visualizer import plot_species_distribution, prepare_plot_data

__all__ = [
    "__version__",
    "PrideApiError",
    "fetch_all_species_counts",
    "load_counts_cache",
    "save_counts_cache",
    "prepare_plot_data",
    "plot_species_distribution",
]
