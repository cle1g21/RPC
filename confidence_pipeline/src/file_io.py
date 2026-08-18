"""Filesystem scanning and delimiter-aware table I/O utilities."""

from __future__ import annotations

# Import csv so delimiter sniffing can inspect the first file line
import csv
# Import Path for directory scanning and output path construction
from pathlib import Path
# Import Sequence so glob patterns can be accepted as a tuple or list
from typing import Sequence

# Import pandas for reading and writing CSV and TSV tables
import pandas as pd


def ensure_dir(path: Path) -> Path:
    """Create a directory and all missing parent folders when needed."""
    # Resolve the incoming path as a concrete Path object
    directory = Path(path)
    # Create the directory tree recursively without raising if it already exists
    directory.mkdir(parents=True, exist_ok=True)
    # Return the resolved directory path for downstream use
    return directory


def scan_input_files(directory: Path, patterns: Sequence[str]) -> list[Path]:
    """Discover input files in a directory using one or more glob patterns."""
    # Resolve the scan directory as a Path object
    root = Path(directory)
    # Raise a clear error when the input directory does not exist
    if not root.is_dir():
        # Stop because there is nothing to scan
        raise FileNotFoundError(f"Input directory not found: {root}")

    # Initialize the list that will collect discovered file paths
    discovered: list[Path] = []
    # Apply every configured glob pattern against the input directory
    for pattern in patterns:
        # Extend the discovered list with paths matching the current pattern
        discovered.extend(sorted(root.glob(pattern)))

    # Remove duplicate paths while preserving sorted order
    unique_paths = sorted(set(discovered))
    # Return the sorted list of unique discovered files
    return unique_paths


def sniff_delimiter(path: Path) -> str:
    """Detect whether a table file is comma- or tab-delimited."""
    # Open the file in text mode for delimiter sniffing
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        # Read the first non-empty line from the file header
        sample = handle.readline()
    # Default to comma when the header line is unexpectedly empty
    if not sample:
        # Return comma because an empty sample cannot be sniffed reliably
        return ","
    # Use Python's csv.Sniffer to guess the delimiter from the header line
    try:
        # Return the sniffed delimiter when detection succeeds
        return csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
    except csv.Error:
        # Fall back to comma when delimiter sniffing fails
        return ","


def read_table(path: Path) -> pd.DataFrame:
    """Read a CSV or TSV file using delimiter sniffing."""
    # Resolve the table path as a Path object
    table_path = Path(path)
    # Detect the delimiter so comma-in-.tsv files parse correctly
    delimiter = sniff_delimiter(table_path)
    # Read the delimited table into a pandas DataFrame
    dataframe = pd.read_csv(table_path, sep=delimiter)
    # Return the parsed DataFrame to the caller
    return dataframe


def write_table(dataframe: pd.DataFrame, path: Path) -> Path:
    """Write a DataFrame to CSV using comma separation."""
    # Resolve the output path as a Path object
    output_path = Path(path)
    # Ensure the parent directory exists before writing the output file
    ensure_dir(output_path.parent)
    # Write the DataFrame to CSV without the pandas index column
    dataframe.to_csv(output_path, index=False)
    # Return the path to the written output file
    return output_path


def derive_output_name(input_path: Path, suffix: str, extension: str = ".csv") -> str:
    """Build an output filename from an input stem plus a suffix."""
    # Take the input stem without any directory components
    stem = Path(input_path).stem
    # Return the stem plus suffix plus extension as the output filename
    return f"{stem}{suffix}{extension}"


def orf_run_name(orf_filename: str) -> str:
    """Map an ORF source filename to a short run directory label."""
    # Import the ORF run label mapping from the configuration module
    from config.config import ORF_RUN_LABELS

    # Return the configured run label when the filename is known
    return ORF_RUN_LABELS.get(orf_filename, Path(orf_filename).stem)
