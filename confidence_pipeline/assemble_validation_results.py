#!/usr/bin/env python3
"""Post-processing assembly of dual-ORF confidence tier outputs."""

from __future__ import annotations

# Import argparse so command-line flags can control assembly behavior
import argparse
# Import logging so assembly progress is visible in the console
import logging
# Import sys so the project root can be added to the module search path
import sys
# Import Path for filesystem path handling throughout the assembly utility
from pathlib import Path

# Resolve the directory that contains this assembly script
_ROOT = Path(__file__).resolve().parent
# Insert the project root at the front of sys.path for local imports
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

# Import ORF source tier labels and default configuration values
from config.config import ORF_FILE_ALL, ORF_FILE_GE30, ORF_SOURCE_TIER_LABELS, default_config
# Import filesystem scanning and table I/O helpers
from src.file_io import ensure_dir, orf_run_name, read_table, scan_input_files, write_table
# Import sequence normalization helpers used for deduplication keys
from src.string_normalizer import normalize_sequence


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the assembly utility."""
    # Create the top-level argument parser for the assembly CLI
    parser = argparse.ArgumentParser(
        description="Merge dual-ORF confidence tier outputs into master files."
    )
    # Add a flag to enable verbose debug logging
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging.")
    # Parse command-line arguments and return the resulting namespace
    return parser.parse_args(argv)


def sequence_key(
    row_value: object,
    sequence_column: str,
    *,
    strip_mods: bool,
) -> str:
    """Build a normalized deduplication key for one sequence value."""
    # Normalize the row sequence using the configured modification stripping rule
    return normalize_sequence(str(row_value), strip_mods=strip_mods, il_equivalent=False)


def label_rows(
    dataframe,
    sequence_column: str,
    label: str,
    *,
    strip_mods: bool,
):
    """Add normalized sequence key and ORF source tier label columns to a DataFrame."""
    # Copy the input DataFrame so the original table is not modified in place
    output_df = dataframe.copy()
    # Add a normalized sequence key column used for deduplication
    output_df["_sequence_key"] = output_df[sequence_column].apply(
        lambda value: sequence_key(value, sequence_column, strip_mods=strip_mods)
    )
    # Add the human-readable ORF source tier label column requested by the specification
    output_df["orf_source_tier"] = label
    # Return the labeled DataFrame copy
    return output_df


def assemble_tier_low_or_high(
    tier_name: str,
    tier_dir: Path,
    sequence_column: str,
    output_path: Path,
    *,
    strip_mods: bool,
) -> dict[str, int]:
    """Merge Low or High tier outputs from both ORF runs into one master file."""
    # Initialize the list that will collect labeled DataFrames from each ORF run
    labeled_frames = []
    # Initialize the set of sequence keys already contributed by the conservative run
    conservative_keys: set[str] = set()

    # Process the conservative ge30aa ORF run first when its directory exists
    ge30_run_dir = tier_dir / orf_run_name(ORF_FILE_GE30)
    # Check whether the ge30aa ORF-run subdirectory exists on disk
    if ge30_run_dir.is_dir():
        # Discover all CSV files in the ge30aa ORF-run subdirectory
        ge30_files = scan_input_files(ge30_run_dir, ("*.csv",))
        # Read and label every ge30aa ORF-run file
        for file_path in ge30_files:
            # Read the current tier output table into a DataFrame
            table = read_table(file_path)
            # Label rows as conservative and add a deduplication key column
            labeled = label_rows(
                table,
                sequence_column,
                ORF_SOURCE_TIER_LABELS[ORF_FILE_GE30],
                strip_mods=strip_mods,
            )
            # Extend the conservative key set with keys from the labeled frame
            conservative_keys.update(labeled["_sequence_key"].tolist())
            # Append the labeled frame to the aggregate list
            labeled_frames.append(labeled)

    # Process the full ORF run second when its directory exists
    all_run_dir = tier_dir / orf_run_name(ORF_FILE_ALL)
    # Check whether the full ORF-run subdirectory exists on disk
    if all_run_dir.is_dir():
        # Discover all CSV files in the full ORF-run subdirectory
        all_files = scan_input_files(all_run_dir, ("*.csv",))
        # Read and label every full ORF-run file
        for file_path in all_files:
            # Read the current tier output table into a DataFrame
            table = read_table(file_path)
            # Label rows as short fragment and add a deduplication key column
            labeled = label_rows(
                table,
                sequence_column,
                ORF_SOURCE_TIER_LABELS[ORF_FILE_ALL],
                strip_mods=strip_mods,
            )
            # Keep only rows whose sequence key is not already present in the conservative run
            novel = labeled[~labeled["_sequence_key"].isin(conservative_keys)].copy()
            # Append the novel short-fragment rows to the aggregate list
            labeled_frames.append(novel)

    # Import pandas here so the module can still be imported without pandas in type checkers
    import pandas as pd

    # Handle the case where no tier outputs were found for assembly
    if not labeled_frames:
        # Import pandas here for empty master file creation
        import pandas as pd

        # Write an empty master file when no tier inputs exist
        empty_df = pd.DataFrame()
        # Ensure the assembled output directory exists before writing
        ensure_dir(output_path.parent)
        # Write an empty CSV master file and return zero counts
        write_table(empty_df, output_path)
        # Return zero input and output row counts
        return {"input_rows": 0, "output_rows": 0}

    # Concatenate all labeled frames into one master DataFrame
    merged = pd.concat(labeled_frames, ignore_index=True)
    # Count how many rows were present before deduplication
    input_rows = len(merged)
    # Drop duplicate rows using the normalized sequence key to prevent duplicate insertion
    deduped = merged.drop_duplicates(subset=["_sequence_key"], keep="first")
    # Remove the internal deduplication helper column before writing the master file
    deduped = deduped.drop(columns=["_sequence_key"])
    # Write the deduplicated master file to the assembled output directory
    write_table(deduped, output_path)
    # Return input and output row counts for logging and summary reporting
    return {"input_rows": input_rows, "output_rows": len(deduped)}


def assemble_tier_med_ntv3(
    med_ntv3_dir: Path,
    sequence_column: str,
    output_path: Path,
    *,
    strip_mods: bool,
) -> dict[str, int]:
    """Merge Mid-NTv3 outputs from both ORF runs into one master file."""
    # Initialize the list that will collect labeled DataFrames from each ORF run
    labeled_frames = []
    # Initialize the set of sequence keys already contributed by the conservative run
    conservative_keys: set[str] = set()

    # Build the expected Tier 2 output filenames for each ORF source file
    expected_files = {
        ORF_FILE_GE30: med_ntv3_dir / "merged_orfs_ge30aa_filtered_medNTv3.csv",
        ORF_FILE_ALL: med_ntv3_dir / "merged_orfs_filtered_medNTv3.csv",
    }

    # Process the conservative ge30aa Tier 2 output first when the file exists
    if expected_files[ORF_FILE_GE30].is_file():
        # Read the ge30aa Mid-NTv3 output table into a DataFrame
        table = read_table(expected_files[ORF_FILE_GE30])
        # Label rows as conservative and add a deduplication key column
        labeled = label_rows(
            table,
            sequence_column,
            ORF_SOURCE_TIER_LABELS[ORF_FILE_GE30],
            strip_mods=strip_mods,
        )
        # Extend the conservative key set with keys from the labeled frame
        conservative_keys.update(labeled["_sequence_key"].tolist())
        # Append the labeled frame to the aggregate list
        labeled_frames.append(labeled)

    # Process the full ORF Tier 2 output second when the file exists
    if expected_files[ORF_FILE_ALL].is_file():
        # Read the full ORF Mid-NTv3 output table into a DataFrame
        table = read_table(expected_files[ORF_FILE_ALL])
        # Label rows as short fragment and add a deduplication key column
        labeled = label_rows(
            table,
            sequence_column,
            ORF_SOURCE_TIER_LABELS[ORF_FILE_ALL],
            strip_mods=strip_mods,
        )
        # Keep only rows whose sequence key is not already present in the conservative run
        novel = labeled[~labeled["_sequence_key"].isin(conservative_keys)].copy()
        # Append the novel short-fragment rows to the aggregate list
        labeled_frames.append(novel)

    # Import pandas here so the module can still be imported without pandas in type checkers
    import pandas as pd

    # Handle the case where no Mid-NTv3 outputs were found for assembly
    if not labeled_frames:
        # Write an empty master file when no Tier 2 inputs exist
        empty_df = pd.DataFrame()
        # Ensure the assembled output directory exists before writing
        ensure_dir(output_path.parent)
        # Write an empty CSV master file and return zero counts
        write_table(empty_df, output_path)
        # Return zero input and output row counts
        return {"input_rows": 0, "output_rows": 0}

    # Concatenate all labeled frames into one master DataFrame
    merged = pd.concat(labeled_frames, ignore_index=True)
    # Count how many rows were present before deduplication
    input_rows = len(merged)
    # Drop duplicate rows using the normalized sequence key to prevent duplicate insertion
    deduped = merged.drop_duplicates(subset=["_sequence_key"], keep="first")
    # Remove the internal deduplication helper column before writing the master file
    deduped = deduped.drop(columns=["_sequence_key"])
    # Write the deduplicated master file to the assembled output directory
    write_table(deduped, output_path)
    # Return input and output row counts for logging and summary reporting
    return {"input_rows": input_rows, "output_rows": len(deduped)}


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for post-processing assembly."""
    # Parse command-line arguments into an argparse namespace
    args = parse_args(argv)
    # Configure root logging level based on the verbose flag
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    # Create a module-level logger for assembly messages
    log = logging.getLogger(__name__)
    # Start from a fresh copy of the default pipeline configuration
    config = default_config()
    # Ensure the assembled output directory exists before writing master files
    ensure_dir(config.confidence_assembled_dir)

    # Assemble Low-tier outputs into a single master CSV file
    low_stats = assemble_tier_low_or_high(
        "low",
        config.confidence_low_dir,
        config.instanovo_sequence_column,
        config.confidence_assembled_dir / "master_low.csv",
        strip_mods=config.strip_modifications,
    )
    # Log the Low-tier assembly row counts
    log.info("Assembled master_low.csv: %s", low_stats)

    # Assemble Mid-NTv3 outputs into a single master CSV file
    med_ntv3_stats = assemble_tier_med_ntv3(
        config.confidence_med_ntv3_dir,
        config.ntv3_sequence_column,
        config.confidence_assembled_dir / "master_med_ntv3.csv",
        strip_mods=config.strip_modifications,
    )
    # Log the Mid-NTv3 assembly row counts
    log.info("Assembled master_med_ntv3.csv: %s", med_ntv3_stats)

    # Assemble High-tier outputs into a single master CSV file
    high_stats = assemble_tier_low_or_high(
        "high",
        config.confidence_high_dir,
        config.instanovo_sequence_column,
        config.confidence_assembled_dir / "master_high.csv",
        strip_mods=config.strip_modifications,
    )
    # Log the High-tier assembly row counts
    log.info("Assembled master_high.csv: %s", high_stats)

    # Return exit code zero to indicate successful completion
    return 0


# Run the CLI entry point only when this file is executed as a script
if __name__ == "__main__":
    # Raise SystemExit with the return code from main()
    raise SystemExit(main())
