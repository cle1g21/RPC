#!/usr/bin/env python3
"""Main orchestrator for proteome filtering and 4-tier confidence validation."""

from __future__ import annotations

# Import argparse so command-line flags can override configuration defaults
import argparse
# Import json so run summaries can be written to disk
import json
# Import logging so pipeline progress is visible in the console
import logging
# Import sys so the project root can be added to the module search path
import sys
# Import Path for filesystem path handling throughout the orchestrator
from pathlib import Path
# Import Any because run summaries use mixed Python value types
from typing import Any

# Resolve the directory that contains this orchestrator script
_ROOT = Path(__file__).resolve().parent
# Insert the project root at the front of sys.path for local imports
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

# Import the default pipeline configuration object
from config.config import PipelineConfig, default_config
# Import proteome retrieval helpers from the local source package
from src.fetch_proteome import fetch_uniprot_proteome, load_proteome_sequences
# Import filesystem scanning and table I/O helpers
from src.file_io import (
    derive_output_name,
    ensure_dir,
    orf_run_name,
    read_table,
    scan_input_files,
    write_table,
)
# Import matching helpers from the local file matcher module
from src.file_matcher import (
    build_protein_substring_index,
    build_sequence_set,
    filter_table_against_proteome,
    match_table,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the confidence pipeline."""
    # Create the top-level argument parser for the pipeline CLI
    parser = argparse.ArgumentParser(
        description=(
            "Filter NTv3 ORFs against the canonical proteome and run "
            "4-tier confidence validation."
        )
    )
    # Add a flag to skip the UniProt proteome download step
    parser.add_argument(
        "--skip-proteome-fetch",
        action="store_true",
        help="Skip UniProt download and require an existing proteome FASTA cache.",
    )
    # Add a flag to force re-download of the UniProt proteome FASTA cache
    parser.add_argument(
        "--refetch-proteome",
        action="store_true",
        help="Force a fresh UniProt proteome download even if cache exists.",
    )
    # Add a flag to skip the post-processing assembly step
    parser.add_argument(
        "--skip-assembly",
        action="store_true",
        help="Do not run assemble_validation_results.py after tier processing.",
    )
    # Add a flag to enable verbose debug logging
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    # Parse command-line arguments and return the resulting namespace
    return parser.parse_args(argv)


def ensure_output_directories(config: PipelineConfig) -> None:
    """Create all pipeline output directories if they do not already exist."""
    # Build the list of directories that must exist before writing outputs
    directories = [
        config.proteome_cache_path.parent,
        config.ntv3_filtered_output_dir,
        config.confidence_low_dir,
        config.confidence_med_ntv3_dir,
        config.confidence_med_instanovo_dir,
        config.confidence_high_dir,
        config.confidence_assembled_dir,
    ]
    # Create each required directory recursively
    for directory in directories:
        # Make the directory and any missing parent folders
        ensure_dir(directory)
    # Create per-ORF-run subdirectories under the Low and High tier folders
    for orf_file in config.orf_files:
        # Resolve the short run label for the current ORF source file
        run_label = orf_run_name(orf_file)
        # Create the Low-tier subdirectory for this ORF run
        ensure_dir(config.confidence_low_dir / run_label)
        # Create the High-tier subdirectory for this ORF run
        ensure_dir(config.confidence_high_dir / run_label)


def filter_ntv3_handoff_files(
    config: PipelineConfig,
    protein_index: Any,
) -> dict[str, dict[str, Any]]:
    """Apply proteome substring exclusion to each configured NTv3 ORF file."""
    # Initialize the aggregate summary dictionary for NTv3 filtering
    summaries: dict[str, dict[str, Any]] = {}
    # Process every configured ORF handoff file independently
    for orf_file in config.orf_files:
        # Resolve the full path to the current ORF handoff source file
        source_path = config.ntv3_handoff_dir / orf_file
        # Read the ORF handoff table using delimiter sniffing
        orf_df = read_table(source_path)
        # Apply proteome exclusion to remove known canonical protein ORF sequences
        filtered_df, stats = filter_table_against_proteome(
            orf_df,
            config.ntv3_sequence_column,
            protein_index,
            config,
        )
        # Build the filtered output filename for the current ORF source file
        output_name = derive_output_name(source_path, "_filtered")
        # Resolve the full output path inside the NTv3 filtered handoff directory
        output_path = config.ntv3_filtered_output_dir / output_name
        # Write the filtered ORF table to disk as CSV
        write_table(filtered_df, output_path)
        # Store per-file statistics in the aggregate summary dictionary
        summaries[orf_file] = {
            **stats,
            "output_path": str(output_path),
        }
    # Return per-file NTv3 filtering summaries
    return summaries


def run_tier_med_instanovo(
    config: PipelineConfig,
    validation_sequences: set[str],
) -> dict[str, dict[str, int]]:
    """Run Tier 3: keep InstaNovo filtered rows matching the validation anchor."""
    # Discover all InstaNovo filtered prediction files in the input directory
    instanovo_files = scan_input_files(
        config.instanovo_filtered_dir,
        (config.instanovo_input_glob,),
    )
    # Initialize the per-file summary dictionary for Tier 3
    summaries: dict[str, dict[str, int]] = {}
    # Process every discovered InstaNovo filtered file independently
    for input_path in instanovo_files:
        # Read the filtered InstaNovo prediction table into a DataFrame
        instanovo_df = read_table(input_path)
        # Count how many rows are present before validation matching
        input_rows = len(instanovo_df)
        # Keep rows whose normalized peptide sequence exactly matches the validation set
        matched_df = match_table(
            instanovo_df,
            config.instanovo_sequence_column,
            validation_sequences,
            mode=config.tier_match_modes["med_instanovo"],
            strip_mods=config.strip_modifications,
            il_equivalent=False,
        )
        # Build the Tier 3 output filename with the configured suffix
        output_name = derive_output_name(input_path, config.tier_suffixes["med_instanovo"])
        # Resolve the full output path inside the Medium-InstaNovo directory
        output_path = config.confidence_med_instanovo_dir / output_name
        # Write the matched InstaNovo rows to the Tier 3 output CSV
        write_table(matched_df, output_path)
        # Store per-file counts for the run summary JSON
        summaries[input_path.name] = {
            "input_rows": input_rows,
            "matched_rows": len(matched_df),
        }
    # Return per-file Tier 3 summaries
    return summaries


def run_tier_low(
    config: PipelineConfig,
    orf_file: str,
    ntv3_filtered_path: Path,
) -> dict[str, dict[str, int]]:
    """Run Tier 1: keep InstaNovo filtered rows whose peptides appear in NTv3 ORFs."""
    # Resolve the short run label for the current ORF source file
    run_label = orf_run_name(orf_file)
    # Read the filtered NTv3 ORF table for this ORF run
    ntv3_df = read_table(ntv3_filtered_path)
    # Extract the list of filtered NTv3 ORF amino acid sequences
    ntv3_sequences = ntv3_df[config.ntv3_sequence_column].dropna().astype(str).tolist()
    # Discover all InstaNovo filtered prediction files in the input directory
    instanovo_files = scan_input_files(
        config.instanovo_filtered_dir,
        (config.instanovo_input_glob,),
    )
    # Initialize the per-file summary dictionary for Tier 1
    summaries: dict[str, dict[str, int]] = {}
    # Process every discovered InstaNovo filtered file independently
    for input_path in instanovo_files:
        # Read the filtered InstaNovo prediction table into a DataFrame
        instanovo_df = read_table(input_path)
        # Count how many rows are present before NTv3 overlap matching
        input_rows = len(instanovo_df)
        # Keep InstaNovo rows whose peptide is a substring of any filtered NTv3 ORF
        matched_df = match_table(
            instanovo_df,
            config.instanovo_sequence_column,
            ntv3_sequences,
            mode=config.tier_match_modes["low"],
            direction="query_in_reference",
            strip_mods=config.strip_modifications,
            il_equivalent=False,
        )
        # Build the Tier 1 output filename retaining the original InstaNovo base name
        output_name = f"{input_path.stem}.csv"
        # Resolve the full output path inside the Low-tier ORF-run subdirectory
        output_path = config.confidence_low_dir / run_label / output_name
        # Write the matched InstaNovo rows to the Tier 1 output CSV
        write_table(matched_df, output_path)
        # Store per-file counts for the run summary JSON
        summaries[input_path.name] = {
            "input_rows": input_rows,
            "matched_rows": len(matched_df),
        }
    # Return per-file Tier 1 summaries for this ORF run
    return summaries


def run_tier_med_ntv3(
    config: PipelineConfig,
    orf_file: str,
    ntv3_filtered_path: Path,
    validation_sequences: set[str],
) -> dict[str, int]:
    """Run Tier 2: keep NTv3 ORF rows that contain a validation peptide."""
    # Read the filtered NTv3 ORF table for this ORF run
    ntv3_df = read_table(ntv3_filtered_path)
    # Count how many ORF rows are present before validation matching
    input_rows = len(ntv3_df)
    # Keep ORF rows that contain any validation peptide as an exact substring
    matched_df = match_table(
        ntv3_df,
        config.ntv3_sequence_column,
        list(validation_sequences),
        mode=config.tier_match_modes["med_ntv3"],
        direction="reference_in_query",
        strip_mods=config.strip_modifications,
        il_equivalent=False,
    )
    # Build the Tier 2 output filename with the configured suffix
    output_name = derive_output_name(ntv3_filtered_path, config.tier_suffixes["med_ntv3"])
    # Resolve the full output path inside the Mid-NTv3 directory
    output_path = config.confidence_med_ntv3_dir / output_name
    # Write the matched NTv3 ORF rows to the Tier 2 output CSV
    write_table(matched_df, output_path)
    # Return a single summary entry for this ORF run Tier 2 output
    return {
        "input_rows": input_rows,
        "matched_rows": len(matched_df),
        "output_path": str(output_path),
    }


def run_tier_high(
    config: PipelineConfig,
    orf_file: str,
    validation_sequences: set[str],
) -> dict[str, dict[str, int]]:
    """Run Tier 4: keep Low-tier rows whose peptides exactly match the validation anchor."""
    # Resolve the short run label for the current ORF source file
    run_label = orf_run_name(orf_file)
    # Resolve the Low-tier ORF-run subdirectory containing Tier 1 outputs
    low_run_dir = config.confidence_low_dir / run_label
    # Discover all Low-tier output files for this ORF run
    low_files = scan_input_files(low_run_dir, ("*.csv",))
    # Initialize the per-file summary dictionary for Tier 4
    summaries: dict[str, dict[str, int]] = {}
    # Process every discovered Low-tier output file independently
    for input_path in low_files:
        # Read the Low-tier InstaNovo prediction table into a DataFrame
        low_df = read_table(input_path)
        # Count how many rows are present before validation matching
        input_rows = len(low_df)
        # Keep Low-tier rows whose normalized peptide exactly matches the validation set
        matched_df = match_table(
            low_df,
            config.instanovo_sequence_column,
            validation_sequences,
            mode=config.tier_match_modes["high"],
            strip_mods=config.strip_modifications,
            il_equivalent=False,
        )
        # Build the Tier 4 output filename with the configured suffix
        output_name = derive_output_name(input_path, config.tier_suffixes["high"])
        # Resolve the full output path inside the High-tier ORF-run subdirectory
        output_path = config.confidence_high_dir / run_label / output_name
        # Write the matched rows to the Tier 4 output CSV
        write_table(matched_df, output_path)
        # Store per-file counts for the run summary JSON
        summaries[input_path.name] = {
            "input_rows": input_rows,
            "matched_rows": len(matched_df),
        }
    # Return per-file Tier 4 summaries for this ORF run
    return summaries


def run_pipeline(
    config: PipelineConfig,
    *,
    skip_proteome_fetch: bool = False,
) -> dict[str, Any]:
    """Execute proteome filtering and all four confidence tiers."""
    # Create a module-level logger for orchestration messages
    log = logging.getLogger(__name__)
    # Initialize the aggregate run summary dictionary
    run_summary: dict[str, Any] = {
        "config": config.as_dict(),
        "ntv3_proteome_filter": {},
        "tiers": {
            "med_instanovo": {},
            "by_orf_run": {},
        },
    }

    # Create required output directories before any file-writing step
    ensure_output_directories(config)

    # Resolve the local FASTA cache path from configuration
    proteome_path = config.proteome_cache_path
    # Download or reuse the canonical proteome FASTA unless fetch is skipped
    if not skip_proteome_fetch:
        # Stream the UniProt FASTA into the configured cache path
        fetch_uniprot_proteome(
            config.uniprot_stream_url,
            proteome_path,
            force=config.refetch_proteome,
        )
    else:
        # Verify that a cached FASTA exists when fetch is skipped by CLI
        if not proteome_path.is_file():
            # Stop because proteome filtering cannot proceed without FASTA input
            raise FileNotFoundError(
                f"--skip-proteome-fetch was set but FASTA cache not found: {proteome_path}"
            )
        # Log that the cached FASTA will be used without downloading
        log.info("Using existing proteome FASTA at %s", proteome_path)

    # Load all canonical protein sequences from the FASTA cache
    protein_sequences = load_proteome_sequences(proteome_path)
    # Store the number of loaded protein sequences in the run summary
    run_summary["proteome_protein_count"] = len(protein_sequences)
    # Build the in-memory protein substring index used for proteome exclusion
    protein_index = build_protein_substring_index(
        protein_sequences,
        il_equivalent=config.treat_leucine_isoleucine_as_identical,
    )

    # Phase 1: filter each NTv3 ORF handoff file against the canonical proteome
    run_summary["ntv3_proteome_filter"] = filter_ntv3_handoff_files(config, protein_index)

    # Read the master validation anchor table into a DataFrame
    validation_df = read_table(config.validation_csv_path)
    # Build the normalized validation peptide set used by multiple tiers
    validation_sequences = build_sequence_set(
        validation_df,
        config.validation_sequence_column,
        strip_mods=config.strip_modifications,
        il_equivalent=False,
    )
    # Log how many unique validation sequences will be used for matching
    log.info("Loaded %d unique validation sequences", len(validation_sequences))

    # Tier 3 is ORF-independent and can run once across all InstaNovo filtered files
    log.info("=== Tier 3: Medium-InstaNovo ===")
    # Run the Medium-InstaNovo validation tier and capture per-file counts
    run_summary["tiers"]["med_instanovo"] = run_tier_med_instanovo(
        config,
        validation_sequences,
    )

    # Process each ORF source file as a separate confidence run
    for orf_file in config.orf_files:
        # Resolve the short run label for the current ORF source file
        run_label = orf_run_name(orf_file)
        # Resolve the filtered NTv3 output path produced in Phase 1
        ntv3_filtered_path = config.ntv3_filtered_output_dir / derive_output_name(
            config.ntv3_handoff_dir / orf_file,
            "_filtered",
        )
        # Log the start of the current ORF-run processing block
        log.info("=== ORF run: %s (%s) ===", run_label, orf_file)

        # Tier 1: InstaNovo filtered peptides overlapping filtered NTv3 ORFs
        log.info("--- Tier 1: Low ---")
        # Run the Low confidence tier for the current ORF run
        low_summary = run_tier_low(config, orf_file, ntv3_filtered_path)

        # Tier 2: filtered NTv3 ORFs overlapping validation peptides
        log.info("--- Tier 2: Mid-NTv3 ---")
        # Run the Mid-NTv3 confidence tier for the current ORF run
        med_ntv3_summary = run_tier_med_ntv3(
            config,
            orf_file,
            ntv3_filtered_path,
            validation_sequences,
        )

        # Tier 4: Low-tier outputs overlapping validation peptides
        log.info("--- Tier 4: High ---")
        # Run the High confidence tier for the current ORF run
        high_summary = run_tier_high(config, orf_file, validation_sequences)

        # Store all tier summaries for the current ORF run
        run_summary["tiers"]["by_orf_run"][run_label] = {
            "orf_file": orf_file,
            "low": low_summary,
            "med_ntv3": med_ntv3_summary,
            "high": high_summary,
        }

    # Serialize the run summary dictionary to formatted JSON on disk
    config.run_summary_path.write_text(
        json.dumps(run_summary, indent=2),
        encoding="utf-8",
    )
    # Log the final location of the run summary JSON file
    log.info("Run summary written to %s", config.run_summary_path)
    # Return the in-memory run summary dictionary
    return run_summary


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for the confidence pipeline."""
    # Parse command-line arguments into an argparse namespace
    args = parse_args(argv)
    # Configure root logging level based on the verbose flag
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    # Start from a fresh copy of the default pipeline configuration
    config = default_config()
    # Force proteome re-download when the CLI refetch flag is set
    if args.refetch_proteome:
        # Set the refetch flag so fetch_proteome ignores an existing cache
        config.refetch_proteome = True
    # Execute the full pipeline using the resolved configuration
    run_pipeline(config, skip_proteome_fetch=args.skip_proteome_fetch)

    # Run the assembly utility unless the caller asked to skip it
    if not args.skip_assembly:
        # Import the assembly main function only when assembly is requested
        from assemble_validation_results import main as assemble_main

        # Execute the post-processing assembly step without forwarding CLI flags
        assemble_main([])

    # Return exit code zero to indicate successful completion
    return 0


# Run the CLI entry point only when this file is executed as a script
if __name__ == "__main__":
    # Raise SystemExit with the return code from main()
    raise SystemExit(main())
