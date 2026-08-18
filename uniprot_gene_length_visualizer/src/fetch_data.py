"""UniProt API streaming, Biopython parsing, and AA-to-bp length conversion."""

# Enable postponed evaluation of annotations used in function signatures
from __future__ import annotations

# Import logging so download and parse progress can be reported to the console
import logging
# Import Path for typed filesystem destinations of cached FASTA files
from pathlib import Path
# Import Mapping and Any for flexible species metadata dictionaries
from typing import Any, Mapping, Sequence

# Import requests to stream UniProt FASTA payloads over HTTP
import requests
# Import pandas to assemble tidy length tables for downstream KDE plotting
import pandas as pd
# Import SeqIO to iterate FASTA sequence records with BioPython
from Bio import SeqIO

# Import centralized cache thresholds and URL templates from config
from config.config import (
    CHUNK_SIZE,
    DATA_DIR,
    FORCE_REFETCH,
    MIN_CACHE_BYTES,
    PROTEOMES,
    REQUEST_TIMEOUT,
    UNIPROT_STREAM_URL_TEMPLATE,
)


# Create a module-level logger for fetch and parse status messages
logger = logging.getLogger(__name__)


def build_stream_url(proteome_id: str) -> str:
    """Return the UniProt REST stream URL for a reviewed proteome FASTA."""
    # Substitute the proteome identifier into the configured URL template
    return UNIPROT_STREAM_URL_TEMPLATE.format(proteome_id=proteome_id)


def cache_is_valid(path: Path, *, min_bytes: int = MIN_CACHE_BYTES) -> bool:
    """Return whether an on-disk FASTA cache looks large enough to reuse."""
    # Normalize the incoming path argument into a Path instance
    cache_path = Path(path)
    # Reject the cache immediately when the destination file is missing
    if not cache_path.is_file():
        # Signal that a network download is required
        return False
    # Measure the existing cache file size in bytes
    file_size = cache_path.stat().st_size
    # Accept the cache only when it meets or exceeds the minimum size threshold
    return file_size >= min_bytes


def download_proteome(
    proteome_id: str,
    dest_path: Path,
    *,
    force: bool = FORCE_REFETCH,
    min_cache_bytes: int = MIN_CACHE_BYTES,
) -> Path:
    """Stream-download a reviewed UniProt proteome FASTA, reusing cache when valid."""
    # Resolve the destination path as a concrete Path object
    output_path = Path(dest_path)
    # Ensure the parent data directory exists before writing the FASTA file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Skip the network call when a usable cache already exists and force is false
    if not force and cache_is_valid(output_path, min_bytes=min_cache_bytes):
        # Log that the pipeline will reuse the local FASTA cache
        logger.info("Reusing cached proteome FASTA at %s", output_path)
        # Return the existing cache path without contacting UniProt
        return output_path
    # Build the UniProt streaming endpoint for this proteome identifier
    url = build_stream_url(proteome_id)
    # Log the start of the UniProt streaming download
    logger.info("Downloading reviewed proteome %s from %s", proteome_id, url)
    # Issue an HTTP GET with streaming enabled for large FASTA payloads
    response = requests.get(url, stream=True, timeout=REQUEST_TIMEOUT)
    # Raise an HTTPError when UniProt returns a non-success status code
    response.raise_for_status()
    # Open the destination FASTA file in binary write mode
    with output_path.open("wb") as handle:
        # Iterate over the response body in fixed-size chunks
        for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
            # Write each non-empty chunk directly to the FASTA cache file
            if chunk:
                # Persist the current network chunk to disk
                handle.write(chunk)
    # Log successful completion of the proteome download
    logger.info("Saved proteome FASTA to %s", output_path)
    # Return the path to the freshly downloaded FASTA cache
    return output_path


def parse_coding_bp_lengths(fasta_path: Path) -> pd.DataFrame:
    """Parse FASTA records into amino-acid and coding base-pair lengths."""
    # Initialize an empty list that will accumulate per-entry length records
    records: list[dict[str, Any]] = []
    # Iterate through every sequence record in the cached UniProt FASTA
    for seq_record in SeqIO.parse(str(fasta_path), "fasta"):
        # Measure the amino-acid length of the current protein sequence
        aa_length = len(seq_record.seq)
        # Convert amino-acid length to coding base-pair length (1 AA = 3 bp)
        bp_length = aa_length * 3
        # Append a tidy row capturing accession, AA length, and bp length
        records.append(
            {
                # Store the UniProt entry identifier from the FASTA header
                "entry_id": seq_record.id,
                # Store the amino-acid sequence length
                "aa_length": aa_length,
                # Store the coding base-pair length after AA x 3 conversion
                "bp_length": bp_length,
            }
        )
    # Construct a DataFrame from the accumulated per-entry dictionaries
    length_df = pd.DataFrame.from_records(records)
    # Return the tidy length table for downstream KDE plotting
    return length_df


def build_species_dataframe(
    species_meta: Mapping[str, str],
    *,
    force: bool = FORCE_REFETCH,
) -> pd.DataFrame:
    """Download and parse one species proteome into a tidy length DataFrame."""
    # Read the UniProt proteome accession from the species metadata mapping
    proteome_id = species_meta["proteome_id"]
    # Read the species display name used on plot legends
    display_name = species_meta["display_name"]
    # Build the absolute FASTA cache path under the configured data directory
    dest_path = DATA_DIR / species_meta["fasta_filename"]
    # Stream or reuse the reviewed proteome FASTA for this species
    fasta_path = download_proteome(proteome_id, dest_path, force=force)
    # Parse amino-acid and coding base-pair lengths from the cached FASTA
    length_df = parse_coding_bp_lengths(fasta_path)
    # Attach the human-readable species name as a column for hue mapping
    length_df["species"] = display_name
    # Attach the UniProt proteome identifier for traceability
    length_df["proteome_id"] = proteome_id
    # Reorder columns into a stable tidy schema for concatenation
    length_df = length_df[
        [
            # Keep species first for readability in exported summaries
            "species",
            # Keep proteome_id next to identify the UniProt source
            "proteome_id",
            # Keep the UniProt entry accession from the FASTA header
            "entry_id",
            # Keep amino-acid length prior to the bp conversion
            "aa_length",
            # Keep coding base-pair length after the AA x 3 conversion
            "bp_length",
        ]
    ]
    # Return the fully annotated species length table
    return length_df


def load_all_proteomes(
    proteomes: Sequence[Mapping[str, str]] | None = None,
    *,
    force: bool = FORCE_REFETCH,
) -> pd.DataFrame:
    """Download and parse all configured proteomes into one concatenated DataFrame."""
    # Default to the centralized PROTEOMES catalog when no override is supplied
    proteome_list = list(proteomes) if proteomes is not None else list(PROTEOMES)
    # Ensure the shared data cache directory exists before downloads begin
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    # Initialize a list that will collect per-species DataFrames
    frames: list[pd.DataFrame] = []
    # Iterate through each model organism proteome definition in order
    for species_meta in proteome_list:
        # Log which species is about to be ingested
        logger.info("Processing proteome for %s", species_meta["display_name"])
        # Build the continuous length table for the current species
        species_df = build_species_dataframe(species_meta, force=force)
        # Append the species table to the accumulation list
        frames.append(species_df)
    # Concatenate all species tables into a single tidy DataFrame
    combined_df = pd.concat(frames, ignore_index=True)
    # Log the total number of protein entries retained across all species
    logger.info("Loaded %d protein entries across %d species", len(combined_df), len(frames))
    # Return the multi-species continuous length table for KDE plotting
    return combined_df
