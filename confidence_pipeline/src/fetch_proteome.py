"""UniProt canonical proteome retrieval and FASTA parsing."""

from __future__ import annotations

# Import logging so download progress can be reported to the console
import logging
# Import Path for filesystem operations on FASTA cache files
from pathlib import Path

# Import requests to stream the UniProt FASTA download over HTTP
import requests
# Import SeqIO to parse downloaded FASTA records with BioPython
from Bio import SeqIO

# Create a module-level logger for proteome fetch messages
logger = logging.getLogger(__name__)

# Default minimum cache size used to decide whether an existing FASTA is usable
DEFAULT_MIN_CACHE_BYTES = 1_000_000


def proteome_cache_is_valid(
    path: Path,
    *,
    min_bytes: int = DEFAULT_MIN_CACHE_BYTES,
) -> bool:
    """Return whether an on-disk proteome cache looks usable."""
    # Convert the incoming path object to a concrete Path instance
    cache_path = Path(path)
    # Reject missing cache files immediately
    if not cache_path.is_file():
        # Return false because there is no cached FASTA to reuse
        return False
    # Measure the cached file size in bytes
    file_size = cache_path.stat().st_size
    # Accept the cache only when it is larger than the configured minimum
    return file_size >= min_bytes


def fetch_uniprot_proteome(
    url: str,
    dest_path: Path,
    *,
    force: bool = False,
    min_cache_bytes: int = DEFAULT_MIN_CACHE_BYTES,
) -> Path:
    """Stream-download the human canonical proteome from UniProt."""
    # Resolve the destination path as a Path object
    output_path = Path(dest_path)
    # Ensure the parent directory exists before writing the FASTA file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Skip the network download when a valid cache already exists and force is false
    if not force and proteome_cache_is_valid(output_path, min_bytes=min_cache_bytes):
        # Log that the existing cache will be reused
        logger.info("Reusing cached proteome FASTA at %s", output_path)
        # Return the existing cache path without downloading again
        return output_path

    # Log the start of the UniProt streaming download
    logger.info("Downloading canonical proteome from %s", url)
    # Issue an HTTP GET request with streaming enabled for large FASTA payloads
    response = requests.get(url, stream=True, timeout=300)
    # Raise an exception when UniProt returns a non-success status code
    response.raise_for_status()

    # Open the destination FASTA file in binary write mode
    with output_path.open("wb") as handle:
        # Iterate over the response body in fixed-size chunks
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            # Write each non-empty chunk directly to the FASTA file
            if chunk:
                handle.write(chunk)

    # Log successful completion of the proteome download
    logger.info("Saved proteome FASTA to %s", output_path)
    # Return the path to the freshly downloaded FASTA cache
    return output_path


def load_proteome_sequences(fasta_path: Path) -> list[str]:
    """Parse a FASTA file and return uppercased protein sequences."""
    # Resolve the FASTA path as a Path object
    path = Path(fasta_path)
    # Raise a clear error when the FASTA file is missing
    if not path.is_file():
        # Stop early because proteome filtering cannot proceed without FASTA input
        raise FileNotFoundError(f"Proteome FASTA not found: {path}")

    # Initialize the list that will store parsed protein sequences
    sequences: list[str] = []
    # Parse every FASTA record using BioPython SeqIO
    for record in SeqIO.parse(path, "fasta"):
        # Convert the record sequence to an uppercase string and append it
        sequences.append(str(record.seq).upper())

    # Raise an error when the FASTA file parsed to zero sequences
    if not sequences:
        # Stop because an empty proteome would invalidate downstream filtering
        raise ValueError(f"No protein sequences parsed from FASTA: {path}")

    # Log how many protein sequences were loaded from the FASTA file
    logger.info("Loaded %d protein sequences from %s", len(sequences), path)
    # Return the full list of canonical protein sequences
    return sequences
