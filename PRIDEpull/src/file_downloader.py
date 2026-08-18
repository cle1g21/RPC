"""Chunked streaming download manager for PRIDE public files."""

# Enable postponed evaluation of type annotations
from __future__ import annotations

# Import ftplib for PRIDE FTP URLs (requests does not support ftp://)
from ftplib import FTP

# Import logging for download progress messages
import logging

# Import Path for destination path handling
from pathlib import Path

# Import urlparse to split ftp:// host and remote path components
from urllib.parse import urlparse

# Import requests for streaming HTTP/HTTPS GET
import requests

# Import config for chunk size and landing directory
from config import config as cfg

# Create a module-level logger
logger = logging.getLogger(__name__)


def _ensure_parent_directory(dest_path: Path) -> None:
    """Create the parent directory for a destination file if it does not exist."""

    # Call mkdir with parents=True so nested accession subdirs are created
    dest_path.parent.mkdir(parents=True, exist_ok=True)


def _download_via_ftp(
    download_url: str,
    partial_path: Path,
    *,
    mode: str,
    resume_offset: int,
) -> int:
    """Download a file from an ftp:// URL using ftplib with chunked writes."""

    # Parse the FTP URL into hostname and remote file path components
    parsed = urlparse(download_url)

    # Read the FTP server hostname from the parsed URL
    host = parsed.hostname

    # Raise when the URL does not contain a hostname
    if not host:
        raise ValueError(f"Invalid FTP URL (no host): {download_url}")

    # Read the remote path on the FTP server (leading slash included)
    remote_path = parsed.path

    # Connect to the PRIDE FTP server with configured read timeout
    ftp = FTP(host, timeout=cfg.PRIDE_TIMEOUT[1])

    # Log in anonymously (PRIDE public FTP requires no credentials)
    ftp.login()

    # Track total bytes written during this download session
    bytes_written = resume_offset

    # Open the partial destination file for binary writing
    with partial_path.open(mode) as out_file:
        # When resuming, tell the FTP server to restart at resume_offset
        if resume_offset > 0:
            ftp.sendcmd(f"REST {resume_offset}")

        # Define a callback that writes each FTP binary chunk to disk
        def _write_chunk(chunk: bytes) -> None:
            nonlocal bytes_written
            out_file.write(chunk)
            bytes_written += len(chunk)

        # Retrieve the remote file in binary mode with configured block size
        ftp.retrbinary(
            f"RETR {remote_path}",
            _write_chunk,
            blocksize=cfg.DOWNLOAD_CHUNK_BYTES,
        )

    # Close the FTP connection cleanly
    ftp.quit()

    # Return total bytes written (including any resume offset)
    return bytes_written


def _download_via_http(
    download_url: str,
    partial_path: Path,
    *,
    mode: str,
    resume_offset: int,
    session: requests.Session,
) -> int:
    """Download a file from http(s):// URL using requests streaming."""

    # Build optional Range header when resuming a partial HTTP download
    headers: dict[str, str] = {}
    if resume_offset > 0:
        headers["Range"] = f"bytes={resume_offset}-"

    # Track total bytes written during this download session
    bytes_written = resume_offset

    # Issue a streaming GET to the PRIDE HTTPS URL
    with session.get(
        download_url,
        stream=True,
        timeout=cfg.PRIDE_TIMEOUT,
        headers=headers,
    ) as response:
        # Raise an exception for HTTP error status codes
        response.raise_for_status()

        # Open the partial destination file for binary writing
        with partial_path.open(mode) as out_file:
            # Iterate over response content in DOWNLOAD_CHUNK_BYTES chunks
            for chunk in response.iter_content(chunk_size=cfg.DOWNLOAD_CHUNK_BYTES):
                # Skip empty keep-alive chunks
                if not chunk:
                    continue
                out_file.write(chunk)
                bytes_written += len(chunk)

    # Return total bytes written
    return bytes_written


def download_pride_file(
    download_url: str,
    dest_path: Path | str,
    *,
    expected_size: int | None = None,
    session: requests.Session | None = None,
) -> Path:
    """
    Stream-download a PRIDE file in chunks to dest_path.

    Supports ftp:// via ftplib and http(s):// via requests.
    Returns the absolute Path to the completed file.
    """

    # Normalize dest_path to a Path object
    dest = Path(dest_path)

    # Ensure the parent directory exists before writing
    _ensure_parent_directory(dest)

    # Use a persistent session when provided, otherwise create a one-off session
    http = session or requests.Session()

    # Path for in-progress partial download (appended during streaming)
    partial_path = dest.with_suffix(dest.suffix + ".partial")

    # Track how many bytes were already downloaded for resume
    resume_offset = 0

    # Check whether a partial file exists from a previous interrupted download
    if partial_path.is_file():
        resume_offset = partial_path.stat().st_size
        logger.info("Resuming download at byte %d for %s", resume_offset, dest.name)

    # Open the partial file in append-binary mode when resuming, else write-new
    mode = "ab" if resume_offset > 0 else "wb"

    # Route ftp:// URLs to ftplib (requests has no FTP adapter)
    if download_url.lower().startswith("ftp://"):
        logger.info("Downloading via FTP: %s", download_url)
        _download_via_ftp(
            download_url,
            partial_path,
            mode=mode,
            resume_offset=resume_offset,
        )
    else:
        logger.info("Downloading via HTTP(S): %s", download_url)
        _download_via_http(
            download_url,
            partial_path,
            mode=mode,
            resume_offset=resume_offset,
            session=http,
        )

    # Rename the completed partial file to the final destination name
    partial_path.rename(dest)

    # Read the final file size from the filesystem
    actual_size = dest.stat().st_size

    # Warn when expected_size from PRIDE metadata differs significantly from actual
    if expected_size is not None and expected_size > 0:
        if abs(actual_size - expected_size) > expected_size * 0.01:
            logger.warning(
                "Download size mismatch for %s: expected %d, got %d",
                dest.name,
                expected_size,
                actual_size,
            )

    # Log successful completion with final byte count
    logger.info("Downloaded %s (%d bytes)", dest, actual_size)

    # Return the absolute path to the completed download
    return dest.resolve()


def project_download_path(accession: str, file_name: str) -> Path:
    """Build the standard landing path under DATA_LANDING_DIR/{accession}/."""

    # Join accession subdirectory and original PRIDE file name
    return cfg.DATA_LANDING_DIR / accession / file_name
