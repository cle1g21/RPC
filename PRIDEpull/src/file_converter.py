"""msconvert subprocess wrapper and streaming MGF subset slicer."""

# Enable postponed evaluation of type annotations
from __future__ import annotations

# Import logging for conversion progress messages
import logging

# Import os for file size checks and path operations
import os

# Import subprocess to invoke msconvert synchronously
import subprocess

# Import sys so the same Python interpreter runs convert_ms.py
import sys

# Import threading for optional conversion lock (used by pipeline_runner)
import threading

# Import Path for MGF path handling
from pathlib import Path

# Import config for msconvert paths and subset parameters
from config import config as cfg

# Create a module-level logger
logger = logging.getLogger(__name__)

# Module-level lock ensuring only one msconvert runs at a time when enabled
_conversion_lock = threading.Lock()


def get_input_stem(input_path: Path) -> str:
    """Return the filename stem used by msconvert for output naming."""

    # Take basename and strip the last extension (.raw, .d, etc.)
    return input_path.stem


def expected_mgf_output_path(input_path: Path, output_dir: Path) -> Path:
    """Build the path where msconvert should write the converted MGF file."""

    # msconvert names output {input_stem}.mgf in the output directory
    return output_dir / f"{get_input_stem(input_path)}.mgf"


def convert_raw_to_mgf(
    raw_path: Path | str,
    output_dir: Path | str | None = None,
) -> Path:
    """
    Invoke /RPC/msconvert/convert_ms.py synchronously to produce an MGF file.

    When SEQUENTIAL_CONVERSION is True, acquires a module-level lock so only
    one conversion runs at a time across the pipeline process.
    """

    # Normalize raw_path to a Path object
    raw = Path(raw_path)

    # Default output directory to the same folder as the raw input file
    if output_dir is None:
        out_dir = raw.parent
    else:
        out_dir = Path(output_dir)

    # Ensure the output directory exists on disk
    out_dir.mkdir(parents=True, exist_ok=True)

    # Compute the expected MGF output path before running conversion
    expected_mgf = expected_mgf_output_path(raw, out_dir)

    # Define the inner conversion function that runs the subprocess
    def _run_conversion() -> Path:
        # Build the subprocess command list pointing at convert_ms.py
        command = [
            sys.executable,
            str(cfg.MSCONVERT_SCRIPT),
            "--input",
            str(raw.resolve()),
            "--output-dir",
            str(out_dir.resolve()),
            "--format",
            "mgf",
            "--image",
            str(cfg.MSCONVERT_IMAGE),
        ]

        # Log the exact command for operator reproducibility
        logger.info("Running msconvert: %s", " ".join(command))

        # Execute msconvert synchronously and capture stdout/stderr
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
        )

        # Log stdout from msconvert when non-empty
        if completed.stdout:
            logger.info("msconvert stdout:\n%s", completed.stdout)

        # Log stderr from msconvert when non-empty
        if completed.stderr:
            logger.warning("msconvert stderr:\n%s", completed.stderr)

        # Raise when msconvert returned a non-zero exit code
        if completed.returncode != 0:
            raise RuntimeError(
                f"msconvert failed with exit code {completed.returncode}"
            )

        # Verify the expected MGF file appeared on disk
        if not expected_mgf.is_file():
            raise FileNotFoundError(
                f"msconvert completed but MGF not found: {expected_mgf}"
            )

        # Log the output MGF path and file size
        logger.info(
            "Conversion complete: %s (%d bytes)",
            expected_mgf,
            expected_mgf.stat().st_size,
        )

        # Return the path to the converted MGF file
        return expected_mgf

    # Acquire the conversion lock when sequential conversion is enabled
    if cfg.SEQUENTIAL_CONVERSION:
        with _conversion_lock:
            return _run_conversion()

    # Run conversion without locking when sequential mode is disabled
    return _run_conversion()


def slice_mgf(
    input_path: Path | str,
    output_path: Path | str,
    max_spectra: int | None = None,
) -> Path:
    """
    Stream-read an MGF file and write the first max_spectra spectrum blocks.

    Never loads the full file into memory. Each spectrum block spans
    BEGIN IONS through END IONS inclusive.
    """

    # Default max_spectra from config when not specified
    if max_spectra is None:
        max_spectra = cfg.SUBSET_SPECTRUM_COUNT

    # Normalize input and output paths
    src = Path(input_path)
    dst = Path(output_path)

    # Ensure the output parent directory exists
    dst.parent.mkdir(parents=True, exist_ok=True)

    # Counter for completed spectrum blocks written to the subset file
    spectra_written = 0

    # Flag indicating we are inside a spectrum block between BEGIN and END IONS
    inside_spectrum = False

    # Buffer lines for the current spectrum block being collected
    current_block: list[str] = []

    # Open source MGF for reading and destination for writing (streaming)
    with src.open(encoding="utf-8", errors="replace") as in_file, dst.open(
        "w", encoding="utf-8"
    ) as out_file:
        # Read the source file one line at a time to bound memory usage
        for line in in_file:
            # Strip trailing newline for consistent block handling
            stripped = line.rstrip("\n\r")

            # Detect start of a new spectrum block
            if stripped.strip() == "BEGIN IONS":
                # Mark that we are collecting lines for one spectrum
                inside_spectrum = True

                # Start a fresh block buffer with the BEGIN IONS line
                current_block = [stripped]

                # Continue to next source line without writing yet
                continue

            # Accumulate lines while inside an active spectrum block
            if inside_spectrum:
                # Append this line to the current block buffer
                current_block.append(stripped)

                # Detect end of the current spectrum block
                if stripped.strip() == "END IONS":
                    # Write the complete block to the subset output file
                    for block_line in current_block:
                        out_file.write(block_line + "\n")

                    # Increment the count of spectra written so far
                    spectra_written += 1

                    # Reset state after closing this spectrum block
                    inside_spectrum = False
                    current_block = []

                    # Stop reading when we have reached the requested spectrum count
                    if spectra_written >= max_spectra:
                        break

                # Continue to next source line
                continue

            # Copy header lines (before first BEGIN IONS) only when no spectra yet
            if spectra_written == 0 and not inside_spectrum:
                out_file.write(stripped + "\n")

    # Log how many spectra were written to the subset file
    logger.info(
        "Wrote subset MGF %s with %d spectra from %s",
        dst,
        spectra_written,
        src.name,
    )

    # Return the path to the subset MGF file
    return dst.resolve()


def count_mgf_spectra(mgf_path: Path | str) -> int:
    """Count BEGIN IONS blocks in an MGF file without loading it fully."""

    # Normalize path
    path = Path(mgf_path)

    # Counter for spectrum blocks
    count = 0

    # Stream-read line by line
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip() == "BEGIN IONS":
                count += 1

    # Return total spectrum count
    return count
