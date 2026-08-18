"""InstaNovo launcher: GPU Slurm jobs (default) or local direct subprocess fallback."""

# Enable postponed evaluation of type annotations
from __future__ import annotations

# Import logging for prediction progress messages
import logging

# Import subprocess to invoke instanovo predict or the Slurm launcher script
import subprocess

# Import shutil to copy Slurm launcher output into canonical predictions dir
import shutil

# Import Path for MGF and CSV path handling
from pathlib import Path

# Import config for InstaNovo paths and routing flags
from config import config as cfg

# Create a module-level logger
logger = logging.getLogger(__name__)


def mgf_size_gb(mgf_path: Path) -> float:
    """Return the MGF file size in gigabytes."""

    # Read byte size from filesystem and convert to GB
    return mgf_path.stat().st_size / (1024**3)


def build_output_csv_path(accession: str, mgf_path: Path) -> Path:
    """Build prediction CSV path aligned with PRIDE accession and MGF stem."""

    # Use MGF stem for disambiguation (includes _subset500 when applicable)
    stem = mgf_path.stem

    # Compose filename: {accession}_{stem}_predictions.csv
    csv_name = f"{accession}_{stem}_predictions.csv"

    # Return full path under PREDICTIONS_OUTPUT_DIR
    return cfg.PREDICTIONS_OUTPUT_DIR / csv_name


def slurm_default_output_csv(mgf_path: Path) -> Path:
    """
    Return the CSV path that run_until_complete.sh writes by convention.

    The launcher writes to instanovo_predictions/{stem}_predictions.csv
    (not the predictions/ subfolder and without PXD accession prefix).
    """

    # Extract MGF stem for default output naming
    stem = mgf_path.stem

    # Build path matching run_until_complete.sh FINAL_OUTPUT variable
    return cfg.PREDICTIONS_OUTPUT_DIR.parent / f"{stem}_predictions.csv"


def run_direct_predict(mgf_path: Path, output_csv: Path) -> Path:
    """
    Run instanovo predict as a local subprocess (login-node CPU fallback).

    Only used when INSTANOVO_USE_SLURM=False. Not recommended on Iridis.
    """

    # Ensure the predictions output directory exists
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    # Build the instanovo predict command with Hydra-style overrides
    command = [
        str(cfg.INSTANOVO_BIN),
        "predict",
        "--data-path",
        str(mgf_path.resolve()),
        "--output-path",
        str(output_csv.resolve()),
        f"batch_size={cfg.INSTANOVO_BATCH_SIZE}",
        f"num_workers={cfg.INSTANOVO_NUM_WORKERS}",
        "stream_predictions=true",
    ]

    # Log the exact command for operator reproducibility
    logger.info("Running local InstaNovo predict (no Slurm): %s", " ".join(command))

    # Execute instanovo predict with working directory set to InstaNovo root
    completed = subprocess.run(
        command,
        cwd=str(cfg.INSTANOVO_ROOT),
        capture_output=True,
        text=True,
        check=False,
    )

    # Log stdout when present
    if completed.stdout:
        logger.info("instanovo stdout:\n%s", completed.stdout)

    # Log stderr when present
    if completed.stderr:
        logger.warning("instanovo stderr:\n%s", completed.stderr)

    # Raise when instanovo returned a non-zero exit code
    if completed.returncode != 0:
        raise RuntimeError(
            f"instanovo predict failed with exit code {completed.returncode}"
        )

    # Verify the output CSV was created
    if not output_csv.is_file():
        raise FileNotFoundError(f"Prediction CSV not found: {output_csv}")

    # Log successful completion
    logger.info("Local prediction complete: %s", output_csv)

    # Return path to the output CSV
    return output_csv


def run_slurm_predict(mgf_path: Path, output_csv: Path) -> Path:
    """
    Submit InstaNovo via run_until_complete.sh and wait for GPU job completion.

    For MGF < 1 GB (including subset files): one direct GPU sbatch job.
    For MGF >= 1 GB: chunked parallel GPU jobs with compile step.
    """

    # Verify the Slurm launcher script exists on the cluster filesystem
    if not cfg.INSTANOVO_SLURM_SCRIPT.is_file():
        raise FileNotFoundError(
            f"Slurm launcher not found: {cfg.INSTANOVO_SLURM_SCRIPT}"
        )

    # Build command: bash run_until_complete.sh /path/to/file.mgf
    command = [
        "bash",
        str(cfg.INSTANOVO_SLURM_SCRIPT),
        str(mgf_path.resolve()),
    ]

    # Log the Slurm launcher invocation and expected routing
    size_gb = mgf_size_gb(mgf_path)
    route = "chunked GPU jobs" if size_gb >= cfg.INSTANOVO_LARGE_FILE_GB else "single GPU job"
    logger.info(
        "Submitting InstaNovo via Slurm (%s, %.3f GB): %s",
        route,
        size_gb,
        " ".join(command),
    )

    # Execute the Slurm-aware launcher; it blocks until predictions complete
    completed = subprocess.run(
        command,
        cwd=str(cfg.INSTANOVO_ROOT),
        capture_output=True,
        text=True,
        check=False,
    )

    # Log stdout when present (includes sbatch job IDs and progress cycles)
    if completed.stdout:
        logger.info("Slurm launcher stdout:\n%s", completed.stdout)

    # Log stderr when present
    if completed.stderr:
        logger.warning("Slurm launcher stderr:\n%s", completed.stderr)

    # Raise when the launcher returned a non-zero exit code
    if completed.returncode != 0:
        raise RuntimeError(
            f"Slurm launcher failed with exit code {completed.returncode}"
        )

    # Resolve where run_until_complete.sh wrote the predictions CSV
    slurm_output = slurm_default_output_csv(mgf_path)

    # Prefer the Slurm default path, then any pre-existing canonical path
    if slurm_output.is_file():
        final_csv = slurm_output
    elif output_csv.is_file():
        final_csv = output_csv
    else:
        raise FileNotFoundError(
            f"Prediction CSV not found after Slurm run: {slurm_output} or {output_csv}"
        )

    # Copy into canonical predictions/ dir with PXD accession prefix when needed
    if final_csv.resolve() != output_csv.resolve():
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(final_csv, output_csv)
        logger.info("Copied Slurm output %s → %s", final_csv, output_csv)
        final_csv = output_csv

    # Log successful completion
    logger.info("Slurm prediction complete: %s", final_csv)

    # Return path to the canonical output CSV
    return final_csv


def run_prediction(
    accession: str,
    mgf_path: Path | str,
    *,
    run_subset_only: bool | None = None,
) -> Path:
    """
    Run InstaNovo predictions on the given MGF file.

    Default (INSTANOVO_USE_SLURM=True): submits GPU Slurm jobs via
    run_until_complete.sh for both subset and full runs.

    Returns the path to the written predictions CSV.
    """

    # Normalize mgf_path to Path
    mgf = Path(mgf_path)

    # Resolve run_subset_only from config when not overridden (used for logging only)
    if run_subset_only is None:
        run_subset_only = cfg.RUN_SUBSET_ONLY

    # Build the canonical output CSV path for this accession and MGF
    output_csv = build_output_csv_path(accession, mgf)

    # Route all runs through Slurm GPU jobs when configured (recommended on Iridis)
    if cfg.INSTANOVO_USE_SLURM:
        mode_label = "subset" if run_subset_only else "full"
        logger.info("InstaNovo mode: %s via Slurm GPU", mode_label)
        return run_slurm_predict(mgf, output_csv)

    # Fallback: local subprocess on the current machine (CPU on login nodes)
    logger.warning(
        "INSTANOVO_USE_SLURM=False: running instanovo predict locally (not recommended on Iridis)"
    )
    return run_direct_predict(mgf, output_csv)
