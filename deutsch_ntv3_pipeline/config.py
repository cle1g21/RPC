"""Global absolute paths and constants for the Deutsch NTv3 validation pipeline."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

# --- IridisX data paths (Deutsch/Kok et al. 2024) ---
DEFAULT_DATA_ROOT = Path("/home/cle1g21/RPC/deutsch_kok_et_al_2024")

# --- DGX / OpenClaw reference paths (Task 2 runs inside sandbox; not invoked by Python) ---
DGX_SANDBOX_NAME = "openclaw-research"
DGX_PROJECT_DIR = "/sandbox/Projects/deltaTE"
DGX_FASTA_STAGING = f"{DGX_PROJECT_DIR}/data/can_nonc_seq.fasta"
DGX_NTV3_OUTPUT_DIR = f"{DGX_PROJECT_DIR}/output"
DGX_NTV3_OUTPUT_FILE = f"{DGX_NTV3_OUTPUT_DIR}/ntv3_predictions.csv"

# --- Default artifact names on IridisX ---
DEFAULT_NTV3_GLOB = "ntv3_predictions*.csv"
MATCHED_CSV_NAME = "matched_peptides.csv"
MATCH_SUMMARY_NAME = "match_summary.json"
SUMMARY_PNG_NAME = "protein_validation_summary.png"


@dataclass(frozen=True)
class PipelineConfig:
    """Resolved file paths for Stages 1, 3, and 4 on IridisX."""

    data_root: Path
    rds_path: Path
    fasta_path: Path
    validation_csv_path: Path
    output_dir: Path
    matched_csv_path: Path
    match_summary_path: Path
    ntv3_output_glob: str
    summary_png_path: Path

    @classmethod
    def from_env(cls) -> PipelineConfig:
        """Build configuration from environment variables and defaults.

        Environment variables:
            PIPELINE_DATA_ROOT: Root folder for input data and FASTA output.
            PIPELINE_OUTPUT_DIR: Folder for NTv3 downloads, matches, and plots.
            NTV3_OUTPUT_GLOB: Glob pattern to find DGX prediction CSV on IridisX.

        Returns:
            Fully resolved pipeline configuration.
        """
        data_root = Path(
            os.environ.get("PIPELINE_DATA_ROOT", str(DEFAULT_DATA_ROOT))
        )
        output_dir = Path(
            os.environ.get("PIPELINE_OUTPUT_DIR", str(data_root / "output"))
        )
        return cls(
            data_root=data_root,
            rds_path=data_root / "can_nonc_seq.RDS",
            fasta_path=data_root / "can_nonc_seq.fasta",
            validation_csv_path=data_root / "41586_2026_10459_MOESM4_ESM.csv",
            output_dir=output_dir,
            matched_csv_path=output_dir / MATCHED_CSV_NAME,
            match_summary_path=output_dir / MATCH_SUMMARY_NAME,
            ntv3_output_glob=os.environ.get("NTV3_OUTPUT_GLOB", DEFAULT_NTV3_GLOB),
            summary_png_path=output_dir / SUMMARY_PNG_NAME,
        )
