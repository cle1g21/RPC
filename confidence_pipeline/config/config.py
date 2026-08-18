"""Centralized configuration for the confidence pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


BASE_DIR = Path("/home/cle1g21/RPC")

NTV3_HANDOFF_DIR = Path("/iridisfs/ddnb/kitsune_labs/predictions/canonical/handoff")

ORF_FILE_ALL = "merged_orfs.tsv"
ORF_FILE_GE30 = "merged_orfs_ge30aa.tsv"

ORF_RUN_LABELS = {
    ORF_FILE_ALL: "all_orfs",
    ORF_FILE_GE30: "ge30aa_orfs",
}

ORF_SOURCE_TIER_LABELS = {
    ORF_FILE_GE30: "Conservative (>=30aa)",
    ORF_FILE_ALL: "Short Fragment (8-29aa)",
}


@dataclass
class PipelineConfig:
    """Resolved paths, column names, and matching rules for the pipeline."""

    base_dir: Path = BASE_DIR
    proteome_cache_path: Path = BASE_DIR / "databases" / "UP000005640_canonical.fasta"
    uniprot_stream_url: str = (
        "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000005640%29"
    )
    refetch_proteome: bool = False
    treat_leucine_isoleucine_as_identical: bool = True
    strip_modifications: bool = True

    ntv3_handoff_dir: Path = NTV3_HANDOFF_DIR
    orf_files: list[str] = field(
        default_factory=lambda: [ORF_FILE_GE30, ORF_FILE_ALL]
    )
    ntv3_filtered_output_dir: Path = BASE_DIR / "NTv3" / "filtered_handoff"
    ntv3_sequence_column: str = "orf_aa_seq"

    instanovo_filtered_dir: Path = BASE_DIR / "instanovo_predictions" / "filtered_predictions"
    instanovo_sequence_column: str = "predictions"
    instanovo_input_glob: str = "*_filtered.csv"

    validation_csv_path: Path = BASE_DIR / "validation" / "41586_2026_10459_MOESM4_ESM.csv"
    validation_sequence_column: str = "sequence"

    confidence_low_dir: Path = BASE_DIR / "confidence_levels" / "low"
    confidence_med_ntv3_dir: Path = BASE_DIR / "confidence_levels" / "med_ntv3"
    confidence_med_instanovo_dir: Path = BASE_DIR / "confidence_levels" / "med_instanovo"
    confidence_high_dir: Path = BASE_DIR / "confidence_levels" / "high"
    confidence_assembled_dir: Path = BASE_DIR / "confidence_levels" / "assembled"
    run_summary_path: Path = BASE_DIR / "confidence_levels" / "run_summary.json"

    tier_match_modes: dict[str, str] = field(
        default_factory=lambda: {
            "low": "substring",
            "med_ntv3": "substring",
            "med_instanovo": "exact",
            "high": "exact",
        }
    )

    tier_suffixes: dict[str, str] = field(
        default_factory=lambda: {
            "low": "",
            "med_ntv3": "_medNTv3",
            "med_instanovo": "_medInstanovo",
            "high": "_high",
        }
    )

    input_globs: tuple[str, ...] = ("*.csv", "*.tsv")

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable summary of key configuration values."""
        return {
            "base_dir": str(self.base_dir),
            "proteome_cache_path": str(self.proteome_cache_path),
            "orf_files": list(self.orf_files),
            "ntv3_handoff_dir": str(self.ntv3_handoff_dir),
            "ntv3_filtered_output_dir": str(self.ntv3_filtered_output_dir),
            "instanovo_filtered_dir": str(self.instanovo_filtered_dir),
            "validation_csv_path": str(self.validation_csv_path),
            "tier_match_modes": dict(self.tier_match_modes),
        }


def default_config() -> PipelineConfig:
    """Return a fresh default pipeline configuration object."""
    return PipelineConfig()
