"""Stage 3: match NTv3 DGX predictions to Nature validation CSV by coordinates."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd

from config import PipelineConfig

logger = logging.getLogger(__name__)

COORD_ALIASES: dict[str, tuple[str, ...]] = {
    "chr": ("chr", "chrm", "chrom", "chromosome"),
    "start": ("start", "starts", "start_pos"),
    "end": ("end", "ends", "end_pos"),
}


class MatcherError(RuntimeError):
    """Raised when coordinate matching fails."""


def load_ntv3_output(path: Path) -> pd.DataFrame:
    """Load NTv3 prediction file from IridisX output directory.

    Args:
        path: CSV or TSV from DGX (e.g. ``ntv3_predictions.csv``).

    Returns:
        DataFrame with normalized ``chr``, ``start``, ``end`` columns.
    """
    if not path.is_file():
        raise FileNotFoundError(f"NTv3 output not found: {path}")
    df = _read_table(path)
    return normalize_coordinates(df, source_label="ntv3")


def load_validation_csv(path: Path) -> pd.DataFrame:
    """Load Nature supplementary validation table.

    Args:
        path: ``41586_2026_10459_MOESM4_ESM.csv``.

    Returns:
        DataFrame with normalized coordinates.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Validation CSV not found: {path}")
    return normalize_coordinates(pd.read_csv(path), source_label="validation")


def normalize_coordinates(df: pd.DataFrame, *, source_label: str) -> pd.DataFrame:
    """Map alias column names to canonical chr/start/end strings."""
    out = df.copy()
    lower_cols = {str(c).lower(): c for c in out.columns}
    rename: dict[str, str] = {}
    for canonical, aliases in COORD_ALIASES.items():
        found = next((lower_cols[a] for a in aliases if a in lower_cols), None)
        if found is None:
            raise MatcherError(
                f"{source_label}: missing '{canonical}' (columns: {list(out.columns)})"
            )
        rename[found] = canonical
    out = out.rename(columns=rename)
    for col in ("chr", "start", "end"):
        out[col] = out[col].map(_coord_to_str)
    return out


def match_peptides(ntv3_df: pd.DataFrame, validation_df: pd.DataFrame) -> pd.DataFrame:
    """Inner-join on identical chr, start, end (high-confidence coordinate match)."""
    return pd.merge(
        ntv3_df,
        validation_df,
        on=["chr", "start", "end"],
        how="inner",
        suffixes=("_ntv3", "_validation"),
    )


def count_coordinate_overlap(ntv3_df: pd.DataFrame, validation_df: pd.DataFrame) -> int:
    """Count unique coordinate triplets present in both datasets."""
    ntv3_keys = set(zip(ntv3_df["chr"], ntv3_df["start"], ntv3_df["end"], strict=True))
    val_keys = set(zip(validation_df["chr"], validation_df["start"], validation_df["end"], strict=True))
    return len(ntv3_keys & val_keys)


def run_stage3(config: PipelineConfig, ntv3_path: Path | None = None) -> dict[str, int]:
    """Run Stage 3 and write matched_peptides.csv plus match_summary.json.

    Args:
        config: Pipeline configuration.
        ntv3_path: Optional explicit NTv3 CSV; auto-discovered from output_dir if None.

    Returns:
        Summary counts for Stage 4 plotting.
    """
    config.output_dir.mkdir(parents=True, exist_ok=True)
    if ntv3_path is None:
        ntv3_path = _resolve_ntv3_output(config)

    ntv3_df = load_ntv3_output(ntv3_path)
    validation_df = load_validation_csv(config.validation_csv_path)
    matched = match_peptides(ntv3_df, validation_df)
    matched.to_csv(config.matched_csv_path, index=False)

    summary = {
        "ntv3_total": int(ntv3_df.drop_duplicates(subset=["chr", "start", "end"]).shape[0]),
        "validation_total": int(
            validation_df.drop_duplicates(subset=["chr", "start", "end"]).shape[0]
        ),
        "matched_total": count_coordinate_overlap(ntv3_df, validation_df),
        "matched_rows": int(len(matched)),
    }
    config.match_summary_path.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    logger.info(
        "Stage 3: %d coordinate overlaps (%d merged rows) → %s",
        summary["matched_total"],
        summary["matched_rows"],
        config.matched_csv_path,
    )
    return summary


def _resolve_ntv3_output(config: PipelineConfig) -> Path:
    matches = sorted(config.output_dir.glob(config.ntv3_output_glob))
    if not matches:
        raise FileNotFoundError(
            f"No file matching {config.ntv3_output_glob} in {config.output_dir}. "
            "Complete Task 2 on DGX (see openclaw_dgx_blueprint.md) and copy output here."
        )
    return matches[0]


def _read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    try:
        return pd.read_csv(path)
    except pd.errors.ParserError:
        return pd.read_csv(path, sep="\t")


def _coord_to_str(value: object) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()
