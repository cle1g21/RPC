"""Stage 4: publication-ready Seaborn validation bar chart."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from config import PipelineConfig

logger = logging.getLogger(__name__)

BAR_A_LABEL = "NTv3 predicted microproteins"
BAR_B_LABEL = "Verified overlap (NTv3 ∩ Nature ESM)"
DEFAULT_DPI = 300
DEFAULT_FIGSIZE = (8.0, 6.0)


def build_summary_counts(ntv3_total: int, matched_total: int) -> pd.DataFrame:
    """Build two-row table for the comparison bar chart.

    Args:
        ntv3_total: Bar A — global NTv3 prediction count (unique coordinates).
        matched_total: Bar B — items in both prediction matrix and validation CSV.

    Returns:
        DataFrame with ``category`` and ``count`` columns.
    """
    return pd.DataFrame(
        {
            "category": [BAR_A_LABEL, BAR_B_LABEL],
            "count": [int(ntv3_total), int(matched_total)],
        }
    )


def plot_validation_bars(
    counts: pd.DataFrame,
    output_path: str | Path,
    *,
    dpi: int = DEFAULT_DPI,
    figsize: tuple[float, float] = DEFAULT_FIGSIZE,
) -> Path:
    """Create whitegrid bar chart with 45° rotated x-labels at 300 DPI.

    Args:
        counts: DataFrame with ``category`` and ``count``.
        output_path: PNG destination.
        dpi: Resolution (default 300 for publication).
        figsize: Figure size in inches.

    Returns:
        Path to saved PNG.
    """
    required = {"category", "count"}
    if not required.issubset(counts.columns) or counts.empty:
        raise ValueError(f"counts must have non-empty {required}")

    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(
        data=counts,
        x="category",
        y="count",
        hue="category",
        palette=["#4C72B0", "#55A868"],
        legend=False,
        edgecolor="0.2",
        linewidth=0.6,
        ax=ax,
    )
    ax.set_title("NTv3 microprotein validation summary", fontsize=16, pad=12)
    ax.set_xlabel("Comparison group", fontsize=13)
    ax.set_ylabel("Number of identified microproteins / peptides", fontsize=13)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    y_max = float(counts["count"].max())
    ax.set_ylim(0, max(1.0, y_max * 1.15))
    for patch, count in zip(ax.patches, counts["count"], strict=True):
        x = patch.get_x() + patch.get_width() / 2
        ax.annotate(
            f"{int(count):,}",
            xy=(x, patch.get_height()),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=11,
        )
    fig.tight_layout()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out


def run_stage4(config: PipelineConfig, summary: dict[str, int] | None = None) -> Path:
    """Run Stage 4 using match_summary.json from Stage 3.

    Args:
        config: Pipeline configuration.
        summary: Optional pre-loaded summary; reads JSON if None.

    Returns:
        Path to ``protein_validation_summary.png``.
    """
    if summary is None:
        if not config.match_summary_path.is_file():
            raise FileNotFoundError(f"Run Stage 3 first: {config.match_summary_path}")
        summary = json.loads(config.match_summary_path.read_text(encoding="utf-8"))

    counts = build_summary_counts(
        ntv3_total=summary["ntv3_total"],
        matched_total=summary["matched_total"],
    )
    out = plot_validation_bars(counts, config.summary_png_path)
    logger.info("Stage 4 complete: %s", out)
    return out
