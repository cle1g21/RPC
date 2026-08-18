"""Publication-ready visualizations of PRIDE species dataset counts."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import LogLocator, ScalarFormatter

OTHER_SPECIES_LABEL = "Other"
DEFAULT_DPI = 300
DEFAULT_BAR_WIDTH_INCHES = 0.45
MIN_FIGURE_WIDTH_INCHES = 12.0
DEFAULT_FIGURE_HEIGHT_INCHES = 8.27

X_LABEL_ROTATION_DEGREES = 45
BAR_VALUE_ROTATION_DEGREES = 0

FONT_TITLE = 18
FONT_AXIS_LABEL = 14
FONT_TICK = 11
FONT_BAR_VALUE = 10

Y_TOP_PADDING_FACTOR = 2.2
BAR_WIDTH = 0.80
BAR_LABEL_OFFSET_PT = 8

def prepare_plot_data(
    df: pd.DataFrame,
    top: int = 30,
    *,
    plot_all: bool = False,
) -> pd.DataFrame:
    """Subset species data for plotting with an aggregated ``Other`` bar.

    ``Other`` is placed first (far left); the top species follow in ascending
    order of dataset count.
    """
    if df.empty:
        return df.copy()

    sorted_df = df.sort_values("dataset_count", ascending=False).reset_index(drop=True)

    if plot_all:
        return sorted_df.sort_values("dataset_count", ascending=True).reset_index(drop=True)

    if top < 1:
        raise ValueError("--top must be at least 1 when not using --plot-all")

    head = sorted_df.head(top).sort_values("dataset_count", ascending=True)
    remainder = sorted_df.iloc[top:]

    if remainder.empty:
        return head.reset_index(drop=True)

    other_count = int(remainder["dataset_count"].sum())
    other_row = pd.DataFrame(
        [{"species": OTHER_SPECIES_LABEL, "dataset_count": other_count}]
    )
    return pd.concat([other_row, head], ignore_index=True)


def _figure_size(n_bars: int, figsize: tuple[float, float] | None) -> tuple[float, float]:
    """Compute figure dimensions from the number of vertical bars."""
    if figsize is not None:
        return figsize
    width = max(MIN_FIGURE_WIDTH_INCHES, DEFAULT_BAR_WIDTH_INCHES * n_bars)
    return (width, DEFAULT_FIGURE_HEIGHT_INCHES)



def plot_species_distribution(
    df: pd.DataFrame,
    output_path: str | Path,
    *,
    top: int = 30,
    plot_all: bool = False,
    figsize: tuple[float, float] | None = None,
    dpi: int = DEFAULT_DPI,
    label_rotation: int = X_LABEL_ROTATION_DEGREES,
) -> Path:
    """Create a vertical log-scaled bar chart of datasets per species.

    Horizontal bar value labels are placed uniformly above each bar.
    """
    plot_df = prepare_plot_data(df, top=top, plot_all=plot_all)
    if plot_df.empty:
        raise ValueError("No species data available to plot")

    n_bars = len(plot_df)
    width, height = _figure_size(n_bars, figsize)
    species_order = plot_df["species"].tolist()

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.0)
    fig, ax = plt.subplots(figsize=(width, height))

    sns.barplot(
        data=plot_df,
        x="species",
        y="dataset_count",
        ax=ax,
        order=species_order,
        width=BAR_WIDTH,
        color="steelblue",
        edgecolor="0.2",
        linewidth=0.6,
    )

    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.get_major_formatter().set_scientific(False)

    ax.set_title("PRIDE Archive datasets per species", fontsize=FONT_TITLE, pad=12)
    ax.set_ylabel(
        "Number of datasets (log10 scale)",
        fontsize=FONT_AXIS_LABEL,
        labelpad=10,
    )
    ax.set_xlabel("Species", fontsize=FONT_AXIS_LABEL, labelpad=12)

    ax.tick_params(axis="y", labelsize=FONT_TICK)
    plt.setp(
        ax.get_xticklabels(),
        rotation=label_rotation,
        ha="right",
        rotation_mode="anchor",
        fontsize=FONT_TICK,
    )

    y_min = max(1, plot_df["dataset_count"].min() * 0.6)
    y_max = float(plot_df["dataset_count"].max())
    y_top = y_max * Y_TOP_PADDING_FACTOR

    counts = plot_df["dataset_count"].tolist()
    for patch, count in zip(ax.patches, counts, strict=True):
        bar_top = patch.get_height()
        x_center = patch.get_x() + patch.get_width() / 2
        ax.annotate(
            f"{int(count):,}",
            xy=(x_center, bar_top),
            xytext=(0, BAR_LABEL_OFFSET_PT),
            textcoords="offset points",
            va="bottom",
            ha="center",
            rotation=BAR_VALUE_ROTATION_DEGREES,
            fontsize=FONT_BAR_VALUE,
            color="0.15",
            clip_on=False,
        )

    ax.set_ylim(bottom=y_min, top=y_top)
    ax.margins(x=0.02, y=0.06)

    fig.tight_layout()

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)

    return out
