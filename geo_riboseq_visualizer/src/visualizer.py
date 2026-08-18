"""Publication-ready plot of GEO Ribo-seq dataset counts by species."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

OTHER_LABEL = "Other"

FONT_TITLE = 18
FONT_AXIS_LABEL = 14
FONT_TICK = 11
FONT_VALUE = 10

DEFAULT_DPI = 300
DEFAULT_FIGSIZE = (13.0, 8.27)  # A4 landscape-ish
BAR_WIDTH = 0.75

LABEL_OFFSET_LOW_PT = 4
LABEL_OFFSET_HIGH_PT = 14
Y_TOP_PADDING_FACTOR = 1.35


def prepare_plot_data(df: pd.DataFrame, top: int = 30) -> pd.DataFrame:
    """Prepare plot data with top-N species plus an aggregated ``Other`` bar.

    Ordering rules:
    - ``Other`` is the first bar (far-left).
    - The top N species follow, ordered from lowest to highest dataset_count.

    Args:
        df: DataFrame with columns ``species`` and ``dataset_count`` (full set).
        top: Number of species to display after ``Other``.

    Returns:
        DataFrame with exactly N+1 rows (unless df has <= N rows), ordered for plotting.

    Raises:
        ValueError: If required columns are missing or ``top < 1``.
    """

    required = {"species", "dataset_count"}
    if not required.issubset(df.columns):
        raise ValueError(f"Input df must contain columns {required}, got {set(df.columns)}")
    if top < 1:
        raise ValueError("--top must be >= 1")

    if df.empty:
        return df.copy()

    sorted_df = df.sort_values("dataset_count", ascending=False).reset_index(drop=True)
    if len(sorted_df) <= top:
        return sorted_df.sort_values("dataset_count", ascending=True).reset_index(drop=True)

    head = sorted_df.head(top).sort_values("dataset_count", ascending=True).reset_index(drop=True)
    other_count = int(sorted_df.iloc[top:]["dataset_count"].sum())
    other_row = pd.DataFrame([{"species": OTHER_LABEL, "dataset_count": other_count}])
    return pd.concat([other_row, head], ignore_index=True)


def _value_label_offset_points(index: int) -> int:
    """Return staggered vertical offset for value labels to reduce collisions."""

    return LABEL_OFFSET_HIGH_PT if index % 2 == 0 else LABEL_OFFSET_LOW_PT


def plot_species_distribution(
    df: pd.DataFrame,
    output_path: str | Path,
    *,
    top: int = 30,
    dpi: int = DEFAULT_DPI,
    figsize: tuple[float, float] = DEFAULT_FIGSIZE,
) -> Path:
    """Create a publication-ready vertical log-scale bar plot.

    Plot rules:
    - X: species (rotation=45, ha='right')
    - Y: dataset_count on log10 scale
    - ``Other`` first (leftmost), then top 30 ascending by count
    - Value labels are horizontal and staggered vertically (odd/even bars)

    Args:
        df: Full counts DataFrame.
        output_path: Output path (png/pdf/svg).
        top: Top species to display after ``Other``.
        dpi: Output resolution for raster formats.
        figsize: Figure size in inches.

    Returns:
        Path to the saved figure file.
    """

    plot_df = prepare_plot_data(df, top=top)
    if plot_df.empty:
        raise ValueError("No data to plot.")

    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=figsize)

    order = plot_df["species"].tolist()
    sns.barplot(
        data=plot_df,
        x="species",
        y="dataset_count",
        order=order,
        width=BAR_WIDTH,
        color="steelblue",
        edgecolor="0.2",
        linewidth=0.6,
        ax=ax,
    )

    ax.set_yscale("log")
    ax.set_title("GEO Ribo-seq datasets per species", fontsize=FONT_TITLE, pad=12)
    ax.set_xlabel("Species", fontsize=FONT_AXIS_LABEL, labelpad=10)
    ax.set_ylabel("Number of datasets (log10 scale)", fontsize=FONT_AXIS_LABEL, labelpad=10)
    ax.tick_params(axis="both", labelsize=FONT_TICK)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # Place value labels above bars with staggered vertical offsets.
    counts = plot_df["dataset_count"].tolist()
    y_max = float(plot_df["dataset_count"].max())
    ax.set_ylim(top=max(1.0, y_max * Y_TOP_PADDING_FACTOR))

    for i, (patch, count) in enumerate(zip(ax.patches, counts, strict=True)):
        x_center = patch.get_x() + patch.get_width() / 2
        bar_top = patch.get_height()
        ax.annotate(
            f"{int(count):,}",
            xy=(x_center, bar_top),
            xytext=(0, _value_label_offset_points(i)),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=FONT_VALUE,
            rotation=0,
            color="0.15",
            clip_on=False,
        )

    fig.tight_layout()

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out

