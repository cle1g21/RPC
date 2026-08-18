"""Publication-ready Seaborn KDE density plotting with linear and log panels."""

# Enable postponed evaluation of annotations used in plot helper signatures
from __future__ import annotations

# Import logging so figure export paths can be reported to the console
import logging
# Import Path for typed figure output destinations
from pathlib import Path

# Import numpy to locate density heights near the annotation base-pair mark
import numpy as np
# Import matplotlib pyplot for figure handles, spines, and annotation
import matplotlib.pyplot as plt
# Import pandas for the tidy multi-species length table
import pandas as pd
# Import seaborn as the primary statistical plotting interface
import seaborn as sns

# Import centralized plot style, annotation, path, and label constants
from config.config import (
    ANNOTATION_BP,
    ANNOTATION_LABEL,
    ANNOTATION_Y_TEXT_FACTOR_LINEAR,
    ANNOTATION_Y_TEXT_FACTOR_LOG,
    BW_ADJUST,
    DPI,
    DUAL_KDE_FIGSIZE,
    DUAL_KDE_STEM,
    KDE_FIGSIZE,
    KDE_STEM,
    KDE_TITLE,
    LOG_X_LIM_MAX,
    LOG_X_LIM_MIN,
    OUTPUT_DIR,
    PANEL_A_TITLE,
    PANEL_B_TITLE,
    SNS_PALETTE,
    SNS_STYLE,
    SPECIES_PALETTE,
    SPINE_COLOR,
    SPINE_LINEWIDTH,
    X_AXIS_LABEL,
    X_LIM_MAX,
    X_LIM_MIN,
    Y_AXIS_LABEL,
)


# Create a module-level logger for plotting status messages
logger = logging.getLogger(__name__)


def apply_publication_style() -> None:
    """Apply Seaborn ticks theme defaults used by the KDE figures."""
    # Enforce the ticks theme and muted palette for a minimal schematic look
    sns.set_theme(style=SNS_STYLE, palette=SNS_PALETTE)


def save_figure(fig: plt.Figure, stem: str, output_dir: Path | None = None) -> list[Path]:
    """Export a Matplotlib figure as both high-resolution PNG and vector SVG."""
    # Default to the configured outputs directory when no override is provided
    target_dir = Path(output_dir) if output_dir is not None else OUTPUT_DIR
    # Ensure the output directory exists before writing figure files
    target_dir.mkdir(parents=True, exist_ok=True)
    # Build the PNG destination path from the shared output stem
    png_path = target_dir / f"{stem}.png"
    # Build the SVG destination path from the shared output stem
    svg_path = target_dir / f"{stem}.svg"
    # Write a 300 DPI raster PNG suitable for manuscript submission
    fig.savefig(png_path, dpi=DPI, bbox_inches="tight")
    # Write a vector SVG suitable for camera-ready editing
    fig.savefig(svg_path, bbox_inches="tight")
    # Close the figure to release Matplotlib memory after export
    plt.close(fig)
    # Log both export destinations for the caller
    logger.info("Saved figures to %s and %s", png_path, svg_path)
    # Return the list of written file paths for the entry-point summary
    return [png_path, svg_path]


def _style_minimal_spines(ax: plt.Axes) -> None:
    """Remove top/right spines and thicken the remaining black axis spines."""
    # Drop the top and right borders for a clean schematic publication look
    sns.despine(ax=ax, top=True, right=True)
    # Iterate over the remaining left and bottom spines to thicken them
    for spine_name in ("left", "bottom"):
        # Select the named spine artist from the axes
        spine = ax.spines[spine_name]
        # Set the spine stroke width to the configured publication thickness
        spine.set_linewidth(SPINE_LINEWIDTH)
        # Force the spine color to solid black for hand-drawn style contrast
        spine.set_color(SPINE_COLOR)
    # Remove background gridlines for a minimal graphic layout
    ax.grid(False)


def _annotation_y_position(ax: plt.Axes, annotation_bp: float) -> float:
    """Estimate a y position near the KDE curves at the annotation x-coordinate."""
    # Initialize the peak density observed near the annotation base-pair mark
    peak_density = 0.0
    # Inspect each plotted line artist contributed by seaborn's kdeplot
    for line in ax.get_lines():
        # Read the x coordinates sampled along the current KDE curve
        x_data = np.asarray(line.get_xdata(), dtype=float)
        # Read the y density values sampled along the current KDE curve
        y_data = np.asarray(line.get_ydata(), dtype=float)
        # Skip empty artists that do not encode a density curve
        if x_data.size == 0 or y_data.size == 0:
            # Continue scanning the remaining line artists
            continue
        # Locate the sampled index closest to the annotation base-pair value
        nearest_index = int(np.argmin(np.abs(x_data - annotation_bp)))
        # Track the maximum density across species at that nearby sample
        peak_density = max(peak_density, float(y_data[nearest_index]))
    # Fall back to a fraction of the y-axis maximum when no line data exists
    if peak_density <= 0.0:
        # Read the current y-axis upper limit as a fallback height reference
        _, y_max = ax.get_ylim()
        # Place the arrow at 60% of the axis height when curve samples are missing
        return float(y_max) * 0.6
    # Return the highest nearby density so the arrow points onto the curve stack
    return peak_density


def annotate_inflection(
    ax: plt.Axes,
    *,
    annotation_bp: float = ANNOTATION_BP,
    label: str = ANNOTATION_LABEL,
    y_text_factor: float = ANNOTATION_Y_TEXT_FACTOR_LINEAR,
) -> None:
    """Draw a vertical annotation arrow pointing down to the key bp inflection."""
    # Estimate the y height on the KDE curves nearest the annotation x position
    y_curve = _annotation_y_position(ax, annotation_bp)
    # Read the current y-axis upper limit to place the label above the curves
    _, y_max = ax.get_ylim()
    # Position the text label above the curves using the panel-specific factor
    y_text = max(y_curve * y_text_factor, y_max * 0.85)
    # Draw a downward arrow annotated with the configured 1,500 bps label
    ax.annotate(
        # Place the human-readable inflection label above the arrow
        label,
        # Anchor the arrow tip at the annotation bp on the density curves
        xy=(annotation_bp, y_curve),
        # Place the text directly above the arrow tip at the same x position
        xytext=(annotation_bp, y_text),
        # Center the label horizontally over the arrow shaft
        ha="center",
        # Anchor the text from its bottom edge so it sits above the arrow
        va="bottom",
        # Use a clean sans-serif weight for publication typography
        fontsize=11,
        # Emphasize the inflection label with bold font weight
        fontweight="bold",
        # Configure a simple black downward arrow matching the schematic style
        arrowprops={
            # Use a standard arrow style pointing toward the curve
            "arrowstyle": "->",
            # Color the arrow shaft and head black for contrast
            "color": "black",
            # Set a visible but not overpowering arrow stroke width
            "lw": 1.5,
        },
    )


def _place_legend(ax: plt.Axes) -> None:
    """Reposition a seaborn legend at the top-right, overlaid on the plot."""
    # Retrieve any legend created automatically by seaborn's hue mapping
    legend = ax.get_legend()
    # Exit early when this axes has no legend to reposition
    if legend is None:
        # Leave the axes unchanged when no legend exists
        return
    # Capture the existing species line handles from the seaborn legend
    handles = legend.legend_handles
    # Capture the existing species display-name labels from the seaborn legend
    labels = [text.get_text() for text in legend.texts]
    # Replace the legend with an overlaid frame in the top-right of the axes
    legend = ax.legend(
        # Reuse the existing species line handles
        handles=handles,
        # Reuse the existing species display-name labels
        labels=labels,
        # Anchor the legend to the upper-right corner of the axes
        loc="upper right",
        # Keep one species per row for a clean vertical legend
        ncol=1,
        # Draw a frame so a translucent grey background can be applied
        frameon=True,
        # Keep legend typography at the requested 12 pt size
        fontsize=12,
        # Drop any residual legend title text
        title=None,
        # Slightly reduce padding so the frame sits tightly around the labels
        borderpad=0.6,
    )
    # Access the legend frame patch to style its background
    frame = legend.get_frame()
    # Set a neutral grey face color for the legend background box
    frame.set_facecolor("0.85")
    # Apply partial transparency so underlying curves remain visible (alpha 0.2)
    frame.set_alpha(0.2)
    # Match the grey edge to the face for a soft, low-contrast border
    frame.set_edgecolor("0.7")
    # Keep the border stroke thin so it does not dominate the plot
    frame.set_linewidth(0.8)


def _draw_kde_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    *,
    title: str,
    x_lim: tuple[float, float],
    log_scale: bool = False,
    legend: bool = True,
    annotate: bool = True,
    y_text_factor: float = ANNOTATION_Y_TEXT_FACTOR_LINEAR,
) -> None:
    """Draw one KDE panel with spines, labels, optional annotation, and legend."""
    # Copy the frame so optional log10 transforms do not mutate caller data
    plot_df = df.copy()
    # Default to plotting raw coding base-pair lengths on a linear axis
    x_column = "bp_length"
    # Default annotation x-position remains the configured 1,500 bp mark
    annotation_x = float(ANNOTATION_BP)
    # Default clip window matches the caller-provided bp limits
    clip_window = x_lim
    # Default axis limits remain in raw base-pair units
    axis_limits = x_lim
    # Switch into log10(bp) space when Panel B requests logarithmic resolution
    if log_scale:
        # Transform coding lengths into log10 base-pair units for true log-space KDE
        plot_df["log10_bp"] = np.log10(plot_df["bp_length"].astype(float))
        # Point seaborn at the transformed log10 length column
        x_column = "log10_bp"
        # Convert the annotation mark into log10 space (log10(1500) ≈ 3.18)
        annotation_x = float(np.log10(ANNOTATION_BP))
        # Convert the visible bp window into matching log10 clip bounds
        clip_window = (float(np.log10(x_lim[0])), float(np.log10(x_lim[1])))
        # Convert the visible bp window into matching log10 axis limits
        axis_limits = clip_window
    # Draw continuous KDE curves for every species onto the supplied axes
    sns.kdeplot(
        # Supply the (optionally log-transformed) multi-species length table
        data=plot_df,
        # Place either raw bp or log10(bp) values on the x-axis
        x=x_column,
        # Distinguish the five species with separate hue-colored curves
        hue="species",
        # Map each species display name to its configured hex color
        palette=SPECIES_PALETTE,
        # Normalize each species independently so shapes remain comparable
        common_norm=False,
        # Tighten bandwidth so secondary modes and inflections remain visible
        bw_adjust=BW_ADJUST,
        # Use a filled-off line style matching schematic density overlays
        fill=False,
        # Set a readable stroke width for each species density curve
        linewidth=2.0,
        # Clip evaluation to the panel-specific visible window in plot units
        clip=clip_window,
        # Bind this KDE rendering to the caller-provided subplot axes
        ax=ax,
        # Suppress the automatic legend when a shared dual-panel legend is used
        legend=legend,
    )
    # Restrict the visible x-axis to the panel-specific window in plot units
    ax.set_xlim(*axis_limits)
    # Ensure the density axis starts at zero for a clean abundance baseline
    ax.set_ylim(bottom=0)
    # Replace log10 tick positions with human-readable bp labels on log panels
    if log_scale:
        # Choose representative ORF-length tick marks spanning the log window
        tick_bps = [100, 300, 500, 1000, 1500, 3000, 8000]
        # Keep only ticks that fall inside the configured positive log window
        tick_bps = [bp for bp in tick_bps if x_lim[0] <= bp <= x_lim[1]]
        # Place ticks at the corresponding log10(bp) coordinates
        ax.set_xticks([np.log10(bp) for bp in tick_bps])
        # Label ticks with the original base-pair values for interpretability
        ax.set_xticklabels([f"{bp:,}" for bp in tick_bps])
    # Label the x-axis as continuous ORF length in base pairs
    ax.set_xlabel(X_AXIS_LABEL)
    # Label the y-axis as kernel density abundance
    ax.set_ylabel(Y_AXIS_LABEL)
    # Set the panel title describing linear or logarithmic scaling
    ax.set_title(title)
    # Apply minimalist thick black spines and remove residual gridlines
    _style_minimal_spines(ax)
    # Optionally draw the vertical annotation arrow at the 1,500 bp mark
    if annotate:
        # Draw the inflection annotation using panel-specific x units and offsets
        annotate_inflection(
            # Target the current panel axes for the annotation overlay
            ax,
            # Pass the panel-specific annotation x in the same units as the curves
            annotation_bp=annotation_x,
            # Keep the human-readable 1,500 bps label text unchanged
            label=ANNOTATION_LABEL,
            # Use the panel-specific vertical text offset factor
            y_text_factor=y_text_factor,
        )
    # Reposition the species legend outside the panel when legend is enabled
    if legend:
        # Move the seaborn hue legend to the exterior of this axes
        _place_legend(ax)


def plot_kde_density(
    df: pd.DataFrame,
    output_stem: str = KDE_STEM,
    output_dir: Path | None = None,
) -> list[Path]:
    """Render a single-panel linear multi-species KDE density plot."""
    # Apply the shared Seaborn ticks theme before drawing axes
    apply_publication_style()
    # Create a Matplotlib figure and axes with the configured single-panel size
    fig, ax = plt.subplots(figsize=KDE_FIGSIZE)
    # Draw the linear-scale KDE panel with legend and without the 1,500 bp arrow
    _draw_kde_panel(
        # Bind the single-panel axes handle for the linear density plot
        ax,
        # Pass the continuous multi-species length table
        df,
        # Leave the title empty so the single-panel figure has no title text
        title="",
        # Restrict the linear axis to the configured 0–8000 bp window
        x_lim=(X_LIM_MIN, X_LIM_MAX),
        # Keep a linear (non-log) x-axis for raw abundance decay
        log_scale=False,
        # Show the species legend on the single-panel figure
        legend=True,
        # Omit the 1,500 bps text and arrow from this single-panel figure
        annotate=False,
        # Retain the linear annotation factor for API compatibility
        y_text_factor=ANNOTATION_Y_TEXT_FACTOR_LINEAR,
    )
    # Remove any residual title artist so the canvas stays title-free
    ax.set_title("")
    # Enlarge axis label text to the requested 16 pt size
    ax.xaxis.label.set_size(16)
    # Enlarge y-axis label text to the requested 16 pt size
    ax.yaxis.label.set_size(16)
    # Enlarge x-axis tick labels to the requested 12 pt size
    ax.tick_params(axis="x", labelsize=12)
    # Enlarge y-axis tick labels to the requested 12 pt size
    ax.tick_params(axis="y", labelsize=12)
    # Tighten layout so labels and the external legend fit within the canvas
    fig.tight_layout()
    # Export PNG and SVG figures and return their written paths
    return save_figure(fig, output_stem, output_dir=output_dir)


def plot_kde_linear_vs_log(
    df: pd.DataFrame,
    output_stem: str = DUAL_KDE_STEM,
    output_dir: Path | None = None,
) -> list[Path]:
    """Render a dual-panel linear-versus-log comparative KDE figure."""
    # Apply the shared Seaborn ticks theme before drawing subplots
    apply_publication_style()
    # Create a 1x2 side-by-side subplot figure for Panel A and Panel B
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=DUAL_KDE_FIGSIZE)
    # Draw Panel A: linear ORF-length KDE showing raw length decay
    _draw_kde_panel(
        # Bind Panel A to the left subplot axes
        ax1,
        # Pass the continuous multi-species length table
        df,
        # Label Panel A as the linear 0–8,000 bp abundance view
        title=PANEL_A_TITLE,
        # Restrict Panel A to the linear 0–8000 bp window
        x_lim=(X_LIM_MIN, X_LIM_MAX),
        # Keep Panel A on a linear x-axis scale
        log_scale=False,
        # Suppress Panel A legend so Panel B owns the shared legend
        legend=False,
        # Use the linear-panel annotation text height factor
        y_text_factor=ANNOTATION_Y_TEXT_FACTOR_LINEAR,
    )
    # Draw Panel B: log10 ORF-length KDE resolving short ORFs and log-normal shape
    _draw_kde_panel(
        # Bind Panel B to the right subplot axes
        ax2,
        # Pass the same continuous multi-species length table
        df,
        # Label Panel B as the logarithmic log10 bps resolution view
        title=PANEL_B_TITLE,
        # Restrict Panel B to a strictly positive log-compatible window
        x_lim=(LOG_X_LIM_MIN, LOG_X_LIM_MAX),
        # Enable seaborn log_scale so the x-axis becomes log10(bps)
        log_scale=True,
        # Keep a single shared species legend on Panel B
        legend=True,
        # Use the log-panel annotation text height factor for denser peaks
        y_text_factor=ANNOTATION_Y_TEXT_FACTOR_LOG,
    )
    # Tighten layout so both panels and the shared legend fit the canvas
    fig.tight_layout()
    # Export PNG and SVG dual-panel figures and return their written paths
    return save_figure(fig, output_stem, output_dir=output_dir)
