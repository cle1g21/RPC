#!/usr/bin/env python3
"""Parse master_low.csv, summarize AA length and confidence, and plot a 2-panel figure."""

# Enable postponed evaluation of annotations for cleaner type hints
from __future__ import annotations

# Import argparse to expose CLI flags for input path, output directory, and DPI
import argparse
# Import math helpers only indirectly via numpy; keep re for UNIMOD token stripping
import re
# Import sys so a non-zero exit code can be returned on fatal errors
import sys
# Import Path for portable filesystem path construction
from pathlib import Path
# Import Sequence for typing ordered collections of candidate column names
from typing import Sequence

# Import numpy for exponential conversion of log-probability scores
import numpy as np
# Import pandas for CSV loading, binning, and summary statistics
import pandas as pd
# Import matplotlib figure API for multi-panel layout and export
import matplotlib.pyplot as plt
# Import seaborn for publication-ready bar and box plots
import seaborn as sns


# Resolve this script file to an absolute Path for locating the package root
_SCRIPT_PATH = Path(__file__).resolve()
# Treat the visualizer package root as the parent of the scripts/ directory
_PACKAGE_ROOT = _SCRIPT_PATH.parent.parent
# Default assembled Low-tier master CSV produced by the confidence pipeline
_DEFAULT_INPUT = Path(
    "/home/cle1g21/RPC/confidence_levels/assembled/master_low.csv"
)
# Default directory for stats text and publication figures
_DEFAULT_OUTPUT_DIR = _PACKAGE_ROOT / "outputs"
# Default raster export resolution in dots per inch
_DEFAULT_DPI = 300
# Ordered preference list for the peptide sequence column
_SEQUENCE_CANDIDATES: Sequence[str] = ("predictions", "sequence")
# Ordered preference list for the model confidence / score column
_SCORE_CANDIDATES: Sequence[str] = (
    "confidence",
    "preds_score",
    "score",
    "log_probs",
    "instanovoplus_prediction_log_probability",
)
# Regular expression that matches bracketed modification tokens such as [UNIMOD:4]
_MOD_TOKEN_RE = re.compile(r"\[[^\]]+\]")
# Numeric edges for pd.cut length bins (open-ended lower/upper via -inf / inf)
_BIN_EDGES = [-np.inf, 7, 11, 15, 20, 30, np.inf]
# Human-readable labels aligned one-to-one with the intervals defined by _BIN_EDGES
_BIN_LABELS = [
    "<8 AA",
    "8-11 AA [HLA Class I]",
    "12-15 AA [HLA Class II]",
    "16-20 AA",
    "21-30 AA",
    ">30 AA",
]
# Stem used for PNG and SVG figure filenames without extension
_FIGURE_STEM = "master_low_length_and_accuracy"
# Filename for the human-readable summary statistics log
_STATS_FILENAME = "master_low_summary_stats.txt"
# Bar fill color for Panel A (muted scientific blue-grey)
_BAR_COLOR = "#4C72B0"
# Boxplot face color for Panel B
_BOX_COLOR = "#55A868"
# Marker color for per-bin mean overlays on Panel B
_MEAN_COLOR = "#C44E52"
# Marker color for per-bin median overlays on Panel B
_MEDIAN_COLOR = "#1A1A1A"
# Spine line width applied after sns.despine for publication axes
_SPINE_LINEWIDTH = 2.0


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for input CSV, output directory, and DPI."""
    # Construct the top-level argument parser with a short program description
    parser = argparse.ArgumentParser(
        description=(
            "Summarize master_low peptide lengths and prediction confidence, "
            "then export a 2-panel Seaborn figure."
        )
    )
    # Register --input pointing at the assembled Low-tier CSV by default
    parser.add_argument(
        "--input",
        type=Path,
        default=_DEFAULT_INPUT,
        help="Path to master_low.csv (default: assembled Low-tier master).",
    )
    # Register --output-dir for stats text and figure exports
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_DEFAULT_OUTPUT_DIR,
        help="Directory for summary text and PNG/SVG figures.",
    )
    # Register --dpi controlling raster PNG resolution
    parser.add_argument(
        "--dpi",
        type=int,
        default=_DEFAULT_DPI,
        help="Raster DPI for PNG export (default: 300).",
    )
    # Parse argv (or sys.argv when argv is None) into a Namespace
    return parser.parse_args(argv)


def strip_mod_tokens(sequence: str) -> str:
    """Remove bracketed modification tokens from a peptide string."""
    # Coerce the cell value to string in case pandas read a non-string dtype
    text = str(sequence)
    # Replace every [UNIMOD:n]-style token with an empty string
    return _MOD_TOKEN_RE.sub("", text)


def resolve_column(columns: Sequence[str], candidates: Sequence[str], kind: str) -> str:
    """Return the first candidate column present in columns or raise ValueError."""
    # Build a set for O(1) membership tests against the DataFrame columns
    available = set(columns)
    # Walk the preference-ordered candidate list
    for name in candidates:
        # Return immediately when the preferred column exists
        if name in available:
            # Successful auto-detection of the requested column kind
            return name
    # Raise when none of the expected columns are present in the CSV
    raise ValueError(
        f"Could not auto-detect a {kind} column among {list(candidates)}; "
        f"available columns: {list(columns)}"
    )


def compute_summary_metrics(series: pd.Series) -> dict[str, float]:
    """Compute minimum, maximum, mean, and median for a numeric Series."""
    # Drop missing values so summary metrics are defined only on observed data
    clean = series.dropna()
    # Return NaN metrics when the series is empty after dropping nulls
    if clean.empty:
        # Empty-bin friendly dictionary of NaNs
        return {"min": float("nan"), "max": float("nan"), "mean": float("nan"), "median": float("nan")}
    # Compute the four requested summary statistics as native floats
    return {
        "min": float(clean.min()),
        "max": float(clean.max()),
        "mean": float(clean.mean()),
        "median": float(clean.median()),
    }


def format_metrics(metrics: dict[str, float], precision: int = 6) -> str:
    """Format a metrics dictionary as a single human-readable line."""
    # Render each metric with fixed floating-point precision for the log file
    return (
        f"min={metrics['min']:.{precision}g}, "
        f"max={metrics['max']:.{precision}g}, "
        f"mean={metrics['mean']:.{precision}g}, "
        f"median={metrics['median']:.{precision}g}"
    )


def load_and_enrich(input_path: Path) -> tuple[pd.DataFrame, str, str, bool]:
    """Load master_low.csv and attach aa_length, length_bin, raw_score, and confidence."""
    # Fail fast when the configured input CSV does not exist on disk
    if not input_path.is_file():
        # Raise FileNotFoundError with the absolute path for easier debugging
        raise FileNotFoundError(f"Input CSV not found: {input_path}")
    # Read the entire Low-tier master table into a DataFrame (71 rows; small)
    df = pd.read_csv(input_path)
    # Auto-detect the peptide sequence column (predictions preferred)
    seq_col = resolve_column(df.columns, _SEQUENCE_CANDIDATES, "sequence")
    # Auto-detect the model score / confidence column (log_probs preferred here)
    score_col = resolve_column(df.columns, _SCORE_CANDIDATES, "score")
    # Strip UNIMOD bracket tokens before measuring amino-acid length
    plain_seq = df[seq_col].map(strip_mod_tokens)
    # Compute amino-acid length as the character length of the stripped sequence
    df["aa_length"] = plain_seq.str.len().astype(int)
    # Copy the raw numeric score into a stable column name for logging
    df["raw_score"] = pd.to_numeric(df[score_col], errors="coerce")
    # Decide whether the detected score is a log-probability needing exp()
    is_log_prob = "log_prob" in score_col.lower()
    # Convert log-probabilities to (0, 1] confidence; otherwise use the score as-is
    if is_log_prob:
        # Exponentiate log_probs to obtain a probability-scale confidence
        df["confidence"] = np.exp(df["raw_score"].to_numpy(dtype=float))
    else:
        # Treat an already-probability or free-form score as the display confidence
        df["confidence"] = df["raw_score"]
    # Bin amino-acid lengths into immunopeptidomics-oriented categorical intervals
    df["length_bin"] = pd.cut(
        df["aa_length"],
        bins=_BIN_EDGES,
        labels=_BIN_LABELS,
        right=True,
        include_lowest=True,
    )
    # Ensure the categorical preserves the biologically ordered bin sequence
    df["length_bin"] = pd.Categorical(df["length_bin"], categories=_BIN_LABELS, ordered=True)
    # Return the enriched frame plus metadata used by the stats writer
    return df, seq_col, score_col, is_log_prob


def build_stats_report(
    df: pd.DataFrame,
    input_path: Path,
    seq_col: str,
    score_col: str,
    is_log_prob: bool,
) -> str:
    """Build the full multi-section summary statistics text for console and file."""
    # Accumulate report lines in a list for efficient joining at the end
    lines: list[str] = []
    # Record the absolute input path for reproducibility
    lines.append(f"Input file: {input_path.resolve()}")
    # Record the total number of peptide rows analyzed
    lines.append(f"Total peptides (rows): {len(df)}")
    # Record which sequence column was auto-detected
    lines.append(f"Sequence column: {seq_col}")
    # Record which score column was auto-detected
    lines.append(f"Score column: {score_col}")
    # Clarify whether confidence was derived via exp(log_probs)
    if is_log_prob:
        # Document the log-probability to confidence transform
        lines.append("Confidence transform: confidence = exp(raw_score) from log-probability")
    else:
        # Document that the raw score was used directly as confidence
        lines.append("Confidence transform: confidence = raw_score (no exp transform)")
    # Insert a blank line before the global length section
    lines.append("")
    # Section header for global peptide length metrics
    lines.append("=== Global peptide length (AA) ===")
    # Compute min/max/mean/median amino-acid length across the full dataset
    length_metrics = compute_summary_metrics(df["aa_length"])
    # Append the formatted global length metrics line
    lines.append(format_metrics(length_metrics, precision=4))
    # Insert a blank line before the global confidence section
    lines.append("")
    # Section header for probability-scale confidence metrics
    lines.append("=== Global prediction confidence (display scale) ===")
    # Compute min/max/mean/median of the display confidence column
    conf_metrics = compute_summary_metrics(df["confidence"])
    # Append the formatted global confidence metrics line
    lines.append(format_metrics(conf_metrics, precision=6))
    # Always also report raw score metrics for log-prob interpretability
    lines.append("")
    # Section header for the unmodified model score column
    lines.append(f"=== Global raw score ({score_col}) ===")
    # Compute min/max/mean/median of the raw score values
    raw_metrics = compute_summary_metrics(df["raw_score"])
    # Append the formatted raw-score metrics line
    lines.append(format_metrics(raw_metrics, precision=6))
    # Insert a blank line before the per-bin breakdown
    lines.append("")
    # Section header for length-bin stratified summaries
    lines.append("=== Per length-bin summaries ===")
    # Count total rows once for percentage denominators
    total_n = len(df)
    # Iterate bins in categorical order so empty HLA bins still appear
    for bin_label in _BIN_LABELS:
        # Boolean mask selecting rows that fall in the current length bin
        mask = df["length_bin"] == bin_label
        # Subset the DataFrame to the current bin
        subset = df.loc[mask]
        # Count peptides in this bin
        n = int(len(subset))
        # Convert the bin count to a percentage of the global total
        pct = (100.0 * n / total_n) if total_n else 0.0
        # Header line naming the bin and reporting n / percent
        lines.append(f"-- {bin_label}: n={n} ({pct:.2f}%)")
        # When the bin is empty, note that metrics are undefined and continue
        if n == 0:
            # Explicit empty-bin message for the log file
            lines.append("   (empty bin; metrics undefined)")
            # Skip metric computation for empty bins
            continue
        # Compute length metrics restricted to this bin
        bin_length = compute_summary_metrics(subset["aa_length"])
        # Append the bin-level length metrics
        lines.append(f"   length AA: {format_metrics(bin_length, precision=4)}")
        # Compute confidence metrics restricted to this bin
        bin_conf = compute_summary_metrics(subset["confidence"])
        # Append the bin-level confidence metrics
        lines.append(f"   confidence: {format_metrics(bin_conf, precision=6)}")
        # Compute raw-score metrics restricted to this bin
        bin_raw = compute_summary_metrics(subset["raw_score"])
        # Append the bin-level raw-score metrics
        lines.append(f"   raw_score:  {format_metrics(bin_raw, precision=6)}")
    # Join all lines with newlines into a single report string
    return "\n".join(lines) + "\n"


def style_spines(ax: plt.Axes) -> None:
    """Apply sns.despine and thicken the remaining left/bottom spines."""
    # Remove the top and right spines for a clean publication look
    sns.despine(ax=ax)
    # Iterate the two retained spines and thicken them in black
    for side in ("left", "bottom"):
        # Access the named spine artist on the axes
        spine = ax.spines[side]
        # Set the spine color to solid black
        spine.set_color("black")
        # Set the spine linewidth to the publication constant
        spine.set_linewidth(_SPINE_LINEWIDTH)


def plot_two_panel_figure(df: pd.DataFrame, output_dir: Path, dpi: int) -> list[Path]:
    """Create Panel A (counts) and Panel B (confidence boxplot) and save PNG+SVG."""
    # Apply a ticks-style Seaborn theme suitable for paper figures
    sns.set_theme(style="ticks", context="paper")
    # Create a side-by-side two-panel figure with a wide landscape canvas
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.5), constrained_layout=True)
    # Alias the left axes as Panel A for abundance bars
    ax_a = axes[0]
    # Alias the right axes as Panel B for confidence distributions
    ax_b = axes[1]
    # Draw a count bar plot of peptides per length bin on Panel A
    sns.countplot(
        data=df,
        x="length_bin",
        order=_BIN_LABELS,
        color=_BAR_COLOR,
        ax=ax_a,
    )
    # Total peptide count used to convert bar heights into percentages
    total_n = max(len(df), 1)
    # Annotate each bar with exact count and percentage of the dataset
    for patch in ax_a.patches:
        # Read the bar height as the integer peptide count
        count = patch.get_height()
        # Skip annotation when a categorical bar is somehow missing
        if count is None or np.isnan(count):
            # Continue to the next bar artist
            continue
        # Convert the count to an integer for display
        count_i = int(count)
        # Compute the percentage of total peptides represented by this bar
        pct = 100.0 * count_i / total_n
        # Place the annotation slightly above the bar top
        ax_a.annotate(
            f"{count_i}\n({pct:.1f}%)",
            xy=(patch.get_x() + patch.get_width() / 2.0, count),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    # Set the Panel A title describing abundance by length bin
    ax_a.set_title("A. Peptide abundance by AA length bin", loc="left", fontweight="bold")
    # Label the Panel A x-axis with the length-bin categories
    ax_a.set_xlabel("Amino acid length bin")
    # Label the Panel A y-axis as peptide count / abundance
    ax_a.set_ylabel("Peptide abundance (count)")
    # Rotate x tick labels so long HLA bin names remain readable
    ax_a.tick_params(axis="x", rotation=25)
    # Apply despine and thick black spines to Panel A
    style_spines(ax_a)
    # Draw a box plot of display confidence stratified by length bin on Panel B
    sns.boxplot(
        data=df,
        x="length_bin",
        y="confidence",
        order=_BIN_LABELS,
        color=_BOX_COLOR,
        ax=ax_b,
        showfliers=True,
    )
    # Compute per-bin mean and median confidence for overlay markers
    grouped = df.groupby("length_bin", observed=False)["confidence"]
    # Mean confidence per ordered length bin (NaN for empty bins)
    means = grouped.mean()
    # Median confidence per ordered length bin (NaN for empty bins)
    medians = grouped.median()
    # Map each bin label to its categorical x-position index
    x_positions = {label: idx for idx, label in enumerate(_BIN_LABELS)}
    # Overlay mean markers (diamonds) for non-empty bins
    for label, value in means.items():
        # Skip empty bins where mean is undefined
        if pd.isna(value):
            # Continue to the next bin label
            continue
        # Scatter a diamond at the mean confidence for this bin
        ax_b.scatter(
            x_positions[str(label)],
            value,
            marker="D",
            s=45,
            color=_MEAN_COLOR,
            zorder=5,
            label="_nolegend_",
        )
        # Annotate the numeric mean slightly above the diamond
        ax_b.annotate(
            f"μ={value:.3f}",
            xy=(x_positions[str(label)], value),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7,
            color=_MEAN_COLOR,
        )
    # Overlay median markers (circles) for non-empty bins
    for label, value in medians.items():
        # Skip empty bins where median is undefined
        if pd.isna(value):
            # Continue to the next bin label
            continue
        # Scatter a circle at the median confidence for this bin
        ax_b.scatter(
            x_positions[str(label)],
            value,
            marker="o",
            s=35,
            facecolors="white",
            edgecolors=_MEDIAN_COLOR,
            linewidths=1.5,
            zorder=6,
            label="_nolegend_",
        )
        # Annotate the numeric median slightly below the circle
        ax_b.annotate(
            f"M={value:.3f}",
            xy=(x_positions[str(label)], value),
            xytext=(0, -12),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=7,
            color=_MEDIAN_COLOR,
        )
    # Add a compact legend explaining mean vs median marker shapes
    ax_b.scatter([], [], marker="D", s=45, color=_MEAN_COLOR, label="Mean")
    # Register an empty median legend handle matching the overlay style
    ax_b.scatter(
        [],
        [],
        marker="o",
        s=35,
        facecolors="white",
        edgecolors=_MEDIAN_COLOR,
        linewidths=1.5,
        label="Median",
    )
    # Show the mean/median legend in the upper right of Panel B
    ax_b.legend(frameon=False, loc="upper right", fontsize=8)
    # Set the Panel B title describing confidence by length bin
    ax_b.set_title(
        "B. Prediction confidence by AA length bin",
        loc="left",
        fontweight="bold",
    )
    # Label the Panel B x-axis with the same length-bin categories
    ax_b.set_xlabel("Amino acid length bin")
    # Label the Panel B y-axis as probability-scale confidence
    ax_b.set_ylabel("Prediction confidence (exp(log_probs))")
    # Rotate x tick labels on Panel B to match Panel A readability
    ax_b.tick_params(axis="x", rotation=25)
    # Apply despine and thick black spines to Panel B
    style_spines(ax_b)
    # Ensure the output directory exists before writing figure files
    output_dir.mkdir(parents=True, exist_ok=True)
    # Build the PNG path from the figure stem
    png_path = output_dir / f"{_FIGURE_STEM}.png"
    # Build the SVG path from the figure stem
    svg_path = output_dir / f"{_FIGURE_STEM}.svg"
    # Export the high-resolution raster PNG at the requested DPI
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    # Export the camera-ready vector SVG
    fig.savefig(svg_path, bbox_inches="tight")
    # Close the figure to free matplotlib memory
    plt.close(fig)
    # Return the list of written figure paths for logging
    return [png_path, svg_path]


def main(argv: Sequence[str] | None = None) -> int:
    """Run the full parse → stats → plot pipeline and return a process exit code."""
    # Parse CLI arguments into a Namespace
    args = parse_args(argv)
    # Resolve the input CSV path to an absolute location
    input_path = args.input.expanduser().resolve()
    # Resolve the output directory to an absolute location
    output_dir = args.output_dir.expanduser().resolve()
    # Load the CSV and attach length, bins, and confidence columns
    df, seq_col, score_col, is_log_prob = load_and_enrich(input_path)
    # Build the multi-section summary statistics report string
    report = build_stats_report(df, input_path, seq_col, score_col, is_log_prob)
    # Ensure the output directory exists before writing the stats file
    output_dir.mkdir(parents=True, exist_ok=True)
    # Construct the path to the summary statistics text log
    stats_path = output_dir / _STATS_FILENAME
    # Write the report to disk using UTF-8 encoding
    stats_path.write_text(report, encoding="utf-8")
    # Mirror the same report to the terminal console
    print(report, end="")
    # Inform the user where the stats file was written
    print(f"Wrote summary stats: {stats_path}")
    # Create and export the two-panel Seaborn figure
    figure_paths = plot_two_panel_figure(df, output_dir, dpi=int(args.dpi))
    # Print each exported figure path for operator confirmation
    for path in figure_paths:
        # Announce a single written figure artifact
        print(f"Wrote figure: {path}")
    # Return success exit code
    return 0


# Execute main only when this file is run as a script, not when imported
if __name__ == "__main__":
    # Propagate the process exit code from main into sys.exit
    sys.exit(main())
