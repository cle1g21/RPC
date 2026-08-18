#!/usr/bin/env python3
"""
Univariate property distribution plots (histogram + KDE) for PeptiVerse scores.

Environment:
    conda activate pep_vis_env
    python properties_histogram.py
"""

from __future__ import annotations

# Import argparse so runtime flags can override CONFIG without editing the script.
import argparse
# Import glob to discover all sharded input CSV files from the results directory.
import glob
# Import os for joining paths when building the input file search pattern.
import os
# Import Path for validating that the output directory already exists on disk.
from pathlib import Path
# Import typing helpers for clear function signatures throughout the module.
from typing import Any

# Import matplotlib.pyplot to create figures, axes, and save finished plots.
import matplotlib.pyplot as plt
# Import Line2D for custom legend handles on distribution plots.
from matplotlib.lines import Line2D
# Import numpy for percentile and summary-statistic calculations on score arrays.
import numpy as np
# Import pandas to load, merge, and filter the peptide property dataset.
import pandas as pd
# Import seaborn for histogram and KDE distribution plotting.
import seaborn as sns


# CONFIG centralizes every path, column name, and visual parameter for instant retuning.
CONFIG: dict[str, Any] = {
    # Directory containing sharded PeptiVerse prediction CSV shard files.
    "input_dir": "/iridisfs/ddnb/shared_files/pepdiverse/results",
    # Glob pattern that selects all property-prediction shard CSVs in input_dir.
    "input_glob": "predictions.Ribo-seq_ORFs.comprehensive.*.csv",
    # Pre-existing directory where all distribution PNG figures will be written.
    "output_dir": "/home/cle1g21/RPC/peptiverse/plots",
    # Column holding the unique alphanumeric peptide identifier for each ORF row.
    "id_column": "id",
    # All five continuous pharmacological score columns to plot as distributions.
    "score_columns": [
        "hemolysis_score",
        "nf_score",
        "solubility_score",
        "permeability_penetrance_score",
        "halflife_score",
    ],
    # Output stem names for the five standalone distribution figures.
    "output_stems": {
        "hemolysis_score": "haemolysis",
        "nf_score": "non_fouling",
        "solubility_score": "solubility",
        "permeability_penetrance_score": "permeability",
        "halflife_score": "half-life",
    },
    # Width and height in inches for every standalone distribution figure.
    "figsize_inches": (7, 7),
    # Resolution in dots-per-inch for crisp publication-quality PNG output.
    "dpi": 300,
    # Number of histogram bins for univariate distribution plots.
    "hist_bins": 50,
    # Random seed for reproducible styling behavior where randomness is involved.
    "random_seed": 42,
    # Axis label font size.
    "axis_label_fontsize": 16,
    # Tick label font size.
    "tick_label_fontsize": 12,
    # Legend text font size.
    "legend_fontsize": 12,
    # Title font size.
    "title_fontsize": 16,
    # Legend face color with transparency (alpha applied via framealpha).
    "legend_facecolor": "gray",
    # Legend frame alpha (0.2 = translucent grey background).
    "legend_framealpha": 0.2,
    # Human-readable labels for distribution plot titles and x-axis text.
    "metric_labels": {
        "hemolysis_score": "Haemolysis score",
        "nf_score": "Non-fouling score",
        "solubility_score": "Solubility score",
        "permeability_penetrance_score": "Permeability (penetrance) score",
        "halflife_score": "Half-life (hours)",
    },
}


def parse_args() -> dict[str, Any]:
    """Parse CLI arguments and return a merged configuration dictionary."""
    # Create the argument parser with a short description of this script's purpose.
    parser = argparse.ArgumentParser(
        description="Generate univariate histogram + KDE distributions for peptide properties."
    )
    # Allow the input directory to be overridden from the command line.
    parser.add_argument("--input-dir", default=CONFIG["input_dir"])
    # Allow the input glob pattern to be overridden from the command line.
    parser.add_argument("--input-glob", default=CONFIG["input_glob"])
    # Allow the output directory to be overridden from the command line.
    parser.add_argument("--output-dir", default=CONFIG["output_dir"])
    # Allow the peptide ID column name to be overridden from the command line.
    parser.add_argument("--id-column", default=CONFIG["id_column"])
    # Allow output DPI to be overridden from the command line.
    parser.add_argument("--dpi", type=int, default=CONFIG["dpi"])
    # Parse argv and store the resulting namespace object.
    args = parser.parse_args()
    # Copy CONFIG so CLI values can be merged without mutating module defaults.
    config = dict(CONFIG)
    # Write the CLI input directory into the runtime config copy.
    config["input_dir"] = args.input_dir
    # Write the CLI input glob into the runtime config copy.
    config["input_glob"] = args.input_glob
    # Write the CLI output directory into the runtime config copy.
    config["output_dir"] = args.output_dir
    # Write the CLI ID column name into the runtime config copy.
    config["id_column"] = args.id_column
    # Write the CLI DPI value into the runtime config copy.
    config["dpi"] = args.dpi
    # Return the fully merged runtime configuration to the caller.
    return config


def validate_output_dir(config: dict[str, Any]) -> None:
    """Confirm the output directory exists; the script never creates it."""
    # Resolve the configured output directory as a Path object.
    output_dir = Path(config["output_dir"])
    # Raise a clear error if the plots directory is missing on disk.
    if not output_dir.is_dir():
        raise FileNotFoundError(
            f"Output directory does not exist (will not be created): {output_dir}"
        )


def load_property_data(config: dict[str, Any]) -> pd.DataFrame:
    """Glob, merge, validate, and clean sharded property prediction CSV files."""
    # Join input_dir and input_glob into one filesystem search pattern.
    pattern = os.path.join(config["input_dir"], config["input_glob"])
    # Collect all matching CSV paths in sorted order for deterministic loading.
    csv_files = sorted(glob.glob(pattern))
    # Stop early with a helpful error when no input shards are found.
    if not csv_files:
        raise FileNotFoundError(f"No input files matched pattern: {pattern}")
    # Read every shard file into a list of DataFrames.
    frames = [pd.read_csv(path) for path in csv_files]
    # Stack all shard tables into one continuous DataFrame with a fresh index.
    df = pd.concat(frames, ignore_index=True)
    # Read the configured peptide identifier column name.
    id_col = config["id_column"]
    # Read the list of required score columns for plotting.
    score_columns = config["score_columns"]
    # Build the full list of columns that must be present after merging.
    required = [id_col] + score_columns
    # Identify any required columns that are absent from the merged table.
    missing = [col for col in required if col not in df.columns]
    # Raise an error listing missing columns so the user can fix input data.
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    # Drop rows with NaN in any score column to keep plots mathematically valid.
    df = df.dropna(subset=score_columns).copy()
    # Remove unconstrained negative half-life predictions before ranking extremes.
    if "halflife_score" in df.columns:
        # Count rows before the half-life filter for a transparent log message.
        n_before = len(df)
        # Keep only peptides with non-negative predicted half-life.
        df = df[df["halflife_score"] >= 0].copy()
        # Report how many negative half-life rows were removed.
        removed = n_before - len(df)
        # Print the filter summary when any rows were dropped.
        if removed > 0:
            print(f"Filtered {removed:,} rows with halflife_score < 0.")
    # Parse ORF identifier, gene symbol, and ORF type from the compound id string.
    df = add_identity_columns(df, config)
    # Return the cleaned merged peptide property table.
    return df


def add_identity_columns(df: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Parse orf_id, gene, and orf_type from the colon-separated id column."""
    # Read the peptide identifier column name from config.
    id_col = config["id_column"]
    # Split each id into colon-separated tokens for field extraction.
    parts = df[id_col].astype(str).str.split(":")
    # ORF id is the first token (e.g. c1riboseqorf661).
    df["orf_id"] = parts.str[0]
    # Gene symbol is the middle token when present (e.g. LINC01409).
    df["gene"] = parts.str[1]
    # ORF type is the last token (e.g. lncRNA, uORF).
    df["orf_type"] = parts.str[-1]
    # Return the dataframe with identity helper columns attached.
    return df


def format_extreme_row(row: pd.Series, col: str) -> str:
    """Build a human-readable multi-line summary for an extreme peptide row."""
    # Prefer the parsed gene symbol when available.
    gene = row.get("gene", "NA")
    # Prefer the amino acid sequence column from the prediction table.
    sequence = row.get("sequence", "NA")
    # Start with the compact gene + AA sequence header the user asked for.
    lines = [f"{gene} (AA sequence {sequence})"]
    # Append every column from the full source row for complete logging.
    for key, value in row.items():
        # Skip duplicate display of gene/sequence already shown in the header.
        if key in {"gene", "orf_id", "orf_type"}:
            continue
        # Format floating-point scores with fixed precision for readability.
        if isinstance(value, float):
            lines.append(f"  {key}: {value:.6f}")
        else:
            lines.append(f"  {key}: {value}")
    # Highlight the property being ranked as the extreme.
    lines.append(f"  [ranked_by]: {col} = {row[col]:.6f}")
    # Join into one printable block.
    return "\n".join(lines)


def print_extreme_genes(df: pd.DataFrame, col: str) -> tuple[pd.Series, pd.Series]:
    """Print and return the highest and lowest scoring peptides for one property."""
    # Find the row index of the maximum score for this property column.
    idx_high = df[col].idxmax()
    # Find the row index of the minimum score for this property column.
    idx_low = df[col].idxmin()
    # Extract the full highest-scoring peptide row.
    high_row = df.loc[idx_high]
    # Extract the full lowest-scoring peptide row.
    low_row = df.loc[idx_low]
    # Print a clear section header for this property.
    print(f"\n=== {col} extremes ===")
    # Print the highest gene with full row information including AA sequence.
    print("HIGHEST:")
    print(format_extreme_row(high_row, col))
    # Print the lowest gene with full row information including AA sequence.
    print("LOWEST:")
    print(format_extreme_row(low_row, col))
    # Return both extreme rows so the plot can annotate them.
    return high_row, low_row


def apply_plot_style(config: dict[str, Any]) -> None:
    """Set global seaborn/matplotlib theme for consistent publication styling."""
    # Fix the random seed so any stochastic drawing steps are reproducible.
    np.random.seed(config["random_seed"])
    # Apply a clean tick-based theme suited to print and slide figures.
    sns.set_theme(style="ticks", context="paper", font_scale=1.05)


def compute_summary_stats(values: np.ndarray) -> dict[str, float]:
    """Compute n, mean, std, min, max, and quartiles for one numeric score array."""
    # Count the number of non-NaN observations in the values array.
    n = int(values.size)
    # Return a dictionary of summary statistics used by distribution overlays.
    return {
        "n": n,
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
        "min": float(np.min(values)),
        "q25": float(np.percentile(values, 25)),
        "median": float(np.percentile(values, 50)),
        "q75": float(np.percentile(values, 75)),
        "max": float(np.max(values)),
    }


def draw_distribution_on_ax(
    ax: plt.Axes,
    df: pd.DataFrame,
    col: str,
    config: dict[str, Any],
    high_row: pd.Series,
    low_row: pd.Series,
) -> None:
    """Draw one histogram + KDE distribution with a single Q1/Med/Q3/min/max legend."""
    # Extract the numeric values for this property as a numpy array.
    values = df[col].to_numpy()
    # Compute summary statistics and quartiles for overlay annotations.
    stats = compute_summary_stats(values)
    # Look up the human-readable metric name for titles and axis labels.
    metric_label = config["metric_labels"].get(col, col)
    # Draw a density-scaled histogram for the property score distribution.
    sns.histplot(
        x=values,
        bins=config["hist_bins"],
        stat="density",
        alpha=0.45,
        color="#4C72B0",
        edgecolor="white",
        ax=ax,
    )
    # Overlay a smooth kernel density estimate curve on the same axes.
    sns.kdeplot(x=values, color="#DD8452", linewidth=1.8, ax=ax)
    # Draw a dashed vertical line at the first quartile (Q1).
    ax.axvline(stats["q25"], color="gray", linestyle="--", linewidth=1.2)
    # Draw a solid vertical line at the median.
    ax.axvline(stats["median"], color="black", linestyle="-", linewidth=1.4)
    # Draw a dashed vertical line at the third quartile (Q3).
    ax.axvline(stats["q75"], color="gray", linestyle="--", linewidth=1.2)
    # Pin non-negative scores so the axis does not start below zero from KDE padding.
    if float(values.min()) >= 0:
        # Read the current upper x limit after hist/KDE drawing.
        x_right = ax.get_xlim()[1]
        # Force the left edge of the axis to start at zero.
        ax.set_xlim(0, x_right)
    # Set the x-axis label to the human-readable property name.
    ax.set_xlabel(metric_label, fontsize=config["axis_label_fontsize"])
    # Set the y-axis label to indicate density scaling on the histogram.
    ax.set_ylabel("Density", fontsize=config["axis_label_fontsize"])
    # Apply configured tick label sizes.
    ax.tick_params(axis="both", labelsize=config["tick_label_fontsize"])
    # Build one legend: Q1/Med/Q3 with line styles; min/max as text-only (no graph lines).
    legend_handles = [
        Line2D(
            [0],
            [0],
            color="gray",
            linestyle="--",
            linewidth=1.2,
            label=f"Q1 = {stats['q25']:.3f}",
        ),
        Line2D(
            [0],
            [0],
            color="black",
            linestyle="-",
            linewidth=1.4,
            label=f"Median = {stats['median']:.3f}",
        ),
        Line2D(
            [0],
            [0],
            color="gray",
            linestyle="--",
            linewidth=1.2,
            label=f"Q3 = {stats['q75']:.3f}",
        ),
        Line2D(
            [0],
            [0],
            color="none",
            marker="",
            linestyle="None",
            label=f"Min = {stats['min']:.3f} ({low_row['gene']})",
        ),
        Line2D(
            [0],
            [0],
            color="none",
            marker="",
            linestyle="None",
            label=f"Max = {stats['max']:.3f} ({high_row['gene']})",
        ),
    ]
    # Show a single legend with translucent grey background.
    ax.legend(
        handles=legend_handles,
        loc="upper right",
        fontsize=config["legend_fontsize"],
        facecolor=config["legend_facecolor"],
        framealpha=config["legend_framealpha"],
        edgecolor="gray",
    )
    # Remove top and right spines for a cleaner distribution plot frame.
    sns.despine(ax=ax)


def plot_univariate_distribution(
    df: pd.DataFrame,
    col: str,
    config: dict[str, Any],
    high_row: pd.Series,
    low_row: pd.Series,
) -> plt.Figure:
    """Create one standalone histogram + KDE distribution figure."""
    # Create a new square figure with the configured dimensions in inches.
    fig, ax = plt.subplots(figsize=config["figsize_inches"])
    # Draw the distribution content onto the single axes.
    draw_distribution_on_ax(ax, df, col, config, high_row, low_row)
    # Tighten layout so labels and legend fit within the saved PNG bounds.
    fig.tight_layout()
    # Return the completed distribution figure to the caller.
    return fig


def save_figure(fig: plt.Figure, output_path: Path, dpi: int) -> None:
    """Save a matplotlib figure as a high-resolution PNG and release its memory."""
    # Write the figure with tight bounding-box cropping at the requested DPI.
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    # Close the figure so matplotlib frees memory before generating the next plot.
    plt.close(fig)


def main() -> None:
    """Orchestrate loading, distribution curves, and saving."""
    # Merge CLI overrides into a runtime copy of CONFIG.
    config = parse_args()
    # Verify the output directory exists before any plotting work begins.
    validate_output_dir(config)
    # Apply the global seaborn/matplotlib visual theme.
    apply_plot_style(config)
    # Load and clean the merged peptide property dataset from sharded CSV files.
    df = load_property_data(config)
    # Resolve the output directory as a Path object for clean path joining.
    output_dir = Path(config["output_dir"])
    # Track saved file paths so they can be printed in a summary at the end.
    saved_paths: list[str] = []
    # Generate one standalone distribution figure for every score column.
    for col in config["score_columns"]:
        # Print and capture the highest and lowest gene rows for this property.
        high_row, low_row = print_extreme_genes(df, col)
        # Build the univariate histogram + KDE figure for this property column.
        fig = plot_univariate_distribution(df, col, config, high_row, low_row)
        # Resolve the user-facing output stem (e.g. haemolysis, non_fouling).
        stem = config["output_stems"][col]
        # Build the output path for this standalone distribution figure.
        out_path = output_dir / f"{stem}.png"
        # Save the distribution figure to disk at the configured DPI.
        save_figure(fig, out_path, config["dpi"])
        # Record the saved path for the final summary printout.
        saved_paths.append(str(out_path))
    # Print how many peptides were loaded so the user can sanity-check the run.
    print(f"\nLoaded {len(df):,} peptides.")
    # Print each saved output path so the user knows where figures were written.
    for path in saved_paths:
        print(f"Saved: {path}")


# Execute main() only when this file is run directly as a script, not when imported.
if __name__ == "__main__":
    main()
