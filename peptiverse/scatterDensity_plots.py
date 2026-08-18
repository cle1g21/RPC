#!/usr/bin/env python3
"""
High-clarity PeptiVerse trade-off scatter plots with KDE density coloring.

For univariate distributions, run properties_histogram.py.

Environment:
    conda activate pep_vis_env
    python scatterDensity_plots.py
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

# Import LinearSegmentedColormap to build the custom density rainbow gradient.
from matplotlib.colors import LinearSegmentedColormap
# Import adjust_text to reposition peptide labels so they do not overlap.
from adjustText import adjust_text
# Import matplotlib.pyplot to create figures, axes, and save finished plots.
import matplotlib.pyplot as plt
# Import numpy for array math used in KDE evaluation and axis padding.
import numpy as np
# Import pandas to load, merge, filter, and annotate the peptide property table.
import pandas as pd
# Import gaussian_kde to estimate localized 2D point density for every scatter point.
from scipy.stats import gaussian_kde
# Import seaborn for global styling helpers such as despine.
import seaborn as sns


# CONFIG centralizes every path, column name, and visual parameter for instant retuning.
CONFIG: dict[str, Any] = {
    # Directory containing sharded PeptiVerse prediction CSV shard files.
    "input_dir": "/iridisfs/ddnb/shared_files/pepdiverse/results",
    # Glob pattern that selects all property-prediction shard CSVs in input_dir.
    "input_glob": "predictions.Ribo-seq_ORFs.comprehensive.*.csv",
    # Pre-existing directory where all output PNG figures will be written.
    "output_dir": "/home/cle1g21/RPC/peptiverse/plots",
    # Column holding the unique alphanumeric peptide identifier for each ORF row.
    "id_column": "id",
    # Optional dedicated ORF-type column; when absent, type is parsed from id suffix.
    "orf_type_column": None,
    # All five continuous pharmacological score columns available in the dataset.
    "score_columns": [
        "hemolysis_score",
        "nf_score",
        "solubility_score",
        "permeability_penetrance_score",
        "halflife_score",
    ],
    # How many peptides to label per individual selection rule on each figure.
    "label_n": 5,
    # Figure size in inches: A4 landscape (297 mm x 210 mm) for print-ready output.
    "figsize_inches": (11.69, 8.27),
    # Resolution in dots-per-inch for crisp publication-quality PNG output.
    "dpi": 300,
    # Marker size for KDE-colored background points on scatter plots.
    "background_size": 10,
    # Marker size for circle rings around highlighted outlier points (scaled for A4).
    "highlight_size": 200,
    # Bold line width for circle rings around highlighted outlier points.
    "highlight_linewidth": 3.5,
    # Edge color for circle rings around highlighted outlier points.
    "highlight_edgecolor": "#000000",
    # Title font size scaled for A4 landscape figures.
    "title_fontsize": 24,
    # X/Y axis label font size scaled for A4 landscape figures.
    "axis_label_fontsize": 22,
    # Tick label font size scaled for A4 landscape figures.
    "tick_label_fontsize": 18,
    # Colorbar axis label font size (legend-equivalent text).
    "colorbar_label_fontsize": 18,
    # Colorbar tick label font size scaled for A4 landscape figures.
    "colorbar_tick_fontsize": 18,
    # Seaborn context font scale multiplier for global text sizing.
    "font_scale": 1.45,
    # Font size for bold peptide point labels (A4-scaled).
    "point_label_fontsize": 18,
    # Line width for connector lines from labels to circled points.
    "connector_linewidth": 2.0,
    # adjustText repulsion strength to separate overlapping point labels.
    "adjust_force_text": (1.2, 1.6),
    # adjustText bounding-box expansion around each label.
    "adjust_expand_text": (1.4, 1.6),
    # Initial label offset as a fraction of axis span from each highlighted point.
    "label_offset_fraction": 0.035,
    # Fractional padding above data max on each axis (lower bound is always 0).
    "axis_margin_fraction": 0.14,
    # Random seed for reproducible styling behavior where randomness is involved.
    "random_seed": 42,
    # Custom colormap stops: purple/blue = sparse (low density), red = dense (high density).
    "density_color_stops": [
        "#9400D3",
        "#0000FF",
        "#00AA00",
        "#FFFF00",
        "#FF8C00",
        "#FF0000",
    ],
    # Trade-off scatter specifications with output paths and labeling modes.
    "tradeoff_plots": [
        {
            "output_path": (
                "/home/cle1g21/RPC/peptiverse/plots/"
                "tradeoff_bioavailability_halflife_vs_solubility2NEW.png"
            ),
            "title": "Bioavailability Window: Half-Life vs. Solubility",
            "x_col": "halflife_score",
            "y_col": "solubility_score",
            "x_label": "Predicted half-life (hours)",
            "y_label": "Solubility score",
            "label_mode": "spatial_extremes",
        },
        {
            "output_path": (
                "/home/cle1g21/RPC/peptiverse/plots/"
                "tradeoff_bioavailability_halflife_vs_solubilityFINAL.png"
            ),
            "title": "",
            "x_col": "halflife_score",
            "y_col": "solubility_score",
            "x_label": "Predicted half-life (hours)",
            "y_label": "Solubility score (0–1 scale)",
            "label_mode": "high_x_high_y",
            "draw_connectors": False,
            "y_lim_max": 1.0,
            "outer_label_margin": True,
        },
        {
            "output_path": (
                "/home/cle1g21/RPC/peptiverse/plots/"
                "tradeoff_solubility_vs_permeabilityFINAL.png"
            ),
            "title": "",
            "x_col": "solubility_score",
            "y_col": "permeability_penetrance_score",
            "x_label": "Solubility score (0–1 scale)",
            "y_label": "Permeability (penetrance) score (0–1 scale)",
            "label_mode": "high_x_high_y",
            "draw_connectors": False,
            "x_lim_max": 1.0,
            "y_lim_max": 1.0,
            "outer_label_margin": True,
        },
        {
            "output_path": (
                "/home/cle1g21/RPC/peptiverse/plots/"
                "tradeoff_hemolysis_vs_nfFINAL.png"
            ),
            "title": "",
            "x_col": "hemolysis_score",
            "y_col": "nf_score",
            "x_label": "Haemolysis score (0–1 scale)",
            "y_label": "Non-fouling score (0–1 scale)",
            "label_mode": "high_y_low_x",
            "draw_connectors": False,
            "x_lim_max": 1.0,
            "y_lim_max": 1.0,
            "outer_label_margin": True,
        },
        {
            "output_path": (
                "/home/cle1g21/RPC/peptiverse/plots/"
                "tradeoff_circulation_halflife_vs_nf2NEW.png"
            ),
            "title": "Clean Circulation Stability: Half-Life vs. Non-Fouling",
            "x_col": "halflife_score",
            "y_col": "nf_score",
            "x_label": "Predicted half-life (hours)",
            "y_label": "Non-fouling score",
            "label_mode": "spatial_extremes",
        },
        {
            "output_path": (
                "/home/cle1g21/RPC/peptiverse/plots/"
                "tradeoff_safe_transport_hemolysis_vs_permeability2NEW.png"
            ),
            "title": "Safe Cellular Transport: Haemolysis vs. Permeability",
            "x_col": "hemolysis_score",
            "y_col": "permeability_penetrance_score",
            "x_label": "Haemolysis score (0–1 scale)",
            "y_label": "Permeability (penetrance) score",
            "label_mode": "spatial_extremes",
        },
    ],
}


def parse_args() -> dict[str, Any]:
    """Parse CLI arguments and return a merged configuration dictionary."""
    # Create the argument parser with a short description of this script's purpose.
    parser = argparse.ArgumentParser(
        description="Generate KDE-colored trade-off scatter plots for peptide properties."
    )
    # Allow the input directory to be overridden from the command line.
    parser.add_argument("--input-dir", default=CONFIG["input_dir"])
    # Allow the input glob pattern to be overridden from the command line.
    parser.add_argument("--input-glob", default=CONFIG["input_glob"])
    # Allow the output directory to be overridden from the command line.
    parser.add_argument("--output-dir", default=CONFIG["output_dir"])
    # Allow the peptide ID column name to be overridden from the command line.
    parser.add_argument("--id-column", default=CONFIG["id_column"])
    # Allow the per-rule label count to be overridden from the command line.
    parser.add_argument("--label-n", type=int, default=CONFIG["label_n"])
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
    # Write the CLI label-N count into the runtime config copy.
    config["label_n"] = args.label_n
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
    """Glob, merge, validate, clean, and annotate sharded property prediction CSVs."""
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
    # Record how many rows exist before removing negative half-life artifacts.
    n_before_halflife_filter = len(df)
    # Remove rows where half-life is negative (unconstrained ML artifact).
    df = df[df["halflife_score"] >= 0].copy()
    # Print how many negative half-life rows were removed for transparency.
    removed = n_before_halflife_filter - len(df)
    if removed > 0:
        print(f"Filtered {removed:,} rows with halflife_score < 0.")
    # Add ORF type and compound label columns used by adjustText annotations.
    df = add_label_columns(df, config)
    # Return the cleaned merged peptide property table.
    return df


def add_label_columns(df: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Add orf_id, gene, orf_type, and compound_label columns for annotation text."""
    # Read the peptide identifier column name from config.
    id_col = config["id_column"]
    # Read the optional dedicated ORF-type column name from config.
    orf_type_col = config.get("orf_type_column")
    # Split each id into colon-separated tokens for field extraction.
    parts = df[id_col].astype(str).str.split(":")
    # Parse the short ORF identifier as the first colon-separated segment.
    df["orf_id"] = parts.str[0]
    # Parse the gene symbol as the middle colon-separated segment.
    df["gene"] = parts.str[1]
    # Use the dedicated type column when it exists in the dataframe.
    if orf_type_col and orf_type_col in df.columns:
        df["orf_type"] = df[orf_type_col].astype(str)
    # Otherwise parse ORF type as the last colon-separated token.
    else:
        df["orf_type"] = parts.str[-1]
    # Build compound labels with gene and ORF type on separate lines.
    df["compound_label"] = (
        df["gene"].astype(str) + "\n" + df["orf_type"].astype(str)
    )
    # Return the dataframe with new annotation helper columns attached.
    return df


def format_compound_label(label: str, config: dict[str, Any]) -> str:
    """Return a compound label string, optionally truncated for readability."""
    # Read the optional maximum label character count from config.
    max_chars = config.get("max_label_chars")
    # Return the full compound label when no truncation limit is configured.
    if max_chars is None:
        return str(label)
    # Convert the label to a plain string for safe slicing operations.
    text = str(label)
    # Return unchanged when the label already fits within the character budget.
    if len(text) <= max_chars:
        return text
    # Truncate and append an ellipsis when the label exceeds the budget.
    return text[: max_chars - 1] + "…"


def select_spatial_extreme_ids(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    id_col: str,
) -> set[str]:
    """Select top-left, bottom-right, furthest-top, and furthest-right peptides."""
    # Copy only the columns needed for spatial ranking to avoid mutating the source table.
    ranked = df[[id_col, x_col, y_col]].copy()
    # Normalize x values to [0, 1] so corner scores are scale-invariant.
    x_min = float(ranked[x_col].min())
    # Read the maximum x value for normalization.
    x_max = float(ranked[x_col].max())
    # Normalize y values to [0, 1] so corner scores are scale-invariant.
    y_min = float(ranked[y_col].min())
    # Read the maximum y value for normalization.
    y_max = float(ranked[y_col].max())
    # Compute a safe x span to avoid division by zero on constant columns.
    x_span = x_max - x_min if x_max > x_min else 1.0
    # Compute a safe y span to avoid division by zero on constant columns.
    y_span = y_max - y_min if y_max > y_min else 1.0
    # Build normalized x coordinates in [0, 1].
    ranked["_x_norm"] = (ranked[x_col] - x_min) / x_span
    # Build normalized y coordinates in [0, 1].
    ranked["_y_norm"] = (ranked[y_col] - y_min) / y_span
    # Top-left score: high y and low x.
    ranked["_top_left"] = ranked["_y_norm"] - ranked["_x_norm"]
    # Bottom-right score: high x and low y.
    ranked["_bottom_right"] = ranked["_x_norm"] - ranked["_y_norm"]
    # Collect unique peptide IDs for the four spatial extremes.
    extreme_ids: set[str] = set()
    # Furthest top: maximum y value.
    extreme_ids.add(str(ranked.loc[ranked[y_col].idxmax(), id_col]))
    # Furthest right: maximum x value.
    extreme_ids.add(str(ranked.loc[ranked[x_col].idxmax(), id_col]))
    # Top-left corner extreme.
    extreme_ids.add(str(ranked.loc[ranked["_top_left"].idxmax(), id_col]))
    # Bottom-right corner extreme.
    extreme_ids.add(str(ranked.loc[ranked["_bottom_right"].idxmax(), id_col]))
    # Return the unique set of spatial extreme IDs (up to 4).
    return extreme_ids


def select_max_axes_and_joint_ids(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    id_col: str,
    include_joint: bool = True,
) -> dict[str, str]:
    """Select max-x, max-y, and optionally joint high-x/high-y peptides keyed by role."""
    # Copy only the columns needed for ranking to avoid mutating the source table.
    ranked = df[[id_col, x_col, y_col]].copy()
    # Normalize x and y to [0, 1] so the joint score is scale-invariant.
    x_min = float(ranked[x_col].min())
    x_max = float(ranked[x_col].max())
    y_min = float(ranked[y_col].min())
    y_max = float(ranked[y_col].max())
    x_span = x_max - x_min if x_max > x_min else 1.0
    y_span = y_max - y_min if y_max > y_min else 1.0
    ranked["_x_norm"] = (ranked[x_col] - x_min) / x_span
    ranked["_y_norm"] = (ranked[y_col] - y_min) / y_span
    # Joint optimum: high on both axes (product rewards balance).
    ranked["_joint"] = ranked["_x_norm"] * ranked["_y_norm"]
    # Return role → peptide id so annotation can place labels without overlap.
    roles = {
        "max_x": str(ranked.loc[ranked[x_col].idxmax(), id_col]),
        "max_y": str(ranked.loc[ranked[y_col].idxmax(), id_col]),
    }
    if include_joint:
        roles["joint_high"] = str(ranked.loc[ranked["_joint"].idxmax(), id_col])
    return roles


def select_high_y_low_x_id(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    id_col: str,
) -> dict[str, str]:
    """Select the single peptide with highest y and lowest x (normalized product)."""
    ranked = df[[id_col, x_col, y_col]].copy()
    x_min = float(ranked[x_col].min())
    x_max = float(ranked[x_col].max())
    y_min = float(ranked[y_col].min())
    y_max = float(ranked[y_col].max())
    x_span = x_max - x_min if x_max > x_min else 1.0
    y_span = y_max - y_min if y_max > y_min else 1.0
    ranked["_x_norm"] = (ranked[x_col] - x_min) / x_span
    ranked["_y_norm"] = (ranked[y_col] - y_min) / y_span
    # High y and low x: (1 - x_norm) * y_norm.
    ranked["_high_y_low_x"] = (1.0 - ranked["_x_norm"]) * ranked["_y_norm"]
    best_id = str(ranked.loc[ranked["_high_y_low_x"].idxmax(), id_col])
    return {"high_y_low_x": best_id}


def select_high_x_high_y_id(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    id_col: str,
) -> dict[str, str]:
    """Select the single peptide with highest x and highest y (normalized product)."""
    ranked = df[[id_col, x_col, y_col]].copy()
    x_min = float(ranked[x_col].min())
    x_max = float(ranked[x_col].max())
    y_min = float(ranked[y_col].min())
    y_max = float(ranked[y_col].max())
    x_span = x_max - x_min if x_max > x_min else 1.0
    y_span = y_max - y_min if y_max > y_min else 1.0
    ranked["_x_norm"] = (ranked[x_col] - x_min) / x_span
    ranked["_y_norm"] = (ranked[y_col] - y_min) / y_span
    # High on both axes: x_norm * y_norm.
    ranked["_high_x_high_y"] = ranked["_x_norm"] * ranked["_y_norm"]
    best_id = str(ranked.loc[ranked["_high_x_high_y"].idxmax(), id_col])
    return {"high_x_high_y": best_id}


def select_high_x_low_y_id(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    id_col: str,
) -> dict[str, str]:
    """Select the single peptide with highest x and lowest y (normalized product)."""
    ranked = df[[id_col, x_col, y_col]].copy()
    x_min = float(ranked[x_col].min())
    x_max = float(ranked[x_col].max())
    y_min = float(ranked[y_col].min())
    y_max = float(ranked[y_col].max())
    x_span = x_max - x_min if x_max > x_min else 1.0
    y_span = y_max - y_min if y_max > y_min else 1.0
    ranked["_x_norm"] = (ranked[x_col] - x_min) / x_span
    ranked["_y_norm"] = (ranked[y_col] - y_min) / y_span
    # High x and low y: x_norm * (1 - y_norm).
    ranked["_high_x_low_y"] = ranked["_x_norm"] * (1.0 - ranked["_y_norm"])
    best_id = str(ranked.loc[ranked["_high_x_low_y"].idxmax(), id_col])
    return {"high_x_low_y": best_id}


def select_ids_by_genes(
    df: pd.DataFrame,
    genes: list[str],
    id_col: str,
    x_col: str | None = None,
    y_col: str | None = None,
) -> dict[str, str]:
    """Map role names to peptide IDs for requested gene symbols.

    When a gene has multiple ORFs, pick the highest joint x*y score if score
    columns are provided; otherwise take the first row.
    """
    roles: dict[str, str] = {}
    for gene in genes:
        matches = df[df["gene"].astype(str) == str(gene)].copy()
        if matches.empty:
            raise ValueError(f"Gene not found in dataset: {gene}")
        if len(matches) > 1 and x_col and y_col:
            matches = matches.assign(
                _joint=matches[x_col].astype(float) * matches[y_col].astype(float)
            ).sort_values("_joint", ascending=False)
        roles[f"gene_{gene}"] = str(matches.iloc[0][id_col])
    return roles


def resolve_label_ids(
    df: pd.DataFrame,
    spec: dict[str, Any],
    config: dict[str, Any],
) -> set[str]:
    """Dispatch to the correct labeling strategy for one trade-off figure spec."""
    # Read which labeling mode this trade-off figure uses.
    mode = spec.get("label_mode", "spatial_extremes")
    # Read the peptide identifier column name from config.
    id_col = config["id_column"]
    # Spatial extremes: top-left, bottom-right, furthest top, furthest right.
    if mode == "spatial_extremes":
        return select_spatial_extreme_ids(df, spec["x_col"], spec["y_col"], id_col)
    # Max x, max y, and optionally joint high-x/high-y product.
    if mode == "max_axes_and_joint":
        return set(
            select_max_axes_and_joint_ids(
                df,
                spec["x_col"],
                spec["y_col"],
                id_col,
                include_joint=spec.get("include_joint", True),
            ).values()
        )
    # Single peptide: highest y and lowest x.
    if mode == "high_y_low_x":
        return set(select_high_y_low_x_id(df, spec["x_col"], spec["y_col"], id_col).values())
    # Single peptide: highest x and highest y.
    if mode == "high_x_high_y":
        return set(select_high_x_high_y_id(df, spec["x_col"], spec["y_col"], id_col).values())
    # Single peptide: highest x and lowest y.
    if mode == "high_x_low_y":
        return set(select_high_x_low_y_id(df, spec["x_col"], spec["y_col"], id_col).values())
    # High x / low y plus one or more named genes.
    if mode == "high_x_low_y_plus_genes":
        roles = select_high_x_low_y_id(df, spec["x_col"], spec["y_col"], id_col)
        roles.update(
            select_ids_by_genes(
                df,
                list(spec.get("extra_genes", [])),
                id_col,
                x_col=spec["x_col"],
                y_col=spec["y_col"],
            )
        )
        return set(roles.values())
    # No labels or circles.
    if mode == "none":
        return set()
    # Raise an error when an unknown labeling mode is encountered in CONFIG.
    raise ValueError(f"Unknown label_mode: {mode}")


def resolve_highlight_ids(
    df: pd.DataFrame,
    spec: dict[str, Any],
    config: dict[str, Any],
) -> set[str]:
    """Return peptide IDs that receive bold circle highlights (same as labels)."""
    # Circle the same spatial extreme peptides that receive text labels.
    return resolve_label_ids(df, spec, config)


def build_density_colormap(config: dict[str, Any]) -> LinearSegmentedColormap:
    """Build colormap where low density = purple and high density = red."""
    # Read the ordered list of hex color stops from CONFIG (sparse → dense).
    colors = config["density_color_stops"]
    # Create a matplotlib colormap: vmin maps to purple, vmax maps to red.
    return LinearSegmentedColormap.from_list("density_rainbow", colors)


def compute_point_densities(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Estimate localized 2D Gaussian KDE density at each (x, y) data point."""
    # Stack x and y into a 2×N array shape required by gaussian_kde.
    stacked = np.vstack([x, y])
    # Fit a Gaussian kernel density estimator to the 2D point cloud.
    kde = gaussian_kde(stacked)
    # Evaluate the KDE at every original data point to get per-point density values.
    densities = kde(stacked)
    # Return the raw density value array aligned with each input point.
    return densities


def normalize_densities(densities: np.ndarray) -> np.ndarray:
    """Scale density values into the 0–1 range for colormap mapping."""
    # Read the minimum density value across all points.
    d_min = float(densities.min())
    # Read the maximum density value across all points.
    d_max = float(densities.max())
    # Return zeros when all densities are identical to avoid division by zero.
    if d_max == d_min:
        return np.zeros_like(densities)
    # Linearly rescale densities so the densest point maps to 1.0 and sparsest to 0.0.
    return (densities - d_min) / (d_max - d_min)


def set_axis_limits_from_origin(
    ax: plt.Axes,
    x_values: np.ndarray,
    y_values: np.ndarray,
    margin_fraction: float,
    y_lim_max: float | None = None,
    x_lim_max: float | None = None,
) -> None:
    """Set axes to start at (0, 0) and pad only the upper limits for label space."""
    # Compute the maximum x value in the plotted data.
    x_max = float(np.max(x_values))
    # Compute the maximum y value in the plotted data.
    y_max = float(np.max(y_values))
    # Use x_max as span reference when padding the upper x limit.
    x_span = x_max if x_max > 0 else 1.0
    # Use y_max as span reference when padding the upper y limit.
    y_span = y_max if y_max > 0 else 1.0
    # Optional fixed x upper bound (e.g. score axes with headroom for labels).
    if x_lim_max is not None:
        ax.set_xlim(0, float(x_lim_max))
    else:
        # Pin x-axis lower bound at zero and add top margin for label whitespace.
        ax.set_xlim(0, x_max + x_span * margin_fraction)
    # Optional fixed y upper bound (e.g. solubility score with headroom for labels).
    if y_lim_max is not None:
        ax.set_ylim(0, float(y_lim_max))
    else:
        # Pin y-axis lower bound at zero and add top margin for label whitespace.
        ax.set_ylim(0, y_max + y_span * margin_fraction)


def apply_plot_style(config: dict[str, Any]) -> None:
    """Set global seaborn/matplotlib theme for consistent publication styling."""
    # Fix the random seed so any stochastic drawing steps are reproducible.
    np.random.seed(config["random_seed"])
    # Apply a clean tick-based theme with A4-scaled fonts via font_scale.
    sns.set_theme(style="ticks", context="paper", font_scale=config["font_scale"])


def plot_kde_density_scatter(
    ax: plt.Axes,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    config: dict[str, Any],
    cmap: LinearSegmentedColormap,
) -> None:
    """Draw the full point cloud colored by per-point 2D KDE density."""
    # Extract x-axis values as a numpy float array.
    x = df[x_col].to_numpy(dtype=float)
    # Extract y-axis values as a numpy float array.
    y = df[y_col].to_numpy(dtype=float)
    # Compute localized KDE density for every point in the scatter.
    densities = compute_point_densities(x, y)
    # Normalize densities to 0–1 so they map cleanly onto the custom colormap.
    density_norm = normalize_densities(densities)
    # Scatter all points with KDE-mapped rainbow colors and small marker size.
    ax.scatter(
        x,
        y,
        c=density_norm,
        cmap=cmap,
        s=config["background_size"],
        linewidths=0,
        rasterized=True,
        zorder=1,
        vmin=0.0,
        vmax=1.0,
    )


def highlight_outlier_points(
    ax: plt.Axes,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    ids_to_highlight: set[str],
    config: dict[str, Any],
) -> None:
    """Draw bold circle rings around selected outlier peptides."""
    # Read the peptide identifier column name from config.
    id_col = config["id_column"]
    # Filter the table to only the rows whose IDs are in the highlight set.
    subset = df[df[id_col].astype(str).isin(ids_to_highlight)]
    # Exit immediately when there are no highlighted points for this axes.
    if subset.empty:
        return
    # Draw bold black circle rings around each highlighted peptide on top of the KDE cloud.
    ax.scatter(
        subset[x_col],
        subset[y_col],
        s=config["highlight_size"],
        facecolors="none",
        edgecolors=config["highlight_edgecolor"],
        linewidths=config["highlight_linewidth"],
        zorder=4,
        clip_on=False,
    )


def annotate_point_labels(
    ax: plt.Axes,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    ids_to_label: set[str],
    config: dict[str, Any],
    label_roles: dict[str, str] | None = None,
    draw_connectors: bool = True,
) -> None:
    """Place bold gene / ORF-type labels next to circled points."""
    # Read the peptide identifier column name from config.
    id_col = config["id_column"]
    # Filter the table to rows whose IDs should receive text labels.
    subset = df[df[id_col].astype(str).isin(ids_to_label)].copy()
    # Exit when there are no points to label on this axes.
    if subset.empty:
        return
    # Map peptide id → role name when role-aware placement is requested.
    id_to_role: dict[str, str] = {}
    if label_roles:
        for role, peptide_id in label_roles.items():
            id_to_role[str(peptide_id)] = role
    # Read axis limits to scale small fixed offsets beside each circle.
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_span = xlim[1] - xlim[0]
    y_span = ylim[1] - ylim[0]
    dx = x_span * config["label_offset_fraction"]
    dy = y_span * config["label_offset_fraction"]
    # Small fixed offsets so labels sit above each circle (may extend past axis=1).
    role_offsets = {
        "max_x": (0.0 * x_span, 0.08 * y_span),
        "max_y": (0.0 * x_span, 0.08 * y_span),
        "joint_high": (0.0 * x_span, 0.08 * y_span),
        "high_y_low_x": (0.0 * x_span, 0.08 * y_span),
        "high_x_high_y": (0.0 * x_span, 0.08 * y_span),
        "high_x_low_y": (0.0 * x_span, 0.08 * y_span),
    }
    # Gene-specific roles use the same above-circle offset.
    for role in list(id_to_role.values()):
        if role not in role_offsets:
            role_offsets[role] = (0.0 * x_span, 0.055 * y_span)
    texts: list[Any] = []
    point_xy: list[tuple[float, float]] = []
    for idx, (_, row) in enumerate(subset.iterrows()):
        x_val = float(row[x_col])
        y_val = float(row[y_col])
        # Label text: gene name and ORF type only.
        label = format_compound_label(row["compound_label"], config)
        role = id_to_role.get(str(row[id_col]))
        if role in role_offsets:
            offset_x, offset_y = role_offsets[role]
        else:
            offset_x = 0.0
            offset_y = dy * 1.2
        texts.append(
            ax.text(
                x_val + offset_x,
                y_val + offset_y,
                label,
                fontsize=config["point_label_fontsize"],
                fontweight="bold",
                color="black",
                zorder=5,
                clip_on=False,
                ha="center",
                va="bottom",
                linespacing=1.1,
            )
        )
        point_xy.append((x_val, y_val))
    # Optionally draw connector lines from labels to circled points.
    if draw_connectors:
        draw_label_connectors(
            ax,
            texts,
            np.array([xy[0] for xy in point_xy]),
            np.array([xy[1] for xy in point_xy]),
            config,
        )


def draw_label_connectors(
    ax: plt.Axes,
    texts: list[Any],
    x_points: np.ndarray,
    y_points: np.ndarray,
    config: dict[str, Any],
) -> None:
    """Draw lines from final label positions to their corresponding circled points."""
    # Loop over each text label paired with its target point coordinate.
    for text, x_pt, y_pt in zip(texts, x_points, y_points):
        # Read the final label position after adjustText repositioning.
        x_txt, y_txt = text.get_position()
        # Draw a line connecting the label to the circled peptide point.
        ax.plot(
            [x_txt, x_pt],
            [y_txt, y_pt],
            color="black",
            linewidth=config["connector_linewidth"],
            solid_capstyle="round",
            zorder=3,
            clip_on=True,
        )


def apply_a4_typography(
    ax: plt.Axes,
    cbar: Any,
    colorbar_label: str,
    config: dict[str, Any],
) -> None:
    """Set tick, axis, title, and colorbar font sizes for A4 landscape readability."""
    # Set tick label font size on both x and y axes.
    ax.tick_params(axis="both", labelsize=config["tick_label_fontsize"])
    # Set the x-axis label font size.
    ax.xaxis.label.set_fontsize(config["axis_label_fontsize"])
    # Set the y-axis label font size.
    ax.yaxis.label.set_fontsize(config["axis_label_fontsize"])
    # Set colorbar tick label font size.
    cbar.ax.tick_params(labelsize=config["colorbar_tick_fontsize"])
    # Set the colorbar axis label with A4-scaled font size.
    cbar.set_label(colorbar_label, fontsize=config["colorbar_label_fontsize"])


def plot_tradeoff_scatter(
    df: pd.DataFrame,
    spec: dict[str, Any],
    config: dict[str, Any],
    cmap: LinearSegmentedColormap,
) -> plt.Figure:
    """Create one KDE-colored trade-off scatter figure from a tradeoff_plots spec entry."""
    # Read the x-axis score column name from this trade-off specification.
    x_col = spec["x_col"]
    # Read the y-axis score column name from this trade-off specification.
    y_col = spec["y_col"]
    # Resolve peptide IDs that receive bold circle highlights on this figure.
    ids_to_highlight = resolve_highlight_ids(df, spec, config)
    # Resolve peptide IDs that receive text labels on this figure.
    ids_to_label = resolve_label_ids(df, spec, config)
    # Role map for FINAL-style placement.
    label_roles = None
    if spec.get("label_mode") == "max_axes_and_joint":
        label_roles = select_max_axes_and_joint_ids(
            df,
            x_col,
            y_col,
            config["id_column"],
            include_joint=spec.get("include_joint", True),
        )
    elif spec.get("label_mode") == "high_y_low_x":
        label_roles = select_high_y_low_x_id(
            df, x_col, y_col, config["id_column"]
        )
    elif spec.get("label_mode") == "high_x_high_y":
        label_roles = select_high_x_high_y_id(
            df, x_col, y_col, config["id_column"]
        )
    elif spec.get("label_mode") == "high_x_low_y":
        label_roles = select_high_x_low_y_id(
            df, x_col, y_col, config["id_column"]
        )
    elif spec.get("label_mode") == "high_x_low_y_plus_genes":
        label_roles = select_high_x_low_y_id(
            df, x_col, y_col, config["id_column"]
        )
        label_roles.update(
            select_ids_by_genes(
                df,
                list(spec.get("extra_genes", [])),
                config["id_column"],
                x_col=x_col,
                y_col=y_col,
            )
        )
    # Create a new square figure with the configured width and height in inches.
    fig, ax = plt.subplots(figsize=config["figsize_inches"])
    # Draw the KDE density-colored background cloud for all peptides.
    plot_kde_density_scatter(ax, df, x_col, y_col, config, cmap)
    # Use configured margin; keep modest padding when y is capped at 1.0.
    margin = config["axis_margin_fraction"]
    y_lim_max = spec.get("y_lim_max")
    x_lim_max = spec.get("x_lim_max")
    if label_roles is not None and y_lim_max is None:
        margin = max(margin, 0.22)
    # Set axis limits from origin (0, 0) with optional fixed upper bounds.
    set_axis_limits_from_origin(
        ax,
        df[x_col].to_numpy(dtype=float),
        df[y_col].to_numpy(dtype=float),
        margin,
        y_lim_max=y_lim_max,
        x_lim_max=x_lim_max,
    )
    # Draw bold circle rings on all highlighted peptides (including unlabeled toxic outliers).
    highlight_outlier_points(ax, df, x_col, y_col, ids_to_highlight, config)
    # Add bold text labels beside circled points (connectors optional per spec).
    annotate_point_labels(
        ax,
        df,
        x_col,
        y_col,
        ids_to_label,
        config,
        label_roles=label_roles,
        draw_connectors=spec.get("draw_connectors", True),
    )
    # Set the x-axis label using the human-readable text from the spec.
    ax.set_xlabel(spec["x_label"], fontsize=config["axis_label_fontsize"])
    # Set the y-axis label using the human-readable text from the spec.
    ax.set_ylabel(spec["y_label"], fontsize=config["axis_label_fontsize"])
    # Set the figure title only when a non-empty title is provided in the spec.
    if str(spec.get("title", "")).strip():
        ax.set_title(spec["title"], pad=14, fontsize=config["title_fontsize"])
    # Add a colorbar explaining that red = densest regions and purple = sparsest.
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    # Ensure the ScalarMappable has a valid cmap for the colorbar draw call.
    sm.set_array([])
    # Draw a vertical colorbar along the right side of the axes.
    cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    # Set the colorbar label describing the KDE density color encoding.
    colorbar_label = "Local point density (purple = sparse, red = dense)"
    # Apply A4-scaled typography to ticks, axes, title, and colorbar.
    apply_a4_typography(ax, cbar, colorbar_label, config)
    # Enable a light y-axis grid to help read values without adding visual noise.
    ax.grid(True, axis="y", linestyle=":", alpha=0.3)
    # Remove top and right spines for a cleaner publication-style frame.
    sns.despine(ax=ax)
    # Leave white figure margin beyond the axes so labels above y=1 / past x=1 remain visible.
    if spec.get("outer_label_margin"):
        fig.subplots_adjust(left=0.10, right=0.88, bottom=0.12, top=0.90)
    else:
        # Tighten subplot margins so labels and colorbar fit within the saved figure bounds.
        fig.tight_layout()
    # Return the completed matplotlib Figure object to the caller.
    return fig


def save_figure(fig: plt.Figure, output_path: Path, dpi: int) -> None:
    """Save a matplotlib figure as a high-resolution PNG and release its memory."""
    # Write the figure to disk with padding so outer labels are not cropped.
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", pad_inches=0.35)
    # Close the figure so matplotlib frees memory before generating the next plot.
    plt.close(fig)


def main() -> None:
    """Orchestrate loading, KDE trade-off scatters, and saving."""
    # Merge CLI overrides into a runtime copy of CONFIG.
    config = parse_args()
    # Verify the output directory exists before any plotting work begins.
    validate_output_dir(config)
    # Apply the global seaborn/matplotlib visual theme.
    apply_plot_style(config)
    # Load, clean, and annotate the merged peptide property dataset.
    df = load_property_data(config)
    # Build the custom density rainbow colormap once for all three figures.
    cmap = build_density_colormap(config)
    # Track saved file paths so they can be printed in a summary at the end.
    saved_paths: list[str] = []
    # Generate each configured trade-off scatter figure in sequence.
    for spec in config["tradeoff_plots"]:
        # Build the KDE-colored trade-off scatter figure from the current spec entry.
        fig = plot_tradeoff_scatter(df, spec, config, cmap)
        # Resolve the absolute output path declared in the trade-off spec.
        out_path = Path(spec["output_path"])
        # Save the trade-off figure to disk at the configured DPI.
        save_figure(fig, out_path, config["dpi"])
        # Record the saved path for the final summary printout.
        saved_paths.append(str(out_path))
    # Print how many peptides remain after cleaning so the user can sanity-check the run.
    print(f"Plotting {len(df):,} peptides after cleaning.")
    # Print each saved output path so the user knows where figures were written.
    for path in saved_paths:
        print(f"Saved: {path}")


# Execute main() only when this file is run directly as a script, not when imported.
if __name__ == "__main__":
    main()
