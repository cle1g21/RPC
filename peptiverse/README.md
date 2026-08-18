# PeptiVerse Peptide Property Visualization Pipeline

Generate publication-ready visualizations of PeptiVerse pharmacological property scores for the Ribo-seq ORF comprehensive dataset.

## Overview

| Script | Outputs |
| --- | --- |
| `scatterDensity_plots.py` | 3 trade-off scatter plots (7×7 in, 300 dpi) |
| `properties_histogram.py` | 5 univariate distribution plots (`dist_*.png`) |

**Input directory:** `/iridisfs/ddnb/shared_files/pepdiverse/results`  
**Input pattern:** `predictions.Ribo-seq_ORFs.comprehensive.*.csv` (15 shards, 28,359 ORFs)  
**Output directory:** `/home/cle1g21/RPC/peptiverse/plots/` (must pre-exist; scripts do not create it)

PeptiVerse input shards live outside this repo (`/iridisfs/ddnb/shared_files/pepdiverse/results`) and are not versioned. `__pycache__/` is gitignored. Generated PNGs under `plots/` can be committed if you want figures in GitHub.

### Trade-off scatter outputs (`scatterDensity_plots.py`)

**v2 (current):** KDE density-colored scatters with compound `{id} ({orf_type})` labels.

| File | Axes | Labeling |
| --- | --- | --- |
| `tradeoff_bioavailability_halflife_vs_solubility2.png` | half-life × solubility | top 5 longest half-lives ∪ top 5 highest solubilities |
| `tradeoff_circulation_halflife_vs_nf2.png` | half-life × non-fouling | top 5 longest half-lives ∪ top 5 highest nf scores |
| `tradeoff_safe_transport_hemolysis_vs_permeability2.png` | hemolysis × permeability | top 5 highest hemolysis (toxic) ∪ top 5 sweet-spot (high permeability − low hemolysis) |

Data cleaning: rows with `halflife_score < 0` are dropped before plotting. ORF type is parsed from the last segment of `id` (e.g. `lncRNA`, `uORF`).

Background points are colored by 2D Gaussian KDE density (red = densest, purple = sparsest).

**v1 (legacy):** `tradeoff_*_*.png` (without `2`) — uniform gray background, top-10 labels.

### Distribution outputs (`properties_histogram.py`)

| File | Property |
| --- | --- |
| `dist_hemolysis_score.png` | Hemolysis |
| `dist_nf_score.png` | Non-fouling |
| `dist_solubility_score.png` | Solubility |
| `dist_permeability_penetrance_score.png` | Permeability (penetrance) |
| `dist_halflife_score.png` | Half-life (hours) |

Each distribution figure includes histogram + KDE, Q1/median/Q3 lines, and a summary box (n, mean, std, min, max).

---

## Quick Start

```bash
conda activate pep_vis_env
cd /home/cle1g21/RPC/peptiverse

# Trade-off scatter plots (3 figures)
python scatterDensity_plots.py

# Univariate distributions (5 figures)
python properties_histogram.py
```

---

## Readable scatter pipeline (`scatterDensity_plots.py`)

### Styling

- Background cloud: 28,359 points at `s=10`, `alpha=0.15`, slate gray
- Outliers: crimson markers (`s=70`) with `adjustText` labels and pointer arrows
- Figure size: 7×7 inches, 300 dpi PNG
- `ensure_inside_axes=True` keeps labels inside plot bounds

### Outlier union filtering

For each trade-off figure, `union_outlier_ids()` applies directional rules from `CONFIG["tradeoff_plots"]`:

1. Sort the full DataFrame by each rule column (ascending or descending).
2. Deduplicate by `id`, keeping the best score per peptide.
3. Take the top `top_n` (default 10) IDs per rule.
4. Union all rule ID sets; annotate that combined set on the scatter.

### Pointer arrow mapping

`annotate_outliers()` highlights union IDs, creates `ax.text()` labels, then calls `adjust_text()` with `arrowprops` drawing thin gray lines from each label to its point.

---

## Distribution pipeline (`properties_histogram.py`)

### CONFIG keys

| Key | Default | Description |
| --- | --- | --- |
| `input_dir` | PeptiVerse results path | Shard CSV directory |
| `output_dir` | `plots/` | Pre-existing output folder |
| `score_columns` | 5 pharmacological scores | One figure per column |
| `figsize_inches` | `(7, 7)` | Figure dimensions |
| `hist_bins` | `50` | Histogram bin count |
| `metric_labels` | dict | Human-readable axis titles |

### Function reference

| Function | Purpose |
| --- | --- |
| `load_property_data(config)` | Merge sharded CSVs |
| `compute_summary_stats(values)` | n, mean, std, min, max, quartiles |
| `plot_univariate_distribution(df, col, config)` | Histogram + KDE + reference lines |
| `main()` | Generate all five `dist_*.png` files |

---

## Environment Setup

Create the isolated Conda environment from the pinned spec:

```bash
cd /home/cle1g21/RPC/peptiverse
conda env create -f environment.yml
conda activate pep_vis_env
```

To update an existing environment:

```bash
conda env update -f environment.yml --prune
```

Pinned packages in `environment.yml`:

| Package | Version |
| --- | --- |
| Python | 3.11.9 |
| pandas | 2.2.3 |
| numpy | 2.1.3 |
| seaborn | 0.13.2 |
| matplotlib | 3.9.2 |
| adjustText | 1.3.0 (pip) |

---

## Data Input Schema

Each shard CSV must contain:

| Column | Type | Required | Description |
| --- | --- | --- | --- |
| `id` | `str` | Yes | Unique ORF identifier (configurable via `id_column`) |
| `hemolysis_score` | `float` | Yes | Hemolysis probability (0–1) |
| `nf_score` | `float` | Yes | Non-fouling probability (0–1) |
| `solubility_score` | `float` | Yes | Solubility probability (0–1) |
| `permeability_penetrance_score` | `float` | Yes | Membrane penetrance probability (0–1) |
| `halflife_score` | `float` | Yes | Predicted half-life in hours |

Rows with NaN in any score column are dropped before plotting.

---

## Interpreting the Plots

| Property | Score meaning |
| --- | --- |
| `hemolysis_score` | Higher = greater red blood cell lysis risk |
| `nf_score` | Higher = more non-fouling (desirable for circulation) |
| `solubility_score` | Higher = more soluble in aqueous conditions |
| `permeability_penetrance_score` | Higher = more membrane permeable |
| `halflife_score` | Predicted serum half-life in hours; desirability depends on drug design goals |

---

## Troubleshooting

| Issue | Solution |
| --- | --- |
| `FileNotFoundError: Output directory does not exist` | Create `/home/cle1g21/RPC/peptiverse/plots/` before running (the script will not create it) |
| `No input files matched pattern` | Check `input_dir` and `input_glob` in CONFIG |
| `Missing required columns` | Verify CSV shards contain `id` and all five score columns |
| Slow adjustText on dense panels | Reduce `--top-n` or set `max_label_chars` in CONFIG |
| `adjustText` import error | Run `pip install adjustText==1.3.0` inside `pep_vis_env` |

---

## Citation

PeptiVerse models and property definitions: [ChatterjeeLab/PeptiVerse](https://huggingface.co/ChatterjeeLab/PeptiVerse)
