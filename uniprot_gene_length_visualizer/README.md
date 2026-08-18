# UniProt Gene Length Visualizer

Production-grade Python pipeline that downloads reviewed UniProt reference proteomes for five model organisms, converts amino-acid sequence lengths to coding base-pair lengths ($1\ \mathrm{AA} = 3\ \mathrm{bp}$), and exports publication-ready continuous Kernel Density Estimate (KDE) abundance curves with a schematic annotation arrow at **1,500 bps**.

**Scientific note:** Coding / ORF base-pair length is approximated strictly as $\mathrm{length_{bp}} = \mathrm{len}(\mathrm{seq_{AA}}) \times 3$. This does **not** include UTRs, introns, or untranslated regions.

## Target proteomes

| Organism | Common name | UniProt proteome ID |
|---|---|---|
| *Homo sapiens* | Human | `UP000005640` |
| *Mus musculus* | Mouse | `UP000000589` |
| *Arabidopsis thaliana* | Thale cress | `UP000006548` |
| *Drosophila melanogaster* | Fruit fly | `UP000000803` |
| *Saccharomyces cerevisiae* | Baker's yeast | `UP000002311` |

Each proteome is streamed from UniProt with:

```text
https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3A{PROTEOME_ID}%29%20AND%20%28reviewed%3Atrue%29
```

## Repository layout

```text
uniprot_gene_length_visualizer/
├── environment.yml             # Conda environment definition
├── README.md                   # This documentation
├── run_visualization.py        # Master entry point (UniProt KDEs)
├── scripts/
│   └── plot_master_low_summary.py  # Low-tier length + confidence summary
├── data/                       # Cached UniProt FASTA downloads (gitignored)
├── outputs/                    # Exported figures and summary artifacts
├── config/
│   ├── __init__.py
│   └── config.py               # Proteome IDs, annotation BP, KDE styles, DPI
└── src/
    ├── __init__.py
    ├── fetch_data.py           # UniProt streamer and SeqIO AA×3 parser
    └── plot_distributions.py   # Seaborn KDE plotting and annotation
```

## Phase 1: Conda environment setup

```bash
cd /home/cle1g21/RPC/uniprot_gene_length_visualizer
conda env create -f environment.yml
conda activate uniprot_viz_env
```

Equivalent one-liner install:

```bash
conda create -n uniprot_viz_env python=3.11 -y
conda activate uniprot_viz_env
conda install -c conda-forge pandas seaborn matplotlib biopython requests -y
```

## Running the pipeline

```bash
conda activate uniprot_viz_env
cd /home/cle1g21/RPC/uniprot_gene_length_visualizer
python run_visualization.py
```

Force fresh UniProt downloads:

```bash
python run_visualization.py --force-refetch
```

### Expected outputs

| Artifact | Path |
|---|---|
| Cached FASTA files | `data/UP00000*_*.fasta` |
| Length summary CSV | `outputs/length_summary.csv` |
| Single-panel linear KDE (raster) | `outputs/gene_length_kde_density.png` |
| Single-panel linear KDE (vector) | `outputs/gene_length_kde_density.svg` |
| Dual-panel linear vs log KDE (raster) | `outputs/gene_length_kde_linear_vs_log.png` |
| Dual-panel linear vs log KDE (vector) | `outputs/gene_length_kde_linear_vs_log.svg` |

Cached FASTA files under `data/` are gitignored. Re-download them with `python run_visualization.py` (or `--force-refetch`) after cloning. Figure PNG/SVG files under `outputs/` are kept in git when present.

## Plot style notes

- Theme: `sns.set_theme(style="ticks")`
- Borders: `sns.despine()`; left/bottom spines thickened (`linewidth=2.5`, black)
- Curves: `sns.kdeplot(..., hue="species", common_norm=False, bw_adjust=0.5)`
- Dual panel: Panel A linear; Panel B `log_scale=True` ($\log_{10}$ bps)
- Annotation: downward arrow + `"1,500 bps"` label at `ANNOTATION_BP = 1500` on both panels
- Axes: `"ORF length (bps)"` × `"Density Abundance"`; no background gridlines

## Why logarithmic ORF-length scales?

Protein coding lengths are typically **right-skewed** on a linear axis: most mass sits below ~2 kb while a long tail of giant ORFs stretches the x-axis and visually compresses short populations. A $\log_{10}(\mathrm{bps})$ transform:

1. **Resolves short ORFs (&lt;500 bp)** that otherwise crowd near zero on a 0–8,000 linear axis.
2. **Compresses long-tail giants**, so multi-kilobase proteins no longer dominate the visual field.
3. Often reveals near **log-normal symmetry** around the typical coding-length mode, making median peaks and the 1,500 bp inflection easier to compare across species.

Panel A therefore shows raw linear decay; Panel B shows the same data under log scaling for short-ORF resolution and shape comparison.

## Function catalog

| Function | Module | Arguments | Inputs | Outputs |
|---|---|---|---|---|
| `build_stream_url` | `src/fetch_data.py` | `proteome_id: str` | UniProt proteome accession | Stream URL string |
| `cache_is_valid` | `src/fetch_data.py` | `path: Path`, `min_bytes: int = MIN_CACHE_BYTES` | Local FASTA path | `True` if cache is reusable |
| `download_proteome` | `src/fetch_data.py` | `proteome_id: str`, `dest_path: Path`, `force: bool = False`, `min_cache_bytes: int` | UniProt REST API / local cache | Path to cached `.fasta` under `data/` |
| `parse_coding_bp_lengths` | `src/fetch_data.py` | `fasta_path: Path` | Cached FASTA file | DataFrame with `entry_id`, `aa_length`, `bp_length` |
| `build_species_dataframe` | `src/fetch_data.py` | `species_meta: Mapping`, `force: bool = False` | One `PROTEOMES` entry | Tidy species DataFrame (`species`, `proteome_id`, `entry_id`, `aa_length`, `bp_length`) |
| `load_all_proteomes` | `src/fetch_data.py` | `proteomes: Sequence \| None = None`, `force: bool = False` | `PROTEOMES` catalog / UniProt / `data/` | Concatenated multi-species DataFrame |
| `apply_publication_style` | `src/plot_distributions.py` | *(none)* | `SNS_STYLE`, `SNS_PALETTE` from config | Global Seaborn theme applied |
| `save_figure` | `src/plot_distributions.py` | `fig: Figure`, `stem: str`, `output_dir: Path \| None = None` | Matplotlib figure | `[stem.png, stem.svg]` under `outputs/` at `DPI=300` |
| `_style_minimal_spines` | `src/plot_distributions.py` | `ax: Axes` | Matplotlib axes | Thick black left/bottom spines; top/right removed |
| `_annotation_y_position` | `src/plot_distributions.py` | `ax: Axes`, `annotation_bp: float` | Plotted KDE line artists | Y height near the annotation x-position |
| `annotate_inflection` | `src/plot_distributions.py` | `ax: Axes`, `annotation_bp`, `label`, `y_text_factor` | KDE axes | Downward arrow annotation at the inflection |
| `_draw_kde_panel` | `src/plot_distributions.py` | `ax`, `df`, `title`, `x_lim`, `log_scale`, `legend`, `y_text_factor` | Length DataFrame + axes | Styled KDE panel with optional legend |
| `plot_kde_density` | `src/plot_distributions.py` | `df: DataFrame`, `output_stem: str = KDE_STEM`, `output_dir: Path \| None = None` | Continuous multi-species length DataFrame | Single-panel linear KDE PNG + SVG |
| `plot_kde_linear_vs_log` | `src/plot_distributions.py` | `df: DataFrame`, `output_stem: str = DUAL_KDE_STEM`, `output_dir: Path \| None = None` | Continuous multi-species length DataFrame | Dual-panel linear vs log PNG + SVG |
| `parse_args` | `run_visualization.py` | *(CLI)* | `sys.argv` | `argparse.Namespace` (`--force-refetch`) |
| `write_summary_csv` | `run_visualization.py` | `df`, `output_dir: Path` | Length DataFrame | `outputs/length_summary.csv` |
| `main` | `run_visualization.py` | *(none)* | Config + UniProt + caches | Exit code `0`; figures + summary CSV |

## Customizing the 1,500 bp annotation and KDE style

Edit [`config/config.py`](config/config.py):

| Constant | Effect |
|---|---|
| `ANNOTATION_BP` | X-position of the annotation arrow (default `1500`) |
| `ANNOTATION_LABEL` | Text drawn above the arrow (default `"1,500 bps"`) |
| `ANNOTATION_Y_TEXT_FACTOR_LINEAR` / `_LOG` | Vertical label offset multipliers per panel |
| `BW_ADJUST` | Seaborn KDE bandwidth scalar (default `0.5`) |
| `X_LIM_MIN` / `X_LIM_MAX` | Visible ORF-length window for linear panels |
| `LOG_X_LIM_MIN` / `LOG_X_LIM_MAX` | Visible ORF-length window for log panels (`LOG_X_LIM_MIN` must be > 0) |
| `ENABLE_LOG_PANEL` | When `True`, also export the dual-panel comparative figure |
| `SNS_STYLE` | Passed to `sns.set_theme(style=...)` (default `"ticks"`) |
| `SNS_PALETTE` | Passed to `sns.set_theme(palette=...)` (default `"muted"`) |
| `SPECIES_PALETTE` | Hex color map keyed by species display name |
| `SPINE_LINEWIDTH` | Thickness of left/bottom spines (default `2.5`) |
| `SPINE_COLOR` | Color of remaining spines (default `"black"`) |
| `KDE_FIGSIZE` / `DUAL_KDE_FIGSIZE` | Single- and dual-panel figure sizes |
| `DPI` | PNG export resolution (default `300`) |
| `KDE_STEM` / `DUAL_KDE_STEM` | Output filename stems without extension |
| `PANEL_A_TITLE` / `PANEL_B_TITLE` | Dual-panel subplot titles |
| `KDE_TITLE`, `X_AXIS_LABEL`, `Y_AXIS_LABEL` | Title and axis strings |

After changing constants, re-run:

```bash
python run_visualization.py
```

Cached FASTAs in `data/` are reused unless you pass `--force-refetch`.

## Pipeline data flow

1. **Ingest:** stream reviewed FASTA proteomes from UniProt (or reuse `data/` caches).
2. **Convert:** parse with `Bio.SeqIO`; compute `bp_length = aa_length * 3`.
3. **Plot:** single-panel linear KDE, plus dual-panel linear vs $\log_{10}$ KDE when `ENABLE_LOG_PANEL` is true.
4. **Export:** `outputs/gene_length_kde_density.{png,svg}`, `outputs/gene_length_kde_linear_vs_log.{png,svg}`, and `outputs/length_summary.csv`.

## Low-tier `master_low` summary plots

Standalone immunopeptidomics summary for the assembled Low confidence tier (InstaNovo peptide ⊂ NTv3 ORF). Parses peptide lengths, bins them into HLA-oriented AA intervals, summarizes prediction confidence, and exports a 2-panel Seaborn figure.

### Input path and column dependencies

| Role | Default / auto-detect | Notes |
|---|---|---|
| Input CSV | `/home/cle1g21/RPC/confidence_levels/assembled/master_low.csv` | Prefer `assembled/`, not `low/assembled/` |
| Sequence | `predictions`, then `sequence` | Bracketed mods (`[UNIMOD:n]`) are stripped before length |
| Score | `confidence` → `preds_score` → `score` → `log_probs` → `instanovoplus_prediction_log_probability` | Current master resolves to `log_probs` |
| Display confidence | `exp(log_probs)` when the score name contains `log_prob` | Reported on a (0, 1] scale in Panel B and the stats log |

### Length bins

| Bin label | AA length |
|---|---|
| `<8 AA` | ≤ 7 |
| `8-11 AA [HLA Class I]` | 8–11 |
| `12-15 AA [HLA Class II]` | 12–15 |
| `16-20 AA` | 16–20 |
| `21-30 AA` | 21–30 |
| `>30 AA` | ≥ 31 |

### How to run

```bash
conda activate uniprot_viz_env
cd /home/cle1g21/RPC/uniprot_gene_length_visualizer
python scripts/plot_master_low_summary.py \
  --input /home/cle1g21/RPC/confidence_levels/assembled/master_low.csv \
  --output-dir outputs
```

Optional `--dpi` (default `300`) controls PNG resolution.

### Expected outputs

| Artifact | Path |
|---|---|
| Summary stats (min/max/mean/median, global + per bin) | `outputs/master_low_summary_stats.txt` |
| 2-panel figure (raster) | `outputs/master_low_length_and_accuracy.png` |
| 2-panel figure (vector) | `outputs/master_low_length_and_accuracy.svg` |

- **Panel A:** bar / count plot of peptide abundance per length bin with count and percentage annotations.
- **Panel B:** box plot of display confidence per length bin with mean (diamond) and median (circle) overlays.

### Customizing bins, colors, and spines

Edit constants at the top of [`scripts/plot_master_low_summary.py`](scripts/plot_master_low_summary.py):

| Constant | Effect |
|---|---|
| `_BIN_EDGES` / `_BIN_LABELS` | `pd.cut` edges and categorical axis labels |
| `_BAR_COLOR` / `_BOX_COLOR` | Panel A / Panel B fill colors |
| `_MEAN_COLOR` / `_MEDIAN_COLOR` | Overlay marker colors on Panel B |
| `_SPINE_LINEWIDTH` | Left/bottom spine thickness after `sns.despine()` (default `2.0`) |
| `_FIGURE_STEM` / `_STATS_FILENAME` | Output basename stems |
| `_DEFAULT_INPUT` / `_DEFAULT_OUTPUT_DIR` / `_DEFAULT_DPI` | CLI defaults |
