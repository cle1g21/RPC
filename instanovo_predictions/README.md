# InstaNovo Predictions and ROC Validation

This folder holds InstaNovo prediction outputs and the amino-acid ROC validation pipeline.

**GitHub:** prediction CSVs, `chunks/`, `filtered_predictions/`, `predictions/`, `plots/`, and `logs/` are excluded by the repository-root `.gitignore` (they are multi-GB). Push only the ROC scripts and this README. Place local CSVs at the paths in the table below before running `generate_roc_curve.py`.

## Resilient prediction pipeline (large datasets)

Use the resilient launcher for production runs with streaming output, checkpointing, and automatic Slurm resubmission.

```bash
cd /iridisfs/ddnb/Charlotte/RPC/InstaNovo

# Run until complete (auto-resubmits on timeout)
./scripts/slurm/run_until_complete.sh /home/cle1g21/RPC/msconvert/output_data/Control4_Neo_SN_114_HLA-I.mgf
./scripts/slurm/run_until_complete.sh /home/cle1g21/RPC/msconvert/output_data/01625b_GB2-TUM_first_pool_10_01_01-3xHCD-1h-R1.mgf
```

| Dataset | Input MGF | Output CSV |
|---|---|---|
| Control4 HLA-I | `.../Control4_Neo_SN_114_HLA-I.mgf` | `Control4_Neo_SN_114_HLA-I_predictions.csv` |
| GB2-TUM pool | `.../01625b_GB2-TUM_first_pool_10_01_01-3xHCD-1h-R1.mgf` | `01625b_GB2-TUM_first_pool_10_01_01-3xHCD-1h-R1_predictions.csv` |

**Route matrix**

| File size | Route | What happens |
|---|---|---|
| ≥ 1 GB | Chunked | Split into 50k-spectrum chunks on scratch; parallel GPU jobs; compile at end |
| < 1 GB | Direct | One GPU job streams to the final CSV |

**Checkpointing:** reruns skip spectra already present in the output CSV (matched by `spectrum_id`).

**Monitor progress:**

```bash
python /iridisfs/ddnb/Charlotte/RPC/InstaNovo/scripts/run_instanovo_resilient.py --status --input /path/to/file.mgf
squeue -u $USER
```

---

## Amino-acid ROC validation

Validate InstaNovo de novo predictions against MaxQuant ground truth and produce publication-ready ROC curves at the **individual amino acid** level.

## Environment setup

Create an isolated conda environment (run once):

```bash
conda create -y -n instanovo_eval python=3.11
conda activate instanovo_eval
conda install -y pandas numpy scikit-learn seaborn matplotlib
```

If `conda` is unavailable, use micromamba:

```bash
micromamba create -y -n instanovo_eval -c conda-forge python=3.11 pandas numpy scikit-learn seaborn matplotlib
micromamba activate instanovo_eval
```

## File paths

| Role | Path |
|---|---|
| InstaNovo predictions | `/home/cle1g21/RPC/instanovo_predictions/01625b_GB2-TUM_first_pool_10_01_01-3xHCD-1h-R1_predictions.csv` (full GB2 run) |
| MGF reference | `/home/cle1g21/RPC/msconvert/output_data/01625b_GB2-TUM_first_pool_10_01_01-3xHCD-1h-R1.mgf` |
| MaxQuant tables | `/iridisfs/ddnb/Charlotte/RPC/InstaNovo/test_data/` |
| Baseline ROC plot | `/home/cle1g21/RPC/instanovo_predictions/plots/instanovoVal_roc_curve.png` |
| Enhanced ROC plot | `/home/cle1g21/RPC/instanovo_predictions/plots/instanovoVal_roc_curveEnhanced.png` |

## Column mapping (locked for this dataset)

| Column | Role |
|---|---|
| `predictions` | Standard model peptide sequence (alignment target) |
| `predictions_tokenised` | Comma-separated residue tokens (used when char tokenization fails) |
| `instanovo_prediction_token_log_probabilities` | Per-amino-acid log probs → `exp()` for ROC scores |
| `instanovo_prediction_log_probability` | **Fallback** for token log probs when the streaming schema is shifted (see below) |
| `token_log_probs` | **Ignored** — empty (NaN) |
| `instanovoplus_*` | **Ignored** — empty without diffusion refinement |

### GB2 CSV schema note

The full GB2 predictions CSV starts with ~500 CombinedPredictor (`*_subset500`) rows with a correct column layout, then continues with the resilient streaming run whose fields are shifted by one starting at `group`. For those rows, numeric token log-probs sit in `instanovo_prediction_log_probability` instead of `instanovo_prediction_token_log_probabilities`.

`generate_roc_curve.py` handles this automatically: it filters to `experiment_name == raw_file_name` (drops subset500 duplicates) and falls back to the shifted column when the primary field is not a numeric list.

## Workflow

### Step 1 — Explore ground truth

```bash
conda activate instanovo_eval
cd /home/cle1g21/RPC/instanovo_predictions
python explore_maxquant_ground_truth.py
```

### Step 2 — Diagnose token mismatches (optional)

Use this when spectra are skipped due to length mismatch between sequence and token probabilities:

```bash
python generate_roc_curve.py --diagnose-mismatch
```

This prints the first few mismatched spectra with char vs tokenised vs prob lengths, plus summary counts.

### Step 3 — Generate ROC curves (baseline + enhanced)

```bash
python generate_roc_curve.py
```

This writes:

- `plots/instanovoVal_roc_curve.png` — single strict ROC (backward compatible)
- `plots/instanovoVal_roc_curveEnhanced.png` — multi-panel figure with strict + mass-equivalent ROC and stratified sub-panels

Override paths or scoring options:

```bash
python generate_roc_curve.py \
  --predictions-csv /path/to/predictions.csv \
  --output-plot /path/to/plots/my_roc.png \
  --output-enhanced-plot /path/to/plots/my_roc_enhanced.png \
  --sample-label "my_cohort" \
  --mass-equiv-kq
```

## CONFIG dictionary (`generate_roc_curve.py`)

Edit the `CONFIG` block at the top of the script for new cohorts:

| Key | Type | Default | Meaning |
|---|---|---|---|
| `predictions_csv` | str | subset500 CSV path | InstaNovo output CSV |
| `mgf_path` | str | MGF path | Source MGF for scan-index mapping |
| `maxquant_search_dir` | str | test_data path | Folder containing `msms.txt`, `msmsScans.txt` |
| `raw_file_name` | str | raw file stem | MaxQuant `Raw file` column value |
| `sample_label` | str | cohort label | Shown in plot titles |
| `scan_index_column` | str | `scan_number` | InstaNovo spectrum index column |
| `prediction_sequence_column` | str | `predictions` | Sequence column |
| `prediction_tokenised_column` | str | `predictions_tokenised` | Comma-separated token column |
| `token_log_probs_column` | str | `instanovo_prediction_token_log_probabilities` | Token log-prob column |
| `max_spectra` | int | `500` | Number of spectra to evaluate |
| `output_plot_png` | str | baseline PNG path | Single strict ROC output |
| `output_plot_enhanced_png` | str | enhanced PNG path | Multi-panel ROC output |
| `mass_equivalent_li` | bool | `True` | Treat L↔I as a match |
| `mass_equivalent_kq` | bool | `False` | Treat K↔Q as a match (enable for low-res instruments) |
| `min_aa_per_stratum` | int | `20` | Minimum amino acids required to plot a stratum |
| `length_bins` | dict | short &lt;10, medium 10–15, long &gt;15 | Peptide length strata |
| `charge_bins` | dict | +2 vs ≥+3 | Precursor charge strata |

CLI flags: `--mass-equiv-li`, `--no-mass-equiv-li`, `--mass-equiv-kq`, `--diagnose-mismatch`, `--output-enhanced-plot`.

## How scan mapping works

1. InstaNovo `scan_number` = 0-based index into the MGF file.
2. MGF `TITLE` field encodes the instrument scan: `{rawfile}.{scan}.{ms_scan}.{charge}`.
3. Join to MaxQuant `msmsScans.txt` on `Scan number`.
4. Fallback: join `msms.txt` on precursor `m/z` + `Charge` (highest Score).

## Token alignment (length-mismatch fix)

The baseline script tokenized `predictions` as one character per residue. That fails on modifications (e.g. `TVC[UNIMOD:4]VSQLK` has 18 characters but 8 residue tokens).

The enhanced script aligns tokens and probabilities in this order:

1. Parse `predictions_tokenised` (comma-separated residue tokens including modifications).
2. Apply InstaNovo residue regex on `predictions` (`A` or `X[UNIMOD:n]`).
3. Use per-character tokenization only when `predictions` contains no `[` (plain peptides).

If lengths still differ (e.g. refined sequence vs beam-0 token probs), the script uses prefix alignment over `min(len(tokens), len(probs))` positions and logs these as truncated. Spectra with no usable overlap are skipped and counted.

After this fix, more modified peptides contribute to ROC. Strict AUC and accuracy may shift slightly upward because previously skipped spectra are included.

## Mass-equivalence rules

Two label arrays are built from the same confidence scores:

- **Strict** — exact letter match at each position (baseline behaviour).
- **Mass-equivalent** — strict match, or L↔I when `mass_equivalent_li=True` (default), or K↔Q when `mass_equivalent_kq=True` (default off, appropriate for FTMS).

Both strict and mass-equivalent AUC and position accuracy are printed. The enhanced plot overlays both global ROC curves.

## Stratification

Each scored amino acid is tagged with:

| Stratum | Rule | Source |
|---|---|---|
| `length_short` | &lt; 10 residues | tokenized prediction length |
| `length_medium` | 10–15 inclusive | same |
| `length_long` | &gt; 15 | same |
| `charge_2` | precursor charge == 2 | MGF merge |
| `charge_3plus` | precursor charge ≥ 3 | MGF merge |

Per-stratum strict and mass-equivalent AUC are printed when the stratum has at least `min_aa_per_stratum` amino acids (default 20). Sub-panels in the enhanced figure show strict ROC by length and charge; strata below the minimum show "Insufficient data".

## ROC methodology

For each spectrum with ground truth:

1. Align prediction tokens with token log probabilities (see above).
2. Tokenize MaxQuant sequence (regex for modifications, characters otherwise).
3. At each aligned position: score = `exp(log_prob_i)`, strict label = 1 if letters match, mass-equiv label uses equivalence rules.
4. Pool all amino acids; compute `sklearn.metrics.roc_curve` and AUC for strict and mass-equivalent labels.

**Position accuracy** = fraction of positions with a matching label (ignores confidence). **AUC** = how well confidence scores rank correct vs incorrect positions.

## Script reference

### `explore_maxquant_ground_truth.py`

| Function | Purpose |
|---|---|
| `print_table_headers` | Print MaxQuant column names |
| `parse_mgf_scan_table` | Map MGF index → instrument scan |
| `load_msms_scans_ground_truth` | Primary ground truth from msmsScans.txt |
| `load_msms_fallback_ground_truth` | Fallback via msms.txt |
| `report_ground_truth_coverage` | Print join statistics |

### `generate_roc_curve.py`

| Function | Purpose |
|---|---|
| `align_tokens_and_probs` | Tokenize and align tokens with prob array |
| `diagnose_token_mismatches` | Print mismatch diagnostics (`--diagnose-mismatch`) |
| `align_predictions_with_ground_truth` | Merge predictions + MaxQuant |
| `build_observation_table` | Per-aa observations with strict and mass-equiv labels |
| `compute_roc_from_observations` | ROC arrays and AUC |
| `compute_stratified_metrics` | Per-bin AUC table |
| `plot_roc_curve` | Baseline single-panel PNG |
| `plot_enhanced_roc_figure` | Multi-panel enhanced PNG |

## Known limitations

- Prefix truncation may drop trailing residues when refined `predictions` length differs from beam-0 token probs.
- Stratified AUC on small bins can be unstable; `min_aa_per_stratum=20` avoids misleading sub-panel curves.
- Not all MGF spectra have MaxQuant identifications; only matched spectra contribute.
- K/Q mass equivalence is off by default (FTMS); enable with `--mass-equiv-kq` for low-resolution data.
