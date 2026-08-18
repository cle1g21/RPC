# PRIDEpull

Automated pipeline to harvest **human (Homo sapiens, taxonomy 9606) immunopeptidomics** datasets from the [PRIDE Archive](https://www.ebi.ac.uk/pride/) API, route native MGF or vendor RAW/.d files through msconvert when needed, run **InstaNovo** de novo sequencing under strict memory constraints, reclaim disk space via delete-on-success cleanup, and maintain a persistent audit ledger.

`PRIDEpull/data/`, `cache/`, `harvest_cache.json`, and `*.partial` downloads are gitignored. Prediction CSVs written under `instanovo_predictions/` are also gitignored at the repository root.

## Quick start

```bash
cd /home/cle1g21/RPC/PRIDEpull
conda env create -f environment.yml
conda activate pride_pull_env
module load apptainer/1.5.0   # required before msconvert fallback

# 1. Harvest human immunopeptidomics projects (writes cache/harvest_cache.json)
python run_pride_pipeline.py --harvest-only -v

# 2. Subset sanity check on one project (default: 500 spectra, GPU Slurm job)
python run_pride_pipeline.py --max-projects 1 -v

# 3. Full production run on a known project
python run_pride_pipeline.py --full-run --accession PXD077095 -v

# 4. Monitor Slurm jobs and audit trail
squeue -u $USER
tail -f processing_history_ledger.txt
```

## Architecture

```text
PRIDE API  →  pride_harvester  →  manifest cache
                      ↓
              pipeline_runner
                      ↓
         file_downloader → PRIDEpull/data/{PXD}/
                      ↓
         ┌────────────┴────────────┐
    native .mgf              .raw / .d only
         │                        │
         │                   msconvert (serialized)
         └────────────┬────────────┘
                      ↓
              [optional subset 500 spectra]
                      ↓
         run_until_complete.sh (GPU Slurm sbatch)
                      ↓
    instanovo_predictions/predictions/{PXD}_{stem}_predictions.csv
                      ↓
              delete data/ intermediates
                      ↓
         processing_history_ledger.txt
```

### Native MGF vs msconvert handoff

| Condition | Route | Action |
|-----------|-------|--------|
| PRIDE provides `.mgf` | **Native MGF streamed direct** | Download to `data/` → optional subset → InstaNovo (no msconvert) |
| Only `.raw` or `.d` | **Vendor RAW converted via msconvert** | Download → `convert_ms.py` (one at a time) → MGF in `data/` → InstaNovo |
| Neither | **Skip** | Project excluded from manifest |

### InstaNovo execution (GPU Slurm by default)

When `INSTANOVO_USE_SLURM=True` (default), **all** runs — including subset sanity checks — go through `run_until_complete.sh`:

| MGF size | Slurm route | What happens |
|----------|-------------|--------------|
| < 1 GB (subset or small full) | Route II | Single GPU `sbatch` job via `submit_direct_prediction.sh` |
| ≥ 1 GB (large full) | Route I | Chunked parallel GPU jobs + compile |

The pipeline blocks until Slurm jobs finish, then copies output to `predictions/{PXD}_{stem}_predictions.csv`.

Set `INSTANOVO_USE_SLURM=False` in `config/config.py` only for local debugging (runs on login-node CPU).

### Subset-first validation

When `RUN_SUBSET_ONLY=True` (default):

1. Slice the first `SUBSET_SPECTRUM_COUNT` (500) MGF spectrum blocks
2. Submit the subset MGF to InstaNovo via **GPU Slurm** (not login-node CPU)
3. Use `--full-run` to predict the full MGF instead

### Storage lifecycle

| Location | Contents | After success |
|----------|----------|---------------|
| `PRIDEpull/data/{PXD}/` | Downloads, converted MGF, subset MGF | **Deleted** (`DELETE_ON_SUCCESS=True`) |
| `instanovo_predictions/predictions/` | Final CSV predictions | **Retained** |
| `processing_history_ledger.txt` | Audit blocks | **Appended** |

Use `--no-cleanup` to keep `data/` files for debugging.

## Configuration

All paths, keywords, and tuning live in [`config/config.py`](config/config.py). Key parameters:

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `DATA_LANDING_DIR` | `PRIDEpull/data` | Download and conversion workspace |
| `RUN_SUBSET_ONLY` | `True` | Subset-first sanity gate |
| `SUBSET_SPECTRUM_COUNT` | `500` | Spectra in subset MGF |
| `MAX_FILES_PER_PROJECT` | `1` | Cap files per PXD (largest MGF first) |
| `DELETE_ON_SUCCESS` | `True` | Reclaim disk after prediction |
| `INSTANOVO_USE_SLURM` | `True` | Submit all InstaNovo runs as GPU Slurm jobs |
| `INSTANOVO_BATCH_SIZE` | `32` | Batch size for local direct runs only |
| `INSTANOVO_LARGE_FILE_GB` | `1.0` | Chunked vs single-job threshold inside Slurm launcher |

## Function reference

| Module | Function | Inputs | Outputs | Purpose |
|--------|----------|--------|---------|---------|
| `pride_api.py` | `create_session()` | — | `Session` | HTTP session with User-Agent |
| `pride_api.py` | `_request_with_retry()` | session, url, params | `Response` | Exponential backoff on 429/5xx |
| `pride_api.py` | `search_projects_page()` | session, page, keyword, organism_facet | `list[dict]` | Paginated project search |
| `pride_api.py` | `fetch_all_project_files()` | session, accession | `list[dict]` | Complete file inventory |
| `pride_api.py` | `extract_download_url()` | file_record | `str \| None` | Best FTP/HTTPS URL |
| `pride_harvester.py` | `harvest_immunopeptidomics_projects()` | refresh, max_projects | `list[dict]` | Full harvest + cache write |
| `pride_harvester.py` | `match_immuno_keywords()` | project dict | `list[str]` | Keyword matrix matches |
| `pride_harvester.py` | `is_human_project()` | project dict | `bool` | Homo sapiens validation |
| `pride_harvester.py` | `inventory_project_files()` | session, accession | file buckets | Classify mgf/raw/d |
| `pride_harvester.py` | `select_project_files()` | buckets, max_files | routing, files | Largest-first selection |
| `pride_harvester.py` | `load_harvest_cache()` | — | `list[dict]` | Load cached manifest |
| `file_downloader.py` | `download_pride_file()` | url, dest, expected_size | `Path` | Chunked streaming download |
| `file_downloader.py` | `project_download_path()` | accession, fileName | `Path` | Standard `data/` path |
| `file_converter.py` | `convert_raw_to_mgf()` | raw_path, output_dir | `Path` | Synchronous msconvert |
| `file_converter.py` | `slice_mgf()` | mgf_path, output_path, max_spectra | `Path` | First-N-spectra subset |
| `file_converter.py` | `count_mgf_spectra()` | mgf_path | `int` | Spectrum block count |
| `instanovo_predictor.py` | `run_slurm_predict()` | mgf_path, output_csv | `Path` | GPU Slurm via run_until_complete.sh |
| `instanovo_predictor.py` | `run_prediction()` | accession, mgf_path | `Path` | Route Slurm vs local direct |
| `instanovo_predictor.py` | `build_output_csv_path()` | accession, mgf_path | `Path` | `{PXD}_{stem}_predictions.csv` |
| `audit_logger.py` | `append_ledger_entry()` | entry dict | — | Append block to ledger |
| `audit_logger.py` | `compute_prediction_stats()` | csv_path | stats dict | Spectrum/peptide/confidence |
| `audit_logger.py` | `cleanup_data_files()` | file paths | confirmations | `os.remove()` / `rmtree` |
| `pipeline_runner.py` | `process_single_file()` | accession, file_entry, routing | ledger dict | One file end-to-end |
| `pipeline_runner.py` | `run_pipeline()` | manifest, flags | summary dict | Full orchestration loop |

## CLI reference

| Flag | Effect |
|------|--------|
| `--harvest-only` | Query PRIDE, write cache, exit |
| `--refresh-harvest` | Ignore existing harvest cache |
| `--accession PXD…` | Process single project |
| `--max-projects N` | Cap batch size |
| `--full-run` | Disable subset mode |
| `--skip-download` | Use files already in `data/` |
| `--no-cleanup` | Keep intermediates after success |
| `-v` / `--verbose` | DEBUG logging |

## Ledger format

Each successful or failed run appends a block to [`processing_history_ledger.txt`](processing_history_ledger.txt). **PRIDE accession** and **original filename** appear first:

```text
================================================================================
PRIDE_ACCESSION: PXD077095
ORIGINAL_PRIDE_FILENAME: P1037.mgf
TIMESTAMP: 2026-07-03T19:45:00Z
PROCESSING_ROUTE: Native MGF streamed direct
RUN_MODE: subset (500 spectra)
--------------------------------------------------------------------------------
DOWNLOAD_SIZE_BYTES: 126541859
...
SPECTRA_PROCESSED: 500
PEPTIDES_SEQUENCED: 487
MEAN_PREDICTION_LOG_PROB: -1.23
...
CLEANUP_DELETED_FILES:
  - /home/cle1g21/RPC/PRIDEpull/data/PXD077095/P1037.mgf [DELETED]
CLEANUP_STATUS: SUCCESS
STATUS: COMPLETED
================================================================================
```

## Output naming

Predictions are written to:

```text
/home/cle1g21/RPC/instanovo_predictions/predictions/{PXD}_{file_stem}_predictions.csv
```

Example: `PXD077095_P1037_predictions.csv`

## Downstream integration

Filtered prediction CSVs in `instanovo_predictions/filtered_predictions/` feed the confidence pipeline:

```bash
python /home/cle1g21/RPC/confidence_pipeline/run_confidence_pipeline.py --verbose
```

## Safety rules

- **Do not modify** `/home/cle1g21/RPC/pride_species_visualizer` — PRIDE API patterns are vendored into `src/pride_api.py`.
- **No hardcoded paths in source modules** — edit `config/config.py` only.
- Use `/projects/{accession}/files` for file inventory — avoid `/files/all` (returns entire catalog).

## Repository layout

```text
PRIDEpull/
├── environment.yml
├── README.md
├── run_pride_pipeline.py
├── processing_history_ledger.txt
├── data/                    # Active downloads (gitignored)
├── cache/                   # harvest_cache.json (gitignored)
├── config/config.py
└── src/
    ├── pride_api.py
    ├── pride_harvester.py
    ├── file_downloader.py
    ├── file_converter.py
    ├── instanovo_predictor.py
    ├── audit_logger.py
    └── pipeline_runner.py
```
