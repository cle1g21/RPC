# Confidence Pipeline: Proteome Filtering + 4-Tier Validation

Production-grade pipeline for filtering NTv3 ORF predictions against the human canonical UniProt proteome and validating novel de novo InstaNovo peptides across four confidence tiers.

## Overview

| Phase | Purpose |
|-------|---------|
| **Phase 1** | Download/cache UniProt canonical proteome; substring-filter NTv3 ORF handoff files |
| **Tier 1 Low** | InstaNovo filtered peptides that appear in filtered NTv3 ORFs |
| **Tier 2 Mid-NTv3** | Filtered NTv3 ORFs that contain a validation anchor peptide |
| **Tier 3 Mid-InstaNovo** | InstaNovo filtered peptides that exactly match the validation anchor |
| **Tier 4 High** | Low-tier peptides that exactly match the validation anchor |
| **Assembly** | Merge dual-ORF runs into master files with `orf_source_tier` labels |

## Quick Start

```bash
conda env create -f /home/cle1g21/RPC/confidence_pipeline/environment.yml
conda activate confidence_pipeline_env
python /home/cle1g21/RPC/confidence_pipeline/run_confidence_pipeline.py --verbose
```

Use `--skip-proteome-fetch` when `/home/cle1g21/RPC/databases/UP000005640_canonical.fasta` already exists.

Run assembly only:

```bash
python /home/cle1g21/RPC/confidence_pipeline/assemble_validation_results.py --verbose
```

Proteome FASTA caches (`databases/`), InstaNovo prediction CSVs, and `NTv3/filtered_handoff/` are listed in the repository-root `.gitignore` and are **not** pushed to GitHub. Download or regenerate them locally before running the pipeline. Proteome exclusion of InstaNovo files uses `filter_root_instanovo_vs_proteome.py` (Slurm wrapper: `run_proteome_filter_slurm.sh`); confidence tiers use `update_confidence_tiers.py`.

## Directory Map

```text
/home/cle1g21/RPC/
├── confidence_pipeline/               # Pipeline source code
├── databases/
│   └── UP000005640_canonical.fasta    # Cached UniProt proteome
├── NTv3/filtered_handoff/             # Phase 1 filtered NTv3 ORF outputs
├── instanovo_predictions/filtered_predictions/  # InstaNovo input (pre-filtered)
├── validation/41586_2026_10459_MOESM4_ESM.csv   # Permanent validation anchor
└── confidence_levels/
    ├── low/<orf_run>/                 # Tier 1 outputs
    ├── med_ntv3/                      # Tier 2 outputs
    ├── med_instanovo/                 # Tier 3 outputs
    ├── high/<orf_run>/                # Tier 4 outputs
    ├── assembled/                     # Master merged outputs
    └── run_summary.json
```

NTv3 raw handoff inputs:

`/iridisfs/ddnb/kitsune_labs/predictions/canonical/handoff/merged_orfs.tsv` (475 ORFs, 8-340 aa)

`/iridisfs/ddnb/kitsune_labs/predictions/canonical/handoff/merged_orfs_ge30aa.tsv` (409 ORFs, >=30 aa)

## Dual ORF Strategy

Both ORF files are processed in separate runs:

| File | ORF run label | Assembly label |
|------|---------------|----------------|
| `merged_orfs_ge30aa.tsv` | `ge30aa_orfs` | `Conservative (>=30aa)` |
| `merged_orfs.tsv` | `all_orfs` | `Short Fragment (8-29aa)` |

Tier 1 and Tier 4 outputs are written under per-run subfolders (`low/ge30aa_orfs/`, `low/all_orfs/`, etc.) to prevent filename collisions.

## Matching Modes

| Tier | Mode | Direction | Description |
|------|------|-----------|-------------|
| Low | `substring` | `query_in_reference` | InstaNovo peptide contained in any NTv3 ORF |
| Mid-NTv3 | `substring` | `reference_in_query` | Validation peptide contained in NTv3 ORF |
| Mid-InstaNovo | `exact` | — | InstaNovo peptide exactly in validation set |
| High | `exact` | — | Low-tier peptide exactly in validation set |

Proteome exclusion uses substring matching with optional L/I isobaric equivalence.

## Delimiter Note

NTv3 handoff files use the `.tsv` extension but are **comma-delimited**. The pipeline sniffs delimiters automatically.

## Function Catalog

| Function | Module | Purpose | Inputs | Outputs |
|----------|--------|---------|--------|---------|
| `default_config` | `config/config.py` | Return default `PipelineConfig` | None | `PipelineConfig` |
| `sanitize_raw_sequence` | `string_normalizer.py` | Strip quotes/whitespace/CR/LF | `sequence: str` | Clean `str` |
| `strip_modifications` | `string_normalizer.py` | Remove `[UNIMOD:n]` tags | `sequence: str` | Mod-free `str` |
| `unify_il` | `string_normalizer.py` | Map I to L | `sequence: str` | I/L-normalized `str` |
| `normalize_sequence` | `string_normalizer.py` | Full normalization pipeline | `sequence`, flags | Normalized `str` |
| `is_valid_peptide_sequence` | `string_normalizer.py` | Validate AA alphabet | `sequence: str` | `bool` |
| `proteome_cache_is_valid` | `fetch_proteome.py` | Check FASTA cache usability | `path`, `min_bytes` | `bool` |
| `fetch_uniprot_proteome` | `fetch_proteome.py` | Stream UniProt FASTA download | `url`, `dest_path` | `Path` |
| `load_proteome_sequences` | `fetch_proteome.py` | Parse FASTA to sequence list | `fasta_path` | `list[str]` |
| `ensure_dir` | `file_io.py` | Create directory tree | `path` | `Path` |
| `scan_input_files` | `file_io.py` | Glob-discover input files | `dir`, `patterns` | `list[Path]` |
| `sniff_delimiter` | `file_io.py` | Detect comma vs tab | `path` | `str` |
| `read_table` | `file_io.py` | Delimiter-aware table read | `path` | `DataFrame` |
| `write_table` | `file_io.py` | Write CSV output | `df`, `path` | `Path` |
| `derive_output_name` | `file_io.py` | Build suffixed filename | `input_path`, `suffix` | `str` |
| `orf_run_name` | `file_io.py` | Map ORF filename to run label | `orf_filename` | `str` |
| `build_protein_substring_index` | `file_matcher.py` | Build proteome index | `proteins`, `il` | `ProteinSubstringIndex` |
| `is_known_fragment` | `file_matcher.py` | Test peptide in proteome | `peptide`, `index` | `bool` |
| `build_sequence_set` | `file_matcher.py` | Normalized sequence set | `df`, `col` | `set[str]` |
| `filter_rows_by_membership` | `file_matcher.py` | Exact `.isin()` filter | `df`, `col`, `ref_set` | `DataFrame` |
| `filter_rows_by_containment` | `file_matcher.py` | Substring containment filter | `df`, `col`, `refs`, `direction` | `DataFrame` |
| `filter_table_against_proteome` | `file_matcher.py` | Proteome exclusion filter | `df`, `col`, `index` | `(df, stats)` |
| `match_table` | `file_matcher.py` | Unified exact/substring matcher | `df`, `col`, `refs`, `mode` | `DataFrame` |
| `run_pipeline` | `run_confidence_pipeline.py` | Full orchestration | `config` | summary `dict` |
| `filter_file` | `filter_root_instanovo_vs_proteome.py` | Chunked proteome exclusion of root InstaNovo CSVs | prediction CSV | `filtered_predictions/*_filtered.csv` |
| `assemble_tier_low_or_high` | `assemble_validation_results.py` | Merge Low/High ORF runs | tier dir, col | stats `dict` |
| `assemble_tier_med_ntv3` | `assemble_validation_results.py` | Merge Mid-NTv3 ORF runs | tier dir, col | stats `dict` |

## Walkthrough: Proteome Exclusion

1. Download or reuse `UP000005640_canonical.fasta` from UniProt.
2. Build an in-memory index of ~147k canonical protein sequences.
3. For each NTv3 ORF in `orf_aa_seq`, test whether the full ORF is a substring of any canonical protein.
4. Drop matching ORFs; write survivors to `NTv3/filtered_handoff/<name>_filtered.csv`.

Optional L/I equivalence maps both isoleucine and leucine to `L` before comparison.

## Walkthrough: Cross-Tier Filtering

**Tier 1 (Low):** For each InstaNovo filtered file and each ORF run, keep prediction rows whose peptide is found inside any filtered NTv3 ORF sequence.

**Tier 2 (Mid-NTv3):** For each ORF run, keep NTv3 ORF rows that contain at least one validation anchor peptide as a substring.

**Tier 3 (Mid-InstaNovo):** For each InstaNovo filtered file, keep rows whose normalized peptide exactly matches a validation anchor sequence. ORF-independent.

**Tier 4 (High):** For each Low-tier output file and ORF run, keep rows whose peptide exactly matches a validation anchor sequence. This is the strongest tier: InstaNovo + NTv3 + validation anchor agreement.

## Assembly Logic

1. Load `ge30aa_orfs` outputs first; label `Conservative (>=30aa)`.
2. Load `all_orfs` outputs; label novel rows `Short Fragment (8-29aa)`.
3. Deduplicate on normalized sequence key.
4. Write `master_low.csv`, `master_med_ntv3.csv`, `master_high.csv`.

## Configuration

Key settings in `config/config.py`:

| Setting | Default |
|---------|---------|
| `orf_files` | `[merged_orfs_ge30aa.tsv, merged_orfs.tsv]` |
| `treat_leucine_isoleucine_as_identical` | `True` (proteome exclusion only) |
| `strip_modifications` | `True` |
| `instanovo_sequence_column` | `predictions` |
| `validation_sequence_column` | `sequence` |
| `ntv3_sequence_column` | `orf_aa_seq` |

## Dependencies

| Package | Purpose |
|---------|---------|
| `pandas` | Table I/O and filtering |
| `requests` | UniProt FASTA streaming |
| `biopython` | FASTA parsing |
| `openpyxl` | Excel support (future extensions) |
