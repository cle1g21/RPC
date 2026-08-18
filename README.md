# RPC

Code for immunopeptidomics *de novo* sequencing (InstaNovo), NTv3-ribo ORF filtering, confidence tiers, and supporting visualizers.

## Packages

| Directory | Purpose |
|---|---|
| `msconvert/` | ProteoWizard MSConvert (Apptainer) → MGF ([separate repo](https://github.com/cle1g21/MSConvert)) |
| `instanovo_predictions/` | ROC scripts; local prediction CSVs are gitignored |
| `confidence_pipeline/` | Proteome exclusion + four confidence tiers |
| `PRIDEpull/` | PRIDE harvest → convert → InstaNovo |
| `peptiverse/` | PeptiVerse property plots |
| `uniprot_gene_length_visualizer/` | UniProt coding-length KDEs + Low-tier summary plots |
| `deutsch_ntv3_pipeline/` | Deutsch/Kok NTv3 coordinate validation |
| `geo_riboseq_visualizer/` | GEO Ribo-seq species counts |
| `pride_species_visualizer/` | PRIDE Archive species counts |

Each package has its own `README.md` and (where used) `environment.yml`.

## Git: large files

The root `.gitignore` excludes multi-GB artefacts so GitHub stays code-only:

- InstaNovo CSVs, `chunks/`, `filtered_predictions/`
- UniProt FASTA caches (`databases/`, `uniprot_gene_length_visualizer/data/`)
- Spectra (`.mgf`, `.mzML`, `.raw`, Bruker `.d`) and `*.sif`
- PRIDEpull `data/` and `cache/`
- `__pycache__/`, `**/logs/`

Regenerate or download those files locally using the per-package READMEs.
