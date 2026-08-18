"""Centralized parameters for UniProt proteome download and KDE plotting."""

# Enable postponed evaluation of type annotations for cleaner forward references
from __future__ import annotations

# Import Path to build filesystem locations relative to this repository
from pathlib import Path


# Resolve the absolute path of this config file on disk
_CONFIG_FILE = Path(__file__).resolve()
# Climb one level from config/ to the repository root directory
REPO_ROOT = _CONFIG_FILE.parents[1]
# Point DATA_DIR at the local FASTA cache folder under the repository root
DATA_DIR = REPO_ROOT / "data"
# Point OUTPUT_DIR at the publication figure export folder under the repository root
OUTPUT_DIR = REPO_ROOT / "outputs"

# Define the UniProt REST streaming URL template with proteome and reviewed filters
UNIPROT_STREAM_URL_TEMPLATE = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?format=fasta"
    "&query=%28proteome%3A{proteome_id}%29%20AND%20%28reviewed%3Atrue%29"
)

# Set the HTTP request timeout in seconds for large proteome downloads
REQUEST_TIMEOUT = 600
# Set the streaming download chunk size to one mebibyte for efficient I/O
CHUNK_SIZE = 1024 * 1024
# Set the minimum acceptable cached FASTA size in bytes before forcing a re-download
MIN_CACHE_BYTES = 10_000
# Default to reusing cached FASTA files unless the caller forces a refetch
FORCE_REFETCH = False

# Catalog the five UniProt reference proteomes targeted by this pipeline
PROTEOMES = [
    # Record Homo sapiens as the first model organism entry
    {
        "species_key": "human",
        "display_name": "Homo sapiens",
        "proteome_id": "UP000005640",
        "fasta_filename": "UP000005640_human.fasta",
    },
    # Record Mus musculus as the second model organism entry
    {
        "species_key": "mouse",
        "display_name": "Mus musculus",
        "proteome_id": "UP000000589",
        "fasta_filename": "UP000000589_mouse.fasta",
    },
    # Record Arabidopsis thaliana as the third model organism entry
    {
        "species_key": "arabidopsis",
        "display_name": "Arabidopsis thaliana",
        "proteome_id": "UP000006548",
        "fasta_filename": "UP000006548_arabidopsis.fasta",
    },
    # Record Drosophila melanogaster as the fourth model organism entry
    {
        "species_key": "drosophila",
        "display_name": "Drosophila melanogaster",
        "proteome_id": "UP000000803",
        "fasta_filename": "UP000000803_drosophila.fasta",
    },
    # Record Saccharomyces cerevisiae as the fifth model organism entry
    {
        "species_key": "yeast",
        "display_name": "Saccharomyces cerevisiae",
        "proteome_id": "UP000002311",
        "fasta_filename": "UP000002311_yeast.fasta",
    },
]

# Mark the key inflection annotation position on the ORF-length axis in base pairs
ANNOTATION_BP = 1500
# Set the visible text label drawn above the annotation arrow
ANNOTATION_LABEL = "1,500 bps"
# Scale the linear-panel annotation text height relative to the curve density
ANNOTATION_Y_TEXT_FACTOR_LINEAR = 1.25
# Scale the log-panel annotation text height relative to denser log-scale peaks
ANNOTATION_Y_TEXT_FACTOR_LOG = 1.35
# Set the Seaborn KDE bandwidth adjustment used to resolve secondary modes
BW_ADJUST = 0.5
# Cap the plotted x-axis upper limit so long tails do not flatten density modes
X_LIM_MAX = 8000
# Set the lower bound of the plotted ORF-length axis in base pairs
X_LIM_MIN = 0
# Set the lower bound for logarithmic panels (must remain strictly positive)
LOG_X_LIM_MIN = 50
# Cap the logarithmic panel upper limit to match the linear long-tail window
LOG_X_LIM_MAX = 8000
# Enable generation of the dual-panel linear-versus-log comparative figure
ENABLE_LOG_PANEL = True

# Select the Seaborn theme style string applied before plotting
SNS_STYLE = "ticks"
# Select the Seaborn default palette name applied before custom species colors
SNS_PALETTE = "muted"
# Assign a distinct hex color to each species display name for consistent legends
SPECIES_PALETTE = {
    "Homo sapiens": "#4C72B0",
    "Mus musculus": "#DD8452",
    "Arabidopsis thaliana": "#55A868",
    "Drosophila melanogaster": "#C44E52",
    "Saccharomyces cerevisiae": "#8172B3",
}
# Set the single-panel KDE figure width and height in inches
KDE_FIGSIZE = (10, 6)
# Set the dual-panel comparative figure width and height in inches
DUAL_KDE_FIGSIZE = (16, 6)
# Set the left and bottom spine linewidth to match schematic publication figures
SPINE_LINEWIDTH = 2.5
# Set the color used for the remaining plot spines
SPINE_COLOR = "black"
# Set raster export resolution for camera-ready PNG figures
DPI = 300
# Name the continuous single-panel KDE density chart output stem
KDE_STEM = "gene_length_kde_density"
# Name the dual-panel linear-versus-log comparative chart output stem
DUAL_KDE_STEM = "gene_length_kde_linear_vs_log"
# Name the optional per-species length summary CSV written under outputs/
SUMMARY_CSV_NAME = "length_summary.csv"
# Set the shared plot title used on the single-panel KDE abundance chart
KDE_TITLE = "ORF length abundance across UniProt model organism proteomes"
# Set Panel A title for the linear raw-abundance KDE
PANEL_A_TITLE = "A) Linear Length Scale (0–8,000 bps)"
# Set Panel B title for the logarithmic log-normal resolution KDE
PANEL_B_TITLE = "B) Logarithmic Scale (Log10 bps)"
# Set the shared x-axis label describing continuous coding ORF length
X_AXIS_LABEL = "CDS length (nts)"
# Set the y-axis label for kernel density abundance on both panels
Y_AXIS_LABEL = "Density"
