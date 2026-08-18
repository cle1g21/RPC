"""Centralized configuration for the PRIDEpull immunopeptidomics pipeline."""

# Import Path so all filesystem locations are defined as pathlib objects
from pathlib import Path

# Resolve the PRIDEpull repository root (parent of the config package)
REPO_ROOT = Path(__file__).resolve().parent.parent

# Load optional environment variables from a .env file if present
try:
    # Import dotenv only when the package is installed in the conda env
    from dotenv import load_dotenv

    # Read .env from the repo root without overriding existing shell variables
    load_dotenv(REPO_ROOT / ".env", override=False)
except ImportError:
    # Continue without dotenv when the dependency is not installed yet
    pass

# ---------------------------------------------------------------------------
# PRIDE Archive REST API
# ---------------------------------------------------------------------------

# Base URL for PRIDE Archive API v3 (no trailing slash)
PRIDE_BASE_URL = "https://www.ebi.ac.uk/pride/ws/archive/v3"

# Number of projects requested per paginated API page
PRIDE_PAGE_SIZE = 100

# HTTP connect and read timeout tuple in seconds (connect, read)
PRIDE_TIMEOUT = (10, 120)

# Maximum retry attempts for transient API failures (429, 5xx)
PRIDE_MAX_RETRIES = 5

# Base wait seconds for exponential backoff between retries
PRIDE_RETRY_BACKOFF_SECONDS = 2.0

# User-Agent string sent with every PRIDE HTTP request
PRIDE_USER_AGENT = "pride-pull/1.0.0"

# ---------------------------------------------------------------------------
# Species / taxonomy
# ---------------------------------------------------------------------------

# NCBI Taxonomy ID for Homo sapiens
HUMAN_TAXONOMY_ID = 9606

# PRIDE organisms_facet value matching taxonomy 9606 (case-sensitive)
HUMAN_ORGANISM_FACET = "Homo sapiens (human)"

# Substring used for client-side human organism validation
HUMAN_ORGANISM_LABEL = "Homo sapiens"

# ---------------------------------------------------------------------------
# Immunopeptidomics keyword matrix (client-side text matching)
# ---------------------------------------------------------------------------

IMMUNO_KEYWORDS: list[str] = [
    "immunopeptidomics",
    "immunopeptidome",
    "peptidome",
    "peptidomics",
    "HLA",
    "HLA-I",
    "HLA-II",
    "human leukocyte antigen",
    "MHC",
    "MHC-I",
    "MHC-II",
    "major histocompatibility complex",
    "neoantigen",
    "neo-antigen",
    "tumor antigen",
    "cancer testis antigen",
    "antigen presentation",
    "antigen processing",
    "eluted peptides",
    "proteogenomics",
    "altORF",
    "non-canonical peptides",
]

# Reduced keyword list used for server-side PRIDE search queries (fewer API calls)
PRIDE_SERVER_KEYWORDS: list[str] = [
    "immunopeptidomics",
    "immunopeptidome",
    "HLA",
    "MHC",
    "neoantigen",
    "eluted peptides",
    "proteogenomics",
]

# ---------------------------------------------------------------------------
# Storage paths
# ---------------------------------------------------------------------------

# Isolated landing directory for all active downloads and conversions
DATA_LANDING_DIR = REPO_ROOT / "data"

# JSON cache written by the harvester to avoid re-querying PRIDE
HARVEST_CACHE_PATH = REPO_ROOT / "cache" / "harvest_cache.json"

# Append-only audit ledger for processing history and statistics
PROCESSING_LEDGER_PATH = REPO_ROOT / "processing_history_ledger.txt"

# Final InstaNovo prediction CSV output directory
PREDICTIONS_OUTPUT_DIR = Path(
    "/home/cle1g21/RPC/instanovo_predictions/predictions"
)

# ---------------------------------------------------------------------------
# msconvert (Apptainer / ProteoWizard)
# ---------------------------------------------------------------------------

# Python wrapper script that invokes msconvert inside Apptainer
MSCONVERT_SCRIPT = Path("/home/cle1g21/RPC/msconvert/convert_ms.py")

# Local Apptainer SIF image containing ProteoWizard msconvert
MSCONVERT_IMAGE = Path("/home/cle1g21/RPC/msconvert/proteowizard.sif")

# Cluster module that must be loaded before apptainer is available
APPTAINER_MODULE = "apptainer/1.5.0"

# ---------------------------------------------------------------------------
# InstaNovo
# ---------------------------------------------------------------------------

# InstaNovo installation root on the Iridis cluster
INSTANOVO_ROOT = Path("/iridisfs/ddnb/Charlotte/RPC/InstaNovo")

# InstaNovo CLI binary inside the micromamba conda environment
INSTANOVO_BIN = Path(
    "/iridisfs/scratch/cle1g21/micromamba/envs/instanovo/bin/instanovo"
)

# Slurm launcher: submits GPU jobs and blocks until predictions are complete
INSTANOVO_SLURM_SCRIPT = (
    INSTANOVO_ROOT / "scripts" / "slurm" / "run_until_complete.sh"
)

# Backwards-compatible alias used by older references
INSTANOVO_RESILIENT_SCRIPT = INSTANOVO_SLURM_SCRIPT

# ---------------------------------------------------------------------------
# Pipeline behaviour flags
# ---------------------------------------------------------------------------

# When True, only the first SUBSET_SPECTRUM_COUNT spectra are predicted (sanity check)
RUN_SUBSET_ONLY = True

# Number of MGF spectrum blocks to slice when RUN_SUBSET_ONLY is True
SUBSET_SPECTRUM_COUNT = 500

# Maximum files to process per PRIDE project (largest native MGF preferred)
MAX_FILES_PER_PROJECT = 1

# Maximum projects to process in one pipeline run (None = no limit)
MAX_PROJECTS: int | None = None

# Process msconvert conversions one at a time to avoid disk I/O thrashing
SEQUENTIAL_CONVERSION = True

# Delete intermediate files in data/ after successful prediction write
DELETE_ON_SUCCESS = True

# ---------------------------------------------------------------------------
# InstaNovo execution routing
# ---------------------------------------------------------------------------

# When True, all InstaNovo runs go through run_until_complete.sh (GPU Slurm jobs).
# When False, instanovo predict runs as a local subprocess (login-node CPU only).
INSTANOVO_USE_SLURM = True

# Batch size passed to instanovo predict for local direct runs (INSTANOVO_USE_SLURM=False)
INSTANOVO_BATCH_SIZE = 32

# Worker count passed to instanovo predict for local direct runs
INSTANOVO_NUM_WORKERS = 4

# MGF files at or above this size (GB) are chunked by run_until_complete.sh
INSTANOVO_LARGE_FILE_GB = 1.0

# ---------------------------------------------------------------------------
# Download settings
# ---------------------------------------------------------------------------

# Bytes read per chunk when streaming a file from PRIDE FTP/HTTPS
DOWNLOAD_CHUNK_BYTES = 8 * 1024 * 1024

# Preferred publicFileLocations protocol names (first match wins)
DOWNLOAD_PROTOCOL_PREFERENCE: list[str] = [
    "FTP Protocol",
    "HTTPS Protocol",
    "Aspera Protocol",
]

# File extensions treated as native MGF peak lists
MGF_EXTENSIONS: tuple[str, ...] = (".mgf",)

# File extensions treated as Thermo vendor RAW binaries
RAW_EXTENSIONS: tuple[str, ...] = (".raw",)

# File extensions treated as Bruker .d folder archives or directories
BRUKER_D_EXTENSIONS: tuple[str, ...] = (".d", ".d.zip", ".tar", ".zip")
