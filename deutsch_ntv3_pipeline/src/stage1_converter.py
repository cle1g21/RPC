"""Stage 1: convert R RDS sequence data to standard FASTA on IridisX."""

from __future__ import annotations

import logging
import re
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from config import PipelineConfig

logger = logging.getLogger(__name__)

SEQUENCE_COLUMN_CANDIDATES = (
    "protein",
    "sequence",
    "seq",
    "nt_sequence",
    "aa_sequence",
)
ID_COLUMN_CANDIDATES = (
    "protein_accession",
    "accession",
    "id",
    "name",
    "orf_id",
)

_R_CONVERT_SCRIPT = """
obj <- readRDS("{rds_path}")
out <- "{fasta_path}"
write_fasta <- function(ids, seqs, path) {{
  con <- file(path, open = "w")
  on.exit(close(con))
  for (i in seq_along(seqs)) {{
    id <- ids[i]
    if (is.na(id) || id == "") id <- paste0("seq_", i)
    writeLines(paste0(">", id), con)
    writeLines(seqs[i], con)
  }}
}}
if (requireNamespace("Biostrings", quietly = TRUE) &&
    (inherits(obj, "XStringSet") || inherits(obj, "DNAStringSet") ||
     inherits(obj, "RNAStringSet") || inherits(obj, "AAStringSet"))) {{
  Biostrings::writeXStringSet(obj, filepath=out)
}} else if (is.character(obj)) {{
  ids <- if (!is.null(names(obj))) names(obj) else paste0("seq_", seq_along(obj))
  write_fasta(ids, as.character(obj), out)
}} else if (is.data.frame(obj)) {{
  seq_col <- intersect(c("protein", "sequence", "seq"), names(obj))[1]
  if (is.na(seq_col)) stop("No sequence column in data.frame")
  id_col <- intersect(c("protein_accession", "accession", "id", "name"), names(obj))[1]
  seqs <- as.character(obj[[seq_col]])
  ids <- if (!is.na(id_col)) as.character(obj[[id_col]]) else paste0("seq_", seq_along(seqs))
  ids[is.na(ids) | ids == ""] <- paste0("seq_", which(is.na(ids) | ids == ""))
  write_fasta(make.unique(ids), seqs, out)
}} else {{
  write_fasta(paste0("seq_", seq_along(obj)), as.character(obj), out)
}}
cat("OK\\n")
"""


class RdsConversionError(RuntimeError):
    """Raised when RDS probing or FASTA conversion fails."""


def probe_rds_class(rds_path: Path) -> str:
    """Probe the R class and column layout of an RDS file.

    Args:
        rds_path: Path to ``can_nonc_seq.RDS``.

    Returns:
        Human-readable description of object class and columns.
    """
    if not rds_path.is_file():
        raise FileNotFoundError(f"RDS not found: {rds_path}")

    rscript = _find_rscript()
    if rscript:
        script = (
            f'obj <- readRDS("{rds_path}"); '
            f'cat(paste(class(obj), collapse=","), "\\n"); '
            f'if (is.data.frame(obj)) cat(paste(names(obj), collapse=","), "\\n")'
        )
        result = subprocess.run(
            [rscript, "--vanilla", "-e", script],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()

    try:
        import pyreadr  # type: ignore

        result = pyreadr.read_r(str(rds_path))
        key = next(iter(result))
        df = result[key]
        return f"pyreadr:{key}:columns={','.join(df.columns)}:rows={len(df)}"
    except ImportError:
        pass

    raise RdsConversionError("Install R (module load R/4.5.1-mkl) or pyreadr.")


def convert_rds_to_fasta(rds_path: Path, fasta_path: Path) -> int:
    """Convert RDS to clean FASTA format.

    Args:
        rds_path: Input RDS path.
        fasta_path: Output FASTA path.

    Returns:
        Number of records written.
    """
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    rscript = _find_rscript()
    if rscript:
        try:
            count = _convert_with_rscript(rscript, rds_path, fasta_path)
            validate_fasta(fasta_path)
            return count
        except RdsConversionError:
            logger.warning("Rscript failed; trying pyreadr")

    count = _convert_with_pyreadr(rds_path, fasta_path)
    validate_fasta(fasta_path)
    return count


def validate_fasta(fasta_path: Path) -> None:
    """Ensure FASTA is non-empty and sequences contain no whitespace."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise RdsConversionError(f"No FASTA records: {fasta_path}")
    ws = re.compile(r"\s")
    for rec in records:
        if not str(rec.seq) or ws.search(str(rec.seq)):
            raise RdsConversionError(f"Invalid sequence in record {rec.id}")


def run_stage1(config: PipelineConfig) -> int:
    """Execute Stage 1 using paths from configuration.

    Args:
        config: Pipeline configuration.

    Returns:
        Number of FASTA records written.
    """
    info = probe_rds_class(config.rds_path)
    logger.info("RDS probe: %s", info.replace("\n", " | "))
    count = convert_rds_to_fasta(config.rds_path, config.fasta_path)
    logger.info("Wrote %d records to %s", count, config.fasta_path)
    return count


def _find_rscript() -> str | None:
    from shutil import which

    if which("Rscript"):
        return which("Rscript")
    for candidate in (
        "/iridisfs/ixsoftware/modules/applications/R/4.5.1-mkl/bin/Rscript",
        "/iridisfs/ixsoftware/modules/applications/R/4.4.3-mkl/bin/Rscript",
    ):
        if Path(candidate).is_file():
            return candidate
    return None


def _convert_with_rscript(rscript: str, rds_path: Path, fasta_path: Path) -> int:
    with tempfile.TemporaryDirectory() as tmp:
        script_path = Path(tmp) / "convert.R"
        script_path.write_text(
            _R_CONVERT_SCRIPT.format(rds_path=str(rds_path), fasta_path=str(fasta_path)),
            encoding="utf-8",
        )
        result = subprocess.run(
            [rscript, "--vanilla", str(script_path)],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0 or "OK" not in result.stdout:
            raise RdsConversionError(result.stderr or result.stdout)
    return sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))


def _convert_with_pyreadr(rds_path: Path, fasta_path: Path) -> int:
    import pyreadr  # type: ignore

    result = pyreadr.read_r(str(rds_path))
    df = result[next(iter(result))]
    seq_col = _pick_column(df.columns, SEQUENCE_COLUMN_CANDIDATES)
    id_col = _pick_column(df.columns, ID_COLUMN_CANDIDATES)
    records: list[SeqRecord] = []
    seen: dict[str, int] = {}
    for idx, row in df.iterrows():
        seq = str(row[seq_col]).strip()
        if not seq or seq.lower() == "nan":
            continue
        rid = (
            str(row[id_col]).strip()
            if id_col and str(row[id_col]).strip() not in {"", "nan"}
            else f"seq_{idx}"
        )
        rid = _unique_id(rid, seen)
        records.append(SeqRecord(Seq(seq), id=rid, description=""))
    if not records:
        raise RdsConversionError("No sequences extracted from RDS")
    SeqIO.write(records, fasta_path, "fasta")
    return len(records)


def _pick_column(columns: object, candidates: tuple[str, ...]) -> str:
    lower = {str(c).lower(): str(c) for c in columns}
    for c in candidates:
        if c in lower:
            return lower[c]
    raise RdsConversionError(f"Missing column; want one of {candidates}")


def _unique_id(rid: str, seen: dict[str, int]) -> str:
    base = re.sub(r"\s+", "_", rid)
    if base not in seen:
        seen[base] = 1
        return base
    seen[base] += 1
    return f"{base}_{seen[base]}"
