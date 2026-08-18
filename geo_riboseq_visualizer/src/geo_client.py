"""NCBI GEO (Entrez) client for Ribo-seq dataset counts by species.

Uses Biopython ``Bio.Entrez`` against the GEO DataSets database (``db='gds'``),
with multi-tiered search terms and deep metadata validation to recover Ribo-seq
studies misclassified under generic library strategies (e.g. ``Other``).
"""

from __future__ import annotations

import logging
import re
import time
from collections import Counter
from dataclasses import dataclass
from typing import Any, Iterator

import pandas as pd
from Bio import Entrez

logger = logging.getLogger(__name__)

DEFAULT_EMAIL = "cle1g21@soton.ac.uk"
DEFAULT_TOOL = "geo-riboseq-visualizer"

# Tier 1: explicit Ribo-seq synonyms across title/abstract/all fields (user spec).
RIBOSEQ_KEYWORD_QUERY = (
    '("ribosome profiling"[All Fields] OR "Ribo-seq"[All Fields] OR '
    '"Ribo-seq"[Library Strategy] OR "Ribo seq"[All Fields] OR '
    '"ribosomal footprinting"[All Fields] OR "ribosome footprints"[All Fields] OR '
    '"ARTseq"[All Fields] OR "translation profiling"[All Fields])'
)

# Tier 2: technical indicators combined with ambiguous library strategies.
RIBOSEQ_INDICATOR_QUERY = (
    '(cycloheximide[All Fields] OR harringtonine[All Fields] OR '
    'lactimidomycin[All Fields] OR "RNase I"[All Fields] OR '
    '"micrococcal nuclease"[All Fields] OR MNase[All Fields] OR '
    '"ribosome-protected fragments"[All Fields] OR RPFs[All Fields] OR '
    'RPF[All Fields])'
)

AMBIGUOUS_STRATEGY_QUERY = (
    '("Other"[Library Strategy] OR "RNA-Seq"[Library Strategy] OR '
    '"RNA-seq"[Library Strategy])'
)

# Text fields scanned during deep validation.
METADATA_TEXT_KEYS = (
    "title",
    "Title",
    "summary",
    "Summary",
    "description",
    "Description",
    "overall_design",
    "Overall Design",
    "overall design",
    "gds_title",
    "GDS_Title",
)

LIBRARY_STRATEGY_KEYS = (
    "library_strategy",
    "Library Strategy",
    "librarystrategy",
    "gdsType",
    "GTYP",
    "type",
    "entrytype",
)

# Case-insensitive indicator patterns for deep metadata mining.
TECHNICAL_INDICATOR_PATTERNS: tuple[re.Pattern[str], ...] = tuple(
    re.compile(pat, re.IGNORECASE)
    for pat in (
        r"\bcycloheximide\b",
        r"\bharringtonine\b",
        r"\blactimidomycin\b",
        r"\bRNase\s*I\b",
        r"\bmicrococcal\s+nuclease\b",
        r"\bMNase\b",
        r"\bribosome[- ]protected fragments?\b",
        r"\bRPFs?\b",
    )
)

EXPLICIT_RIBOSEQ_PATTERNS: tuple[re.Pattern[str], ...] = tuple(
    re.compile(pat, re.IGNORECASE)
    for pat in (
        r"\bribosome profiling\b",
        r"\bRibo[- ]?seq\b",
        r"\bribosomal footprinting\b",
        r"\bribosome footprints?\b",
        r"\bARTseq\b",
        r"\btranslation profiling\b",
    )
)

AMBIGUOUS_LIBRARY_STRATEGIES = frozenset(
    {"other", "rna-seq", "rna seq", "rna sequencing", "unknown", ""}
)


class GeoEntrezError(RuntimeError):
    """Raised when Entrez querying fails unrecoverably."""


@dataclass(frozen=True)
class EntrezSearchHandle:
    """History server handle returned from Entrez.esearch(usehistory='y')."""

    webenv: str
    query_key: str
    total: int


@dataclass
class RiboseqValidationStats:
    """Counters describing how records were classified during aggregation."""

    total_summaries: int = 0
    included: int = 0
    excluded: int = 0
    rescued_by_indicators: int = 0
    explicit_strategy_or_terms: int = 0


def create_entrez_config(email: str = DEFAULT_EMAIL, api_key: str | None = None) -> None:
    """Configure Biopython Entrez identity and optional API key.

    Args:
        email: Email address to set as ``Entrez.email`` (required by NCBI).
        api_key: Optional NCBI API key. If set, higher rate limits may apply.
    """
    Entrez.email = email
    Entrez.tool = DEFAULT_TOOL
    if api_key:
        Entrez.api_key = api_key


def build_riboseq_query() -> str:
    """Build the multi-tiered Entrez query for GEO Series (GSE) Ribo-seq discovery.

    Tier 1 searches explicit Ribo-seq synonyms across all indexed fields.
    Tier 2 searches technical Ribo-seq indicators within studies whose library
    strategy is ambiguous (``Other`` / ``RNA-Seq``), recovering misclassified sets.

    Returns:
        Entrez query term suitable for ``db='gds'``.
    """
    tier2 = f"({RIBOSEQ_INDICATOR_QUERY} AND {AMBIGUOUS_STRATEGY_QUERY})"
    return f"(({RIBOSEQ_KEYWORD_QUERY}) OR {tier2}) AND gse[Entry Type]"


def _coerce_summary_dict(summary: Any) -> dict[str, Any]:
    """Normalize a Biopython Entrez summary object to a plain dict."""
    if isinstance(summary, dict):
        return summary
    try:
        return dict(summary)
    except Exception:  # noqa: BLE001
        return {}


def _field_value(summary: dict[str, Any], keys: tuple[str, ...]) -> str:
    """Return the first non-empty string value among candidate keys."""
    for key in keys:
        value = summary.get(key)
        if value is None:
            continue
        if isinstance(value, list):
            value = value[0] if value else None
        text = str(value).strip()
        if text:
            return text
    return ""


def extract_metadata_text(summary: dict[str, Any]) -> str:
    """Concatenate title, summary, and overall design text for deep scanning.

    Args:
        summary: GEO DataSets document summary.

    Returns:
        Lowercased, space-joined metadata text.
    """
    parts: list[str] = []
    for key in METADATA_TEXT_KEYS:
        value = summary.get(key)
        if value is None:
            continue
        if isinstance(value, list):
            parts.extend(str(v).strip() for v in value if str(v).strip())
        else:
            text = str(value).strip()
            if text:
                parts.append(text)
    return " ".join(parts)


def extract_library_strategy(summary: dict[str, Any]) -> str:
    """Extract the library strategy / dataset type label from a GEO summary.

    Args:
        summary: GEO DataSets document summary.

    Returns:
        Library strategy string, or empty string if unavailable.
    """
    return _field_value(summary, LIBRARY_STRATEGY_KEYS)


def has_technical_riboseq_indicators(text: str) -> bool:
    """Return True if text contains Ribo-seq technical indicator phrases.

    Indicators include translation inhibitors, footprinting nucleases, and RPF terms.

    Args:
        text: Metadata text to scan (case-insensitive).

    Returns:
        True when at least one indicator pattern matches.
    """
    if not text:
        return False
    return any(pattern.search(text) for pattern in TECHNICAL_INDICATOR_PATTERNS)


def has_explicit_riboseq_terms(text: str) -> bool:
    """Return True if text contains explicit Ribo-seq nomenclature.

    Args:
        text: Metadata text to scan (case-insensitive).

    Returns:
        True when at least one explicit Ribo-seq term pattern matches.
    """
    if not text:
        return False
    return any(pattern.search(text) for pattern in EXPLICIT_RIBOSEQ_PATTERNS)


def strategy_indicates_riboseq(strategy: str) -> bool:
    """Return True when the library strategy field explicitly suggests Ribo-seq.

    Args:
        strategy: Library strategy label from GEO metadata.

    Returns:
        True for non-ambiguous Ribo-seq strategy labels.
    """
    normalized = strategy.strip().lower()
    if not normalized or normalized in AMBIGUOUS_LIBRARY_STRATEGIES:
        return False
    ribo_markers = (
        "ribo",
        "ribosome",
        "footprint",
        "profiling",
        "artseq",
        "translation profiling",
        "rpf",
    )
    return any(marker in normalized for marker in ribo_markers)


def is_riboseq_dataset(summary: dict[str, Any]) -> tuple[bool, str]:
    """Validate whether a GEO Series record should count as Ribo-seq.

    Records with ambiguous library strategies (``Other``, ``RNA-Seq``) are not
    discarded outright. They are retained when deep metadata scanning of
    ``Title``, ``Summary``, and ``Overall Design`` finds explicit Ribo-seq terms
    or technical Ribo-seq indicators (cycloheximide, harringtonine, RNase I, RPFs).

    Args:
        summary: GEO DataSets document summary dict.

    Returns:
        Tuple of (is_riboseq, reason) where reason is a short classification tag.
    """
    text = extract_metadata_text(summary)
    strategy = extract_library_strategy(summary)
    normalized_strategy = strategy.strip().lower()

    if strategy_indicates_riboseq(strategy):
        return True, "explicit_library_strategy"

    if has_explicit_riboseq_terms(text):
        return True, "explicit_metadata_terms"

    if has_technical_riboseq_indicators(text):
        if normalized_strategy in AMBIGUOUS_LIBRARY_STRATEGIES:
            return True, "rescued_by_technical_indicators"
        return True, "technical_indicators"

    if normalized_strategy in AMBIGUOUS_LIBRARY_STRATEGIES:
        return False, "ambiguous_strategy_no_indicators"

    return False, "no_riboseq_evidence"


def esearch_series_uids(term: str) -> EntrezSearchHandle:
    """Search GEO DataSets for matching Series and return History server handle.

    Args:
        term: Entrez search term.

    Returns:
        History server handle including total record count.

    Raises:
        GeoEntrezError: If the request fails or returns malformed results.
    """
    try:
        with Entrez.esearch(db="gds", term=term, usehistory="y", retmax=0) as handle:
            record = Entrez.read(handle)
    except Exception as exc:  # noqa: BLE001
        raise GeoEntrezError(f"Entrez.esearch failed for term: {term!r}") from exc

    try:
        webenv = str(record["WebEnv"])
        query_key = str(record["QueryKey"])
        total = int(record["Count"])
    except Exception as exc:  # noqa: BLE001
        raise GeoEntrezError(f"Unexpected esearch record structure: {record!r}") from exc

    return EntrezSearchHandle(webenv=webenv, query_key=query_key, total=total)


def iter_series_summaries(
    search: EntrezSearchHandle,
    *,
    batch_size: int = 200,
    sleep_seconds: float = 0.34,
) -> Iterator[dict[str, Any]]:
    """Iterate GEO Series summaries using the Entrez History server.

    Args:
        search: History server handle from ``esearch_series_uids``.
        batch_size: Number of records to fetch per ``esummary`` call.
        sleep_seconds: Delay between requests for NCBI rate limiting.

    Yields:
        Parsed summary dicts (one per GEO Series record).

    Raises:
        GeoEntrezError: If an esummary call fails.
    """
    total = search.total
    for retstart in range(0, total, batch_size):
        try:
            with Entrez.esummary(
                db="gds",
                query_key=search.query_key,
                WebEnv=search.webenv,
                retstart=retstart,
                retmax=batch_size,
            ) as handle:
                record = Entrez.read(handle)
        except Exception as exc:  # noqa: BLE001
            raise GeoEntrezError(
                f"Entrez.esummary failed at retstart={retstart}, retmax={batch_size}"
            ) from exc

        if hasattr(record, "get"):
            docs: Any = record.get("DocumentSummarySet", {}).get("DocumentSummary", [])
        else:
            docs = record

        if not isinstance(docs, list):
            docs = [docs]

        for doc in docs:
            yield _coerce_summary_dict(doc)

        time.sleep(sleep_seconds)


def extract_species(summary: dict[str, Any]) -> str | None:
    """Extract organism/species string from a GEO DataSets summary.

    Args:
        summary: One document summary dict from ``iter_series_summaries``.

    Returns:
        Species name if found, else None.
    """
    candidates = (
        "taxon",
        "organism",
        "Organism",
        "taxname",
        "Taxon",
        "TaxName",
    )
    for key in candidates:
        value = summary.get(key)
        if value is None:
            continue
        if isinstance(value, list):
            value = value[0] if value else None
        if value is None:
            continue
        text = str(value).strip()
        if text:
            return text
    return None


def fetch_species_counts(
    term: str | None = None,
    *,
    batch_size: int = 200,
    sleep_seconds: float = 0.34,
    unknown_label: str = "Unknown",
) -> pd.DataFrame:
    """Fetch GEO Series, validate Ribo-seq metadata, and aggregate counts by species.

    Args:
        term: Entrez query term. If None, uses ``build_riboseq_query()``.
        batch_size: esummary records per request.
        sleep_seconds: Delay between requests.
        unknown_label: Label used when organism cannot be extracted.

    Returns:
        DataFrame with columns ``species`` and ``dataset_count`` (descending).
    """
    query = term or build_riboseq_query()
    search = esearch_series_uids(query)
    logger.info("Entrez term: %s", query)
    logger.info("GEO Series hits (pre-validation): %d", search.total)

    counts: Counter[str] = Counter()
    stats = RiboseqValidationStats()

    for i, raw_summary in enumerate(
        iter_series_summaries(search, batch_size=batch_size, sleep_seconds=sleep_seconds),
        start=1,
    ):
        stats.total_summaries += 1
        is_ribo, reason = is_riboseq_dataset(raw_summary)
        if not is_ribo:
            stats.excluded += 1
            logger.debug("Excluded (%s): %s", reason, _field_value(raw_summary, ("title", "Title")))
            continue

        stats.included += 1
        if reason == "rescued_by_technical_indicators":
            stats.rescued_by_indicators += 1
        if reason in {"explicit_library_strategy", "explicit_metadata_terms"}:
            stats.explicit_strategy_or_terms += 1

        species = extract_species(raw_summary) or unknown_label
        counts[species] += 1

        if i % 500 == 0:
            logger.info(
                "Processed %d/%d | included=%d excluded=%d rescued=%d",
                i,
                search.total,
                stats.included,
                stats.excluded,
                stats.rescued_by_indicators,
            )

    logger.info(
        "Validation complete: total=%d included=%d excluded=%d "
        "rescued_by_indicators=%d explicit=%d species=%d",
        stats.total_summaries,
        stats.included,
        stats.excluded,
        stats.rescued_by_indicators,
        stats.explicit_strategy_or_terms,
        len(counts),
    )

    df = pd.DataFrame(
        [{"species": species, "dataset_count": int(count)} for species, count in counts.items()]
    ).sort_values("dataset_count", ascending=False)
    return df.reset_index(drop=True)
