#!/usr/bin/env python3
"""Proteome-substring filter for root-level InstaNovo prediction CSVs.

Follows confidence_pipeline README Phase 1 logic:
1. Reuse UniProt UP000005640 canonical FASTA
2. Build in-memory protein index (6-mer accelerated, I/L equivalent)
3. For each prediction peptide, test substring membership in any protein
4. Drop matches; write survivors to filtered_predictions/*_filtered.csv
"""

from __future__ import annotations

import json
import logging
import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/home/cle1g21/RPC/confidence_pipeline")
from src.fetch_proteome import load_proteome_sequences
from src.string_normalizer import is_valid_peptide_sequence, normalize_sequence

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger("proteome_filter_instanovo")

BASE = Path("/home/cle1g21/RPC")
PROTEOME_FASTA = BASE / "databases" / "UP000005640_canonical.fasta"
OUT_DIR = BASE / "instanovo_predictions" / "filtered_predictions"
SUMMARY_PATH = BASE / "confidence_levels" / "root_instanovo_proteome_filter_summary.json"

INPUTS = [
    BASE / "instanovo_predictions" / "01625b_GB2-TUM_first_pool_10_01_01-3xHCD-1h-R1_predictions.csv",
    BASE / "instanovo_predictions" / "Control4_Neo_SN_114_HLA-I_predictions.csv",
]

SEQ_COL = "predictions"
KMER = 6
CHUNKSIZE = 50_000


def build_kmer_index(proteins: list[str], k: int = KMER) -> dict[str, list[int]]:
    index: dict[str, list[int]] = defaultdict(list)
    for protein_idx, protein in enumerate(proteins):
        if len(protein) < k:
            continue
        seen: set[str] = set()
        for start in range(len(protein) - k + 1):
            kmer = protein[start : start + k]
            if kmer in seen:
                continue
            seen.add(kmer)
            index[kmer].append(protein_idx)
    return index


def peptide_in_proteome(
    peptide: str,
    proteins: list[str],
    kmer_index: dict[str, list[int]],
    k: int = KMER,
) -> bool:
    if not is_valid_peptide_sequence(peptide):
        return False
    if len(peptide) < k:
        return any(peptide in protein for protein in proteins)
    candidates = kmer_index.get(peptide[:k])
    if not candidates:
        return False
    return any(peptide in proteins[idx] for idx in candidates)


def filter_file(
    input_csv: Path,
    output_csv: Path,
    proteins: list[str],
    kmer_index: dict[str, list[int]],
) -> dict[str, object]:
    log.info("Filtering %s", input_csv)
    t0 = time.perf_counter()
    unique_cache: dict[str, bool] = {}
    input_rows = 0
    removed_rows = 0
    retained_rows = 0
    first_chunk = True

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    if output_csv.exists():
        output_csv.unlink()

    for chunk_idx, chunk in enumerate(
        pd.read_csv(input_csv, chunksize=CHUNKSIZE, low_memory=False),
        start=1,
    ):
        if SEQ_COL not in chunk.columns:
            raise KeyError(f"Missing column '{SEQ_COL}' in {input_csv}")

        input_rows += len(chunk)
        norms = []
        for value in chunk[SEQ_COL].astype(str):
            norm = normalize_sequence(value, strip_mods=True, il_equivalent=True)
            norms.append(norm)
            if norm not in unique_cache:
                unique_cache[norm] = peptide_in_proteome(norm, proteins, kmer_index)

        known_mask = [unique_cache[norm] for norm in norms]
        kept = chunk.loc[[not known for known in known_mask]].copy()
        removed_rows += int(sum(known_mask))
        retained_rows += len(kept)

        kept.to_csv(
            output_csv,
            mode="w" if first_chunk else "a",
            header=first_chunk,
            index=False,
        )
        first_chunk = False

        if chunk_idx % 10 == 0 or chunk_idx == 1:
            log.info(
                "  %s chunk %d: cumulative input=%d removed=%d retained=%d unique=%d",
                input_csv.name,
                chunk_idx,
                input_rows,
                removed_rows,
                retained_rows,
                len(unique_cache),
            )

    stats = {
        "input_path": str(input_csv),
        "output_path": str(output_csv),
        "input_rows": input_rows,
        "unique_peptides": len(unique_cache),
        "removed_rows": removed_rows,
        "retained_rows": retained_rows,
        "elapsed_seconds": round(time.perf_counter() - t0, 1),
    }
    log.info(
        "Done %s: input=%d removed=%d retained=%d (%.1fs)",
        input_csv.name,
        input_rows,
        removed_rows,
        retained_rows,
        stats["elapsed_seconds"],
    )
    return stats


def main() -> int:
    t0 = time.perf_counter()
    if not PROTEOME_FASTA.is_file():
        raise FileNotFoundError(f"Proteome FASTA missing: {PROTEOME_FASTA}")

    log.info("Loading proteome from %s", PROTEOME_FASTA)
    proteins = [seq.replace("I", "L") for seq in load_proteome_sequences(PROTEOME_FASTA)]
    log.info("Building %d-mer index over %d proteins", KMER, len(proteins))
    kmer_index = build_kmer_index(proteins, KMER)
    log.info("Index ready: %d unique %d-mers", len(kmer_index), KMER)

    summary: dict[str, object] = {
        "proteome_protein_count": len(proteins),
        "proteome_cache_path": str(PROTEOME_FASTA),
        "files": {},
    }

    for input_csv in INPUTS:
        if not input_csv.is_file():
            raise FileNotFoundError(input_csv)
        output_csv = OUT_DIR / f"{input_csv.stem}_filtered.csv"
        summary["files"][input_csv.name] = filter_file(
            input_csv, output_csv, proteins, kmer_index
        )

    summary["elapsed_seconds"] = round(time.perf_counter() - t0, 1)
    SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    log.info("Summary written to %s", SUMMARY_PATH)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
