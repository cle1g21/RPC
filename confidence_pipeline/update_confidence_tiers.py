#!/usr/bin/env python3
"""Update confidence tiers + assembled masters from current filtered InstaNovo files.

Uses chunked I/O so Control4 (~4.8M rows / 11 GB) does not need a full in-memory load.
"""

from __future__ import annotations

import json
import logging
import sys
import time
from pathlib import Path

import pandas as pd

_ROOT = Path("/home/cle1g21/RPC/confidence_pipeline")
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from config.config import PipelineConfig, default_config
from src.file_io import (
    derive_output_name,
    ensure_dir,
    orf_run_name,
    read_table,
    scan_input_files,
    write_table,
)
from src.file_matcher import build_sequence_set, match_table
from src.string_normalizer import normalize_sequence

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger("update_tiers")

CHUNKSIZE = 50_000
SEQ_COL = "predictions"


def ensure_dirs(config: PipelineConfig) -> None:
    ensure_dir(config.confidence_low_dir)
    ensure_dir(config.confidence_med_ntv3_dir)
    ensure_dir(config.confidence_med_instanovo_dir)
    ensure_dir(config.confidence_high_dir)
    ensure_dir(config.confidence_assembled_dir)
    for orf_file in config.orf_files:
        run_label = orf_run_name(orf_file)
        ensure_dir(config.confidence_low_dir / run_label)
        ensure_dir(config.confidence_high_dir / run_label)


def build_orf_peptide_lookup(orf_sequences: list[str], *, strip_mods: bool) -> set[str]:
    """Build a set of all ORF substrings in the typical peptide length range."""
    orfs = [
        normalize_sequence(seq, strip_mods=strip_mods, il_equivalent=False)
        for seq in orf_sequences
    ]
    orfs = [seq for seq in orfs if seq]
    # Cover observed InstaNovo peptide lengths (short fragments through long ORFs).
    min_len = 3
    max_len = max((len(orf) for orf in orfs), default=3)
    lookup: set[str] = set()
    for orf in orfs:
        upper = min(max_len, len(orf))
        for length in range(min_len, upper + 1):
            for start in range(0, len(orf) - length + 1):
                lookup.add(orf[start : start + length])
    return lookup


def filter_instanovo_chunked_substring(
    input_csv: Path,
    output_csv: Path,
    orf_sequences: list[str],
    *,
    strip_mods: bool,
) -> dict[str, int]:
    """Keep rows whose peptide is a substring of any ORF (Tier 1 Low)."""
    lookup = build_orf_peptide_lookup(orf_sequences, strip_mods=strip_mods)
    log.info("ORF substring lookup size: %d", len(lookup))
    cache: dict[str, bool] = {}
    input_rows = 0
    matched_rows = 0
    first = True
    if output_csv.exists():
        output_csv.unlink()

    for chunk in pd.read_csv(input_csv, chunksize=CHUNKSIZE, low_memory=False):
        input_rows += len(chunk)
        norms = []
        for value in chunk[SEQ_COL].astype(str):
            norm = normalize_sequence(value, strip_mods=strip_mods, il_equivalent=False)
            norms.append(norm)
            if norm not in cache:
                cache[norm] = norm in lookup
        keep_mask = [cache[n] for n in norms]
        kept = chunk.loc[keep_mask].copy()
        matched_rows += len(kept)
        kept.to_csv(output_csv, mode="w" if first else "a", header=first, index=False)
        first = False

    if first:
        # No chunks / empty input — write empty placeholder file.
        pd.DataFrame().to_csv(output_csv, index=False)

    return {"input_rows": input_rows, "matched_rows": matched_rows}


def filter_instanovo_chunked_exact(
    input_csv: Path,
    output_csv: Path,
    validation_set: set[str],
    *,
    strip_mods: bool,
) -> dict[str, int]:
    """Keep rows whose peptide exactly matches the validation set."""
    input_rows = 0
    matched_rows = 0
    first = True
    if output_csv.exists():
        output_csv.unlink()

    for chunk in pd.read_csv(input_csv, chunksize=CHUNKSIZE, low_memory=False):
        input_rows += len(chunk)
        norms = chunk[SEQ_COL].astype(str).apply(
            lambda value: normalize_sequence(value, strip_mods=strip_mods, il_equivalent=False)
        )
        keep_mask = norms.isin(validation_set)
        kept = chunk.loc[keep_mask].copy()
        matched_rows += len(kept)
        kept.to_csv(output_csv, mode="w" if first else "a", header=first, index=False)
        first = False

    if first:
        pd.DataFrame().to_csv(output_csv, index=False)

    return {"input_rows": input_rows, "matched_rows": matched_rows}


def run_update(config: PipelineConfig) -> dict:
    t0 = time.perf_counter()
    ensure_dirs(config)

    validation_df = read_table(config.validation_csv_path)
    validation_set = build_sequence_set(
        validation_df,
        config.validation_sequence_column,
        strip_mods=config.strip_modifications,
        il_equivalent=False,
    )
    log.info("Validation sequences: %d", len(validation_set))

    instanovo_files = scan_input_files(
        config.instanovo_filtered_dir,
        (config.instanovo_input_glob,),
    )
    log.info("InstaNovo filtered files: %s", [p.name for p in instanovo_files])

    summary: dict = {
        "instanovo_files": [p.name for p in instanovo_files],
        "tiers": {"med_instanovo": {}, "by_orf_run": {}},
    }

    # Tier 3 Medium-InstaNovo (ORF-independent)
    log.info("=== Tier 3: Medium-InstaNovo ===")
    for input_path in instanovo_files:
        out_name = derive_output_name(input_path, config.tier_suffixes["med_instanovo"])
        out_path = config.confidence_med_instanovo_dir / out_name
        stats = filter_instanovo_chunked_exact(
            input_path,
            out_path,
            validation_set,
            strip_mods=config.strip_modifications,
        )
        summary["tiers"]["med_instanovo"][input_path.name] = stats
        log.info("Tier3 %s: %s", input_path.name, stats)

    for orf_file in config.orf_files:
        run_label = orf_run_name(orf_file)
        ntv3_filtered_path = config.ntv3_filtered_output_dir / derive_output_name(
            config.ntv3_handoff_dir / orf_file,
            "_filtered",
        )
        if not ntv3_filtered_path.is_file():
            raise FileNotFoundError(
                f"Missing NTv3 filtered handoff: {ntv3_filtered_path}. "
                "Run Phase 1 proteome filter on NTv3 ORFs first."
            )

        log.info("=== ORF run: %s ===", run_label)
        ntv3_df = read_table(ntv3_filtered_path)
        ntv3_seqs = ntv3_df[config.ntv3_sequence_column].dropna().astype(str).tolist()

        # Tier 1 Low
        log.info("--- Tier 1: Low ---")
        low_summary: dict[str, dict[str, int]] = {}
        for input_path in instanovo_files:
            out_path = config.confidence_low_dir / run_label / f"{input_path.stem}.csv"
            stats = filter_instanovo_chunked_substring(
                input_path,
                out_path,
                ntv3_seqs,
                strip_mods=config.strip_modifications,
            )
            low_summary[input_path.name] = stats
            log.info("Low %s/%s: %s", run_label, input_path.name, stats)

        # Tier 2 Mid-NTv3
        log.info("--- Tier 2: Mid-NTv3 ---")
        med_matched = match_table(
            ntv3_df,
            config.ntv3_sequence_column,
            list(validation_set),
            mode="substring",
            direction="reference_in_query",
            strip_mods=config.strip_modifications,
            il_equivalent=False,
        )
        med_out = config.confidence_med_ntv3_dir / derive_output_name(
            ntv3_filtered_path, config.tier_suffixes["med_ntv3"]
        )
        write_table(med_matched, med_out)
        med_summary = {
            "input_rows": len(ntv3_df),
            "matched_rows": len(med_matched),
            "output_path": str(med_out),
        }
        log.info("Mid-NTv3 %s: %s", run_label, med_summary)

        # Tier 4 High (from Low outputs)
        log.info("--- Tier 4: High ---")
        high_summary: dict[str, dict[str, int]] = {}
        low_dir = config.confidence_low_dir / run_label
        for low_path in scan_input_files(low_dir, ("*.csv",)):
            # Skip empty files with no columns
            try:
                peek = pd.read_csv(low_path, nrows=0)
            except pd.errors.EmptyDataError:
                high_out = config.confidence_high_dir / run_label / derive_output_name(
                    low_path, config.tier_suffixes["high"]
                )
                write_table(pd.DataFrame(), high_out)
                high_summary[low_path.name] = {"input_rows": 0, "matched_rows": 0}
                continue

            if SEQ_COL not in peek.columns:
                # Empty placeholder CSV
                high_out = config.confidence_high_dir / run_label / derive_output_name(
                    low_path, config.tier_suffixes["high"]
                )
                write_table(pd.DataFrame(), high_out)
                high_summary[low_path.name] = {"input_rows": 0, "matched_rows": 0}
                continue

            high_out = config.confidence_high_dir / run_label / derive_output_name(
                low_path, config.tier_suffixes["high"]
            )
            stats = filter_instanovo_chunked_exact(
                low_path,
                high_out,
                validation_set,
                strip_mods=config.strip_modifications,
            )
            high_summary[low_path.name] = stats
            log.info("High %s/%s: %s", run_label, low_path.name, stats)

        summary["tiers"]["by_orf_run"][run_label] = {
            "orf_file": orf_file,
            "low": low_summary,
            "med_ntv3": med_summary,
            "high": high_summary,
        }

    summary["elapsed_seconds"] = round(time.perf_counter() - t0, 1)
    config.run_summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    log.info("Run summary: %s", config.run_summary_path)

    # Assembly
    log.info("=== Assembly ===")
    from assemble_validation_results import main as assemble_main

    assemble_main([])
    return summary


def main() -> int:
    summary = run_update(default_config())
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
