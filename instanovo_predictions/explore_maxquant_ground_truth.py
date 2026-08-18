#!/usr/bin/env python3
"""
Diagnostic script: explore MaxQuant ground-truth tables and verify scan mapping.

Run this once before generate_roc_curve.py to confirm paths, column headers,
and how InstaNovo scan indices link to MaxQuant sequences.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

CONFIG = {
    "maxquant_search_dir": "/iridisfs/ddnb/Charlotte/RPC/InstaNovo/test_data",
    "mgf_path": "/home/cle1g21/RPC/msconvert/output_data/01625b_GB2-TUM_first_pool_10_01_01-3xHCD-1h-R1.mgf",
    "predictions_csv": "/home/cle1g21/RPC/instanovo_predictions/01625b_GB2-TUM_first_pool_10_01_01-3xHCD-1h-R1_subset500.csv",
    "raw_file_name": "01625b_GB2-TUM_first_pool_10_01_01-3xHCD-1h-R1",
    "max_spectra": 500,
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Explore MaxQuant ground-truth files.")
    parser.add_argument("--maxquant-dir", default=None)
    parser.add_argument("--mgf-path", default=None)
    parser.add_argument("--predictions-csv", default=None)
    parser.add_argument("--raw-file-name", default=None)
    return parser.parse_args()


def print_table_headers(maxquant_dir: Path) -> None:
    for table_name in ("msms.txt", "peptides.txt", "msmsScans.txt"):
        table_path = maxquant_dir / table_name
        if not table_path.is_file():
            print(f"MISSING: {table_path}")
            continue
        columns = table_path.read_text().splitlines()[0].split("\t")
        print(f"\n=== {table_name} ({len(columns)} columns) ===")
        print(", ".join(columns[:20]), "..." if len(columns) > 20 else "")


def parse_mgf_scan_table(mgf_path: Path, max_spectra: int) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    title: str | None = None
    precursor_mz: float | None = None
    charge: int | None = None
    spectrum_index = 0

    with mgf_path.open() as mgf_file:
        for line in mgf_file:
            stripped = line.strip()
            if stripped.startswith("TITLE="):
                title = stripped.split("=", 1)[1]
            elif stripped.startswith("PEPMASS="):
                precursor_mz = float(stripped.split()[1])
            elif stripped.startswith("CHARGE="):
                charge = int(re.sub(r"[^0-9]", "", stripped))
            elif stripped == "END IONS":
                parts = title.split(".")  # type: ignore[union-attr]
                records.append(
                    {
                        "instanovo_scan_number": spectrum_index,
                        "instrument_scan_number": int(parts[1]) if len(parts) > 1 else -1,
                        "precursor_mz": precursor_mz,
                        "precursor_charge": charge,
                        "mgf_title": title,
                    }
                )
                spectrum_index += 1
                title = precursor_mz = charge = None
                if spectrum_index >= max_spectra:
                    break
    return pd.DataFrame(records)


def load_msms_scans_ground_truth(maxquant_dir: Path, raw_file_name: str) -> pd.DataFrame:
    scans = pd.read_csv(maxquant_dir / "msmsScans.txt", sep="\t")
    scans = scans[scans["Raw file"] == raw_file_name].copy()
    scans["Sequence"] = scans["Sequence"].fillna("").astype(str).str.strip()
    scans = scans[(scans["Identified"] == "+") & (scans["Sequence"] != "")]
    return scans.rename(
        columns={"Scan number": "instrument_scan_number", "Sequence": "ground_truth_sequence"}
    )[["instrument_scan_number", "ground_truth_sequence", "m/z", "Charge"]].drop_duplicates(
        "instrument_scan_number", keep="first"
    )


def load_msms_fallback_ground_truth(maxquant_dir: Path, raw_file_name: str) -> pd.DataFrame:
    msms = pd.read_csv(maxquant_dir / "msms.txt", sep="\t", low_memory=False)
    msms = msms[msms["Raw file"] == raw_file_name].sort_values("Score", ascending=False)
    msms = msms.drop_duplicates(["m/z", "Charge"], keep="first")
    return msms.rename(
        columns={"Sequence": "ground_truth_sequence_fb", "m/z": "precursor_mz", "Charge": "precursor_charge"}
    )[["precursor_mz", "precursor_charge", "ground_truth_sequence_fb", "Score"]]


def report_ground_truth_coverage(
    mgf_table: pd.DataFrame,
    scans_truth: pd.DataFrame,
    msms_fallback: pd.DataFrame,
) -> None:
    merged = mgf_table.merge(scans_truth, on="instrument_scan_number", how="left")
    primary_hits = int(merged["ground_truth_sequence"].notna().sum())

    missing = merged[merged["ground_truth_sequence"].isna()].copy()
    fallback = missing.merge(msms_fallback, on=["precursor_mz", "precursor_charge"], how="left")
    fallback_hits = int(fallback["ground_truth_sequence_fb"].notna().sum())
    total_hits = primary_hits + fallback_hits

    print(f"\n=== Ground-truth coverage (first {len(mgf_table)} MGF spectra) ===")
    print(f"Primary (msmsScans):  {primary_hits}")
    print(f"Fallback (msms.txt):  {fallback_hits}")
    print(f"Total with sequence:  {total_hits} / {len(mgf_table)} ({100 * total_hits / len(mgf_table):.1f}%)")

    example = merged[merged["ground_truth_sequence"].notna()].head(3)
    for _, row in example.iterrows():
        print(
            f"  scan {row['instanovo_scan_number']} → instrument {row['instrument_scan_number']} "
            f"→ {row['ground_truth_sequence']}"
        )


def main() -> int:
    args = parse_arguments()
    maxquant_dir = Path(args.maxquant_dir or CONFIG["maxquant_search_dir"])
    mgf_path = Path(args.mgf_path or CONFIG["mgf_path"])
    predictions_csv = Path(args.predictions_csv or CONFIG["predictions_csv"])
    raw_file_name = args.raw_file_name or CONFIG["raw_file_name"]

    print("=== MaxQuant ground-truth exploration ===")
    print_table_headers(maxquant_dir)
    mgf_table = parse_mgf_scan_table(mgf_path, CONFIG["max_spectra"])
    print(mgf_table.head().to_string(index=False))
    report_ground_truth_coverage(
        mgf_table,
        load_msms_scans_ground_truth(maxquant_dir, raw_file_name),
        load_msms_fallback_ground_truth(maxquant_dir, raw_file_name),
    )

    if predictions_csv.is_file():
        pred = pd.read_csv(predictions_csv, nrows=1)
        col = "instanovo_prediction_token_log_probabilities"
        print(f"\n{col} populated:", pred[col].notna().all())
        print("predictions example:", pred.iloc[0]["predictions"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
