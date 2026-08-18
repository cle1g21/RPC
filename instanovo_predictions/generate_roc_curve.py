#!/usr/bin/env python3
"""
Amino-acid-level ROC validation: InstaNovo predictions vs MaxQuant ground truth.

Enhanced pipeline with token alignment fix, mass-equivalent scoring,
stratified metrics, and multi-panel publication figures.
"""

from __future__ import annotations

import argparse
import ast
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec
from sklearn.metrics import auc, roc_curve

# =============================================================================
# CONFIG — edit paths and behaviour here for new cohorts
# =============================================================================

CONFIG: dict[str, Any] = {
    "predictions_csv": "/home/cle1g21/RPC/instanovo_predictions/01625b_GB2-TUM_first_pool_10_01_01-3xHCD-1h-R1_predictions.csv",
    "mgf_path": "/home/cle1g21/RPC/msconvert/output_data/01625b_GB2-TUM_first_pool_10_01_01-3xHCD-1h-R1.mgf",
    "maxquant_search_dir": "/iridisfs/ddnb/Charlotte/RPC/InstaNovo/test_data",
    "raw_file_name": "01625b_GB2-TUM_first_pool_10_01_01-3xHCD-1h-R1",
    "sample_label": "01625b_GB2-TUM_full",
    "scan_index_column": "scan_number",
    "prediction_sequence_column": "predictions",
    "prediction_tokenised_column": "predictions_tokenised",
    "token_log_probs_column": "instanovo_prediction_token_log_probabilities",
    # When streaming resumed onto a CombinedPredictor header, token log-probs may
    # land in this shifted column instead (see README).
    "token_log_probs_fallback_column": "instanovo_prediction_log_probability",
    # Drop subset500 seed rows so scans 0–499 are not double-counted vs the full run.
    "prefer_experiment_name": "01625b_GB2-TUM_first_pool_10_01_01-3xHCD-1h-R1",
    "max_spectra": None,  # None = use all spectra in the predictions CSV
    "output_plot_png": "/home/cle1g21/RPC/instanovo_predictions/plots/instanovoVal_roc_curve.png",
    "output_plot_enhanced_png": "/home/cle1g21/RPC/instanovo_predictions/plots/instanovoVal_roc_curveEnhanced.png",
    "mass_equivalent_li": True,
    "mass_equivalent_kq": False,
    "min_aa_per_stratum": 20,
    "length_bins": {"short": (0, 9), "medium": (10, 15), "long": (16, 999)},
    "charge_bins": {"plus2": [2], "plus3plus": [3, 4, 5, 6, 7, 8, 9, 10]},
}

IGNORED_COLUMNS = ("token_log_probs", "instanovoplus_prediction_token_log_probabilities")

INSTANOVO_RESIDUE_PATTERN = re.compile(r"([A-Z]\[UNIMOD:\d+\]|[A-Z])")


@dataclass
class AlignmentResult:
    """Tokens, probabilities, and how they were aligned for one spectrum."""

    tokens: list[str]
    token_log_probs: list[float]
    alignment_mode: str
    tokenization_path: str


@dataclass
class Observation:
    """One scored amino acid for ROC and stratification."""

    label_strict: int
    label_mass_equiv: int
    score: float
    length_bin: str
    charge_bin: str
    scan_number: object = None


@dataclass
class PipelineStats:
    """Summary counters printed at the end of a run."""

    spectra_with_ground_truth: int = 0
    spectra_used: int = 0
    spectra_skipped_probs: int = 0
    spectra_exact_match: int = 0
    spectra_recovered_via_tokenised: int = 0
    spectra_truncated_alignment: int = 0
    spectra_still_unusable: int = 0
    amino_acids_scored: int = 0
    strict_accuracy: float = 0.0
    mass_equiv_accuracy: float = 0.0
    strict_auc: float = 0.0
    mass_equiv_auc: float = 0.0
    stratified: dict[str, dict[str, float]] = field(default_factory=dict)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate amino-acid ROC curves for InstaNovo.")
    parser.add_argument("--predictions-csv", default=None)
    parser.add_argument("--mgf-path", default=None)
    parser.add_argument("--maxquant-dir", default=None)
    parser.add_argument("--output-plot", default=None)
    parser.add_argument("--output-enhanced-plot", default=None)
    parser.add_argument("--sample-label", default=None)
    parser.add_argument("--raw-file-name", default=None)
    parser.add_argument("--mass-equiv-li", action="store_true", default=None)
    parser.add_argument("--no-mass-equiv-li", action="store_true")
    parser.add_argument("--mass-equiv-kq", action="store_true", default=None)
    parser.add_argument("--diagnose-mismatch", action="store_true")
    return parser.parse_args()


def resolve_config(args: argparse.Namespace) -> dict[str, Any]:
    cfg = dict(CONFIG)
    if args.predictions_csv:
        cfg["predictions_csv"] = args.predictions_csv
    if args.mgf_path:
        cfg["mgf_path"] = args.mgf_path
    if args.maxquant_dir:
        cfg["maxquant_search_dir"] = args.maxquant_dir
    if args.output_plot:
        cfg["output_plot_png"] = args.output_plot
    if args.output_enhanced_plot:
        cfg["output_plot_enhanced_png"] = args.output_enhanced_plot
    if args.sample_label:
        cfg["sample_label"] = args.sample_label
    if args.raw_file_name:
        cfg["raw_file_name"] = args.raw_file_name
    if args.mass_equiv_li is True:
        cfg["mass_equivalent_li"] = True
    if args.no_mass_equiv_li:
        cfg["mass_equivalent_li"] = False
    if args.mass_equiv_kq is True:
        cfg["mass_equivalent_kq"] = True
    return cfg


def parse_token_log_probs(raw_value: object) -> list[float]:
    if raw_value is None or (isinstance(raw_value, float) and np.isnan(raw_value)):
        raise ValueError("empty token log probs")
    text = str(raw_value).strip()
    if not text or text.lower() == "nan":
        raise ValueError("empty token log probs")
    # Must look like a Python list of numbers; peptide strings parse as Names.
    if not text.startswith("["):
        raise ValueError(f"expected list literal, got {text[:40]!r}")
    parsed = ast.literal_eval(text)
    if not isinstance(parsed, list) or not parsed:
        raise ValueError(f"expected non-empty list, got {type(parsed)}")
    if not all(isinstance(x, (int, float)) for x in parsed):
        raise ValueError("list does not contain numeric log probs")
    return [float(x) for x in parsed]


def extract_token_log_probs(row: pd.Series, cfg: dict[str, Any]) -> list[float]:
    """Read token log-probs, recovering from the GB2 streaming column shift."""
    primary = cfg["token_log_probs_column"]
    fallback = cfg.get("token_log_probs_fallback_column")
    errors: list[str] = []
    for column in (primary, fallback):
        if not column or column not in row.index:
            continue
        try:
            return parse_token_log_probs(row[column])
        except (ValueError, SyntaxError, TypeError) as exc:
            errors.append(f"{column}: {exc}")
    raise ValueError("; ".join(errors) if errors else "no token log-prob column available")


def tokenize_from_csv_field(raw_value: object) -> list[str]:
    if raw_value is None or (isinstance(raw_value, float) and np.isnan(raw_value)):
        return []
    return [token.strip() for token in str(raw_value).split(",") if token.strip()]


def tokenize_with_instanovo_regex(sequence: str) -> list[str]:
    return INSTANOVO_RESIDUE_PATTERN.findall(str(sequence).strip().upper())


def tokenize_as_characters(sequence: str) -> list[str]:
    return list(str(sequence).strip().upper())


def tokenize_ground_truth(sequence: str) -> list[str]:
    if "[" in sequence:
        tokens = tokenize_with_instanovo_regex(sequence)
        if tokens:
            return tokens
    return tokenize_as_characters(sequence)


def parse_mgf_scan_table(mgf_path: Path, max_spectra: int | None) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    title: str | None = None
    precursor_mz: float | None = None
    charge: int | None = None
    spectrum_index = 0

    with mgf_path.open() as handle:
        for line in handle:
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
                    }
                )
                spectrum_index += 1
                title = precursor_mz = charge = None
                if max_spectra is not None and spectrum_index >= max_spectra:
                    break
    return pd.DataFrame(records)


def load_ground_truth(maxquant_dir: Path, raw_file_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    scans = pd.read_csv(maxquant_dir / "msmsScans.txt", sep="\t")
    scans = scans[scans["Raw file"] == raw_file_name].copy()
    scans["Sequence"] = scans["Sequence"].fillna("").astype(str).str.strip()
    scans = scans[(scans["Identified"] == "+") & (scans["Sequence"] != "")]
    primary = scans.rename(
        columns={"Scan number": "instrument_scan_number", "Sequence": "ground_truth_sequence"}
    )[["instrument_scan_number", "ground_truth_sequence"]].drop_duplicates("instrument_scan_number")

    msms = pd.read_csv(maxquant_dir / "msms.txt", sep="\t", low_memory=False)
    msms = msms[msms["Raw file"] == raw_file_name].sort_values("Score", ascending=False)
    msms = msms.drop_duplicates(["m/z", "Charge"], keep="first")
    fallback = msms.rename(
        columns={"Sequence": "ground_truth_sequence_fb", "m/z": "precursor_mz", "Charge": "precursor_charge"}
    )[["precursor_mz", "precursor_charge", "ground_truth_sequence_fb"]]
    return primary, fallback


def align_predictions_with_ground_truth(cfg: dict[str, Any]) -> pd.DataFrame:
    predictions = pd.read_csv(cfg["predictions_csv"], low_memory=False)
    prefer = cfg.get("prefer_experiment_name")
    if prefer and "experiment_name" in predictions.columns:
        preferred = predictions["experiment_name"].astype(str) == str(prefer)
        if preferred.any():
            n_drop = int((~preferred).sum())
            predictions = predictions.loc[preferred].copy()
            print(f"Filtered to experiment_name={prefer!r} ({len(predictions)} rows; dropped {n_drop})")
        else:
            predictions = predictions.copy()
            print(f"Warning: prefer_experiment_name={prefer!r} matched 0 rows; keeping all {len(predictions)}")
    else:
        predictions = predictions.copy()
    if cfg.get("max_spectra") is not None:
        predictions = predictions.head(int(cfg["max_spectra"]))
    mgf_map = parse_mgf_scan_table(Path(cfg["mgf_path"]), cfg.get("max_spectra"))
    primary_truth, fallback_truth = load_ground_truth(Path(cfg["maxquant_search_dir"]), cfg["raw_file_name"])

    predictions = predictions.drop(columns=["precursor_mz", "precursor_charge"], errors="ignore")
    merged = predictions.merge(
        mgf_map, left_on=cfg["scan_index_column"], right_on="instanovo_scan_number", how="left"
    )
    merged = merged.merge(primary_truth, on="instrument_scan_number", how="left")
    merged = merged.merge(fallback_truth, on=["precursor_mz", "precursor_charge"], how="left")
    merged["ground_truth_sequence"] = merged["ground_truth_sequence"].fillna(merged["ground_truth_sequence_fb"])
    return merged[merged["ground_truth_sequence"].notna()].copy()


def validate_token_prob_column(predictions: pd.DataFrame, cfg: dict[str, Any]) -> None:
    column = cfg["token_log_probs_column"]
    if column not in predictions.columns:
        raise ValueError(f"Missing column: {column}")
    fallback = cfg.get("token_log_probs_fallback_column")
    primary_ok = predictions[column].notna().any()
    fallback_ok = bool(fallback) and fallback in predictions.columns and predictions[fallback].notna().any()
    if not primary_ok and not fallback_ok:
        raise ValueError(f"Column '{column}' is empty — cannot compute ROC.")


def align_tokens_and_probs(row: pd.Series, cfg: dict[str, Any]) -> AlignmentResult | None:
    sequence_text = str(row[cfg["prediction_sequence_column"]])
    tokenised_field = row.get(cfg["prediction_tokenised_column"])

    try:
        token_log_probs = extract_token_log_probs(row, cfg)
    except (ValueError, SyntaxError):
        return None

    candidate_paths: list[tuple[str, list[str]]] = []

    tokenised_tokens = tokenize_from_csv_field(tokenised_field)
    if tokenised_tokens:
        candidate_paths.append(("predictions_tokenised", tokenised_tokens))

    regex_tokens = tokenize_with_instanovo_regex(sequence_text)
    if regex_tokens:
        candidate_paths.append(("instanovo_regex", regex_tokens))

    if "[" not in sequence_text:
        char_tokens = tokenize_as_characters(sequence_text)
        if char_tokens:
            candidate_paths.append(("characters", char_tokens))

    for path_name, tokens in candidate_paths:
        if len(tokens) == len(token_log_probs):
            return AlignmentResult(tokens, token_log_probs, "exact", path_name)

    if not candidate_paths:
        return None

    path_name, tokens = candidate_paths[0]
    if not tokens or not token_log_probs:
        return None

    aligned_length = min(len(tokens), len(token_log_probs))
    if aligned_length == 0:
        return None

    return AlignmentResult(
        tokens=tokens[:aligned_length],
        token_log_probs=token_log_probs[:aligned_length],
        alignment_mode="truncated",
        tokenization_path=path_name,
    )


def diagnose_token_mismatches(aligned: pd.DataFrame, cfg: dict[str, Any], limit: int = 5) -> dict[str, int]:
    counts = {
        "exact": 0,
        "truncated": 0,
        "unusable": 0,
        "recovered_via_tokenised": 0,
    }
    printed = 0

    print("\n=== Token mismatch diagnostics (first mismatches) ===")
    for _, row in aligned.iterrows():
        sequence_text = str(row[cfg["prediction_sequence_column"]])
        char_tokens = tokenize_as_characters(sequence_text)
        try:
            prob_len = len(extract_token_log_probs(row, cfg))
        except (ValueError, SyntaxError):
            counts["unusable"] += 1
            continue

        alignment = align_tokens_and_probs(row, cfg)
        if alignment is None:
            counts["unusable"] += 1
            if printed < limit:
                printed += 1
                print(f"\nScan {row.get(cfg['scan_index_column'])}: UNUSABLE")
                print(f"  predictions          : {sequence_text}")
                print(f"  predictions_tokenised: {row.get(cfg['prediction_tokenised_column'])}")
                print(f"  prob length          : {prob_len}")
                print(f"  char length          : {len(char_tokens)}")
                print(f"  beam-0               : {row.get('instanovo_predictions_beam_0', '')}")
            continue

        counts[alignment.alignment_mode] += 1
        if alignment.alignment_mode == "truncated":
            if printed < limit:
                printed += 1
                print(f"\nScan {row.get(cfg['scan_index_column'])}: TRUNCATED via {alignment.tokenization_path}")
                print(f"  predictions          : {sequence_text}")
                print(f"  predictions_tokenised: {row.get(cfg['prediction_tokenised_column'])}")
                print(f"  aligned tokens       : {len(alignment.tokens)}")
                print(f"  prob length          : {prob_len}")
        elif (
            len(char_tokens) != prob_len
            and alignment.tokenization_path == "predictions_tokenised"
        ):
            counts["recovered_via_tokenised"] += 1
            if printed < limit:
                printed += 1
                print(f"\nScan {row.get(cfg['scan_index_column'])}: RECOVERED via predictions_tokenised")
                print(f"  predictions          : {sequence_text}")
                print(f"  char length          : {len(char_tokens)}")
                print(f"  tokenised length     : {len(alignment.tokens)}")
                print(f"  prob length          : {prob_len}")

    print("\nSummary:")
    print(f"  exact alignments           : {counts['exact']}")
    print(f"  truncated alignments       : {counts['truncated']}")
    print(f"  recovered via tokenised    : {counts['recovered_via_tokenised']}")
    print(f"  unusable                   : {counts['unusable']}")
    return counts


def residues_match_strict(predicted: str, true: str) -> bool:
    return predicted == true


def residues_match_mass_equivalent(predicted: str, true: str, cfg: dict[str, Any]) -> bool:
    if residues_match_strict(predicted, true):
        return True
    if cfg["mass_equivalent_li"] and {predicted, true} == {"L", "I"}:
        return True
    if cfg["mass_equivalent_kq"] and {predicted, true} == {"K", "Q"}:
        return True
    return False


def assign_length_bin(length: int, cfg: dict[str, Any]) -> str:
    bins: dict[str, tuple[int, int]] = cfg["length_bins"]
    if bins["short"][0] <= length <= bins["short"][1]:
        return "short"
    if bins["medium"][0] <= length <= bins["medium"][1]:
        return "medium"
    return "long"


def assign_charge_bin(charge: object, cfg: dict[str, Any]) -> str:
    charge_value = int(charge)
    if charge_value in cfg["charge_bins"]["plus2"]:
        return "charge_2"
    if charge_value in cfg["charge_bins"]["plus3plus"]:
        return "charge_3plus"
    return "other_charge"


def build_observation_table(aligned: pd.DataFrame, cfg: dict[str, Any]) -> tuple[list[Observation], PipelineStats]:
    observations: list[Observation] = []
    stats = PipelineStats(spectra_with_ground_truth=len(aligned))

    for _, row in aligned.iterrows():
        alignment = align_tokens_and_probs(row, cfg)
        if alignment is None:
            stats.spectra_skipped_probs += 1
            stats.spectra_still_unusable += 1
            continue

        if alignment.alignment_mode == "exact":
            stats.spectra_exact_match += 1
            if alignment.tokenization_path == "predictions_tokenised":
                sequence_text = str(row[cfg["prediction_sequence_column"]])
                if len(tokenize_as_characters(sequence_text)) != len(alignment.tokens):
                    stats.spectra_recovered_via_tokenised += 1
        else:
            stats.spectra_truncated_alignment += 1

        stats.spectra_used += 1
        true_tokens = tokenize_ground_truth(str(row["ground_truth_sequence"]))
        peptide_length = len(alignment.tokens)
        length_bin = assign_length_bin(peptide_length, cfg)
        charge_bin = assign_charge_bin(row.get("precursor_charge"), cfg)

        for position, predicted_token in enumerate(alignment.tokens):
            token_probability = float(np.exp(alignment.token_log_probs[position]))
            true_token = true_tokens[position] if position < len(true_tokens) else ""
            label_strict = int(residues_match_strict(predicted_token, true_token))
            label_mass_equiv = int(
                residues_match_mass_equivalent(predicted_token, true_token, cfg)
            )
            observations.append(
                Observation(
                    label_strict=label_strict,
                    label_mass_equiv=label_mass_equiv,
                    score=token_probability,
                    length_bin=length_bin,
                    charge_bin=charge_bin,
                    scan_number=row.get(cfg["scan_index_column"]),
                )
            )

    if not observations:
        raise RuntimeError("No amino-acid observations collected.")

    stats.amino_acids_scored = len(observations)
    stats.strict_accuracy = float(np.mean([obs.label_strict for obs in observations]))
    stats.mass_equiv_accuracy = float(np.mean([obs.label_mass_equiv for obs in observations]))
    return observations, stats


def compute_roc_from_observations(
    observations: list[Observation],
    label_field: str,
) -> tuple[np.ndarray, np.ndarray, float]:
    labels = np.array([getattr(obs, label_field) for obs in observations])
    scores = np.array([obs.score for obs in observations])
    if len(np.unique(labels)) < 2:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0]), float("nan")
    fpr, tpr, _ = roc_curve(labels, scores)
    return fpr, tpr, float(auc(fpr, tpr))


def compute_stratified_metrics(
    observations: list[Observation],
    cfg: dict[str, Any],
) -> dict[str, dict[str, float]]:
    strata_definitions = {
        "length_short": lambda obs: obs.length_bin == "short",
        "length_medium": lambda obs: obs.length_bin == "medium",
        "length_long": lambda obs: obs.length_bin == "long",
        "charge_2": lambda obs: obs.charge_bin == "charge_2",
        "charge_3plus": lambda obs: obs.charge_bin == "charge_3plus",
    }

    results: dict[str, dict[str, float]] = {}
    min_aa = int(cfg["min_aa_per_stratum"])

    print("\n=== Stratified metrics ===")
    print(f"{'Stratum':<16} {'N_AA':>6} {'Strict_AUC':>12} {'MassEq_AUC':>12} {'Strict_Acc':>12}")
    print("-" * 62)

    for stratum_name, predicate in strata_definitions.items():
        subset = [obs for obs in observations if predicate(obs)]
        if len(subset) < min_aa:
            print(f"{stratum_name:<16} {len(subset):>6} {'(insufficient data)':>26}")
            continue

        _, _, strict_auc = compute_roc_from_observations(subset, "label_strict")
        _, _, mass_auc = compute_roc_from_observations(subset, "label_mass_equiv")
        strict_acc = float(np.mean([obs.label_strict for obs in subset]))

        results[stratum_name] = {
            "n_aa": len(subset),
            "strict_auc": strict_auc,
            "mass_equiv_auc": mass_auc,
            "strict_accuracy": strict_acc,
        }
        strict_auc_text = f"{strict_auc:.3f}" if not np.isnan(strict_auc) else "n/a"
        mass_auc_text = f"{mass_auc:.3f}" if not np.isnan(mass_auc) else "n/a"
        print(
            f"{stratum_name:<16} {len(subset):>6} {strict_auc_text:>12} {mass_auc_text:>12} {strict_acc:>12.3f}"
        )

    return results


def plot_roc_curve(
    fpr: np.ndarray,
    tpr: np.ndarray,
    auc_score: float,
    sample_label: str,
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
    ax.plot(fpr, tpr, color="#1f77b4", linewidth=2, label=f"Strict (AUC = {auc_score:.3f})")
    ax.plot([0, 1], [0, 1], linestyle="--", color="#888888", linewidth=1, label="Random (y = x)")
    ax.set_xlabel("False Positive Rate (1 - Specificity)")
    ax.set_ylabel("True Positive Rate (Sensitivity)")
    ax.set_title(f"InstaNovo Amino-Acid ROC — {sample_label}")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(output_path, format="png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def _plot_stratified_panel(
    ax: plt.Axes,
    observations: list[Observation],
    strata: list[tuple[str, str]],
    title: str,
    min_aa: int,
) -> None:
    colors = ["#4c72b0", "#dd8452", "#55a868", "#c44e52"]
    plotted = False

    for index, (stratum_name, label_field_filter) in enumerate(strata):
        if "length" in stratum_name:
            subset = [obs for obs in observations if obs.length_bin == stratum_name.replace("length_", "")]
        else:
            subset = [obs for obs in observations if obs.charge_bin == stratum_name]

        if len(subset) < min_aa:
            continue

        fpr, tpr, auc_score = compute_roc_from_observations(subset, "label_strict")
        if np.isnan(auc_score):
            continue

        ax.plot(
            fpr,
            tpr,
            color=colors[index % len(colors)],
            linewidth=1.8,
            label=f"{stratum_name.replace('_', ' ')} (AUC={auc_score:.2f})",
        )
        plotted = True

    ax.plot([0, 1], [0, 1], linestyle="--", color="#bbbbbb", linewidth=1)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    ax.set_title(title)
    if plotted:
        ax.legend(loc="lower right", fontsize=8)
    else:
        ax.text(0.5, 0.5, "Insufficient data", ha="center", va="center", transform=ax.transAxes)


def plot_enhanced_roc_figure(
    strict_curve: tuple[np.ndarray, np.ndarray, float],
    mass_curve: tuple[np.ndarray, np.ndarray, float],
    observations: list[Observation],
    cfg: dict[str, Any],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sns.set_style("whitegrid")

    fig = plt.figure(figsize=(7, 7), dpi=300)
    grid = GridSpec(2, 2, figure=fig, height_ratios=[1.4, 1.0], hspace=0.35, wspace=0.28)
    ax_main = fig.add_subplot(grid[0, :])
    ax_length = fig.add_subplot(grid[1, 0])
    ax_charge = fig.add_subplot(grid[1, 1])

    strict_fpr, strict_tpr, strict_auc = strict_curve
    mass_fpr, mass_tpr, mass_auc = mass_curve

    ax_main.plot(
        strict_fpr,
        strict_tpr,
        color="#1f77b4",
        linewidth=2.2,
        label=f"Strict (AUC = {strict_auc:.3f})",
    )
    ax_main.plot(
        mass_fpr,
        mass_tpr,
        color="#ff7f0e",
        linewidth=2.2,
        label=f"Mass-equivalent (AUC = {mass_auc:.3f})",
    )
    ax_main.plot([0, 1], [0, 1], linestyle="--", color="#888888", linewidth=1, label="Random (y = x)")
    ax_main.set_xlabel("False Positive Rate (1 - Specificity)")
    ax_main.set_ylabel("True Positive Rate (Sensitivity)")
    ax_main.set_title(f"InstaNovo Amino-Acid ROC — {cfg['sample_label']}")
    ax_main.set_xlim(0.0, 1.0)
    ax_main.set_ylim(0.0, 1.0)
    ax_main.legend(loc="lower right")

    min_aa = int(cfg["min_aa_per_stratum"])
    _plot_stratified_panel(
        ax_length,
        observations,
        [("length_short", ""), ("length_medium", ""), ("length_long", "")],
        "By peptide length (strict)",
        min_aa,
    )
    _plot_stratified_panel(
        ax_charge,
        observations,
        [("charge_2", ""), ("charge_3plus", "")],
        "By precursor charge (strict)",
        min_aa,
    )

    fig.savefig(output_path, format="png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main() -> int:
    args = parse_arguments()
    cfg = resolve_config(args)

    print("=== InstaNovo amino-acid ROC validation (enhanced) ===")
    print(f"Sequence column : {cfg['prediction_sequence_column']}")
    print(f"Token prob col  : {cfg['token_log_probs_column']}")
    print(f"Mass-equiv L/I  : {cfg['mass_equivalent_li']}")
    print(f"Mass-equiv K/Q  : {cfg['mass_equivalent_kq']}")
    print(f"Ignored columns : {', '.join(IGNORED_COLUMNS)}")

    raw_predictions = pd.read_csv(cfg["predictions_csv"], low_memory=False)
    validate_token_prob_column(raw_predictions, cfg)

    aligned = align_predictions_with_ground_truth(cfg)
    print(f"Spectra with ground truth: {len(aligned)}")

    if args.diagnose_mismatch:
        diagnose_token_mismatches(aligned, cfg)
        return 0

    observations, stats = build_observation_table(aligned, cfg)
    strict_fpr, strict_tpr, stats.strict_auc = compute_roc_from_observations(observations, "label_strict")
    mass_fpr, mass_tpr, stats.mass_equiv_auc = compute_roc_from_observations(
        observations, "label_mass_equiv"
    )
    stats.stratified = compute_stratified_metrics(observations, cfg)

    print("\n=== Global summary ===")
    print(f"Spectra used                    : {stats.spectra_used}")
    print(f"Exact token/prob alignment      : {stats.spectra_exact_match}")
    print(f"Recovered via tokenised field   : {stats.spectra_recovered_via_tokenised}")
    print(f"Truncated alignments            : {stats.spectra_truncated_alignment}")
    print(f"Skipped / unusable              : {stats.spectra_still_unusable + stats.spectra_skipped_probs}")
    print(f"Amino acids scored              : {stats.amino_acids_scored}")
    print(f"Strict position accuracy        : {stats.strict_accuracy:.4f}")
    print(f"Mass-equivalent accuracy        : {stats.mass_equiv_accuracy:.4f}")
    print(f"Strict AUC-ROC                  : {stats.strict_auc:.4f}")
    print(f"Mass-equivalent AUC-ROC         : {stats.mass_equiv_auc:.4f}")

    baseline_path = Path(cfg["output_plot_png"])
    enhanced_path = Path(cfg["output_plot_enhanced_png"])

    plot_roc_curve(strict_fpr, strict_tpr, stats.strict_auc, cfg["sample_label"], baseline_path)
    plot_enhanced_roc_figure(
        (strict_fpr, strict_tpr, stats.strict_auc),
        (mass_fpr, mass_tpr, stats.mass_equiv_auc),
        observations,
        cfg,
        enhanced_path,
    )

    print(f"\nBaseline plot saved to : {baseline_path}")
    print(f"Enhanced plot saved to : {enhanced_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
