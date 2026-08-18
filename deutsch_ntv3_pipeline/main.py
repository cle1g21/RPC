#!/usr/bin/env python3
"""CLI for Deutsch NTv3 validation pipeline (IridisX Stages 1, 3, 4).

Stage 2 (NTv3 GPU inference on DGX) is executed by the OpenClaw Agent inside the
``openclaw-research`` sandbox. See ``openclaw_dgx_blueprint.md`` for instructions.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Ensure project root is on sys.path for ``config`` and ``src`` imports.
_ROOT = Path(__file__).resolve().parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from config import PipelineConfig
from src.stage1_converter import run_stage1
from src.stage3_matcher import run_stage3
from src.stage4_plotter import run_stage4


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Optional argument list (defaults to ``sys.argv``).

    Returns:
        Parsed namespace with ``stage``, ``ntv3_output``, and ``verbose``.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Deutsch/Kok NTv3 validation pipeline (IridisX). "
            "Stages 1, 3, 4 run here; Stage 2 runs on DGX via OpenClaw (see openclaw_dgx_blueprint.md)."
        )
    )
    parser.add_argument(
        "--stage",
        choices=["1", "3", "4", "all"],
        default="all",
        help="Stage to run: 1=RDS→FASTA, 3=match, 4=plot, all=1+3+4 (default: all).",
    )
    parser.add_argument(
        "--ntv3-output",
        type=Path,
        default=None,
        help="Explicit NTv3 prediction CSV for Stage 3 (auto-discovered if omitted).",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug logging.")
    return parser.parse_args(argv)


def run_pipeline(args: argparse.Namespace, config: PipelineConfig) -> int:
    """Execute selected pipeline stages.

    Args:
        args: Parsed CLI arguments.
        config: Path configuration.

    Returns:
        Exit code 0 on success.
    """
    stages = ["1", "3", "4"] if args.stage == "all" else [args.stage]
    summary: dict[str, int] | None = None
    log = logging.getLogger(__name__)

    for stage in stages:
        if stage == "1":
            log.info("=== Stage 1: RDS → FASTA (IridisX) ===")
            run_stage1(config)
        elif stage == "3":
            log.info("=== Stage 3: coordinate matching (IridisX) ===")
            summary = run_stage3(config, ntv3_path=args.ntv3_output)
        elif stage == "4":
            log.info("=== Stage 4: validation bar chart (IridisX) ===")
            run_stage4(config, summary=summary)

    return 0


def main(argv: list[str] | None = None) -> int:
    """Entry point."""
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    config = PipelineConfig.from_env()
    try:
        return run_pipeline(args, config)
    except Exception:
        logging.getLogger(__name__).exception("Pipeline failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
