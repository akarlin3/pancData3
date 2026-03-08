#!/usr/bin/env python3
"""Orchestrator for the pancData3 analysis suite.

Runs the full post-pipeline analysis workflow:
1. Vision-based graph analysis (batch_graph_analysis.py)
2. Direct log file parsing (parse_log_metrics.py)
3. Direct CSV parsing (parse_csv_results.py)
4. Markdown report generation (generate_report.py)

Usage:
    python run_analysis.py                     # auto-detect latest folder
    python run_analysis.py --folder PATH       # specify folder
    python run_analysis.py --skip-vision       # skip Claude API calls
    python run_analysis.py --report-only       # only generate report
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

from shared import find_latest_saved_folder, setup_utf8_stdout

setup_utf8_stdout()

ANALYSIS_DIR = Path(__file__).resolve().parent


def _run_script(name: str, folder: Path) -> bool:
    """Run a sibling Python script, return True on success."""
    script = ANALYSIS_DIR / name
    if not script.exists():
        print(f"  [SKIP] {name} not found")
        return False

    print(f"\n  Running {name} ...")
    t0 = time.time()
    result = subprocess.run(
        [sys.executable, str(script), str(folder)],
        cwd=str(ANALYSIS_DIR),
        capture_output=False,
    )
    elapsed = time.time() - t0
    if result.returncode == 0:
        print(f"  Done: {name} ({elapsed:.1f}s)")
        return True
    else:
        print(f"  FAILED: {name} (exit code {result.returncode}, {elapsed:.1f}s)")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="pancData3 analysis suite orchestrator"
    )
    parser.add_argument(
        "--folder",
        type=str,
        default=None,
        help="Path to saved_files_* folder (default: auto-detect latest)",
    )
    parser.add_argument(
        "--skip-vision",
        action="store_true",
        help="Skip Claude API vision analysis (use existing CSV if available)",
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="Only generate the Markdown report from existing data",
    )
    args = parser.parse_args()

    # Resolve folder
    if args.folder:
        folder = Path(args.folder)
        if not folder.is_dir():
            sys.exit(f"ERROR: Folder does not exist: {folder}")
    else:
        folder = find_latest_saved_folder()

    print("=" * 70)
    print("  pancData3 Analysis Suite")
    print("=" * 70)
    print(f"  Output folder: {folder}")
    print(f"  Skip vision:   {args.skip_vision}")
    print(f"  Report only:   {args.report_only}")
    print()

    results = {}

    if not args.report_only:
        # Step 1: Vision-based graph analysis
        if not args.skip_vision:
            csv_path = folder / "graph_analysis_results.csv"
            print(f"  Vision CSV exists: {csv_path.exists()}")
            results["vision"] = _run_script("batch_graph_analysis.py", folder)
        else:
            csv_path = folder / "graph_analysis_results.csv"
            if csv_path.exists():
                print(f"  Skipping vision analysis (existing CSV: {csv_path})")
            else:
                print("  Skipping vision analysis (no existing CSV)")
            results["vision"] = "skipped"

        # Step 2: Parse log files
        results["logs"] = _run_script("parse_log_metrics.py", folder)

        # Step 3: Parse CSV exports
        results["csvs"] = _run_script("parse_csv_results.py", folder)

    # Step 4: Generate report
    results["report"] = _run_script("generate_report.py", folder)

    # Summary
    print("\n" + "=" * 70)
    print("  Analysis Complete")
    print("=" * 70)
    for step, status in results.items():
        icon = "OK" if status is True else ("SKIP" if status == "skipped" else "FAIL")
        print(f"  [{icon:>4}] {step}")

    report_path = folder / "analysis_report.md"
    if report_path.exists():
        print(f"\n  Report: {report_path}")
    print()


if __name__ == "__main__":
    main()
