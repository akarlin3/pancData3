#!/usr/bin/env python3
"""Orchestrator for the pancData3 analysis suite.

Runs the full post-pipeline analysis workflow:
1. Vision-based graph analysis (batch_graph_analysis.py)
2. Direct log file parsing (parse_log_metrics.py)
3. Direct CSV parsing (parse_csv_results.py)
4. HTML report generation (generate_report.py)

Each step is executed as a subprocess so that failures in one step do not
prevent subsequent steps from running.  The final summary table indicates
which steps succeeded, failed, or were skipped.

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

# Ensure emoji and special characters print correctly on Windows consoles.
setup_utf8_stdout()

# Absolute path to the analysis/ directory (where this script lives).
ANALYSIS_DIR = Path(__file__).resolve().parent


def _run_script(name: str, folder: Path) -> bool:
    """Run a sibling Python script as a subprocess.

    Parameters
    ----------
    name : str
        Filename of the script to run (e.g. ``"parse_log_metrics.py"``).
        Must reside in the same directory as this orchestrator.
    folder : Path
        Path to the ``saved_files_*`` output folder, passed as the first
        positional argument to the child script.

    Returns
    -------
    bool
        ``True`` if the script exited with code 0, ``False`` otherwise
        (including when the script file is not found).
    """
    script = ANALYSIS_DIR / name
    if not script.exists():
        print(f"  [SKIP] {name} not found")
        return False

    print(f"\n  Running {name} ...")
    t0 = time.time()
    # Run the child script using the same Python interpreter as this process.
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
    """Entry point: parse CLI arguments and execute the analysis pipeline."""
    # ── Argument parsing ──
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
        help="Skip Gemini API vision analysis (use existing CSV if available)",
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="Only generate the HTML report from existing data",
    )
    parser.add_argument(
        "--no-pdf",
        action="store_true",
        help="Skip PDF generation (produce HTML report only)",
    )
    args = parser.parse_args()

    # ── Resolve output folder ──
    if args.folder:
        folder = Path(args.folder)
        if not folder.is_dir():
            sys.exit(f"ERROR: Folder does not exist: {folder}")
    else:
        # Auto-detect the most recent saved_files_* directory.
        folder = find_latest_saved_folder()

    # ── Print banner ──
    print("=" * 70)
    print("  pancData3 Analysis Suite")
    print("=" * 70)
    print(f"  Output folder: {folder}")
    print(f"  Skip vision:   {args.skip_vision}")
    print(f"  Report only:   {args.report_only}")
    print()

    # Track success/failure/skip status for each pipeline step.
    results = {}

    if not args.report_only:
        # Step 1: Vision-based graph analysis (requires GEMINI_API_KEY).
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

        # Step 2: Parse MATLAB diary log files for structured metrics.
        results["logs"] = _run_script("parse_log_metrics.py", folder)

        # Step 3: Parse pipeline-exported CSV files (significance tables).
        results["csvs"] = _run_script("parse_csv_results.py", folder)

        # Step 3.5: Parse MATLAB .mat files (core comparison, dosimetry).
        results["mat"] = _run_script("parse_mat_metrics.py", folder)

    # Step 4: Assemble the final HTML (+PDF) report from all collected data.
    report_script = ANALYSIS_DIR / "generate_report.py"
    if report_script.exists():
        report_args = [sys.executable, str(report_script), str(folder)]
        if args.no_pdf:
            report_args.append("--no-pdf")
        print(f"\n  Running generate_report.py ...")
        t0 = time.time()
        result = subprocess.run(
            report_args, cwd=str(ANALYSIS_DIR), capture_output=False,
        )
        elapsed = time.time() - t0
        if result.returncode == 0:
            print(f"  Done: generate_report.py ({elapsed:.1f}s)")
            results["report"] = True
        else:
            print(f"  FAILED: generate_report.py (exit code {result.returncode}, {elapsed:.1f}s)")
            results["report"] = False
    else:
        results["report"] = False

    # ── Summary table ──
    print("\n" + "=" * 70)
    print("  Analysis Complete")
    print("=" * 70)
    for step, status in results.items():
        icon = "OK" if status is True else ("SKIP" if status == "skipped" else "FAIL")
        print(f"  [{icon:>4}] {step}")

    report_path = folder / "analysis_report.html"
    pdf_path = folder / "analysis_report.pdf"
    if report_path.exists():
        print(f"\n  HTML Report: {report_path}")
    if pdf_path.exists():
        print(f"  PDF Report:  {pdf_path}")
    print()


if __name__ == "__main__":
    main()
