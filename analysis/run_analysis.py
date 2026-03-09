#!/usr/bin/env python3
"""Orchestrator for the pancData3 analysis suite.

Runs the full post-pipeline analysis workflow:
0. Pre-flight checks: verify requirements and run test suite
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
    python run_analysis.py --skip-checks       # skip pre-flight checks
"""

from __future__ import annotations

import argparse
import io
import subprocess
import sys
import time
from pathlib import Path

from tqdm import tqdm

from shared import find_latest_saved_folder, get_config, load_analysis_config, reset_config_cache, setup_utf8_stdout

# Ensure emoji and special characters print correctly on Windows consoles.
setup_utf8_stdout()

# Absolute path to the analysis/ directory (where this script lives).
ANALYSIS_DIR = Path(__file__).resolve().parent


class TeeWriter:
    """Write to both a terminal stream and a log file simultaneously.

    Implements enough of the file-like interface (``write``, ``flush``,
    ``isatty``, ``encoding``, ``fileno``) to be usable as a drop-in
    ``sys.stdout`` / ``sys.stderr`` replacement.
    """

    def __init__(self, terminal: io.TextIOBase, log_file: io.TextIOBase):
        self.terminal = terminal
        self.log_file = log_file
        self.encoding = getattr(terminal, "encoding", "utf-8")

    def write(self, message: str) -> int:
        self.terminal.write(message)
        self.log_file.write(message)
        return len(message)

    def flush(self) -> None:
        self.terminal.flush()
        self.log_file.flush()

    def isatty(self) -> bool:
        """Return the terminal stream's isatty status."""
        return hasattr(self.terminal, "isatty") and self.terminal.isatty()

    def fileno(self) -> int:
        """Delegate fileno to the terminal stream (for subprocess piping)."""
        return self.terminal.fileno()


def _check_requirements() -> bool:
    """Verify that all packages listed in requirements.txt are installed.

    Returns True if every requirement is satisfied, False otherwise.
    Missing packages are printed to stdout.
    """
    req_file = ANALYSIS_DIR / "requirements.txt"
    if not req_file.exists():
        print("  [WARN] requirements.txt not found — skipping requirement check")
        return True

    missing: list[str] = []
    import importlib.metadata as _meta

    for line in req_file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Strip version specifiers to get the package name.
        pkg_name = line.split(">=")[0].split("<=")[0].split("==")[0].split("!=")[0].split("<")[0].split(">")[0].strip()
        try:
            _meta.version(pkg_name)
        except _meta.PackageNotFoundError:
            missing.append(line)

    if missing:
        print("  Missing Python packages:")
        for m in missing:
            print(f"    - {m}")
        print(f"\n  Install with:  pip install -r {req_file}")
        return False

    print("  All requirements satisfied.")
    return True


def _run_tests() -> bool:
    """Run the analysis test suite via pytest.

    Returns True if all tests pass, False otherwise.
    """
    test_dir = ANALYSIS_DIR / "tests"
    if not test_dir.is_dir():
        print("  [WARN] tests/ directory not found — skipping test run")
        return True

    print("  Running analysis test suite (pytest) ...")
    t0 = time.time()
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-v", str(test_dir)],
        cwd=str(ANALYSIS_DIR),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    elapsed = time.time() - t0

    # Show the last portion of output (summary) to keep output manageable.
    if result.stdout:
        # Print full output so failures are visible.
        sys.stdout.write(result.stdout)
    if result.stderr:
        sys.stderr.write(result.stderr)

    if result.returncode == 0:
        print(f"  All tests passed ({elapsed:.1f}s)")
        return True
    else:
        print(f"  Tests FAILED (exit code {result.returncode}, {elapsed:.1f}s)")
        return False


def _run_script(name: str, folder: Path, log_file: io.TextIOBase | None = None) -> bool:
    """Run a sibling Python script as a subprocess.

    Parameters
    ----------
    name : str
        Filename of the script to run (e.g. ``"parse_log_metrics.py"``).
        Must reside in the same directory as this orchestrator.
    folder : Path
        Path to the ``saved_files_*`` output folder, passed as the first
        positional argument to the child script.
    log_file : io.TextIOBase or None
        If provided, subprocess stdout and stderr are captured and written
        to both the terminal and this log file.

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
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    elapsed = time.time() - t0

    # Emit captured output to terminal and log file.
    for stream_output in (result.stdout, result.stderr):
        if stream_output:
            sys.stdout.write(stream_output)
            if log_file is not None:
                log_file.write(stream_output)

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
        help="Only generate the report from existing data (skip parsing steps)",
    )
    parser.add_argument(
        "--no-pdf",
        action="store_true",
        help="Skip PDF generation (produce HTML report only)",
    )
    parser.add_argument(
        "--html",
        action="store_true",
        help="Also write the HTML report to disk (default: PDF only)",
    )
    parser.add_argument(
        "--gemini-model",
        type=str,
        default=None,
        help="Override the Gemini vision model (e.g. gemini-2.0-flash)",
    )
    parser.add_argument(
        "--concurrency",
        type=int,
        default=None,
        help="Max concurrent Gemini API requests (default: from config)",
    )
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to analysis_config.json (default: auto-detect)",
    )
    parser.add_argument(
        "--skip-checks",
        action="store_true",
        help="Skip pre-flight checks (requirements verification and test suite)",
    )
    args = parser.parse_args()

    # ── Apply CLI overrides to the centralised config ──
    # Load the config (from file or defaults), then patch with CLI args.
    import os as _os
    import shared as _shared_mod
    cfg = load_analysis_config(config_path=args.config)
    if args.gemini_model:
        cfg["vision"]["gemini_model"] = args.gemini_model
        # Propagate to subprocess children via environment variable.
        _os.environ["PANCDATA3_GEMINI_MODEL"] = args.gemini_model
    if args.concurrency is not None:
        cfg["vision"]["max_concurrent_requests"] = args.concurrency
        _os.environ["PANCDATA3_GEMINI_CONCURRENCY"] = str(args.concurrency)
    # Install the patched config into the shared module cache so that
    # child scripts (imported as modules) see the CLI overrides.
    _shared_mod._config_cache = cfg

    # ── Resolve output folder ──
    if args.folder:
        folder = Path(args.folder)
        if not folder.is_dir():
            sys.exit(f"ERROR: Folder does not exist: {folder}")
    else:
        # Auto-detect the most recent saved_files_* directory.
        folder = find_latest_saved_folder()

    # ── Open log file and tee all output ──
    log_path = folder / "run_analysis_output.log"
    log_file = open(log_path, "w", encoding="utf-8")
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = TeeWriter(original_stdout, log_file)
    sys.stderr = TeeWriter(original_stderr, log_file)

    try:
        # ── Print banner ──
        print("=" * 70)
        print("  pancData3 Analysis Suite")
        print("=" * 70)
        print(f"  Output folder: {folder}")
        print(f"  Skip vision:   {args.skip_vision}")
        print(f"  Report only:   {args.report_only}")
        print(f"  Skip checks:   {args.skip_checks}")
        print(f"  Gemini model:  {cfg['vision']['gemini_model']}")
        print(f"  Concurrency:   {cfg['vision']['max_concurrent_requests']}")
        print(f"  Log file:      {log_path}")
        print()

        # ── Pre-flight checks ──
        if not args.skip_checks:
            print("── Pre-flight checks ──")
            reqs_ok = _check_requirements()
            if not reqs_ok:
                sys.exit("Pre-flight check failed: missing requirements. "
                         "Install them or re-run with --skip-checks.")
            tests_ok = _run_tests()
            if not tests_ok:
                sys.exit("Pre-flight check failed: test suite did not pass. "
                         "Fix failing tests or re-run with --skip-checks.")
            print("── Pre-flight checks passed ──\n")
        else:
            print("  Skipping pre-flight checks (--skip-checks)\n")

        # Track success/failure/skip status for each pipeline step.
        results = {}

        # Determine total number of steps for the progress bar.
        total_steps = 1  # report generation always runs
        if not args.report_only:
            total_steps += 3  # logs, csvs, mat
            if not args.skip_vision:
                total_steps += 1  # vision

        pipeline_bar = tqdm(
            total=total_steps,
            desc="Analysis pipeline",
            unit="step",
            bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} steps [{elapsed}<{remaining}] {postfix}",
        )

        if not args.report_only:
            # Step 1: Vision-based graph analysis (requires GEMINI_API_KEY).
            if not args.skip_vision:
                csv_path = folder / "graph_analysis_results.csv"
                print(f"  Vision CSV exists: {csv_path.exists()}")
                pipeline_bar.set_postfix_str("vision analysis", refresh=True)
                results["vision"] = _run_script("batch_graph_analysis.py", folder, log_file)
                pipeline_bar.update(1)
            else:
                csv_path = folder / "graph_analysis_results.csv"
                if csv_path.exists():
                    print(f"  Skipping vision analysis (existing CSV: {csv_path})")
                else:
                    print("  Skipping vision analysis (no existing CSV)")
                results["vision"] = "skipped"

            # Step 2: Parse MATLAB diary log files for structured metrics.
            pipeline_bar.set_postfix_str("parsing logs", refresh=True)
            results["logs"] = _run_script("parse_log_metrics.py", folder, log_file)
            pipeline_bar.update(1)

            # Step 3: Parse pipeline-exported CSV files (significance tables).
            pipeline_bar.set_postfix_str("parsing CSVs", refresh=True)
            results["csvs"] = _run_script("parse_csv_results.py", folder, log_file)
            pipeline_bar.update(1)

            # Step 3.5: Parse MATLAB .mat files (core comparison, dosimetry).
            pipeline_bar.set_postfix_str("parsing MAT files", refresh=True)
            results["mat"] = _run_script("parse_mat_metrics.py", folder, log_file)
            pipeline_bar.update(1)

        # Step 4: Assemble the final HTML (+PDF) report from all collected data.
        report_script = ANALYSIS_DIR / "generate_report.py"
        if report_script.exists():
            report_args = [sys.executable, str(report_script), str(folder)]
            if args.no_pdf:
                report_args.append("--no-pdf")
            if not args.html:
                report_args.append("--no-html")
            pipeline_bar.set_postfix_str("generating report", refresh=True)
            print(f"\n  Running generate_report.py ...")
            t0 = time.time()
            result = subprocess.run(
                report_args, cwd=str(ANALYSIS_DIR),
                capture_output=True, text=True,
                encoding="utf-8", errors="replace",
            )
            elapsed = time.time() - t0
            for stream_output in (result.stdout, result.stderr):
                if stream_output:
                    sys.stdout.write(stream_output)
            if result.returncode == 0:
                print(f"  Done: generate_report.py ({elapsed:.1f}s)")
                results["report"] = True
            else:
                print(f"  FAILED: generate_report.py (exit code {result.returncode}, {elapsed:.1f}s)")
                results["report"] = False
        else:
            results["report"] = False

        pipeline_bar.update(1)
        pipeline_bar.set_postfix_str("complete", refresh=True)
        pipeline_bar.close()

        # ── Summary table ──
        print("\n" + "=" * 70)
        print("  Analysis Complete")
        print("=" * 70)
        for step, status in results.items():
            icon = "OK" if status is True else ("SKIP" if status == "skipped" else "FAIL")
            print(f"  [{icon:>4}] {step}")

        report_path = folder / "analysis_report.html"
        pdf_path = folder / "analysis_report.pdf"
        if report_path.exists() and args.html:
            print(f"\n  HTML Report: {report_path}")
        if pdf_path.exists():
            print(f"  PDF Report:  {pdf_path}")
        print(f"  Log file:    {log_path}")
        print()

    finally:
        # Restore original streams and close the log file.
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        log_file.close()


if __name__ == "__main__":
    main()
