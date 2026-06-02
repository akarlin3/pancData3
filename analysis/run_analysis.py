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
    python run_analysis.py --skip-vision       # skip vision API calls
    python run_analysis.py --provider both     # run both Gemini & Claude
    python run_analysis.py --report-only       # only generate report
    python run_analysis.py --skip-checks       # skip pre-flight checks

Subprocess and output helpers (TeeWriter, _run_script, requirement/test
checks, API-key masking, Windows path handling) live in the sibling module
:mod:`run_analysis_helpers` and are re-exported here for backward
compatibility.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

from tqdm import tqdm  # type: ignore

from shared import find_latest_saved_folder, get_api_key, get_config, load_analysis_config, reset_config_cache, setup_utf8_stdout  # type: ignore

# Re-exported from run_analysis_helpers (kept importable from run_analysis for
# backward compatibility — tests import several of these directly).
from run_analysis_helpers import (  # noqa: F401  (re-export)
    ANALYSIS_DIR,
    TeeWriter,
    _check_requirements,
    _handle_windows_path,
    _mask_api_keys,
    _mask_key_value,
    _run_script,
    _run_tests,
)

# Ensure emoji and special characters print correctly on Windows consoles.
setup_utf8_stdout()


def _ensure_config_key(provider: str, display_name: str, config_key: str,
                       get_url: str) -> None:
    """Ensure an API key is available via config or environment.

    Resolution order (via ``get_api_key``):
    1. ``vision.<config_key>`` in ``analysis_config.json``.
    2. The corresponding environment variable.

    If still missing, prompts interactively and saves to
    ``analysis_config.json`` for future runs.
    """
    if get_api_key(provider):
        return

    print("\n" + "!" * 70)
    print(f"  WARNING: {display_name} API key is not set.")
    print(f"  This key is required for {display_name} vision analysis.")
    print(f"  You can get an API key from: {get_url}")
    print(f"  Set vision.{config_key} in analysis_config.json")
    print("!" * 70)

    import getpass
    key = getpass.getpass(f"\n  Enter your {display_name} API key (hidden): ").strip()

    if not key:
        print(f"  [WARN] No API key provided. {display_name} vision analysis will fail if not skipped.\n")
        return

    # Inject into running config so downstream code sees it immediately
    os.environ[{"gemini": "GEMINI_API_KEY", "claude": "ANTHROPIC_API_KEY"}[provider]] = key

    # Persist to analysis_config.json for future runs
    config_path = ANALYSIS_DIR.parent / "analysis_config.json"
    try:
        if config_path.is_file():
            with open(config_path, "r", encoding="utf-8") as f:
                cfg_data = json.load(f)
        else:
            cfg_data = {}
        cfg_data.setdefault("vision", {})[config_key] = key
        with open(config_path, "w", encoding="utf-8", newline="") as f:
            json.dump(cfg_data, f, indent=4)
            f.write("\n")
        print(f"  [INFO] Saved API key to {config_path}\n")
    except Exception as e:
        print(f"  [WARN] Failed to save key to config: {e}\n")


def _ensure_api_keys(provider: str) -> None:
    """Ensure all required API keys are available for the given provider.

    Parameters
    ----------
    provider : str
        One of ``"gemini"``, ``"claude"``, or ``"both"``.
    """
    if provider in ("gemini", "both"):
        _ensure_config_key(
            "gemini", "Gemini", "gemini_api_key",
            "https://aistudio.google.com/")
    if provider in ("claude", "both"):
        _ensure_config_key(
            "claude", "Claude (Anthropic)", "anthropic_api_key",
            "https://console.anthropic.com/")


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
        help="Skip vision API analysis (use existing CSV if available)",
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
        "--provider",
        type=str,
        default=None,
        choices=["gemini", "claude", "both"],
        help="Vision API provider: gemini (default), claude, or both (runs both and compares)",
    )
    parser.add_argument(
        "--gemini-model",
        type=str,
        default=None,
        help="Override the Gemini vision model (e.g. gemini-2.0-flash)",
    )
    parser.add_argument(
        "--claude-model",
        type=str,
        default=None,
        help="Override the Claude vision model (e.g. claude-opus-4-6)",
    )
    parser.add_argument(
        "--concurrency",
        type=int,
        default=None,
        help="Max concurrent API requests (default: from config)",
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
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Also generate an interactive HTML report with client-side filtering",
    )
    args = parser.parse_args()

    # ── Apply CLI overrides to the centralised config ──
    # Load the config (from file or defaults), then patch with CLI args.
    import os as _os
    import shared as _shared_mod  # type: ignore
    cfg = load_analysis_config(config_path=args.config)
    if args.provider:
        cfg["vision"]["provider"] = args.provider
        _os.environ["PANCDATA3_VISION_PROVIDER"] = args.provider
    if args.gemini_model:
        cfg["vision"]["gemini_model"] = args.gemini_model
        _os.environ["PANCDATA3_GEMINI_MODEL"] = args.gemini_model
    if args.claude_model:
        cfg["vision"]["claude_model"] = args.claude_model
        _os.environ["PANCDATA3_CLAUDE_MODEL"] = args.claude_model
    if args.concurrency is not None:
        cfg["vision"]["max_concurrent_requests"] = args.concurrency
        _os.environ["PANCDATA3_GEMINI_CONCURRENCY"] = str(args.concurrency)
    # Install the patched config into the shared module cache so that
    # child scripts (imported as modules) see the CLI overrides.
    _shared_mod._config_cache = cfg

    # ── Resolve output folder ──
    if args.folder:
        folder = Path(args.folder)
        # Apply Windows path handling
        folder = _handle_windows_path(folder)
        if not folder.is_dir():
            sys.exit(f"ERROR: Folder does not exist: {folder}")
    else:
        # Auto-detect the most recent saved_files_* directory.
        folder = find_latest_saved_folder()
        # Apply Windows path handling
        folder = _handle_windows_path(folder)

    # ── Open log file and tee all output ──
    log_path = folder / "run_analysis_output.log"
    log_file = open(log_path, "w", encoding="utf-8")
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = TeeWriter(original_stdout, log_file)  # type: ignore
    sys.stderr = TeeWriter(original_stderr, log_file)  # type: ignore

    try:
        # ── Print banner ──
        print("=" * 70)
        print("  pancData3 Analysis Suite")
        print("=" * 70)
        print(f"  Output folder: {folder}")
        print(f"  Skip vision:   {args.skip_vision}")
        print(f"  Report only:   {args.report_only}")
        print(f"  Skip checks:   {args.skip_checks}")
        print(f"  Provider:      {cfg['vision'].get('provider', 'gemini')}")
        print(f"  Gemini model:  {cfg['vision']['gemini_model']}")
        print(f"  Claude model:  {cfg['vision'].get('claude_model', 'claude-opus-4-6')}")
        print(f"  Concurrency:   {cfg['vision']['max_concurrent_requests']}")
        print(f"  Log file:      {log_path}")
        print()

        # ── Setup Environment ──
        vision_provider = cfg["vision"].get("provider", "gemini")
        if not args.skip_vision:
            _ensure_api_keys(vision_provider)

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
                vision_timeout = cfg["vision"].get("script_timeout_seconds")
                results["vision"] = _run_script(
                    "parsers/batch_graph_analysis.py", folder, log_file,
                    timeout=vision_timeout,
                )
                pipeline_bar.update(1)
            else:
                csv_path = folder / "graph_analysis_results.csv"
                if csv_path.exists():
                    print(f"  Skipping vision analysis (existing CSV: {csv_path})")
                else:
                    print("  Skipping vision analysis (no existing CSV)")
                results["vision"] = "skipped"  # type: ignore

            # Step 2: Parse MATLAB diary log files for structured metrics.
            pipeline_bar.set_postfix_str("parsing logs", refresh=True)
            results["logs"] = _run_script("parsers/parse_log_metrics.py", folder, log_file)
            pipeline_bar.update(1)

            # Step 3: Parse pipeline-exported CSV files (significance tables).
            pipeline_bar.set_postfix_str("parsing CSVs", refresh=True)
            results["csvs"] = _run_script("parsers/parse_csv_results.py", folder, log_file)
            pipeline_bar.update(1)

            # Step 3.5: Parse MATLAB .mat files (core comparison, dosimetry).
            pipeline_bar.set_postfix_str("parsing MAT files", refresh=True)
            results["mat"] = _run_script("parsers/parse_mat_metrics.py", folder, log_file)
            pipeline_bar.update(1)

            # Step 3.75: Cross-DWI agreement analysis (Bland-Altman, CCC, ICC).
            pipeline_bar.set_postfix_str("cross-DWI agreement", refresh=True)
            results["agreement"] = _run_script(
                "cross_reference/cross_dwi_agreement.py", folder, log_file
            )
            pipeline_bar.update(1)

        # Step 4: Assemble the final HTML (+PDF) report from all collected data.
        report_script = ANALYSIS_DIR / "report" / "generate_report.py"
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

        # Step 5 (optional): Generate interactive report.
        if args.interactive:
            interactive_script = ANALYSIS_DIR / "report" / "generate_interactive_report.py"
            if interactive_script.exists():
                pipeline_bar.total += 1
                pipeline_bar.refresh()
                pipeline_bar.set_postfix_str("interactive report", refresh=True)
                print(f"\n  Running generate_interactive_report.py ...")
                t0 = time.time()
                result = subprocess.run(
                    [sys.executable, str(interactive_script), str(folder)],
                    cwd=str(ANALYSIS_DIR),
                    capture_output=True, text=True,
                    encoding="utf-8", errors="replace",
                )
                elapsed = time.time() - t0
                for stream_output in (result.stdout, result.stderr):
                    if stream_output:
                        sys.stdout.write(stream_output)
                if result.returncode == 0:
                    print(f"  Done: generate_interactive_report.py ({elapsed:.1f}s)")
                    results["interactive"] = True
                else:
                    print(f"  FAILED: generate_interactive_report.py (exit code {result.returncode}, {elapsed:.1f}s)")
                    results["interactive"] = False
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

        report_path = folder / "analysis_report.html"  # type: ignore
        pdf_path = folder / "analysis_report.pdf"  # type: ignore
        if report_path.exists() and args.html:
            print(f"\n  HTML Report: {report_path}")
        if pdf_path.exists():
            print(f"  PDF Report:  {pdf_path}")
        interactive_path = folder / "interactive_report.html"  # type: ignore
        if interactive_path.exists():
            print(f"  Interactive: {interactive_path}")
        print(f"  Log file:    {log_path}")
        print()

    finally:
        # Restore original streams and close the log file.
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        log_file.close()


if __name__ == "__main__":
    main()
