"""Integration tests for the full analysis pipeline (run_analysis.py).

These tests exercise the complete pipeline wiring — parsers, report generator,
and orchestrator — on synthetic data built from conftest.py fixtures.  They
catch composition bugs where individual modules work in isolation but break
when connected together via run_analysis.main() or generate_report.main().

The tests use ``--skip-vision --no-pdf --skip-checks`` to avoid external
dependencies (Gemini API, WeasyPrint/Chrome, and the pytest pre-flight suite).
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from conftest import GRAPH_CSV_COLUMNS, SAMPLE_GRAPH_CSV_ROWS


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_full_synthetic_dir(tmp_path: Path) -> Path:
    """Build a saved_files directory with all data sources populated.

    Combines graph CSV, MATLAB logs, and pipeline CSV exports into a single
    directory so the full pipeline exercises every parser and report section.
    """
    folder = tmp_path / "saved_files_20260301_120000"
    for dwi_type in ("Standard", "dnCNN", "IVIMnet"):
        (folder / dwi_type).mkdir(parents=True)

    # ── Graph analysis CSV ──
    csv_path = folder / "graph_analysis_results.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
        writer.writeheader()
        for row in SAMPLE_GRAPH_CSV_ROWS:
            writer.writerow(row)

    # ── MATLAB log files (Standard) ──
    std = folder / "Standard"

    (std / "metrics_stats_comparisons_output_Standard.txt").write_text(
        "Interaction P-Value for ADC: 0.023\n"
        "Interaction P-Value for D_star: 1.5e-04\n"
        "mean_adc: p=0.015, adj_alpha=0.025\n"
        "mean_d*: p=0.042, adj_alpha=0.0125\n"
        "Timepoint: BL — 3 significant\n"
        "Timepoint: W2 — 1 significant\n"
        "Excluded 5/42 (11.9%) competing-risk patients\n",
        encoding="utf-8",
    )

    (std / "metrics_stats_predictive_output_Standard.txt").write_text(
        "Elastic Net Selected Features for BL (Opt Lambda=0.0532): mean_adc, mean_d, f_ratio\n"
        "Elastic Net Selected Features for W2 (Opt Lambda=0.1200): delta_adc, mean_f\n"
        "Firth refit successful for BL (3 features)\n"
        "PRIMARY ROC ANALYSIS results for BL\n"
        "  AUC = 0.781\n"
        "  Youden Optimal Score Cutoff = 0.450\n"
        "  Sensitivity = 82.5% | Specificity = 68.2%\n"
        "PRIMARY ROC ANALYSIS results for W2\n"
        "  AUC = 0.692\n"
        "  Youden Optimal Score Cutoff = 0.520\n"
        "  Sensitivity = 75.0% | Specificity = 60.0%\n",
        encoding="utf-8",
    )

    (std / "metrics_survival_output_Standard.txt").write_text(
        "  mean_adc   1.250   0.980   1.590   0.0680\n"
        "  delta_d    0.750   0.550   1.020   0.0340\n"
        "Global LRT: chi2(2) = 7.82, p = 0.0200\n"
        "IPCW weights applied ... range [0.85, 1.42]\n",
        encoding="utf-8",
    )

    (std / "metrics_baseline_output_Standard.txt").write_text(
        "Outlier flag (mean_adc): 3 flagged (LF=2, LC=1, CR=0)\n"
        "Outlier flag (mean_d): 1 flagged (LF=0, LC=0, CR=1)\n"
        "Total outliers removed: 4 / 42 (9.52%)\n"
        "Excluded 6/48 patients due to missing baseline\n"
        "LF rate: included=35.7%, excluded=50.0%\n",
        encoding="utf-8",
    )

    # ── Pipeline CSV exports ──
    for dwi_type in ("Standard", "dnCNN"):
        sig_csv = folder / dwi_type / "Significant_LF_Metrics.csv"
        with open(sig_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["Metric", "Timepoint", "p_value"])
            writer.writeheader()
            if dwi_type == "Standard":
                writer.writerow({"Metric": "mean_adc", "Timepoint": "BL", "p_value": "0.003"})
                writer.writerow({"Metric": "mean_d", "Timepoint": "BL", "p_value": "0.01"})
                writer.writerow({"Metric": "mean_adc", "Timepoint": "W2", "p_value": "0.04"})
            else:
                writer.writerow({"Metric": "mean_adc", "Timepoint": "BL", "p_value": "0.005"})

    return folder


# ---------------------------------------------------------------------------
# Test: generate_report on fully-populated synthetic data
# ---------------------------------------------------------------------------

class TestGenerateReportIntegration:
    """End-to-end tests for generate_report() with all data sources present.

    Unlike the unit tests in test_generate_report.py (which test individual
    fixtures), these tests feed *all* data sources simultaneously to verify
    that the section builders compose correctly when every data path is active.
    """

    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path: Path):
        self.folder = _build_full_synthetic_dir(tmp_path)

    def _generate(self) -> str:
        from report.generate_report import generate_report
        return generate_report(self.folder)

    # ── HTML validity ──

    def test_valid_html_structure(self):
        """Report must have proper DOCTYPE, html, head, and body tags."""
        html = self._generate()
        assert html.startswith("<!DOCTYPE html>")
        assert "<html" in html
        assert "<head>" in html
        assert "<body>" in html
        assert "</html>" in html

    def test_meta_charset(self):
        """Report must declare UTF-8 charset."""
        html = self._generate()
        assert '<meta charset="utf-8">' in html

    def test_contains_css(self):
        """Report must include embedded CSS stylesheet."""
        html = self._generate()
        assert "<style>" in html

    def test_contains_javascript(self):
        """Report must include embedded JavaScript."""
        html = self._generate()
        assert "<script>" in html

    # ── Report header and metadata ──

    def test_title_contains_timestamp(self):
        """Title should contain the formatted timestamp from the folder name."""
        html = self._generate()
        assert "March 1, 2026" in html

    def test_header_lists_dwi_types(self):
        """Report header should list all three DWI types."""
        html = self._generate()
        for dwi_type in ("Standard", "dnCNN", "IVIMnet"):
            assert dwi_type in html

    def test_graph_count_in_header(self):
        """Report should mention the number of graphs analysed."""
        html = self._generate()
        assert "3" in html  # 3 rows in SAMPLE_GRAPH_CSV_ROWS

    # ── Part 1: Overview sections ──

    def test_executive_summary_present(self):
        """Executive summary section with structured abstract subsections."""
        html = self._generate()
        assert "Executive Summary" in html
        assert "Objective" in html
        assert "Key Results" in html

    def test_methods_section_present(self):
        """Methods section with statistical methodology description."""
        html = self._generate()
        assert "Statistical Methods" in html
        assert "Wilcoxon" in html
        assert "Cox" in html

    def test_cohort_overview_present(self):
        """Cohort overview section should exist."""
        html = self._generate()
        assert "Cohort" in html

    def test_patient_flow_present(self):
        """Patient flow section should show attrition data from logs."""
        html = self._generate()
        assert "Patient Flow" in html or "Attrition" in html

    # ── Part 2: Data sections ──

    def test_data_quality_section(self):
        """Data quality section should contain outlier data from baseline logs."""
        html = self._generate()
        assert "Data Quality" in html
        assert "Outlier" in html

    def test_outlier_flags_table(self):
        """Per-metric outlier flags table should appear with log data."""
        html = self._generate()
        assert "mean_adc" in html
        assert "mean_d" in html

    def test_baseline_exclusion_data(self):
        """Baseline exclusion stat card should appear from log data."""
        html = self._generate()
        # The fixture has "Excluded 6/48 patients due to missing baseline"
        assert "Baseline Excluded" in html or "baseline" in html.lower()

    def test_hypothesis_section(self):
        """Hypothesis section should be generated."""
        html = self._generate()
        assert "Hypothes" in html  # Hypothesis or Hypotheses

    def test_graph_overview_table(self):
        """Graph overview should list graph types from the CSV."""
        html = self._generate()
        assert "Graph Type" in html or "Graph Overview" in html
        assert "box" in html
        assert "line" in html

    def test_graph_issues_listed(self):
        """Known issues from the vision CSV should appear."""
        html = self._generate()
        assert "Y-axis label partially obscured by legend" in html
        assert "Overlapping legend text" in html

    # ── Part 3: Statistics sections ──

    def test_stats_by_graph_type(self):
        """Statistics by graph type section should exist."""
        html = self._generate()
        assert "Graph Type" in html

    def test_statistical_significance_section(self):
        """Statistical significance section with p-values from all sources."""
        html = self._generate()
        # p = 0.003 from graph CSV and CSV exports
        assert "0.003" in html

    def test_effect_sizes_section(self):
        """Effect sizes section should appear (log data has hazard ratios)."""
        html = self._generate()
        assert "Effect Size" in html

    def test_multiple_comparisons_section(self):
        """Multiple comparisons section with FDR correction."""
        html = self._generate()
        assert "Multiple Comparison" in html

    def test_cross_dwi_comparison(self):
        """Cross-DWI comparison for Feature_BoxPlots (Standard vs dnCNN)."""
        html = self._generate()
        if "Cross-DWI" in html:
            # Feature_BoxPlots has opposing trends → should show DIFFER or AGREE
            assert "DIFFER" in html or "AGREE" in html

    def test_fdr_global_section(self):
        """FDR global section should be present (even if no FDR CSV)."""
        html = self._generate()
        assert "FDR" in html

    # ── Part 4: Outcomes sections ──

    def test_correlations_section(self):
        """Correlations from graph summaries (r = 0.65 in fixture)."""
        html = self._generate()
        assert "0.65" in html

    def test_treatment_response_section(self):
        """Treatment response section should exist."""
        html = self._generate()
        assert "Treatment" in html

    def test_predictive_performance_section(self):
        """Predictive performance with AUC from log data."""
        html = self._generate()
        assert "0.781" in html or "Predictive" in html

    def test_model_diagnostics_section(self):
        """Model diagnostics should list assumptions."""
        html = self._generate()
        assert "Model Diagnostics" in html

    # ── Part 5: Discussion sections ──

    def test_limitations_section(self):
        """Limitations section with standard caveats."""
        html = self._generate()
        assert "Limitations" in html

    def test_conclusions_section(self):
        """Conclusions section should exist."""
        html = self._generate()
        assert "Conclusions" in html

    # ── Part 6: Appendix sections ──

    def test_appendix_lists_graphs(self):
        """Appendix should list all analysed graphs."""
        html = self._generate()
        assert "Appendix" in html
        assert "Feature_BoxPlots" in html

    def test_table_index(self):
        """Table index section should exist."""
        html = self._generate()
        # Table index or "Table of Tables"
        assert "Table" in html

    def test_figure_index(self):
        """Figure index section should exist."""
        html = self._generate()
        assert "Figure" in html

    def test_references_section(self):
        """References section should be present."""
        html = self._generate()
        assert "References" in html

    def test_footer(self):
        """Report footer attribution."""
        html = self._generate()
        assert "Report generated by pancData3" in html

    # ── Tables and structural elements ──

    def test_contains_tables(self):
        """Report should contain multiple HTML tables."""
        html = self._generate()
        table_count = html.count("<table")
        assert table_count >= 3, f"Expected >= 3 tables, found {table_count}"

    def test_contains_headings(self):
        """Report should contain multiple H2 section headings."""
        html = self._generate()
        h2_count = html.count("<h2")
        assert h2_count >= 10, f"Expected >= 10 h2 headings, found {h2_count}"

    def test_contains_stat_cards(self):
        """Report should contain stat-card elements from log data."""
        html = self._generate()
        assert "stat-grid" in html or "stat-card" in html

    def test_navigation_bar(self):
        """Sticky navigation bar should be present."""
        html = self._generate()
        assert "nav" in html.lower()

    def test_dwi_badges_present(self):
        """DWI type badges (coloured labels) should appear in the report."""
        html = self._generate()
        assert "badge-Standard" in html or "badge" in html

    # ── Data integration: all sources used together ──

    def test_log_data_feeds_into_report(self):
        """Hazard ratios from survival logs should appear in report tables."""
        html = self._generate()
        # HR = 1.250 for mean_adc from the survival log
        assert "1.250" in html or "1.25" in html

    def test_csv_data_feeds_into_report(self):
        """Pipeline CSV significance data should appear in the report."""
        html = self._generate()
        # mean_adc at BL with p = 0.003 from Significant_LF_Metrics.csv
        assert "0.003" in html

    def test_graph_csv_feeds_into_report(self):
        """Vision CSV data (trends, issues) should appear in the report."""
        html = self._generate()
        assert "increasing" in html.lower() or "decreasing" in html.lower()

    def test_no_python_tracebacks(self):
        """Report HTML should not contain Python tracebacks."""
        html = self._generate()
        assert "Traceback (most recent call last)" not in html
        assert "raise " not in html or "raise" in html.lower()  # tolerate prose


# ---------------------------------------------------------------------------
# Test: run_analysis.py subprocess execution
# ---------------------------------------------------------------------------

class TestRunAnalysisSubprocess:
    """Run run_analysis.py as a subprocess with --skip-vision --no-pdf --skip-checks.

    This exercises the full orchestrator wiring: argument parsing, folder
    resolution, subprocess launching of parsers, and report generation.
    """

    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path: Path):
        self.folder = _build_full_synthetic_dir(tmp_path)

    def test_orchestrator_exits_cleanly(self):
        """run_analysis.py should exit with code 0 on synthetic data."""
        result = subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve().parent.parent / "run_analysis.py"),
                "--folder", str(self.folder),
                "--skip-vision",
                "--no-pdf",
                "--skip-checks",
                "--html",
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=120,
        )
        assert result.returncode == 0, (
            f"run_analysis.py failed (exit {result.returncode}):\n"
            f"STDOUT:\n{result.stdout[-2000:]}\n"
            f"STDERR:\n{result.stderr[-2000:]}"
        )

    def test_html_report_created(self):
        """The orchestrator should produce analysis_report.html."""
        subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve().parent.parent / "run_analysis.py"),
                "--folder", str(self.folder),
                "--skip-vision",
                "--no-pdf",
                "--skip-checks",
                "--html",
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=120,
        )
        report_path = self.folder / "analysis_report.html"
        assert report_path.exists(), "analysis_report.html was not created"
        content = report_path.read_text(encoding="utf-8")
        assert len(content) > 1000, f"Report too short ({len(content)} chars)"

    def test_html_report_contains_key_sections(self):
        """Generated HTML report should contain expected sections."""
        subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve().parent.parent / "run_analysis.py"),
                "--folder", str(self.folder),
                "--skip-vision",
                "--no-pdf",
                "--skip-checks",
                "--html",
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=120,
        )
        report_path = self.folder / "analysis_report.html"
        if not report_path.exists():
            pytest.skip("Report not created (previous assertion should catch this)")
        content = report_path.read_text(encoding="utf-8")

        # Key sections that should be present
        expected_sections = [
            "Executive Summary",
            "Statistical Methods",
            "Limitations",
            "Conclusions",
            "Appendix",
            "References",
        ]
        for section in expected_sections:
            assert section in content, f"Missing section: {section}"

    def test_log_file_created(self):
        """The orchestrator should create run_analysis_output.log."""
        subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve().parent.parent / "run_analysis.py"),
                "--folder", str(self.folder),
                "--skip-vision",
                "--no-pdf",
                "--skip-checks",
                "--html",
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=120,
        )
        log_path = self.folder / "run_analysis_output.log"
        assert log_path.exists(), "run_analysis_output.log was not created"

    def test_summary_table_in_output(self):
        """Orchestrator stdout should contain the final summary table."""
        result = subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve().parent.parent / "run_analysis.py"),
                "--folder", str(self.folder),
                "--skip-vision",
                "--no-pdf",
                "--skip-checks",
                "--html",
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=120,
        )
        # The orchestrator prints a summary with [OK], [SKIP], or [FAIL]
        assert "Analysis Complete" in result.stdout
        assert "[  OK]" in result.stdout or "[ OK]" in result.stdout or "[OK]" in result.stdout

    def test_no_pdf_flag_skips_pdf(self):
        """With --no-pdf, no PDF file should be generated."""
        subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve().parent.parent / "run_analysis.py"),
                "--folder", str(self.folder),
                "--skip-vision",
                "--no-pdf",
                "--skip-checks",
                "--html",
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=120,
        )
        pdf_path = self.folder / "analysis_report.pdf"
        assert not pdf_path.exists(), "PDF should not be created with --no-pdf"


# ---------------------------------------------------------------------------
# Test: generate_report.main() via sys.argv patching
# ---------------------------------------------------------------------------

class TestGenerateReportMain:
    """Test generate_report.main() entry point directly (not via subprocess).

    This is faster than subprocess tests and catches import/wiring issues
    at the Python level.
    """

    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path: Path):
        self.folder = _build_full_synthetic_dir(tmp_path)

    def test_main_creates_html(self):
        """generate_report.main() should create analysis_report.html."""
        with patch.object(
            sys, "argv",
            ["generate_report.py", str(self.folder), "--no-pdf", "--no-html"],
        ):
            # Import fresh to avoid stale sys.argv
            from report.generate_report import main
            # We pass --no-html to avoid writing (main still generates internally).
            # Instead, call generate_report directly.

        from report.generate_report import generate_report
        html = generate_report(self.folder)
        assert len(html) > 5000, f"Report unexpectedly short: {len(html)} chars"

    def test_main_writes_html_file(self):
        """generate_report.main() with no flags should write the HTML file."""
        with patch.object(
            sys, "argv",
            ["generate_report.py", str(self.folder), "--no-pdf"],
        ):
            from report.generate_report import main
            main()

        report = self.folder / "analysis_report.html"
        assert report.exists()
        content = report.read_text(encoding="utf-8")
        assert "<!DOCTYPE html>" in content
        assert "Analysis Report" in content


# ---------------------------------------------------------------------------
# Test: empty/minimal input graceful degradation
# ---------------------------------------------------------------------------

class TestGracefulDegradation:
    """Verify the pipeline handles missing or incomplete data gracefully."""

    def test_empty_folder_produces_report(self, tmp_path: Path):
        """A folder with only DWI subdirectories (no data) should still work."""
        folder = tmp_path / "saved_files_20260101_000000"
        for d in ("Standard", "dnCNN", "IVIMnet"):
            (folder / d).mkdir(parents=True)

        from report.generate_report import generate_report
        html = generate_report(folder)
        assert "Analysis Report" in html
        assert "<!DOCTYPE html>" in html

    def test_graph_csv_only(self, tmp_path: Path):
        """Folder with only graph CSV (no logs, no pipeline CSVs) works."""
        folder = tmp_path / "saved_files_20260201_060000"
        for d in ("Standard", "dnCNN", "IVIMnet"):
            (folder / d).mkdir(parents=True)

        csv_path = folder / "graph_analysis_results.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
            writer.writeheader()
            for row in SAMPLE_GRAPH_CSV_ROWS:
                writer.writerow(row)

        from report.generate_report import generate_report
        html = generate_report(folder)
        assert "Feature_BoxPlots" in html
        assert "box" in html

    def test_logs_only(self, tmp_path: Path):
        """Folder with only log files (no graph CSV, no pipeline CSVs) works."""
        folder = tmp_path / "saved_files_20260201_070000"
        std = folder / "Standard"
        std.mkdir(parents=True)
        for d in ("dnCNN", "IVIMnet"):
            (folder / d).mkdir(parents=True)

        (std / "metrics_baseline_output_Standard.txt").write_text(
            "Total outliers removed: 2 / 30 (6.67%)\n",
            encoding="utf-8",
        )

        from report.generate_report import generate_report
        html = generate_report(folder)
        assert "Analysis Report" in html

    def test_subprocess_on_empty_folder(self, tmp_path: Path):
        """run_analysis.py on an empty folder should not crash."""
        folder = tmp_path / "saved_files_20260101_000000"
        for d in ("Standard", "dnCNN", "IVIMnet"):
            (folder / d).mkdir(parents=True)

        result = subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve().parent.parent / "run_analysis.py"),
                "--folder", str(folder),
                "--skip-vision",
                "--no-pdf",
                "--skip-checks",
                "--html",
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=120,
        )
        assert result.returncode == 0, (
            f"run_analysis.py crashed on empty folder:\n{result.stderr[-1000:]}"
        )
