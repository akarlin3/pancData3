"""Integration tests for the full analysis pipeline.

Runs ``run_analysis.py --skip-vision --no-pdf --skip-checks --html`` on a
synthetic ``saved_files_*`` directory assembled from conftest fixtures, then
verifies that the output HTML report exists and contains the expected
sections, tables, and figures.
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

# Re-use constants and fixtures from the shared conftest.
from conftest import GRAPH_CSV_COLUMNS, SAMPLE_GRAPH_CSV_ROWS

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

ANALYSIS_DIR = Path(__file__).resolve().parent.parent


def _write_graph_csv(folder: Path) -> None:
    """Write the standard 3-row graph_analysis_results.csv."""
    csv_path = folder / "graph_analysis_results.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
        writer.writeheader()
        for row in SAMPLE_GRAPH_CSV_ROWS:
            writer.writerow(row)


def _write_logs(folder: Path) -> None:
    """Write synthetic MATLAB log files under Standard/."""
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


def _write_csvs(folder: Path) -> None:
    """Write Significant_LF_Metrics.csv for Standard and dnCNN."""
    for dwi_type in ("Standard", "dnCNN"):
        csv_path = folder / dwi_type / "Significant_LF_Metrics.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["Metric", "Timepoint", "p_value"])
            writer.writeheader()
            if dwi_type == "Standard":
                writer.writerow({"Metric": "mean_adc", "Timepoint": "BL", "p_value": "0.003"})
                writer.writerow({"Metric": "mean_d", "Timepoint": "BL", "p_value": "0.01"})
                writer.writerow({"Metric": "mean_adc", "Timepoint": "W2", "p_value": "0.04"})
            else:
                writer.writerow({"Metric": "mean_adc", "Timepoint": "BL", "p_value": "0.005"})


# ---------------------------------------------------------------------------
# Fixture: fully-populated synthetic output folder
# ---------------------------------------------------------------------------

@pytest.fixture
def integration_folder(tmp_path: Path) -> Path:
    """Build a complete synthetic saved_files_* directory for integration testing.

    Combines graph CSV, MATLAB log files, and pipeline CSV exports into a
    single directory tree that exercises all data-loading paths in the
    report generator.
    """
    folder = tmp_path / "saved_files_20260301_120000"
    for dwi_type in ("Standard", "dnCNN", "IVIMnet"):
        (folder / dwi_type).mkdir(parents=True)

    _write_graph_csv(folder)
    _write_logs(folder)
    _write_csvs(folder)
    return folder


# ---------------------------------------------------------------------------
# Fixture: run the pipeline once and cache the HTML output
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def pipeline_result(tmp_path_factory) -> dict:
    """Run the full analysis pipeline and return the results.

    Uses module scope so the pipeline only runs once for all tests in this
    file, significantly reducing test execution time.

    Returns a dict with keys:
        - ``folder``: Path to the saved_files_* directory
        - ``returncode``: subprocess exit code
        - ``stdout``: captured stdout
        - ``stderr``: captured stderr
        - ``html_path``: Path to the generated analysis_report.html
        - ``html``: contents of the HTML report (empty string if not generated)
    """
    tmp_path = tmp_path_factory.mktemp("integration")
    folder = tmp_path / "saved_files_20260301_120000"
    for dwi_type in ("Standard", "dnCNN", "IVIMnet"):
        (folder / dwi_type).mkdir(parents=True)

    _write_graph_csv(folder)
    _write_logs(folder)
    _write_csvs(folder)

    result = subprocess.run(
        [
            sys.executable,
            str(ANALYSIS_DIR / "run_analysis.py"),
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

    html_path = folder / "analysis_report.html"
    html = html_path.read_text(encoding="utf-8") if html_path.exists() else ""

    return {
        "folder": folder,
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "html_path": html_path,
        "html": html,
    }


# ===========================================================================
# Test classes
# ===========================================================================

class TestPipelineExecution:
    """Tests that the pipeline subprocess completes successfully."""

    def test_exit_code_zero(self, pipeline_result):
        assert pipeline_result["returncode"] == 0, (
            f"run_analysis.py failed (exit {pipeline_result['returncode']}):\n"
            f"STDOUT:\n{pipeline_result['stdout']}\n"
            f"STDERR:\n{pipeline_result['stderr']}"
        )

    def test_html_report_created(self, pipeline_result):
        assert pipeline_result["html_path"].exists(), "analysis_report.html was not created"

    def test_html_report_nonempty(self, pipeline_result):
        assert len(pipeline_result["html"]) > 1000, (
            f"HTML report too short ({len(pipeline_result['html'])} chars)"
        )

    def test_log_file_created(self, pipeline_result):
        log_path = pipeline_result["folder"] / "run_analysis_output.log"
        assert log_path.exists(), "run_analysis_output.log was not created"

    def test_stdout_contains_summary(self, pipeline_result):
        out = pipeline_result["stdout"]
        assert "Analysis Complete" in out or "Analysis Suite" in out


class TestHTMLStructure:
    """Tests that the HTML report has valid structure."""

    def test_doctype(self, pipeline_result):
        assert pipeline_result["html"].startswith("<!DOCTYPE html>")

    def test_html_head(self, pipeline_result):
        html = pipeline_result["html"]
        assert "<head>" in html
        assert '<meta charset="utf-8">' in html
        assert "<title>" in html

    def test_html_body(self, pipeline_result):
        html = pipeline_result["html"]
        assert "<body>" in html
        assert "</body>" in html
        assert "</html>" in html

    def test_css_included(self, pipeline_result):
        assert "<style>" in pipeline_result["html"]

    def test_javascript_included(self, pipeline_result):
        assert "<script>" in pipeline_result["html"]


class TestReportSections:
    """Tests that major report sections are present in the HTML output."""

    # -- Part 1: Overview --
    def test_executive_summary(self, pipeline_result):
        assert "Executive Summary" in pipeline_result["html"]

    def test_methods_section(self, pipeline_result):
        assert "Methods" in pipeline_result["html"]

    def test_cohort_overview(self, pipeline_result):
        assert "Cohort" in pipeline_result["html"]

    def test_patient_flow(self, pipeline_result):
        assert "Patient Flow" in pipeline_result["html"]

    # -- Part 2: Data --
    def test_data_quality(self, pipeline_result):
        html = pipeline_result["html"]
        # Data quality / outlier section
        assert "Outlier" in html or "Data Quality" in html

    def test_hypothesis(self, pipeline_result):
        assert "Hypothesis" in pipeline_result["html"]

    def test_graph_overview(self, pipeline_result):
        assert "Graph Analysis Overview" in pipeline_result["html"]

    def test_graph_issues(self, pipeline_result):
        html = pipeline_result["html"]
        assert "Graph Issues" in html or "Issues" in html

    # -- Part 3: Statistics --
    def test_statistical_significance(self, pipeline_result):
        assert "Statistical Significance" in pipeline_result["html"]

    def test_effect_sizes(self, pipeline_result):
        assert "Effect Size" in pipeline_result["html"]

    def test_multiple_comparisons(self, pipeline_result):
        assert "Multiple Comparison Correction" in pipeline_result["html"]

    def test_cross_dwi_comparison(self, pipeline_result):
        assert "Cross-DWI" in pipeline_result["html"]

    # -- Part 4: Outcomes --
    def test_correlations(self, pipeline_result):
        assert "Correlation" in pipeline_result["html"]

    def test_treatment_response(self, pipeline_result):
        assert "Treatment Response" in pipeline_result["html"]

    def test_predictive_performance(self, pipeline_result):
        assert "Predictive" in pipeline_result["html"]

    # -- Part 5: Discussion --
    def test_limitations(self, pipeline_result):
        assert "Limitation" in pipeline_result["html"]

    def test_conclusions(self, pipeline_result):
        assert "Conclusion" in pipeline_result["html"]

    # -- Part 6: Appendix --
    def test_appendix(self, pipeline_result):
        assert "Appendix" in pipeline_result["html"]

    def test_reporting_checklist(self, pipeline_result):
        assert "Reporting Checklist" in pipeline_result["html"] or "STROBE" in pipeline_result["html"]


class TestPartBreaks:
    """Tests that the report contains the expected part dividers."""

    @pytest.mark.parametrize("part", [
        "Part 1",
        "Part 2",
        "Part 3",
        "Part 4",
        "Part 5",
        "Part 6",
    ])
    def test_part_break_present(self, pipeline_result, part):
        assert part in pipeline_result["html"]


class TestReportMetadata:
    """Tests that report metadata is correctly embedded."""

    def test_dwi_types_listed(self, pipeline_result):
        html = pipeline_result["html"]
        assert "Standard" in html
        assert "dnCNN" in html
        assert "IVIMnet" in html

    def test_timestamp_present(self, pipeline_result):
        html = pipeline_result["html"]
        # The folder name is saved_files_20260301_120000, which gets
        # formatted as "March 1, 2026 at 12:00:00" in the report.
        assert "March 1, 2026" in html or "20260301" in html

    def test_generation_date(self, pipeline_result):
        html = pipeline_result["html"]
        assert "Generated" in html


class TestTablesPresent:
    """Tests that key tables are present in the HTML output."""

    def test_table_tags_exist(self, pipeline_result):
        html = pipeline_result["html"]
        assert "<table" in html
        assert "</table>" in html

    def test_multiple_tables(self, pipeline_result):
        html = pipeline_result["html"]
        table_count = html.count("<table")
        assert table_count >= 3, f"Expected at least 3 tables, found {table_count}"

    def test_table_headers(self, pipeline_result):
        html = pipeline_result["html"]
        assert "<th" in html

    def test_table_rows(self, pipeline_result):
        html = pipeline_result["html"]
        tr_count = html.count("<tr")
        assert tr_count >= 5, f"Expected at least 5 table rows, found {tr_count}"


class TestGraphData:
    """Tests that graph analysis data from the CSV appears in the report."""

    def test_graph_types_mentioned(self, pipeline_result):
        html = pipeline_result["html"]
        # The fixture includes box and line graph types.
        assert "box" in html.lower()
        assert "line" in html.lower()

    def test_p_values_appear(self, pipeline_result):
        html = pipeline_result["html"]
        assert "0.003" in html

    def test_sample_sizes(self, pipeline_result):
        html = pipeline_result["html"]
        assert "42" in html or "n=42" in html

    def test_trend_directions(self, pipeline_result):
        html = pipeline_result["html"]
        lower = html.lower()
        assert "increasing" in lower or "decreasing" in lower


class TestLogParsedData:
    """Tests that data parsed from MATLAB logs appears in the report."""

    def test_auc_values(self, pipeline_result):
        html = pipeline_result["html"]
        assert "0.781" in html or "0.78" in html

    def test_hazard_ratios(self, pipeline_result):
        html = pipeline_result["html"]
        assert "1.250" in html or "1.25" in html

    def test_elastic_net_features(self, pipeline_result):
        html = pipeline_result["html"]
        assert "mean_adc" in html

    def test_outlier_counts(self, pipeline_result):
        html = pipeline_result["html"]
        # Baseline log has "Total outliers removed: 4 / 42"
        assert "outlier" in html.lower()


class TestCSVParsedData:
    """Tests that pipeline CSV exports appear in the report."""

    def test_significant_metrics(self, pipeline_result):
        html = pipeline_result["html"]
        # Standard has mean_adc significant at BL
        assert "mean_adc" in html

    def test_dwi_badges(self, pipeline_result):
        html = pipeline_result["html"]
        # Report should include DWI type badges
        assert "Standard" in html


class TestNavigationBar:
    """Tests that the sticky navigation bar is present."""

    def test_nav_element(self, pipeline_result):
        html = pipeline_result["html"]
        assert "<nav" in html or "nav-bar" in html or "navbar" in html.lower()


class TestNoVisionAnalysis:
    """Tests that vision analysis was properly skipped."""

    def test_vision_skipped_in_stdout(self, pipeline_result):
        out = pipeline_result["stdout"]
        assert "skip" in out.lower() or "vision" in out.lower()


class TestReportWithMinimalData:
    """Tests that the pipeline handles a minimal (nearly empty) folder gracefully."""

    def test_empty_folder_does_not_crash(self, tmp_path):
        """Pipeline should produce a report even with no data files."""
        folder = tmp_path / "saved_files_20260315_100000"
        for dwi_type in ("Standard", "dnCNN", "IVIMnet"):
            (folder / dwi_type).mkdir(parents=True)

        result = subprocess.run(
            [
                sys.executable,
                str(ANALYSIS_DIR / "run_analysis.py"),
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
            f"Pipeline crashed on empty folder:\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
        html_path = folder / "analysis_report.html"
        assert html_path.exists(), "Report not generated for empty folder"

        html = html_path.read_text(encoding="utf-8")
        assert "<!DOCTYPE html>" in html
        assert "Executive Summary" in html
