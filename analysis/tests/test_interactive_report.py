"""Tests for generate_interactive_report.py — interactive HTML report.

Covers:
- Interactive report generation with empty and populated data
- Sidebar filter rendering (DWI types, patients, core methods)
- Chart.js script inclusion
- Data injection (JSON blob in script tag)
- Section builders: overview, patient explorer, visualisations,
  significance tables, graph explorer, core comparison, dosimetry
- Filtering-related data attributes on table rows
- Tab panel structure
- Helper functions: _extract_patients, _extract_core_methods
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest  # type: ignore

from report.generate_interactive_report import (  # type: ignore
    _build_core_comparison,
    _build_dosimetry_section,
    _build_graph_explorer,
    _build_overview_section,
    _build_patient_explorer,
    _build_significance_table,
    _build_visualizations_section,
    _dwi_badge,
    _esc,
    _extract_core_methods,
    _extract_patients,
    _sig_class,
    _trend_tag,
    generate_interactive_report,
)


# ---------------------------------------------------------------------------
# Helper function tests
# ---------------------------------------------------------------------------


class TestEsc:
    """Verify HTML escaping."""

    def test_escapes_angle_brackets(self):
        assert "&lt;" in _esc("<script>")

    def test_plain_text_unchanged(self):
        assert _esc("hello") == "hello"


class TestDwiBadge:
    """Verify DWI type badge rendering."""

    def test_standard_badge(self):
        badge = _dwi_badge("Standard")
        assert "badge-standard" in badge
        assert "Standard" in badge

    def test_dncnn_badge(self):
        assert "badge-dncnn" in _dwi_badge("dnCNN")

    def test_ivimnet_badge(self):
        assert "badge-ivimnet" in _dwi_badge("IVIMnet")


class TestSigClass:
    """Verify p-value CSS class mapping."""

    def test_highly_significant(self):
        assert _sig_class(0.0001) == "sig-3"

    def test_significant(self):
        assert _sig_class(0.005) == "sig-2"

    def test_noteworthy(self):
        assert _sig_class(0.03) == "sig-1"

    def test_not_significant(self):
        assert _sig_class(0.5) == ""


class TestTrendTag:
    """Verify trend tag rendering."""

    def test_increasing(self):
        tag = _trend_tag("increasing")
        assert "trend-incr" in tag
        assert "\u2191" in tag

    def test_decreasing(self):
        tag = _trend_tag("decreasing")
        assert "trend-decr" in tag

    def test_stable(self):
        assert "trend-flat" in _trend_tag("stable")

    def test_non_monotonic(self):
        assert "trend-nm" in _trend_tag("U-shaped")


# ---------------------------------------------------------------------------
# Data extraction tests
# ---------------------------------------------------------------------------


class TestExtractPatients:
    """Verify patient data extraction from log and MAT data."""

    def test_empty_data_returns_empty(self):
        result = _extract_patients(None, {})
        assert result == []

    def test_extracts_from_mat_longitudinal(self):
        mat_data = {
            "Standard": {
                "longitudinal": {
                    "PT001": [
                        {"timepoint": "Fx1", "adc_mean": 0.0012, "d_mean": 0.001},
                        {"timepoint": "Fx5", "adc_mean": 0.0014},
                    ],
                    "PT002": [
                        {"timepoint": "Fx1", "adc_mean": 0.0011},
                    ],
                },
            },
        }
        patients = _extract_patients(None, mat_data)
        assert len(patients) == 2
        ids = [p["id"] for p in patients]
        assert "PT001" in ids
        assert "PT002" in ids
        pt1 = next(p for p in patients if p["id"] == "PT001")
        assert len(pt1["timepoints"]) == 2
        assert pt1["timepoints"][0]["adc_mean"] == 0.0012

    def test_skips_nan_values(self):
        mat_data = {
            "Standard": {
                "longitudinal": {
                    "PT001": [{"timepoint": "Fx1", "adc_mean": float("nan")}],
                },
            },
        }
        patients = _extract_patients(None, mat_data)
        pt = patients[0]
        tp = pt["timepoints"][0]
        assert "adc_mean" not in tp  # NaN should be skipped

    def test_sorted_by_id(self):
        mat_data = {
            "Standard": {
                "longitudinal": {
                    "Z_patient": [{"timepoint": "Fx1", "adc_mean": 0.001}],
                    "A_patient": [{"timepoint": "Fx1", "adc_mean": 0.002}],
                },
            },
        }
        patients = _extract_patients(None, mat_data)
        assert patients[0]["id"] == "A_patient"
        assert patients[1]["id"] == "Z_patient"


class TestExtractCoreMethods:
    """Verify core method name extraction from MAT data."""

    def test_empty_data(self):
        assert _extract_core_methods({}) == []

    def test_extracts_method_names(self):
        mat_data = {
            "Standard": {
                "core_comparison": {"otsu": {}, "kmeans": {}, "gmm": {}},
            },
        }
        methods = _extract_core_methods(mat_data)
        assert sorted(methods) == ["gmm", "kmeans", "otsu"]


# ---------------------------------------------------------------------------
# Section builder tests
# ---------------------------------------------------------------------------


class TestBuildOverview:
    """Verify the overview section builder."""

    def test_contains_patient_count(self):
        patients = [{"id": "PT001"}, {"id": "PT002"}]
        h = _build_overview_section(None, ["Standard"], [], None, {}, patients)
        html = "\n".join(h)
        assert "2" in html
        assert "Patients" in html

    def test_contains_dwi_types(self):
        h = _build_overview_section(None, ["Standard", "dnCNN"], [], None, {}, [])
        html = "\n".join(h)
        assert "Standard" in html
        assert "dnCNN" in html


class TestBuildPatientExplorer:
    """Verify the patient explorer section."""

    def test_empty_patients(self):
        h = _build_patient_explorer([], ["Standard"])
        html = "\n".join(h)
        assert "No patient data available" in html

    def test_renders_patient_cards(self):
        patients = [
            {
                "id": "PT001",
                "outcome": "",
                "timepoints": [
                    {"dwi_type": "Standard", "label": "Fx1", "adc_mean": 0.0012},
                ],
            },
        ]
        h = _build_patient_explorer(patients, ["Standard"])
        html = "\n".join(h)
        assert "PT001" in html
        assert "patient-card" in html
        assert "adc_mean" in html

    def test_search_bar_present(self):
        h = _build_patient_explorer([], ["Standard"])
        html = "\n".join(h)
        assert "search-input" in html


class TestBuildVisualizationsSection:
    """Verify the visualisations section with tabs and chart canvases."""

    def test_tab_bar_present(self):
        h = _build_visualizations_section(["Standard"])
        html = "\n".join(h)
        assert "tab-bar" in html
        assert "Longitudinal" in html
        assert "Distribution" in html
        assert "DWI Comparison" in html

    def test_chart_canvases_present(self):
        h = _build_visualizations_section(["Standard"])
        html = "\n".join(h)
        assert "chart-longitudinal" in html
        assert "chart-distribution" in html
        assert "chart-dwi-compare" in html


class TestBuildSignificanceTable:
    """Verify the significance table section."""

    def test_no_data(self):
        h = _build_significance_table(None, None, [])
        html = "\n".join(h)
        assert "No significance data available" in html

    def test_fdr_table_rendered(self):
        csv_data = {
            "fdr_global": {
                "Standard": [
                    {"Metric": "ADC_mean", "p_value": "0.003", "q_value": "0.01"},
                ],
            },
        }
        h = _build_significance_table(csv_data, None, ["Standard"])
        html = "\n".join(h)
        assert "FDR Global" in html
        assert "ADC_mean" in html
        assert "sortable" in html

    def test_hazard_ratios_rendered(self):
        log_data = {
            "Standard": {
                "hazard_ratios": [
                    {"feature": "ADC_mean", "hr": 1.5, "ci_lo": 1.1, "ci_hi": 2.0, "p": 0.02},
                ],
            },
        }
        h = _build_significance_table(None, log_data, ["Standard"])
        html = "\n".join(h)
        assert "Hazard Ratios" in html
        assert "ADC_mean" in html


class TestBuildGraphExplorer:
    """Verify the graph analysis explorer section."""

    def test_no_rows(self):
        h = _build_graph_explorer([], {})
        html = "\n".join(h)
        assert "No graph analysis data" in html

    def test_with_rows(self):
        rows = [
            {"graph_name": "Feature_BoxPlots", "dwi_type": "Standard",
             "x_axis": "Feature", "summary": "Box plot"},
        ]
        h = _build_graph_explorer(rows, {})
        html = "\n".join(h)
        assert "Feature_BoxPlots" in html
        assert "Table View" in html
        assert "Cross-DWI" in html


class TestBuildCoreComparison:
    """Verify the core method comparison section."""

    def test_no_data(self):
        h = _build_core_comparison({})
        html = "\n".join(h)
        assert "No core method comparison data" in html

    def test_dice_matrix_rendered(self):
        mat_data = {
            "Standard": {
                "core_comparison": {
                    "dice_matrix": [[1.0, 0.8], [0.8, 1.0]],
                    "method_names": ["otsu", "kmeans"],
                },
            },
        }
        h = _build_core_comparison(mat_data)
        html = "\n".join(h)
        assert "otsu" in html
        assert "kmeans" in html
        assert "1.000" in html
        assert "0.800" in html


class TestBuildDosimetry:
    """Verify the dosimetry section."""

    def test_no_data(self):
        h = _build_dosimetry_section({}, ["Standard"])
        html = "\n".join(h)
        assert "No dosimetry data" in html

    def test_with_data(self):
        mat_data = {"Standard": {"dosimetry": {"D95": 45.2, "V50": 89.1}}}
        h = _build_dosimetry_section(mat_data, ["Standard"])
        html = "\n".join(h)
        assert "45.20" in html
        assert "89.10" in html


# ---------------------------------------------------------------------------
# Integration tests
# ---------------------------------------------------------------------------


class TestGenerateInteractiveReport:
    """Integration tests for the full interactive report."""

    def test_report_contains_header(self, saved_files_with_graph_csv: Path):
        report = generate_interactive_report(saved_files_with_graph_csv)
        assert "Interactive Analysis Report" in report
        assert "March 1, 2026" in report

    def test_report_is_valid_html(self, saved_files_with_graph_csv: Path):
        report = generate_interactive_report(saved_files_with_graph_csv)
        assert report.startswith("<!DOCTYPE html>")
        assert "</html>" in report
        assert "<body>" in report

    def test_report_contains_chartjs(self, saved_files_with_graph_csv: Path):
        report = generate_interactive_report(saved_files_with_graph_csv)
        assert "chart.js" in report.lower() or "Chart" in report

    def test_report_contains_json_data_blob(self, saved_files_with_graph_csv: Path):
        report = generate_interactive_report(saved_files_with_graph_csv)
        assert "REPORT_DATA" in report

    def test_report_contains_sidebar(self, saved_files_with_graph_csv: Path):
        report = generate_interactive_report(saved_files_with_graph_csv)
        assert "sidebar" in report
        assert "dwi-filter-checks" in report
        assert "patient-select" in report

    def test_report_contains_tabs(self, saved_files_with_graph_csv: Path):
        report = generate_interactive_report(saved_files_with_graph_csv)
        assert "tab-bar" in report
        assert "tab-btn" in report

    def test_report_contains_sections(self, saved_files_with_graph_csv: Path):
        report = generate_interactive_report(saved_files_with_graph_csv)
        assert 'id="overview"' in report
        assert 'id="visualizations"' in report
        assert 'id="patients"' in report
        assert 'id="significance"' in report
        assert 'id="graphs"' in report

    def test_report_footer(self, saved_files_with_graph_csv: Path):
        report = generate_interactive_report(saved_files_with_graph_csv)
        assert "Interactive report generated by pancData3" in report

    def test_empty_folder(self, saved_files_dir: Path):
        """An empty folder should still produce a valid report."""
        report = generate_interactive_report(saved_files_dir)
        assert "Interactive Analysis Report" in report
        assert "REPORT_DATA" in report

    def test_report_contains_dwi_types(self, saved_files_with_graph_csv: Path):
        report = generate_interactive_report(saved_files_with_graph_csv)
        assert "Standard" in report
        assert "dnCNN" in report
        assert "IVIMnet" in report

    def test_data_blob_is_valid_json(self, saved_files_with_graph_csv: Path):
        """The embedded JSON data blob should be parseable."""
        report = generate_interactive_report(saved_files_with_graph_csv)
        # Extract JSON between "window.REPORT_DATA = " and ";"
        marker = "window.REPORT_DATA = "
        start = report.index(marker) + len(marker)
        end = report.index(";", start)
        json_str = report[start:end]
        data = json.loads(json_str)
        assert "patients" in data
        assert "dwi_types" in data
        assert isinstance(data["dwi_types"], list)

    def test_sortable_tables(self, saved_files_with_graph_csv: Path):
        """Tables should have the sortable class for interactive sorting."""
        report = generate_interactive_report(saved_files_with_graph_csv)
        assert "sortable" in report

    def test_graph_rows_in_table(self, saved_files_with_graph_csv: Path):
        """Graph analysis rows should be rendered in the table."""
        report = generate_interactive_report(saved_files_with_graph_csv)
        assert "Feature_BoxPlots" in report or "Feature BoxPlots" in report
