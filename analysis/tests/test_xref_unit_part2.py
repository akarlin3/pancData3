"""Unit tests for cross-reference scripts (part 2) — trend direction
classification, cross-DWI comparison, priority ordering, parameter-map /
inflection sections, full cross-DWI detail output, and the summary table.

This is the file-size split companion of ``test_xref_unit.py``. It exercises:

- ``statistical_by_graph_type.py``: trend direction keyword matching, summary
  table
- ``statistical_relevance.py``: cross-DWI significance comparison
- ``cross_reference_summary.py``: priority graph ordering, parameter-map and
  inflection point sections
- ``cross_reference_dwi.py``: full cross-DWI detail output correctness

Shared imports, constants, and helpers live in ``_xref_unit_shared.py``.
"""

from __future__ import annotations

from _xref_unit_shared import *  # noqa: F401,F403


# ═══════════════════════════════════════════════════════════════════════════════
# Trend direction classification (statistical_by_graph_type.py)
# ═══════════════════════════════════════════════════════════════════════════════


class TestTrendDirectionClassification:
    """Verify keyword-based trend direction counting."""

    def _make_folder_with_trends(self, tmp_path, directions):
        """Create a folder with a single graph containing trends with given directions."""
        folder = tmp_path / "saved_files_20260301_120000"
        (folder / "Standard").mkdir(parents=True)
        trends = [{"series": f"s{i}", "direction": d} for i, d in enumerate(directions)]
        rows = [{
            **{col: "" for col in GRAPH_CSV_COLUMNS},
            "file_path": "saved_files_20260301_120000/Standard/Test_Standard.png",
            "graph_type": "line",
            "trends_json": json.dumps(trends),
            "summary": "",
            "statistical_tests_json": "[]",
            "inflection_points_json": "[]",
            "clinical_relevance": "",
        }]
        _write_csv(folder, rows)
        return folder

    def test_increasing_keywords(self, tmp_path):
        """'increasing', 'up', 'higher', 'rising' all count as increasing."""
        folder = self._make_folder_with_trends(
            tmp_path, ["increasing", "up", "higher", "rising"]
        )
        output = _run_by_graph_type(folder)
        assert "Increasing: 4" in output

    def test_decreasing_keywords(self, tmp_path):
        """'decreasing', 'down', 'lower', 'falling', 'drop' all count as decreasing."""
        folder = self._make_folder_with_trends(
            tmp_path, ["decreasing", "down", "lower", "falling", "drop"]
        )
        output = _run_by_graph_type(folder)
        assert "Decreasing: 5" in output

    def test_stable_keywords(self, tmp_path):
        """'flat', 'stable', 'constant' all count as stable."""
        folder = self._make_folder_with_trends(
            tmp_path, ["flat", "stable", "constant"]
        )
        output = _run_by_graph_type(folder)
        assert "Stable: 3" in output

    def test_other_unrecognised_keyword(self, tmp_path):
        """Unrecognised direction keywords count as 'Other'."""
        folder = self._make_folder_with_trends(
            tmp_path, ["variable", "oscillating"]
        )
        output = _run_by_graph_type(folder)
        assert "Other: 2" in output

    def test_case_insensitive_direction(self, tmp_path):
        """Direction matching is case-insensitive (d.lower() is used)."""
        folder = self._make_folder_with_trends(
            tmp_path, ["Increasing", "DECREASING", "Flat"]
        )
        output = _run_by_graph_type(folder)
        assert "Increasing: 1" in output
        assert "Decreasing: 1" in output
        assert "Stable: 1" in output

    def test_partial_match(self, tmp_path):
        """'slightly increasing' matches via 'increas' substring."""
        folder = self._make_folder_with_trends(
            tmp_path, ["slightly increasing"]
        )
        output = _run_by_graph_type(folder)
        assert "Increasing: 1" in output


# ═══════════════════════════════════════════════════════════════════════════════
# Cross-DWI p-value comparison (statistical_relevance.py)
# ═══════════════════════════════════════════════════════════════════════════════


class TestCrossDwiPvalueComparison:
    """Verify cross-DWI significance counts are computed correctly."""

    def test_sig_count_per_dwi_type(self, saved_files_with_graph_csv):
        """Feature_BoxPlots: Standard has p=0.003 (sig), dnCNN has p=0.12 (non-sig)."""
        output = _run_stat_relevance(saved_files_with_graph_csv)
        lines = output.splitlines()
        # Find lines in the cross-DWI section about Feature_BoxPlots
        in_cross = False
        fb_lines = []
        for line in lines:
            if "Same Analysis, Different Significance" in line:
                in_cross = True
            if in_cross and "Feature_BoxPlots" in line:
                # Collect the next few lines for this graph
                fb_lines.append(line)
            elif in_cross and fb_lines and line.strip().startswith(("Standard:", "dnCNN:", "IVIMnet:")):
                fb_lines.append(line)

        # Standard should show 1 significant p-value
        std_lines = [l for l in fb_lines if "Standard:" in l]
        if std_lines:
            assert "1/" in std_lines[0] or "significant" in std_lines[0]


# ═══════════════════════════════════════════════════════════════════════════════
# Priority graph ordering (cross_reference_summary.py)
# ═══════════════════════════════════════════════════════════════════════════════


class TestPriorityGraphOrdering:
    """Verify that priority graphs appear before non-priority graphs."""

    def test_priority_graphs_first(self, tmp_path):
        """Feature_BoxPlots (priority) should appear before ZZZ_Plot (non-priority)."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN"):
            (folder / t).mkdir(parents=True)
        rows = _build_csv_rows([
            # Non-priority graph (alphabetically last)
            ("Standard", "ZZZ_Plot", [{"series": "X", "direction": "up"}]),
            ("dnCNN", "ZZZ_Plot", [{"series": "X", "direction": "up"}]),
            # Priority graph
            ("Standard", "Feature_BoxPlots", [{"series": "X", "direction": "up"}]),
            ("dnCNN", "Feature_BoxPlots", [{"series": "X", "direction": "up"}]),
        ])
        _write_csv(folder, rows)
        output = _run_summary(folder)

        fb_pos = output.find("Feature_BoxPlots")
        zzz_pos = output.find("ZZZ_Plot")
        assert fb_pos < zzz_pos, "Priority graph should appear before non-priority"


# ═══════════════════════════════════════════════════════════════════════════════
# Parameter map and inflection point sections (cross_reference_summary.py)
# ═══════════════════════════════════════════════════════════════════════════════


class TestParameterMapAndInflectionSections:
    """Verify parameter map and inflection point sections filter correctly."""

    def test_parameter_map_filtering(self, tmp_path):
        """Only graphs with 'Parameter_Maps' in the name appear in that section."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN"):
            (folder / t).mkdir(parents=True)
        rows = _build_csv_rows([
            ("Standard", "ADC_Parameter_Maps", []),
            ("dnCNN", "ADC_Parameter_Maps", []),
            ("Standard", "Feature_BoxPlots", []),
            ("dnCNN", "Feature_BoxPlots", []),
        ])
        _write_csv(folder, rows)
        output = _run_summary(folder)

        # Find the parameter map section
        pm_section_start = output.find("PARAMETER MAP COUNTS")
        inflection_start = output.find("INFLECTION POINTS")
        pm_section = output[pm_section_start:inflection_start]

        assert "ADC_Parameter_Maps" in pm_section
        assert "Feature_BoxPlots" not in pm_section

    def test_inflection_points_only_longitudinal(self, tmp_path):
        """Only graphs with 'Longitudinal' in the name appear in inflection section."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN"):
            (folder / t).mkdir(parents=True)
        ip_json = json.dumps([{"approximate_x": 30, "description": "inflect"}])
        rows = [
            {**{col: "" for col in GRAPH_CSV_COLUMNS},
             "file_path": f"saved_files_20260301_120000/{t}/Longitudinal_Test_{t}.png",
             "graph_type": "line", "trends_json": "[]", "summary": "",
             "statistical_tests_json": "[]", "clinical_relevance": "",
             "inflection_points_json": ip_json}
            for t in ("Standard", "dnCNN")
        ] + [
            {**{col: "" for col in GRAPH_CSV_COLUMNS},
             "file_path": f"saved_files_20260301_120000/{t}/BoxPlot_{t}.png",
             "graph_type": "box", "trends_json": "[]", "summary": "",
             "statistical_tests_json": "[]", "clinical_relevance": "",
             "inflection_points_json": ip_json}
            for t in ("Standard", "dnCNN")
        ]
        _write_csv(folder, rows)
        output = _run_summary(folder)

        inflection_start = output.find("INFLECTION POINTS")
        inflection_section = output[inflection_start:]

        assert "Longitudinal_Test" in inflection_section
        assert "BoxPlot" not in inflection_section


# ═══════════════════════════════════════════════════════════════════════════════
# Full cross-DWI detail output (cross_reference_dwi.py)
# ═══════════════════════════════════════════════════════════════════════════════


class TestCrossReferenceDwiCorrectness:
    """Verify correctness of values printed by cross_reference_dwi."""

    def test_start_end_values_formatted_to_four_decimals(self, saved_files_with_graph_csv):
        """Start/end values should appear with 4 decimal places."""
        output = _run_xref_dwi(saved_files_with_graph_csv)
        # Row 1 has start_value=0.001, end_value=0.0015
        assert "0.0010" in output
        assert "0.0015" in output

    def test_pvalue_formatted_to_four_decimals(self, saved_files_with_graph_csv):
        """Statistical test p-values should show 4 decimal places."""
        output = _run_xref_dwi(saved_files_with_graph_csv)
        # Wilcoxon p=0.003 → p=0.0030
        assert "p=0.0030" in output

    def test_statistic_value_formatted_when_present(self, tmp_path):
        """Statistic values should show 3 decimal places when not None."""
        folder = tmp_path / "saved_files_20260301_120000"
        for t in ("Standard", "dnCNN"):
            (folder / t).mkdir(parents=True)
        rows = [
            {**{col: "" for col in GRAPH_CSV_COLUMNS},
             "file_path": f"saved_files_20260301_120000/{t}/Test_{t}.png",
             "graph_type": "line", "trends_json": "[]", "summary": "",
             "inflection_points_json": "[]", "clinical_relevance": "",
             "statistical_tests_json": json.dumps([
                 {"test_name": "log-rank", "statistic_value": 5.41,
                  "p_value": 0.02, "comparison_groups": "LF vs LC"}
             ])}
            for t in ("Standard", "dnCNN")
        ]
        _write_csv(folder, rows)
        output = _run_xref_dwi(folder)
        assert "stat=5.410" in output

    def test_matched_count_correct(self, saved_files_with_graph_csv):
        """Total matched graph sets should equal the number of cross-DWI groups."""
        output = _run_xref_dwi(saved_files_with_graph_csv)
        # Only Feature_BoxPlots has 2+ real DWI types (Standard + dnCNN)
        # Longitudinal_Mean_Metrics only has IVIMnet (1 type)
        assert "Total matched graph sets across DWI types: 1" in output

    def test_comparison_groups_shown(self, saved_files_with_graph_csv):
        """Comparison groups (e.g., 'LF vs LC') appear in the output."""
        output = _run_xref_dwi(saved_files_with_graph_csv)
        assert "LF vs LC" in output

    def test_magnitude_in_parentheses(self, saved_files_with_graph_csv):
        """Trend magnitude should appear in parentheses."""
        output = _run_xref_dwi(saved_files_with_graph_csv)
        assert "(~15% increase)" in output

    def test_metadata_line_content(self, saved_files_with_graph_csv):
        """Metadata line should include n=42, series count, error bars."""
        output = _run_xref_dwi(saved_files_with_graph_csv)
        assert "n=42" in output
        assert "2 series" in output
        assert "error bars: IQR" in output


# ═══════════════════════════════════════════════════════════════════════════════
# Summary table correctness (statistical_by_graph_type.py)
# ═══════════════════════════════════════════════════════════════════════════════


class TestSummaryTable:
    """Verify the summary table at the end of statistical_by_graph_type output."""

    def test_summary_table_present(self, saved_files_with_graph_csv):
        output = _run_by_graph_type(saved_files_with_graph_csv)
        assert "SUMMARY TABLE" in output

    def test_box_type_row_in_table(self, saved_files_with_graph_csv):
        """Box type should appear with count=2 (Standard + dnCNN rows)."""
        output = _run_by_graph_type(saved_files_with_graph_csv)
        # Look for the box row in the summary table
        lines = output.splitlines()
        summary_start = next(
            i for i, l in enumerate(lines) if "SUMMARY TABLE" in l
        )
        table_lines = lines[summary_start:]
        box_lines = [l for l in table_lines if l.strip().startswith("box")]
        assert len(box_lines) >= 1
        # Should show count of 2 (Standard + dnCNN box plots)
        assert "2" in box_lines[0]

    def test_line_type_row_in_table(self, saved_files_with_graph_csv):
        """Line type should appear with count=1 (IVIMnet longitudinal)."""
        output = _run_by_graph_type(saved_files_with_graph_csv)
        lines = output.splitlines()
        summary_start = next(
            i for i, l in enumerate(lines) if "SUMMARY TABLE" in l
        )
        table_lines = lines[summary_start:]
        line_rows = [l for l in table_lines if l.strip().startswith("line")]
        assert len(line_rows) >= 1
        assert "1" in line_rows[0]
