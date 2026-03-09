"""Tests for the script-oriented analysis modules that use main() functions.

Covers:
- cross_reference_dwi.main: detailed cross-DWI comparison output
- cross_reference_summary.main: concise cross-DWI summary
- statistical_relevance.main: p-value and correlation extraction
- statistical_by_graph_type.main: graph-type-filtered analysis
- run_analysis._run_script: subprocess launcher helper

These scripts are primarily output-oriented (print to stdout), so we test
them by capturing stdout and verifying key content. We also test that they
exit cleanly when the required CSV is missing.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

import pytest


# ---------------------------------------------------------------------------
# cross_reference_dwi
# ---------------------------------------------------------------------------

class TestCrossReferenceDwi:
    """Test cross_reference_dwi.main stdout output.

    This script compares graph analysis results across DWI types (Standard,
    dnCNN, IVIMnet), looking for graphs that appear in multiple types and
    comparing their trends and summaries.
    """

    def test_outputs_matched_graphs(self, saved_files_with_graph_csv: Path, capsys):
        """Feature_BoxPlots appears in both Standard and dnCNN → should be matched.

        Patches sys.argv to pass the saved_files directory as the CLI argument,
        then captures stdout and checks that the matched graph name and both
        DWI types appear in the output.
        """
        with patch.object(sys, "argv", ["script", str(saved_files_with_graph_csv)]):
            from cross_reference_dwi import main
            main()
        out = capsys.readouterr().out
        assert "Feature_BoxPlots" in out
        assert "Standard" in out
        assert "dnCNN" in out

    def test_exits_on_missing_csv(self, saved_files_dir: Path):
        """Should sys.exit when graph_analysis_results.csv is absent."""
        with patch.object(sys, "argv", ["script", str(saved_files_dir)]):
            from cross_reference_dwi import main
            with pytest.raises(SystemExit):
                main()


# ---------------------------------------------------------------------------
# cross_reference_summary
# ---------------------------------------------------------------------------

class TestCrossReferenceSummary:
    """Test cross_reference_summary.main output.

    The summary script focuses on a curated list of priority graphs and
    produces a concise AGREE/DIFFER verdict for each, rather than the
    full detailed comparison from cross_reference_dwi.
    """

    def test_outputs_priority_graphs(self, saved_files_with_graph_csv: Path, capsys):
        """Feature_BoxPlots is a priority graph; its trend agreement/disagreement is reported."""
        with patch.object(sys, "argv", ["script", str(saved_files_with_graph_csv)]):
            from cross_reference_summary import main
            main()
        out = capsys.readouterr().out
        # Feature_BoxPlots is a priority graph and has 2 DWI types
        assert "Feature_BoxPlots" in out
        assert "AGREE" in out or "DIFFER" in out

    def test_exits_on_missing_csv(self, saved_files_dir: Path):
        """Should sys.exit when graph_analysis_results.csv is absent."""
        with patch.object(sys, "argv", ["script", str(saved_files_dir)]):
            from cross_reference_summary import main
            with pytest.raises(SystemExit):
                main()


# ---------------------------------------------------------------------------
# statistical_relevance
# ---------------------------------------------------------------------------

class TestStatisticalRelevance:
    """Test statistical_relevance.main output.

    This script extracts p-values and correlation coefficients from graph
    summaries and reports significant findings (p < 0.05) and notable
    correlations (|r| > some threshold).
    """

    def test_finds_significant_pvalues(self, saved_files_with_graph_csv: Path, capsys):
        """The fixture has p = 0.003 in Standard's summary, which should be flagged."""
        with patch.object(sys, "argv", ["script", str(saved_files_with_graph_csv)]):
            from statistical_relevance import main
            main()
        out = capsys.readouterr().out
        # p = 0.003 from the Standard Feature_BoxPlots summary
        assert "SIGNIFICANT" in out.upper() or "p=" in out

    def test_finds_correlations(self, saved_files_with_graph_csv: Path, capsys):
        """The fixture has r = 0.65 in Standard's summary, which should be reported."""
        with patch.object(sys, "argv", ["script", str(saved_files_with_graph_csv)]):
            from statistical_relevance import main
            main()
        out = capsys.readouterr().out
        assert "CORRELATION" in out.upper() or "r=" in out

    def test_exits_on_missing_csv(self, saved_files_dir: Path):
        """Should sys.exit when graph_analysis_results.csv is absent."""
        with patch.object(sys, "argv", ["script", str(saved_files_dir)]):
            from statistical_relevance import main
            with pytest.raises(SystemExit):
                main()


# ---------------------------------------------------------------------------
# statistical_by_graph_type
# ---------------------------------------------------------------------------

class TestStatisticalByGraphType:
    """Test statistical_by_graph_type.main output.

    This script groups graph analysis results by graph type (box, line,
    scatter, etc.) and reports per-type statistics: significant p-values,
    trend directions, and a summary table at the end.
    """

    def test_groups_by_graph_type(self, saved_files_with_graph_csv: Path, capsys):
        """Output should contain sections for the 'box' and 'line' types from the fixture."""
        with patch.object(sys, "argv", ["script", str(saved_files_with_graph_csv)]):
            from statistical_by_graph_type import main
            main()
        out = capsys.readouterr().out
        # The fixture has 'box' and 'line' graph types
        assert "BOX" in out.upper()
        assert "LINE" in out.upper()

    def test_summary_table(self, saved_files_with_graph_csv: Path, capsys):
        """The overall summary table should appear at the end."""
        with patch.object(sys, "argv", ["script", str(saved_files_with_graph_csv)]):
            from statistical_by_graph_type import main
            main()
        out = capsys.readouterr().out
        assert "SUMMARY TABLE" in out.upper()

    def test_trend_directions_counted(self, saved_files_with_graph_csv: Path, capsys):
        """Trend direction counts (Increasing/Decreasing) from the fixture appear in output."""
        with patch.object(sys, "argv", ["script", str(saved_files_with_graph_csv)]):
            from statistical_by_graph_type import main
            main()
        out = capsys.readouterr().out
        # The fixture has increasing and decreasing trends
        assert "Increasing" in out or "Decreasing" in out

    def test_exits_on_missing_csv(self, saved_files_dir: Path):
        """Should sys.exit when graph_analysis_results.csv is absent."""
        with patch.object(sys, "argv", ["script", str(saved_files_dir)]):
            from statistical_by_graph_type import main
            with pytest.raises(SystemExit):
                main()


# ---------------------------------------------------------------------------
# run_analysis._run_script
# ---------------------------------------------------------------------------

class TestRunAnalysisHelper:
    """Test the _run_script subprocess helper."""

    def test_run_missing_script(self, saved_files_dir: Path, capsys):
        from run_analysis import _run_script
        result = _run_script("nonexistent_script.py", saved_files_dir)
        assert result is False
        out = capsys.readouterr().out
        assert "SKIP" in out
