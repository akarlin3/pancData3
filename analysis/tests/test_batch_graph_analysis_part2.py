"""Tests for batch_graph_analysis.py — part 2 (split from test_batch_graph_analysis.py).

Covers timeout configuration, rate-limit detection/coordination, CLI flag
parsing, local fallback analysis, shared response parsing, and dual-provider
comparison.

See test_batch_graph_analysis.py for the companion suite (image discovery,
encoding, Pydantic schema, and flatten).
"""

from __future__ import annotations

import asyncio
import json
from pathlib import Path


from parsers.batch_graph_analysis import (
    COMPARISON_COLUMNS,
    CSV_COLUMNS,
    GraphAnalysis,
    InflectionPoint,
    REQUEST_TIMEOUT,
    Trend,
    _RateLimitCoordinator,
    _compare_results,
    _get_provider,
    _is_claude_rate_limit_error,
    _is_rate_limit_error,
    _parse_vision_response,
    _should_use_local,
    analyze_image_local,
    flatten,
)


# ---------------------------------------------------------------------------
# Timeout configuration
# ---------------------------------------------------------------------------

class TestTimeoutConfig:
    """Verify that timeout settings are loaded from config."""

    def test_request_timeout_is_positive(self):
        """REQUEST_TIMEOUT must be a positive number (prevents indefinite hangs)."""
        assert isinstance(REQUEST_TIMEOUT, (int, float))
        assert REQUEST_TIMEOUT > 0

    def test_request_timeout_matches_config(self):
        """REQUEST_TIMEOUT should match the vision config value."""
        from shared import get_config
        expected = get_config()["vision"]["request_timeout_seconds"]
        assert REQUEST_TIMEOUT == expected


# ---------------------------------------------------------------------------
# _is_rate_limit_error
# ---------------------------------------------------------------------------

class TestIsRateLimitError:
    """Verify detection of rate-limit / quota error messages."""

    def test_429_status_code(self):
        """HTTP 429 status code in error message is detected."""
        assert _is_rate_limit_error(Exception("HTTP 429 Too Many Requests"))

    def test_rate_keyword(self):
        """'rate' keyword in error message is detected."""
        assert _is_rate_limit_error(Exception("Rate limit exceeded"))

    def test_resource_exhausted(self):
        """'resource' keyword (Google API resource exhausted) is detected."""
        assert _is_rate_limit_error(Exception("Resource exhausted"))

    def test_quota_exceeded(self):
        """'quota' keyword is detected."""
        assert _is_rate_limit_error(Exception("Quota exceeded for project"))

    def test_non_rate_limit_error(self):
        """Normal errors (network, auth) are not misclassified."""
        assert not _is_rate_limit_error(Exception("Connection refused"))
        assert not _is_rate_limit_error(Exception("Invalid API key"))
        assert not _is_rate_limit_error(Exception("Internal server error 500"))


# ---------------------------------------------------------------------------
# _RateLimitCoordinator
# ---------------------------------------------------------------------------

class TestRateLimitCoordinator:
    """Verify shared rate-limit coordination across workers."""

    def test_no_cooldown_initially(self):
        """A fresh coordinator should not block."""
        async def _run():
            coord = _RateLimitCoordinator()
            # Should return immediately (no cooldown active).
            await coord.wait_if_cooling_down()
        asyncio.run(_run())

    def test_signal_sets_cooldown(self):
        """After signal(), the coordinator should report a cooldown."""
        async def _run():
            coord = _RateLimitCoordinator()
            await coord.signal()
            # _resume_at should be in the future.
            loop = asyncio.get_event_loop()
            now = loop.time()
            assert coord._resume_at > now
        asyncio.run(_run())

    def test_consecutive_hits_increase_cooldown(self):
        """Multiple consecutive signals should increase the cooldown."""
        async def _run():
            coord = _RateLimitCoordinator()
            await coord.signal()
            first_resume = coord._resume_at
            await coord.signal()
            second_resume = coord._resume_at
            # Second cooldown should extend further into the future.
            assert second_resume >= first_resume
        asyncio.run(_run())

    def test_success_resets_consecutive_counter(self):
        """record_success() should reset the consecutive hit counter."""
        async def _run():
            coord = _RateLimitCoordinator()
            await coord.signal()
            await coord.signal()
            assert coord._consecutive_hits == 2
            await coord.record_success()
            assert coord._consecutive_hits == 0
        asyncio.run(_run())


# ---------------------------------------------------------------------------
# _should_use_local
# ---------------------------------------------------------------------------

class TestShouldUseLocal:
    """Verify --local CLI flag detection."""

    def test_local_flag_present(self):
        """--local flag is detected in argv."""
        assert _should_use_local(["script.py", "--local", "/some/path"])

    def test_local_flag_absent(self):
        """Returns False when --local is not in argv."""
        assert not _should_use_local(["script.py", "/some/path"])

    def test_empty_argv(self):
        """Returns False for empty argv."""
        assert not _should_use_local([])


# ---------------------------------------------------------------------------
# analyze_image_local
# ---------------------------------------------------------------------------

class TestAnalyzeImageLocal:
    """Verify filename-based local fallback analysis."""

    def test_line_graph_from_longitudinal_filename(self):
        """Longitudinal filenames are inferred as 'line' graph type."""
        p = Path("saved_files_20240115/Standard/Longitudinal_Mean_Metrics_Standard.png")
        result = analyze_image_local(p)
        assert result.graph_type == "line"
        assert result.file_path == str(p)
        assert "local fallback" in result.summary.lower()

    def test_heatmap_from_dice_filename(self):
        """Dice heatmap filenames are inferred as 'heatmap' graph type."""
        p = Path("saved_files_20240115/Standard/core_method_dice_heatmap.png")
        result = analyze_image_local(p)
        assert result.graph_type == "heatmap"

    def test_bar_chart_from_volume_comparison(self):
        """Volume comparison filenames are inferred as 'bar' graph type."""
        p = Path("saved_files_20240115/Standard/core_method_volume_comparison.png")
        result = analyze_image_local(p)
        assert result.graph_type == "bar"

    def test_scatter_from_correlation_filename(self):
        """Correlation filenames are inferred as 'scatter' graph type."""
        p = Path("saved_files_20240115/Standard/Dose_Correlation_Standard.png")
        result = analyze_image_local(p)
        assert result.graph_type == "scatter"

    def test_histogram_from_hist_filename(self):
        """Histogram filenames are inferred as 'histogram' graph type."""
        p = Path("saved_files_20240115/Standard/Feature_Histograms_Standard.png")
        result = analyze_image_local(p)
        assert result.graph_type == "histogram"

    def test_box_from_boxplot_filename(self):
        """Boxplot filenames are inferred as 'box' graph type."""
        p = Path("saved_files_20240115/Standard/Feature_BoxPlots_Standard.png")
        result = analyze_image_local(p)
        assert result.graph_type == "box"

    def test_unknown_type_for_unrecognised_filename(self):
        """Unrecognised filenames get 'unknown' graph type."""
        p = Path("saved_files_20240115/Standard/mystery_plot.png")
        result = analyze_image_local(p)
        assert result.graph_type == "unknown"

    def test_clinical_relevance_from_dwi_path(self):
        """DWI type in path is reflected in clinical_relevance."""
        p = Path("saved_files_20240115/dnCNN/some_graph.png")
        result = analyze_image_local(p)
        assert result.clinical_relevance is not None
        assert "dnCNN" in result.clinical_relevance

    def test_no_clinical_relevance_for_root_path(self):
        """Files not inside a DWI-type subfolder get no clinical relevance."""
        p = Path("saved_files_20240115/some_root_graph.png")
        result = analyze_image_local(p)
        assert result.clinical_relevance is None

    def test_comparison_type_longitudinal(self):
        """Longitudinal filenames get comparison_type='longitudinal'."""
        p = Path("saved_files_20240115/Standard/Longitudinal_Mean_Metrics.png")
        result = analyze_image_local(p)
        assert result.comparison_type == "longitudinal"

    def test_comparison_type_dose_response(self):
        """Dose_vs filenames get comparison_type='dose-response'."""
        p = Path("saved_files_20240115/Standard/Dose_vs_Diffusion.png")
        result = analyze_image_local(p)
        assert result.comparison_type == "dose-response"

    def test_graph_title_from_filename(self):
        """Graph title is derived from filename with underscores replaced."""
        p = Path("saved_files_20240115/Standard/Feature_BoxPlots_Standard.png")
        result = analyze_image_local(p)
        assert result.graph_title == "Feature BoxPlots Standard"

    def test_x_axis_for_longitudinal(self):
        """Longitudinal graphs get a time-related x-axis."""
        p = Path("saved_files_20240115/Standard/Longitudinal_Mean_Metrics.png")
        result = analyze_image_local(p)
        assert result.x_axis is not None
        assert "time" in result.x_axis.label.lower() or "timepoint" in result.x_axis.label.lower()

    def test_y_axis_for_adc_graph(self):
        """ADC-related filenames get an ADC y-axis hint."""
        p = Path("saved_files_20240115/Standard/ADC_distribution.png")
        result = analyze_image_local(p)
        assert result.y_axis is not None
        assert "adc" in result.y_axis.label.lower()

    def test_figure_quality_is_unknown(self):
        """Local fallback sets figure_quality to 'unknown' (no visual inspection)."""
        p = Path("saved_files_20240115/Standard/some_graph.png")
        result = analyze_image_local(p)
        assert result.figure_quality == "unknown"

    def test_result_is_valid_graph_analysis(self):
        """The returned object is a valid GraphAnalysis that can be flattened."""
        p = Path("saved_files_20240115/Standard/Longitudinal_Mean_Metrics_Standard.png")
        result = analyze_image_local(p)
        row = flatten(result)
        for col in CSV_COLUMNS:
            assert col in row, f"Missing column: {col}"

    def test_survival_inferred_as_line(self):
        """Survival/Kaplan-Meier filenames are inferred as 'line' type."""
        p = Path("saved_files_20240115/Standard/Kaplan_Meier_Curve_Standard.png")
        result = analyze_image_local(p)
        assert result.graph_type == "line"

    def test_parameter_map_inferred(self):
        """Parameter map filenames are inferred as 'parameter_map' type."""
        p = Path("saved_files_20240115/Standard/ADC_parameter_map_Standard.png")
        result = analyze_image_local(p)
        assert result.graph_type == "parameter_map"


# ---------------------------------------------------------------------------
# _is_claude_rate_limit_error
# ---------------------------------------------------------------------------

class TestIsClaudeRateLimitError:
    """Verify detection of Claude rate-limit error messages."""

    def test_429_status_code(self):
        """HTTP 429 in error message is detected."""
        assert _is_claude_rate_limit_error(Exception("Error 429 rate_limit_error"))

    def test_rate_limit_keyword(self):
        """'rate_limit' keyword is detected."""
        assert _is_claude_rate_limit_error(Exception("rate_limit_error: too many requests"))

    def test_rate_limit_with_space(self):
        """'rate limit' with space is detected."""
        assert _is_claude_rate_limit_error(Exception("Rate limit exceeded"))

    def test_overloaded(self):
        """'overloaded' keyword is detected."""
        assert _is_claude_rate_limit_error(Exception("API is overloaded"))

    def test_non_rate_limit_error(self):
        """Normal errors are not misclassified."""
        assert not _is_claude_rate_limit_error(Exception("Invalid API key"))
        assert not _is_claude_rate_limit_error(Exception("Connection refused"))


# ---------------------------------------------------------------------------
# _get_provider
# ---------------------------------------------------------------------------

class TestGetProvider:
    """Verify --provider CLI flag extraction."""

    def test_gemini_provider(self):
        """--provider gemini returns 'gemini'."""
        assert _get_provider(["script.py", "--provider", "gemini"]) == "gemini"

    def test_claude_provider(self):
        """--provider claude returns 'claude'."""
        assert _get_provider(["script.py", "--provider", "claude"]) == "claude"

    def test_both_provider(self):
        """--provider both returns 'both'."""
        assert _get_provider(["script.py", "--provider", "both"]) == "both"

    def test_default_without_flag(self):
        """Without --provider flag, returns the config default."""
        result = _get_provider(["script.py", "/some/path"])
        # Should be the default from config (typically "gemini")
        assert isinstance(result, str)

    def test_mixed_flags(self):
        """--provider works alongside other flags."""
        result = _get_provider(["script.py", "--local", "--provider", "claude", "/path"])
        assert result == "claude"

    def test_provider_case_insensitive(self):
        """--provider BOTH is lowercased to 'both'."""
        assert _get_provider(["script.py", "--provider", "BOTH"]) == "both"


# ---------------------------------------------------------------------------
# _parse_vision_response
# ---------------------------------------------------------------------------

class TestParseVisionResponse:
    """Verify shared JSON response parsing."""

    def test_valid_json_response(self):
        """A well-formed JSON response is parsed correctly."""
        raw = json.dumps({
            "graph_type": "line",
            "summary": "A simple line graph.",
            "trends": [],
            "inflection_points": [],
            "statistical_tests": [],
            "outliers": [],
            "reference_lines": [],
            "issues": [],
            "annotations": [],
            "legend_items": [],
        })
        result = _parse_vision_response(raw, Path("test.png"), "test.png")
        assert result.graph_type == "line"
        assert result.summary == "A simple line graph."

    def test_json_with_code_fences(self):
        """Markdown code fences around JSON are stripped."""
        raw = '```json\n{"graph_type": "bar", "summary": "Bar chart."}\n```'
        result = _parse_vision_response(raw, Path("test.png"), "test.png")
        assert result.graph_type == "bar"

    def test_invalid_json_returns_unknown(self):
        """Unparseable JSON returns graph_type='unknown'."""
        result = _parse_vision_response("not json at all", Path("test.png"), "test.png")
        assert result.graph_type == "unknown"
        assert "parse error" in result.summary.lower()

    def test_file_path_injected(self):
        """The file_path field is set from the rel_path argument."""
        raw = json.dumps({"graph_type": "scatter", "summary": "dots"})
        result = _parse_vision_response(raw, Path("a.png"), "my/path.png")
        assert result.file_path == "my/path.png"


# ---------------------------------------------------------------------------
# _compare_results
# ---------------------------------------------------------------------------

class TestCompareResults:
    """Verify dual-provider comparison logic."""

    def _make_ga(self, **kwargs) -> GraphAnalysis:
        """Helper to create a GraphAnalysis with defaults."""
        defaults = {
            "file_path": "test.png",
            "graph_type": "line",
            "summary": "A graph.",
        }
        defaults.update(kwargs)
        return GraphAnalysis(**defaults)

    def test_identical_results_no_differences(self):
        """Two identical results produce differences='none'."""
        ga = self._make_ga()
        row = _compare_results(ga, ga)
        assert row["differences"] == "none"
        assert row["graph_type_match"] == "yes"

    def test_different_graph_type_noted(self):
        """Different graph types are flagged."""
        gem = self._make_ga(graph_type="line")
        cla = self._make_ga(graph_type="scatter")
        row = _compare_results(gem, cla)
        assert row["graph_type_match"] == "no"
        assert "graph_type" in row["differences"]

    def test_different_trend_counts_noted(self):
        """Different numbers of trends are flagged."""
        gem = self._make_ga(
            trends=[Trend(direction="increasing", description="up")])
        cla = self._make_ga()
        row = _compare_results(gem, cla)
        assert "num_trends" in row["differences"]

    def test_different_trend_directions_noted(self):
        """Different trend directions are flagged."""
        gem = self._make_ga(
            trends=[Trend(direction="increasing", description="up")])
        cla = self._make_ga(
            trends=[Trend(direction="decreasing", description="down")])
        row = _compare_results(gem, cla)
        assert row["trend_directions_match"] == "no"
        assert "trend_directions" in row["differences"]

    def test_comparison_columns_complete(self):
        """All COMPARISON_COLUMNS are present in the output."""
        ga = self._make_ga()
        row = _compare_results(ga, ga)
        for col in COMPARISON_COLUMNS:
            assert col in row, f"Missing comparison column: {col}"

    def test_different_figure_quality(self):
        """Different figure quality is flagged."""
        gem = self._make_ga(figure_quality="high")
        cla = self._make_ga(figure_quality="medium")
        row = _compare_results(gem, cla)
        assert "figure_quality" in row["differences"]

    def test_summaries_truncated(self):
        """Long summaries are truncated to 300 characters."""
        long_summary = "x" * 500
        gem = self._make_ga(summary=long_summary)
        cla = self._make_ga(summary="short")
        row = _compare_results(gem, cla)
        assert len(row["summary_gemini"]) == 300
        assert row["summary_claude"] == "short"
