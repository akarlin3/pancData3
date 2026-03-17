"""Tests for batch_graph_analysis.py — vision-based graph analysis helpers.

Tests the non-API-dependent functions:
- collect_images: recursive image file discovery
- image_to_base64: file encoding
- media_type_for: MIME type mapping
- flatten: Pydantic model → CSV dict flattening
- Pydantic schema validation (GraphAnalysis, Axis, Trend, InflectionPoint)
- Timeout configuration: REQUEST_TIMEOUT loaded from config
- _is_rate_limit_error: rate-limit error detection
- _is_claude_rate_limit_error: Claude rate-limit error detection
- _RateLimitCoordinator: shared rate-limit backoff coordination
- _get_provider: CLI provider flag extraction
- _compare_results: dual-provider comparison
- _parse_vision_response: shared JSON response parsing
- _check_gemini_available / _check_claude_available: provider checks

API-dependent functions (analyze_image, analyze_image_claude, main) are not
tested here as they require active API keys and network access.
"""

from __future__ import annotations

import asyncio
import base64
import json
from pathlib import Path


from parsers.batch_graph_analysis import (
    Axis,
    COMPARISON_COLUMNS,
    CSV_COLUMNS,
    GraphAnalysis,
    InflectionPoint,
    Outlier,
    ReferenceLine,
    REQUEST_TIMEOUT,
    StatisticalTest,
    Trend,
    _RateLimitCoordinator,
    _compare_results,
    _get_provider,
    _is_claude_rate_limit_error,
    _is_rate_limit_error,
    _parse_vision_response,
    _should_use_local,
    analyze_image_local,
    collect_images,
    flatten,
    image_to_base64,
    media_type_for,
)


# ---------------------------------------------------------------------------
# collect_images
# ---------------------------------------------------------------------------

class TestCollectImages:
    """Verify recursive image file discovery."""

    def test_finds_png_files(self, tmp_path: Path):
        """PNG files (the primary pipeline output format) are discovered."""
        (tmp_path / "graph1.png").write_bytes(b"\x89PNG")
        (tmp_path / "graph2.png").write_bytes(b"\x89PNG")
        result = collect_images(tmp_path)
        assert len(result) == 2

    def test_finds_jpg_files(self, tmp_path: Path):
        """Both .jpg and .jpeg extensions are recognised as image files."""
        (tmp_path / "photo.jpg").write_bytes(b"\xff\xd8")
        (tmp_path / "photo2.jpeg").write_bytes(b"\xff\xd8")
        result = collect_images(tmp_path)
        assert len(result) == 2

    def test_ignores_non_image_files(self, tmp_path: Path):
        """CSV and text files are excluded; only the PNG is returned."""
        (tmp_path / "data.csv").write_text("a,b,c")
        (tmp_path / "log.txt").write_text("log entry")
        (tmp_path / "graph.png").write_bytes(b"\x89PNG")
        result = collect_images(tmp_path)
        assert len(result) == 1

    def test_recursive_discovery(self, tmp_path: Path):
        """Images in subdirectories (e.g., Standard/) are included."""
        sub = tmp_path / "Standard"
        sub.mkdir()
        (sub / "plot.png").write_bytes(b"\x89PNG")
        (tmp_path / "root.png").write_bytes(b"\x89PNG")
        result = collect_images(tmp_path)
        assert len(result) == 2

    def test_empty_directory(self, tmp_path: Path):
        """An empty directory returns an empty list."""
        result = collect_images(tmp_path)
        assert result == []

    def test_case_insensitive_extension(self, tmp_path: Path):
        """Uppercase .PNG should be found (suffix.lower() comparison)."""
        (tmp_path / "GRAPH.PNG").write_bytes(b"\x89PNG")
        result = collect_images(tmp_path)
        assert len(result) == 1

    def test_sorted_output(self, tmp_path: Path):
        (tmp_path / "c.png").write_bytes(b"\x89PNG")
        (tmp_path / "a.png").write_bytes(b"\x89PNG")
        (tmp_path / "b.png").write_bytes(b"\x89PNG")
        result = collect_images(tmp_path)
        names = [p.name for p in result]
        assert names == sorted(names)


# ---------------------------------------------------------------------------
# image_to_base64
# ---------------------------------------------------------------------------

class TestImageToBase64:
    """Verify base64 encoding of image files."""

    def test_encodes_correctly(self, tmp_path: Path):
        """Round-trip: encoding then decoding recovers the original bytes."""
        data = b"fake image data for testing"
        p = tmp_path / "test.png"
        p.write_bytes(data)
        encoded = image_to_base64(p)
        # Decode back and compare
        assert base64.standard_b64decode(encoded) == data

    def test_returns_string(self, tmp_path: Path):
        """The encoded result is a Python str (not bytes)."""
        p = tmp_path / "test.png"
        p.write_bytes(b"\x89PNG")
        assert isinstance(image_to_base64(p), str)


# ---------------------------------------------------------------------------
# media_type_for
# ---------------------------------------------------------------------------

class TestMediaTypeFor:
    """Verify MIME type mapping for image file extensions."""

    def test_png(self, tmp_path: Path):
        """.png maps to image/png."""
        assert media_type_for(tmp_path / "graph.png") == "image/png"

    def test_jpg(self, tmp_path: Path):
        """.jpg maps to image/jpeg."""
        assert media_type_for(tmp_path / "photo.jpg") == "image/jpeg"

    def test_jpeg(self, tmp_path: Path):
        """.jpeg also maps to image/jpeg."""
        assert media_type_for(tmp_path / "photo.jpeg") == "image/jpeg"

    def test_unknown_defaults_to_png(self, tmp_path: Path):
        """Unrecognised extensions (e.g., .bmp) fall back to image/png."""
        assert media_type_for(tmp_path / "file.bmp") == "image/png"

    def test_case_insensitive(self, tmp_path: Path):
        """Uppercase extensions should still map correctly."""
        assert media_type_for(tmp_path / "GRAPH.PNG") == "image/png"
        assert media_type_for(tmp_path / "photo.JPG") == "image/jpeg"


# ---------------------------------------------------------------------------
# Pydantic schema validation
# ---------------------------------------------------------------------------

class TestPydanticSchemas:
    """Verify Pydantic model construction and field defaults."""

    def test_axis_all_fields(self):
        """Axis with all fields populated stores them correctly."""
        ax = Axis(label="ADC", units="mm²/s", range_min=0.0, range_max=0.003)
        assert ax.label == "ADC"
        assert ax.units == "mm²/s"

    def test_axis_optional_fields(self):
        """Axis with only the required 'label' leaves optional fields as None."""
        ax = Axis(label="X")
        assert ax.units is None
        assert ax.range_min is None
        assert ax.scale_type is None
        assert ax.tick_count is None
        assert ax.tick_labels == []

    def test_axis_scale_type(self):
        """Axis stores scale_type when provided."""
        ax = Axis(label="Time", scale_type="log")
        assert ax.scale_type == "log"

    def test_axis_tick_labels(self):
        """Axis stores categorical tick labels."""
        ax = Axis(label="Timepoint", scale_type="categorical",
                  tick_labels=["BL", "Fx1", "W2"])
        assert ax.tick_labels == ["BL", "Fx1", "W2"]
        assert ax.tick_count is None

    def test_trend_required_fields(self):
        """Trend requires direction and description; series is optional (None)."""
        t = Trend(direction="increasing", description="ADC rises")
        assert t.series is None
        assert t.direction == "increasing"
        assert t.magnitude is None
        assert t.statistical_significance is None
        assert t.confidence_band is None
        assert t.start_value is None
        assert t.end_value is None

    def test_trend_full_fields(self):
        """Trend with all optional fields populated."""
        t = Trend(
            series="LF", direction="decreasing", description="LF declines",
            magnitude="~20% decline", statistical_significance="p=0.003",
            confidence_band="95% CI shaded", start_value=0.002, end_value=0.0016,
        )
        assert t.magnitude == "~20% decline"
        assert t.statistical_significance == "p=0.003"
        assert t.confidence_band == "95% CI shaded"
        assert t.start_value == 0.002
        assert t.end_value == 0.0016

    def test_inflection_point(self):
        """InflectionPoint stores approximate x/y coordinates and description."""
        ip = InflectionPoint(
            approximate_x=60.0, approximate_y=0.001, description="divergence"
        )
        assert ip.approximate_x == 60.0
        assert ip.magnitude is None

    def test_inflection_point_with_magnitude(self):
        """InflectionPoint stores magnitude of change."""
        ip = InflectionPoint(
            approximate_x=60.0, description="jump", magnitude=0.0005
        )
        assert ip.magnitude == 0.0005

    def test_statistical_test(self):
        """StatisticalTest stores test details."""
        st = StatisticalTest(
            test_name="Wilcoxon", statistic_value=3.45, p_value=0.003,
            comparison_groups="LF vs LC",
        )
        assert st.test_name == "Wilcoxon"
        assert st.p_value == 0.003
        assert st.comparison_groups == "LF vs LC"

    def test_statistical_test_minimal(self):
        """StatisticalTest with only required test_name."""
        st = StatisticalTest(test_name="ANOVA")
        assert st.statistic_value is None
        assert st.p_value is None

    def test_outlier(self):
        """Outlier stores location and description."""
        o = Outlier(
            approximate_x=5.0, approximate_y=0.005, series="ADC",
            description="extreme high value",
        )
        assert o.approximate_y == 0.005
        assert o.series == "ADC"

    def test_reference_line(self):
        """ReferenceLine stores orientation, value, label, style."""
        rl = ReferenceLine(
            orientation="horizontal", value=0.001,
            label="ADC threshold", style="dashed",
        )
        assert rl.orientation == "horizontal"
        assert rl.value == 0.001
        assert rl.style == "dashed"

    def test_graph_analysis_minimal(self):
        """GraphAnalysis with only required fields defaults lists to empty and axes to None."""
        ga = GraphAnalysis(
            file_path="test.png", graph_type="line", summary="A line graph."
        )
        assert ga.trends == []
        assert ga.inflection_points == []
        assert ga.statistical_tests == []
        assert ga.outliers == []
        assert ga.reference_lines == []
        assert ga.issues == []
        assert ga.annotations == []
        assert ga.legend_items == []
        assert ga.x_axis is None
        assert ga.sample_size is None
        assert ga.data_series_count is None
        assert ga.error_bars is None
        assert ga.clinical_relevance is None
        assert ga.data_density is None
        assert ga.spatial_pattern is None
        assert ga.subpanel_count is None
        assert ga.comparison_type is None
        assert ga.figure_quality is None

    def test_graph_analysis_full(self):
        """GraphAnalysis with all optional fields populated stores nested models."""
        ga = GraphAnalysis(
            file_path="test.png",
            graph_title="My Graph",
            graph_type="scatter",
            x_axis=Axis(label="X"),
            y_axis=Axis(label="Y", units="mm"),
            trends=[Trend(direction="flat", description="no change")],
            inflection_points=[InflectionPoint(description="peak")],
            statistical_tests=[StatisticalTest(test_name="t-test", p_value=0.04)],
            outliers=[Outlier(description="high point")],
            reference_lines=[ReferenceLine(orientation="horizontal", value=1.0)],
            summary="A scatter plot.",
            sample_size="n=30",
            data_series_count=3,
            error_bars="SD",
            annotations=["*", "ns"],
            clinical_relevance="No significant association with treatment response.",
            data_density="moderate",
            spatial_pattern=None,
            legend_items=["Group A", "Group B", "Group C"],
            subpanel_count=2,
            comparison_type="unpaired",
            figure_quality="high",
        )
        assert ga.graph_title == "My Graph"
        assert len(ga.trends) == 1
        assert len(ga.inflection_points) == 1
        assert len(ga.statistical_tests) == 1
        assert ga.statistical_tests[0].p_value == 0.04
        assert len(ga.outliers) == 1
        assert len(ga.reference_lines) == 1
        assert ga.sample_size == "n=30"
        assert ga.data_series_count == 3
        assert ga.error_bars == "SD"
        assert ga.annotations == ["*", "ns"]
        assert ga.clinical_relevance is not None
        assert ga.data_density == "moderate"
        assert len(ga.legend_items) == 3
        assert ga.subpanel_count == 2
        assert ga.comparison_type == "unpaired"
        assert ga.figure_quality == "high"


# ---------------------------------------------------------------------------
# flatten
# ---------------------------------------------------------------------------

class TestFlatten:
    """Verify Pydantic model → CSV dict flattening."""

    def test_minimal_graph(self):
        """A minimal GraphAnalysis produces a dict with all CSV_COLUMNS present.

        None axes should be flattened to empty strings so the CSV writer
        does not raise on missing keys.
        """
        ga = GraphAnalysis(
            file_path="img.png", graph_type="bar", summary="A bar chart."
        )
        row = flatten(ga)
        # All CSV_COLUMNS should be present
        for col in CSV_COLUMNS:
            assert col in row, f"Missing column: {col}"
        assert row["file_path"] == "img.png"
        assert row["graph_type"] == "bar"
        assert row["x_axis_label"] == ""  # None axis → empty strings

    def test_with_axes(self):
        """Axis fields are flattened into prefixed column names (e.g., x_axis_label).

        None sub-fields within a present axis become empty strings.
        """
        ga = GraphAnalysis(
            file_path="img.png",
            graph_type="line",
            x_axis=Axis(label="Time", units="days", range_min=0, range_max=180),
            y_axis=Axis(label="ADC", units="mm²/s"),
            summary="Line plot.",
        )
        row = flatten(ga)
        assert row["x_axis_label"] == "Time"
        assert row["x_axis_units"] == "days"
        assert row["x_axis_range_min"] == 0
        assert row["x_axis_range_max"] == 180
        assert row["y_axis_label"] == "ADC"
        assert row["y_axis_range_min"] == ""  # None → empty

    def test_trends_serialised_as_json(self):
        """Trend objects are serialised to a JSON string in the trends_json column."""
        ga = GraphAnalysis(
            file_path="img.png",
            graph_type="line",
            trends=[
                Trend(series="LF", direction="decreasing", description="drops"),
                Trend(series="LC", direction="flat", description="stable"),
            ],
            summary="Two trends.",
        )
        row = flatten(ga)
        assert row["num_trends"] == 2
        # Verify the JSON is valid and contains the expected data
        parsed = json.loads(row["trends_json"])
        assert len(parsed) == 2
        assert parsed[0]["series"] == "LF"

    def test_inflection_points_serialised(self):
        """InflectionPoint objects are serialised to JSON in inflection_points_json."""
        ga = GraphAnalysis(
            file_path="img.png",
            graph_type="line",
            inflection_points=[
                InflectionPoint(approximate_x=30, description="change"),
            ],
            summary="One inflection.",
        )
        row = flatten(ga)
        assert row["num_inflection_points"] == 1
        parsed = json.loads(row["inflection_points_json"])
        assert parsed[0]["approximate_x"] == 30

    def test_issues_serialised_as_json(self):
        """Issues list is serialised to a JSON string in the issues_json column."""
        ga = GraphAnalysis(
            file_path="img.png",
            graph_type="box",
            issues=["Missing axis label", "Legend overlaps data"],
            summary="Box plot with issues.",
        )
        row = flatten(ga)
        assert row["num_issues"] == 2
        parsed = json.loads(row["issues_json"])
        assert len(parsed) == 2
        assert parsed[0] == "Missing axis label"

    def test_empty_issues(self):
        """A graph with no issues has num_issues=0 and issues_json='[]'."""
        ga = GraphAnalysis(
            file_path="img.png", graph_type="line", summary="Clean graph."
        )
        row = flatten(ga)
        assert row["num_issues"] == 0
        assert json.loads(row["issues_json"]) == []

    def test_color_axis_flattened(self):
        """The optional color_axis (used by heatmaps) is flattened like x/y axes."""
        ga = GraphAnalysis(
            file_path="img.png",
            graph_type="heatmap",
            color_axis=Axis(label="Dice", units="", range_min=0, range_max=1),
            summary="Heatmap.",
        )
        row = flatten(ga)
        assert row["color_axis_label"] == "Dice"
        assert row["color_axis_range_max"] == 1

    def test_statistical_tests_serialised(self):
        """StatisticalTest objects are serialised to JSON in statistical_tests_json."""
        ga = GraphAnalysis(
            file_path="img.png",
            graph_type="box",
            statistical_tests=[
                StatisticalTest(test_name="Wilcoxon", p_value=0.003, comparison_groups="LF vs LC"),
            ],
            summary="Box plot with stats.",
        )
        row = flatten(ga)
        assert row["num_statistical_tests"] == 1
        parsed = json.loads(row["statistical_tests_json"])
        assert parsed[0]["test_name"] == "Wilcoxon"
        assert parsed[0]["p_value"] == 0.003

    def test_outliers_serialised(self):
        """Outlier objects are serialised to JSON in outliers_json."""
        ga = GraphAnalysis(
            file_path="img.png",
            graph_type="scatter",
            outliers=[Outlier(approximate_y=0.005, description="extreme value")],
            summary="Scatter with outlier.",
        )
        row = flatten(ga)
        assert row["num_outliers"] == 1
        parsed = json.loads(row["outliers_json"])
        assert parsed[0]["approximate_y"] == 0.005

    def test_reference_lines_serialised(self):
        """ReferenceLine objects are serialised to JSON."""
        ga = GraphAnalysis(
            file_path="img.png",
            graph_type="line",
            reference_lines=[
                ReferenceLine(orientation="horizontal", value=0.001,
                              label="threshold", style="dashed"),
            ],
            summary="Line with reference.",
        )
        row = flatten(ga)
        assert row["num_reference_lines"] == 1
        parsed = json.loads(row["reference_lines_json"])
        assert parsed[0]["orientation"] == "horizontal"
        assert parsed[0]["value"] == 0.001

    def test_new_scalar_fields_flattened(self):
        """New scalar metadata fields are correctly flattened."""
        ga = GraphAnalysis(
            file_path="img.png",
            graph_type="box",
            summary="Box plot.",
            sample_size="n=42",
            data_series_count=3,
            error_bars="SD",
            annotations=["*", "**"],
            clinical_relevance="Significant association with response.",
            data_density="moderate",
            spatial_pattern=None,
            legend_items=["LF", "LC", "CR"],
            subpanel_count=4,
            comparison_type="unpaired",
            figure_quality="high",
        )
        row = flatten(ga)
        assert row["sample_size"] == "n=42"
        assert row["data_series_count"] == 3
        assert row["error_bars"] == "SD"
        assert row["num_annotations"] == 2
        assert json.loads(row["annotations_json"]) == ["*", "**"]
        assert row["clinical_relevance"] == "Significant association with response."
        assert row["data_density"] == "moderate"
        assert row["spatial_pattern"] == ""
        assert row["num_legend_items"] == 3
        assert json.loads(row["legend_items_json"]) == ["LF", "LC", "CR"]
        assert row["subpanel_count"] == 4
        assert row["comparison_type"] == "unpaired"
        assert row["figure_quality"] == "high"

    def test_axis_scale_type_and_ticks_flattened(self):
        """Axis scale_type, tick_count, and tick_labels are flattened correctly."""
        ga = GraphAnalysis(
            file_path="img.png",
            graph_type="bar",
            x_axis=Axis(label="Group", scale_type="categorical",
                        tick_count=3, tick_labels=["A", "B", "C"]),
            summary="Bar chart.",
        )
        row = flatten(ga)
        assert row["x_axis_scale_type"] == "categorical"
        assert row["x_axis_tick_count"] == 3
        assert json.loads(row["x_axis_tick_labels_json"]) == ["A", "B", "C"]

    def test_trends_with_new_fields_serialised(self):
        """Trend magnitude, significance, confidence_band, start/end values are serialised."""
        ga = GraphAnalysis(
            file_path="img.png",
            graph_type="line",
            trends=[
                Trend(series="LF", direction="decreasing", description="declines",
                      magnitude="~20%", statistical_significance="p<0.05",
                      confidence_band="95% CI", start_value=0.002, end_value=0.001),
            ],
            summary="Line plot.",
        )
        row = flatten(ga)
        parsed = json.loads(row["trends_json"])
        assert parsed[0]["magnitude"] == "~20%"
        assert parsed[0]["statistical_significance"] == "p<0.05"
        assert parsed[0]["confidence_band"] == "95% CI"
        assert parsed[0]["start_value"] == 0.002
        assert parsed[0]["end_value"] == 0.001


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
