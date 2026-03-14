"""Tests for batch_graph_analysis.py — vision-based graph analysis helpers.

Tests the non-API-dependent functions:
- collect_images: recursive image file discovery
- image_to_base64: file encoding
- media_type_for: MIME type mapping
- flatten: Pydantic model → CSV dict flattening
- Pydantic schema validation (GraphAnalysis, Axis, Trend, InflectionPoint)
- Timeout configuration: REQUEST_TIMEOUT loaded from config
- _is_rate_limit_error: rate-limit error detection
- _RateLimitCoordinator: shared rate-limit backoff coordination

API-dependent functions (analyze_image, main) are not tested here as they
require an active GEMINI_API_KEY and network access.
"""

from __future__ import annotations

import asyncio
import base64
import json
from pathlib import Path


from parsers.batch_graph_analysis import (
    Axis,
    CSV_COLUMNS,
    GraphAnalysis,
    InflectionPoint,
    Outlier,
    ReferenceLine,
    REQUEST_TIMEOUT,
    StatisticalTest,
    Trend,
    _RateLimitCoordinator,
    _is_rate_limit_error,
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
