"""Tests for batch_graph_analysis.py — vision-based graph analysis helpers.

Tests the non-API-dependent functions:
- collect_images: recursive image file discovery
- image_to_base64: file encoding
- media_type_for: MIME type mapping
- flatten: Pydantic model → CSV dict flattening
- Pydantic schema validation (GraphAnalysis, Axis, Trend, InflectionPoint)

API-dependent functions (analyze_image, main) are not tested here as they
require an active GEMINI_API_KEY and network access.
"""

from __future__ import annotations

import base64
import json
from pathlib import Path

import pytest

from batch_graph_analysis import (
    Axis,
    CSV_COLUMNS,
    GraphAnalysis,
    InflectionPoint,
    Trend,
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

    def test_trend_required_fields(self):
        """Trend requires direction and description; series is optional (None)."""
        t = Trend(direction="increasing", description="ADC rises")
        assert t.series is None
        assert t.direction == "increasing"

    def test_inflection_point(self):
        """InflectionPoint stores approximate x/y coordinates and description."""
        ip = InflectionPoint(
            approximate_x=60.0, approximate_y=0.001, description="divergence"
        )
        assert ip.approximate_x == 60.0

    def test_graph_analysis_minimal(self):
        """GraphAnalysis with only required fields defaults lists to empty and axes to None."""
        ga = GraphAnalysis(
            file_path="test.png", graph_type="line", summary="A line graph."
        )
        assert ga.trends == []
        assert ga.inflection_points == []
        assert ga.issues == []
        assert ga.x_axis is None

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
            summary="A scatter plot.",
        )
        assert ga.graph_title == "My Graph"
        assert len(ga.trends) == 1
        assert len(ga.inflection_points) == 1


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
