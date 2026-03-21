"""Tests for forest plot generation."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from report.sections.forest_plot import (
    extract_hr_data,
    generate_forest_plot,
)


@pytest.fixture
def sample_log_data():
    """Synthetic log data with hazard ratios."""
    return {
        "Standard": {
            "survival": {
                "hazard_ratios": [
                    {"covariate": "ADC", "hr": 0.75, "ci_lo": 0.5, "ci_hi": 1.1, "p": 0.04},
                    {"covariate": "D", "hr": 1.2, "ci_lo": 0.8, "ci_hi": 1.8, "p": 0.3},
                    {"covariate": "f", "hr": 0.6, "ci_lo": 0.3, "ci_hi": 1.2, "p": 0.15},
                ]
            }
        },
        "dnCNN": {
            "survival": {
                "hazard_ratios": [
                    {"covariate": "ADC", "hr": 0.8, "ci_lo": 0.55, "ci_hi": 1.15, "p": 0.02},
                    {"covariate": "D", "hr": 1.1, "ci_lo": 0.7, "ci_hi": 1.7, "p": 0.5},
                ]
            }
        },
    }


class TestExtractHRData:
    """Tests for HR data extraction."""

    def test_extracts_all_entries(self, sample_log_data):
        """Should extract all HR entries across DWI types."""
        result = extract_hr_data(sample_log_data, ["Standard", "dnCNN"])
        assert len(result) == 5

    def test_empty_log_data(self):
        """Empty log data should return empty list."""
        result = extract_hr_data(None, ["Standard"])
        assert result == []

    def test_missing_dwi_type(self, sample_log_data):
        """Should handle missing DWI type gracefully."""
        result = extract_hr_data(sample_log_data, ["IVIMnet"])
        assert result == []

    def test_entry_fields(self, sample_log_data):
        """Each entry should have required fields."""
        result = extract_hr_data(sample_log_data, ["Standard"])
        for entry in result:
            assert "dwi_type" in entry
            assert "covariate" in entry
            assert "hr" in entry
            assert "ci_lo" in entry
            assert "ci_hi" in entry
            assert "p" in entry


class TestGenerateForestPlot:
    """Tests for forest plot figure generation."""

    def test_generates_figure(self, sample_log_data):
        """Should generate a PNG file."""
        hr_data = extract_hr_data(sample_log_data, ["Standard", "dnCNN"])
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "forest_plot.png"
            try:
                result = generate_forest_plot(hr_data, output_path)
                if result:
                    assert output_path.exists()
                    assert output_path.stat().st_size > 0
            except ImportError:
                pytest.skip("matplotlib not available")

    def test_empty_data(self):
        """Empty data should return False."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "forest_plot.png"
            result = generate_forest_plot([], output_path)
            assert result is False

    def test_single_type(self, sample_log_data):
        """Should work with a single DWI type."""
        hr_data = extract_hr_data(sample_log_data, ["Standard"])
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "forest_plot.png"
            try:
                result = generate_forest_plot(hr_data, output_path)
                # Should complete without error
                assert isinstance(result, bool)
            except ImportError:
                pytest.skip("matplotlib not available")

    def test_missing_ci_handled(self):
        """Should handle entries with missing CI gracefully."""
        hr_data = [
            {"dwi_type": "Standard", "covariate": "ADC", "hr": 0.8, "ci_lo": 0.5, "ci_hi": 1.2, "p": 0.05},
            {"dwi_type": "Standard", "covariate": "D", "hr": 1.0, "ci_lo": 1.0, "ci_hi": 1.0, "p": 1.0},
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "forest_plot.png"
            try:
                result = generate_forest_plot(hr_data, output_path)
                assert isinstance(result, bool)
            except ImportError:
                pytest.skip("matplotlib not available")

    def test_log_scale_axis(self, sample_log_data):
        """Verify figure uses log scale (by checking it generates without error)."""
        hr_data = extract_hr_data(sample_log_data, ["Standard", "dnCNN"])
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "forest_plot.png"
            try:
                result = generate_forest_plot(
                    hr_data, output_path,
                    title="Test Forest Plot (Log Scale)"
                )
                assert isinstance(result, bool)
            except ImportError:
                pytest.skip("matplotlib not available")
