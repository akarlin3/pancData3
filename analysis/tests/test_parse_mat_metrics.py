"""Unit tests for parse_mat_metrics.py.

Tests MAT file parsing with mocked scipy.io.loadmat to avoid real .mat file
dependencies, plus CLI entry point and graceful degradation when scipy is
unavailable.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

# Ensure analysis/ is on sys.path
ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

import parse_mat_metrics


# Uses saved_files_dir fixture from conftest.py


# ── parse_mat_files_for_dwi ──────────────────────────────────────────────


class TestParseMatFilesForDwi:
    def test_empty_folder_returns_empty_dicts(self, saved_files_dir):
        """When no MAT files exist, all keys should be empty dicts."""
        result = parse_mat_metrics.parse_mat_files_for_dwi(saved_files_dir, "Standard")
        assert result["core_method"] == {}
        assert result["dosimetry"] == {}
        assert result["longitudinal"] == {}

    def test_returns_three_keys(self, saved_files_dir):
        """Output dict should always have exactly these three keys."""
        result = parse_mat_metrics.parse_mat_files_for_dwi(saved_files_dir, "Standard")
        assert set(result.keys()) == {"core_method", "dosimetry", "longitudinal"}

    @pytest.mark.skipif(
        parse_mat_metrics.scipy_io is None,
        reason="scipy not installed"
    )
    def test_core_method_parsing(self, saved_files_dir):
        """Mock the core comparison MAT file and verify parsing."""
        import numpy as np

        # Create a mock compare_results struct
        methods = np.array(["adc_threshold", "otsu", "gmm"])
        dice_matrix = np.array([[1.0, 0.8, 0.7], [0.8, 1.0, 0.75], [0.7, 0.75, 1.0]])

        compare_results = SimpleNamespace(
            method_names=methods,
            mean_dice_matrix=dice_matrix,
        )
        mat_data = {"compare_results": compare_results}

        with patch.object(parse_mat_metrics.scipy_io, "loadmat", return_value=mat_data):
            # Create the expected file so the existence check passes
            core_file = saved_files_dir / "Standard" / "compare_core_results_Standard.mat"
            core_file.touch()

            result = parse_mat_metrics.parse_mat_files_for_dwi(saved_files_dir, "Standard")

        assert "methods" in result["core_method"]
        assert len(result["core_method"]["methods"]) == 3
        assert "mean_dice_matrix" in result["core_method"]

    @pytest.mark.skipif(
        parse_mat_metrics.scipy_io is None,
        reason="scipy not installed"
    )
    def test_dosimetry_parsing(self, saved_files_dir):
        """Mock dosimetry MAT file and verify mean computation."""
        import numpy as np

        mat_data = {
            "d95_adc_sub": np.array([40.0, 50.0, 60.0]),
            "v50_adc_sub": np.array([0.8, 0.9, 0.85]),
            "d95_d_sub": np.array([35.0, 45.0]),
            "v50_d_sub": np.array([0.7, 0.75]),
        }

        with patch.object(parse_mat_metrics.scipy_io, "loadmat", return_value=mat_data):
            dosimetry_file = saved_files_dir / "Standard" / "metrics_dosimetry_results_Standard.mat"
            dosimetry_file.touch()

            result = parse_mat_metrics.parse_mat_files_for_dwi(saved_files_dir, "Standard")

        assert result["dosimetry"]["d95_adc_mean"] == pytest.approx(50.0)
        assert result["dosimetry"]["v50_adc_mean"] == pytest.approx(0.85)

    @pytest.mark.skipif(
        parse_mat_metrics.scipy_io is None,
        reason="scipy not installed"
    )
    def test_longitudinal_parsing(self, saved_files_dir):
        """Mock summary metrics MAT file to extract cohort dimensions."""
        import numpy as np

        summary = SimpleNamespace(ADC_abs=np.zeros((30, 5)))
        mat_data = {"summary_metrics": summary}

        with patch.object(parse_mat_metrics.scipy_io, "loadmat", return_value=mat_data):
            summary_file = saved_files_dir / "Standard" / "summary_metrics_Standard.mat"
            summary_file.touch()

            result = parse_mat_metrics.parse_mat_files_for_dwi(saved_files_dir, "Standard")

        assert result["longitudinal"]["num_patients"] == 30
        assert result["longitudinal"]["num_timepoints"] == 5

    def test_graceful_on_corrupt_file(self, saved_files_dir):
        """If scipy raises an exception parsing a file, it should not crash."""
        if parse_mat_metrics.scipy_io is None:
            pytest.skip("scipy not installed")

        with patch.object(parse_mat_metrics.scipy_io, "loadmat", side_effect=Exception("corrupt")):
            core_file = saved_files_dir / "Standard" / "compare_core_results_Standard.mat"
            core_file.touch()

            result = parse_mat_metrics.parse_mat_files_for_dwi(saved_files_dir, "Standard")

        # Should return empty dicts, not crash
        assert result["core_method"] == {}


# ── Without scipy ────────────────────────────────────────────────────────


class TestWithoutScipy:
    def test_returns_empty_when_scipy_missing(self, saved_files_dir):
        """When scipy_io is None, function should return empty data."""
        original = parse_mat_metrics.scipy_io
        parse_mat_metrics.scipy_io = None
        try:
            result = parse_mat_metrics.parse_mat_files_for_dwi(saved_files_dir, "Standard")
            assert result["core_method"] == {}
            assert result["dosimetry"] == {}
            assert result["longitudinal"] == {}
        finally:
            parse_mat_metrics.scipy_io = original


# ── CLI main ─────────────────────────────────────────────────────────────


class TestMain:
    def test_nonexistent_folder(self, capsys):
        """main() should print error for nonexistent folder."""
        with patch("sys.argv", ["parse_mat_metrics.py", "/nonexistent/path"]):
            parse_mat_metrics.main()
        captured = capsys.readouterr()
        assert "not found" in captured.out.lower() or "Folder not found" in captured.out

    def test_no_dwi_folders(self, tmp_path, capsys):
        """main() should report no DWI type folders when directory is empty."""
        empty_dir = tmp_path / "empty_saved_files"
        empty_dir.mkdir()
        with patch("sys.argv", ["parse_mat_metrics.py", str(empty_dir)]):
            parse_mat_metrics.main()
        captured = capsys.readouterr()
        assert "No DWI type folders" in captured.out

    @pytest.mark.skipif(
        parse_mat_metrics.scipy_io is None,
        reason="scipy not installed"
    )
    def test_writes_json_output(self, saved_files_dir):
        """main() should write parsed_mat_metrics_{dwi}.json files."""
        with patch("sys.argv", ["parse_mat_metrics.py", str(saved_files_dir)]):
            parse_mat_metrics.main()

        # Should write one JSON per DWI type subfolder
        for dwi in ("Standard", "dnCNN", "IVIMnet"):
            json_file = saved_files_dir / f"parsed_mat_metrics_{dwi}.json"
            assert json_file.exists()
            data = json.loads(json_file.read_text())
            assert "core_method" in data
