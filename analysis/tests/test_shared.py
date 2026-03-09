"""Tests for shared.py — the utility module used by all analysis scripts.

Covers:
- parse_dwi_info: DWI type extraction and filename normalisation
- extract_pvalues: regex-based p-value extraction from free text
- extract_correlations: regex-based correlation coefficient extraction
- load_graph_csv / group_by_graph_name: CSV loading and grouping
- find_latest_saved_folder / resolve_folder: CLI folder resolution
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

from conftest import GRAPH_CSV_COLUMNS, SAMPLE_GRAPH_CSV_ROWS
from shared import (
    DWI_TYPES,
    extract_correlations,
    extract_pvalues,
    find_latest_saved_folder,
    group_by_graph_name,
    load_graph_csv,
    parse_dwi_info,
    resolve_folder,
)


# ---------------------------------------------------------------------------
# parse_dwi_info
# ---------------------------------------------------------------------------

class TestParseDwiInfo:
    """Verify DWI type extraction and base name normalisation from file paths."""

    def test_standard_path(self):
        """Standard DWI type is detected and its suffix stripped from the base name."""
        dwi, name = parse_dwi_info(
            "saved_files_20260301/Standard/Feature_BoxPlots_Standard.png"
        )
        assert dwi == "Standard"
        assert name == "Feature_BoxPlots"

    def test_dncnn_path(self):
        """dnCNN DWI type is detected; the '_dnCNN' suffix is stripped."""
        dwi, name = parse_dwi_info(
            "saved_files_20260301/dnCNN/Longitudinal_dnCNN.png"
        )
        assert dwi == "dnCNN"
        assert name == "Longitudinal"

    def test_ivimnet_path(self):
        """IVIMnet type is detected even with an absolute path prefix."""
        dwi, name = parse_dwi_info(
            "/abs/saved_files_20260301/IVIMnet/Dose_vs_Diffusion_IVIMnet.png"
        )
        assert dwi == "IVIMnet"
        assert name == "Dose_vs_Diffusion"

    def test_root_path_without_dwi_subfolder(self):
        """Files directly under saved_files_* (no DWI subfolder) → 'Root'."""
        dwi, name = parse_dwi_info(
            "saved_files_20260301/some_graph.png"
        )
        assert dwi == "Root"
        assert name == "some_graph"

    def test_windows_backslash_path(self):
        """Backslash separators should be normalised."""
        dwi, name = parse_dwi_info(
            r"C:\data\saved_files_20260301\Standard\plot_Standard.png"
        )
        assert dwi == "Standard"
        assert name == "plot"

    def test_multiple_dwi_suffixes_stripped(self):
        """All DWI suffixes are removed from the base name."""
        _, name = parse_dwi_info(
            "saved_files_20260301/Standard/A_Standard_Standard.png"
        )
        assert "_Standard" not in name

    def test_no_saved_files_in_path(self):
        """Paths without 'saved_files' should still return something sensible."""
        dwi, name = parse_dwi_info("some/other/path/graph.png")
        assert dwi == "Root"
        assert name == "graph"


# ---------------------------------------------------------------------------
# extract_pvalues
# ---------------------------------------------------------------------------

class TestExtractPvalues:
    """Verify regex-based p-value extraction from free text."""

    def test_basic_equality(self):
        """'p = 0.032' with an equals sign is the most common form."""
        results = extract_pvalues("The result was p = 0.032 significant")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(0.032)

    def test_less_than_sign(self):
        """'p < 0.001' using a less-than comparator should also match."""
        results = extract_pvalues("p < 0.001 was observed")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(0.001)

    def test_scientific_notation(self):
        """Scientific notation like '3.5e-04' must be parsed correctly."""
        results = extract_pvalues("p = 3.5e-04 highly significant")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(3.5e-04)

    def test_pvalue_keyword_form(self):
        """The hyphenated form 'p-value = ...' should be recognised."""
        results = extract_pvalues("p-value = 0.045 for the test")
        assert len(results) >= 1
        assert any(abs(v - 0.045) < 1e-6 for v, _ in results)

    def test_multiple_pvalues(self):
        """Multiple p-values in one string should all be captured."""
        text = "First test p = 0.01, second test p = 0.5, third p = 1e-5"
        results = extract_pvalues(text)
        values = [v for v, _ in results]
        assert 0.01 in values or any(abs(v - 0.01) < 1e-6 for v in values)
        assert len(results) >= 3

    def test_no_pvalues(self):
        """Text without any p-value pattern should return an empty list."""
        assert extract_pvalues("No statistical results here") == []

    def test_case_insensitive(self):
        """Uppercase 'P' should match just like lowercase 'p'."""
        results = extract_pvalues("P = 0.05 uppercase match")
        assert len(results) >= 1

    def test_context_snippet_returned(self):
        """Context around the match should be included."""
        results = extract_pvalues("This is the context around p = 0.01 in the text")
        assert len(results) == 1
        _, context = results[0]
        assert "context" in context.lower() or "p" in context.lower()


# ---------------------------------------------------------------------------
# extract_correlations
# ---------------------------------------------------------------------------

class TestExtractCorrelations:
    """Verify regex-based correlation coefficient extraction."""

    def test_pearson_r(self):
        """Standard 'r = 0.85' Pearson correlation is extracted."""
        results = extract_correlations("Pearson r = 0.85 strong positive")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(0.85)

    def test_negative_r(self):
        """Negative correlation coefficients (r = -0.42) are supported."""
        results = extract_correlations("r = -0.42 negative correlation")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(-0.42)

    def test_spearman_rs(self):
        """Spearman rank correlation 'rs = 0.72' is recognised."""
        results = extract_correlations("Spearman rs = 0.72")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(0.72)

    def test_r_squared(self):
        """r² (with Unicode superscript) should be extracted."""
        results = extract_correlations("r\u00b2 = 0.45 fit quality")
        assert len(results) == 1
        assert results[0][0] == pytest.approx(0.45)

    def test_no_correlations(self):
        """Text without correlation patterns returns an empty list."""
        assert extract_correlations("No correlations reported") == []

    def test_multiple_correlations(self):
        """Multiple correlation types (r, rs, r-squared) in one string."""
        text = "r = 0.5, rs = 0.6, r² = 0.7"
        # r\u00b2 is r² — check for at least 2 matches
        results = extract_correlations(text.replace("²", "\u00b2"))
        assert len(results) >= 2


# ---------------------------------------------------------------------------
# load_graph_csv / group_by_graph_name
# ---------------------------------------------------------------------------

class TestLoadGraphCsv:
    """Verify CSV loading and graceful degradation."""

    def test_load_existing_csv(self, saved_files_with_graph_csv: Path):
        """All rows from the fixture CSV are loaded with correct column names."""
        rows = load_graph_csv(saved_files_with_graph_csv)
        assert len(rows) == len(SAMPLE_GRAPH_CSV_ROWS)
        assert "file_path" in rows[0]

    def test_load_missing_csv(self, saved_files_dir: Path):
        """Missing CSV should return an empty list, not raise."""
        rows = load_graph_csv(saved_files_dir)
        assert rows == []


class TestGroupByGraphName:
    """Verify grouping of CSV rows by normalised graph name."""

    def test_grouping_structure(self, saved_files_with_graph_csv: Path):
        """Feature_BoxPlots rows are grouped under the same key, keyed by DWI type."""
        rows = load_graph_csv(saved_files_with_graph_csv)
        groups = group_by_graph_name(rows)

        # Feature_BoxPlots should have Standard and dnCNN
        assert "Feature_BoxPlots" in groups
        fb = groups["Feature_BoxPlots"]
        assert "Standard" in fb
        assert "dnCNN" in fb

    def test_longitudinal_grouped_separately(self, saved_files_with_graph_csv: Path):
        """Longitudinal_Mean_Metrics (a different graph name) gets its own group."""
        rows = load_graph_csv(saved_files_with_graph_csv)
        groups = group_by_graph_name(rows)

        assert "Longitudinal_Mean_Metrics" in groups
        lm = groups["Longitudinal_Mean_Metrics"]
        assert "IVIMnet" in lm

    def test_empty_rows(self):
        """An empty row list produces an empty group dictionary."""
        groups = group_by_graph_name([])
        assert groups == {}


# ---------------------------------------------------------------------------
# find_latest_saved_folder / resolve_folder
# ---------------------------------------------------------------------------

class TestFindLatestSavedFolder:
    """Verify auto-detection of the most recent saved_files_* directory."""

    def test_finds_latest(self, tmp_path: Path):
        """When multiple saved_files_* dirs exist, the lexicographically last is returned."""
        # Create two timestamped folders; the later timestamp should win
        (tmp_path / "saved_files_20260101_100000").mkdir()
        (tmp_path / "saved_files_20260301_120000").mkdir()

        result = find_latest_saved_folder(str(tmp_path))
        # Lexicographic sort, reversed → 20260301 is "latest"
        assert "20260301" in str(result)

    def test_exits_when_none_found(self, tmp_path: Path):
        """sys.exit is called when no saved_files_* directories exist."""
        with pytest.raises(SystemExit):
            find_latest_saved_folder(str(tmp_path))


class TestResolveFolder:
    """Verify CLI argument resolution of the output folder."""

    def test_explicit_folder(self, saved_files_dir: Path):
        """An explicit folder path given as argv[1] is used directly."""
        result = resolve_folder(["script.py", str(saved_files_dir)])
        assert result == saved_files_dir

    def test_nonexistent_folder_exits(self, tmp_path: Path):
        """A nonexistent explicit path triggers sys.exit."""
        fake = str(tmp_path / "nonexistent")
        with pytest.raises(SystemExit):
            resolve_folder(["script.py", fake])

    def test_auto_detect_with_no_args(self, tmp_path: Path, monkeypatch):
        """With no CLI args, falls back to find_latest_saved_folder."""
        folder = tmp_path / "saved_files_20260301_120000"
        folder.mkdir()
        # Patch find_latest_saved_folder so resolve_folder's fallback
        # path returns our temp directory without needing __file__ tricks.
        monkeypatch.setattr(
            "shared.find_latest_saved_folder", lambda base_dir=None: folder
        )
        result = resolve_folder([])
        assert result == folder
