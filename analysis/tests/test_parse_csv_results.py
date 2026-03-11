"""Tests for parse_csv_results.py — pipeline CSV export parsing.

Covers:
- _read_csv: graceful degradation for missing files
- parse_significant_metrics: loading Significant_LF_Metrics.csv per DWI type
- parse_fdr_global: loading FDR_Sig_Global.csv per DWI type
- cross_reference_significance: cross-DWI consistency analysis
- parse_all_csvs: top-level orchestrator
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest  # type: ignore

from parse_csv_results import (  # type: ignore
    _read_csv,
    cross_reference_significance,
    parse_all_csvs,
    parse_fdr_global,
    parse_significant_metrics,
)


# ---------------------------------------------------------------------------
# _read_csv
# ---------------------------------------------------------------------------

class TestReadCsv:
    """Verify CSV reading with graceful missing-file handling."""

    def test_reads_existing_csv(self, tmp_path: Path):
        """A valid CSV file is loaded into a list of dicts with correct values."""
        p = tmp_path / "data.csv"
        with open(p, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["A", "B"])
            w.writeheader()
            w.writerow({"A": "1", "B": "2"})
        rows = _read_csv(p)
        assert len(rows) == 1
        assert rows[0]["A"] == "1"

    def test_returns_empty_for_missing(self, tmp_path: Path):
        """A missing CSV file returns an empty list instead of raising."""
        rows = _read_csv(tmp_path / "nonexistent.csv")
        assert rows == []


# ---------------------------------------------------------------------------
# parse_significant_metrics
# ---------------------------------------------------------------------------

class TestParseSignificantMetrics:
    """Verify loading of Significant_LF_Metrics.csv from DWI subfolders."""

    def test_loads_per_dwi_type(self, saved_files_with_csvs: Path):
        """CSVs from Standard (3 rows) and dnCNN (1 row) are loaded separately."""
        result = parse_significant_metrics(saved_files_with_csvs)
        assert "Standard" in result
        assert len(result["Standard"]) == 3  # 3 rows in the fixture
        assert "dnCNN" in result
        assert len(result["dnCNN"]) == 1

    def test_skips_missing_files(self, saved_files_dir: Path):
        """DWI subfolders without CSV files are not included."""
        result = parse_significant_metrics(saved_files_dir)
        assert result == {}

    def test_ivimnet_not_present(self, saved_files_with_csvs: Path):
        """IVIMnet has no CSV in the fixture → not in results."""
        result = parse_significant_metrics(saved_files_with_csvs)
        assert "IVIMnet" not in result


# ---------------------------------------------------------------------------
# parse_fdr_global
# ---------------------------------------------------------------------------

class TestParseFdrGlobal:
    """Verify loading of FDR_Sig_Global.csv from DWI subfolders."""

    def test_loads_when_present(self, saved_files_dir: Path):
        """An FDR_Sig_Global.csv created on-the-fly is correctly loaded."""
        # Create an FDR CSV in Standard
        csv_path = saved_files_dir / "Standard" / "FDR_Sig_Global.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["Metric", "p_raw", "p_fdr", "Significant"])
            w.writeheader()
            w.writerow({"Metric": "adc", "p_raw": "0.001", "p_fdr": "0.01", "Significant": "yes"})
        result = parse_fdr_global(saved_files_dir)
        assert "Standard" in result
        assert len(result["Standard"]) == 1

    def test_empty_when_no_files(self, saved_files_dir: Path):
        """No FDR CSV files in any DWI subfolder returns an empty dict."""
        result = parse_fdr_global(saved_files_dir)
        assert result == {}


# ---------------------------------------------------------------------------
# cross_reference_significance
# ---------------------------------------------------------------------------

class TestCrossReferenceSignificance:
    """Verify cross-DWI consistency analysis of significant metrics."""

    def test_consistent_metric(self):
        """A metric significant in all DWI types is flagged as consistent."""
        sig = {
            "Standard": [{"Metric": "mean_adc", "Timepoint": "BL"}],
            "dnCNN": [{"Metric": "mean_adc", "Timepoint": "BL"}],
            "IVIMnet": [{"Metric": "mean_adc", "Timepoint": "BL"}],
        }
        result = cross_reference_significance(sig)
        assert len(result) == 1
        assert result[0]["consistent"] is True
        assert result[0]["not_significant_in"] == []

    def test_inconsistent_metric(self):
        """A metric significant in only some DWI types → inconsistent."""
        sig = {
            "Standard": [{"Metric": "mean_d", "Timepoint": "BL"}],
            # dnCNN and IVIMnet do not have this metric
        }
        result = cross_reference_significance(sig)
        assert len(result) == 1
        assert result[0]["consistent"] is False
        assert "dnCNN" in result[0]["not_significant_in"]
        assert "IVIMnet" in result[0]["not_significant_in"]

    def test_multiple_metrics_mixed_consistency(self):
        """Mix of consistent and inconsistent metrics."""
        sig = {
            "Standard": [
                {"Metric": "mean_adc", "Timepoint": "BL"},
                {"Metric": "mean_d", "Timepoint": "W2"},
            ],
            "dnCNN": [
                {"Metric": "mean_adc", "Timepoint": "BL"},
            ],
            "IVIMnet": [
                {"Metric": "mean_adc", "Timepoint": "BL"},
            ],
        }
        result = cross_reference_significance(sig)
        # Timepoints are normalised to lowercase, so "BL" → "bl", "W2" → "w2".
        by_key = {f"{r['metric']}@{r['timepoint']}": r for r in result}

        assert by_key["mean_adc@bl"]["consistent"] is True
        assert by_key["mean_d@w2"]["consistent"] is False

    def test_empty_input(self):
        """No DWI types at all yields an empty cross-reference list."""
        result = cross_reference_significance({})
        assert result == []

    def test_timepoint_aware_keys(self):
        """Same metric at different timepoints should be separate entries."""
        sig = {
            "Standard": [
                {"Metric": "adc", "Timepoint": "BL"},
                {"Metric": "adc", "Timepoint": "W2"},
            ],
        }
        result = cross_reference_significance(sig)
        assert len(result) == 2


# ---------------------------------------------------------------------------
# parse_all_csvs (integration)
# ---------------------------------------------------------------------------

class TestParseAllCsvs:
    """Integration test for the top-level CSV parsing orchestrator."""

    def test_returns_all_keys(self, saved_files_with_csvs: Path):
        """The orchestrator returns all three expected top-level keys."""
        result = parse_all_csvs(saved_files_with_csvs)
        assert "significant_metrics" in result
        assert "fdr_global" in result
        assert "cross_reference" in result

    def test_cross_reference_populated(self, saved_files_with_csvs: Path):
        """Cross-reference identifies at least 3 metric/timepoint combinations.

        Given the fixture data (Standard has 3 entries, dnCNN has 1),
        no metric is fully consistent across all three DWI types (IVIMnet
        has no CSV), so at least one entry should be marked inconsistent.
        """
        result = parse_all_csvs(saved_files_with_csvs)
        cr = result["cross_reference"]
        # Standard has mean_adc@BL, mean_d@BL, mean_adc@W2
        # dnCNN has mean_adc@BL only
        # → mean_adc@BL is in both (but not IVIMnet) → inconsistent
        # → mean_d@BL is only in Standard → inconsistent
        # → mean_adc@W2 is only in Standard → inconsistent
        assert len(cr) >= 3
        inconsistent = [c for c in cr if not c["consistent"]]
        assert len(inconsistent) >= 1

    def test_empty_folder(self, saved_files_dir: Path):
        """A folder with no CSV exports returns empty containers for all keys."""
        result = parse_all_csvs(saved_files_dir)
        assert result["significant_metrics"] == {}
        assert result["cross_reference"] == []
