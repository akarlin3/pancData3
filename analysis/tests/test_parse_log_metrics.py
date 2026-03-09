"""Tests for parse_log_metrics.py — MATLAB log file parsing.

Covers:
- parse_stats_comparisons: GLME interaction p-values, FDR timepoints, exclusions
- parse_stats_predictive: Elastic Net features, ROC/AUC, Firth refit
- parse_survival: hazard ratios, global LRT, IPCW weights
- parse_baseline: outlier flags, total outliers, baseline exclusion
- parse_all_logs: orchestrator that reads from DWI type subfolders
- _read_log: graceful file-missing handling
"""

from __future__ import annotations

from pathlib import Path

import pytest

from parse_log_metrics import (
    _read_log,
    parse_all_logs,
    parse_baseline,
    parse_stats_comparisons,
    parse_stats_predictive,
    parse_survival,
)


# ---------------------------------------------------------------------------
# parse_stats_comparisons
# ---------------------------------------------------------------------------

class TestParseStatsComparisons:
    """Verify extraction of GLME interactions, details, FDR, and exclusions."""

    def test_glme_interactions(self):
        """GLME interaction p-values are extracted from 'Interaction P-Value for ...: <value>' lines.

        Tests both decimal (0.023) and scientific notation (1.5e-04) formats.
        """
        text = (
            "Interaction P-Value for ADC: 0.023\n"
            "Interaction P-Value for D*: 1.5e-04\n"
        )
        result = parse_stats_comparisons(text)
        assert result["glme_interactions"] == pytest.approx([0.023, 1.5e-04])

    def test_glme_details(self):
        """Per-metric Wilcoxon detail lines ('metric: p=..., adj_alpha=...') are parsed."""
        text = "mean_adc: p=0.015, adj_alpha=0.025\n"
        result = parse_stats_comparisons(text)
        assert len(result["glme_details"]) == 1
        d = result["glme_details"][0]
        assert d["metric"] == "mean_adc"
        assert d["p"] == pytest.approx(0.015)
        assert d["adj_alpha"] == pytest.approx(0.025)

    def test_fdr_timepoints(self):
        """FDR timepoint summary lines ('Timepoint: XX -- N significant') are parsed."""
        text = (
            "Timepoint: BL — 3 significant\n"
            "Timepoint: W2 — 1 significant\n"
        )
        result = parse_stats_comparisons(text)
        assert len(result["fdr_timepoints"]) == 2
        assert result["fdr_timepoints"][0]["timepoint"] == "BL"
        assert result["fdr_timepoints"][0]["n_significant"] == 3

    def test_glme_excluded(self):
        """Competing-risk exclusion line is parsed for count, total, and percentage."""
        text = "Excluded 5/42 (11.9%) competing-risk patients\n"
        result = parse_stats_comparisons(text)
        ex = result["glme_excluded"]
        assert ex is not None
        assert ex["n_excluded"] == 5
        assert ex["n_total"] == 42
        assert ex["pct"] == pytest.approx(11.9)

    def test_empty_text(self):
        """Empty input yields empty lists and None for optional fields."""
        result = parse_stats_comparisons("")
        assert result["glme_interactions"] == []
        assert result["glme_details"] == []
        assert result["fdr_timepoints"] == []
        assert result["glme_excluded"] is None


# ---------------------------------------------------------------------------
# parse_stats_predictive
# ---------------------------------------------------------------------------

class TestParseStatsPredictive:
    """Verify extraction of feature selections and ROC analyses."""

    def test_elastic_net_features(self):
        """Single elastic net line: timepoint, lambda, and comma-separated feature names."""
        text = "Elastic Net Selected Features for BL (Opt Lambda=0.0532): mean_adc, mean_d, f_ratio\n"
        result = parse_stats_predictive(text)
        assert len(result["feature_selections"]) == 1
        fs = result["feature_selections"][0]
        assert fs["timepoint"] == "BL"
        assert fs["lambda"] == pytest.approx(0.0532)
        assert fs["features"] == ["mean_adc", "mean_d", "f_ratio"]

    def test_multiple_timepoints(self):
        """Two elastic net lines at different timepoints are both captured."""
        text = (
            "Elastic Net Selected Features for BL (Opt Lambda=0.05): a, b\n"
            "Elastic Net Selected Features for W2 (Opt Lambda=0.10): c\n"
        )
        result = parse_stats_predictive(text)
        assert len(result["feature_selections"]) == 2
        assert result["feature_selections"][1]["timepoint"] == "W2"

    def test_roc_analysis_single_block(self):
        """A complete ROC block with AUC, Youden cutoff, sensitivity, and specificity."""
        text = (
            "PRIMARY ROC ANALYSIS results for BL\n"
            "  AUC = 0.781\n"
            "  Youden Optimal Score Cutoff = 0.450\n"
            "  Sensitivity = 82.5% | Specificity = 68.2%\n"
        )
        result = parse_stats_predictive(text)
        assert len(result["roc_analyses"]) == 1
        roc = result["roc_analyses"][0]
        assert roc["timepoint"] == "BL"
        assert roc["auc"] == pytest.approx(0.781)
        assert roc["youden_cutoff"] == pytest.approx(0.450)
        assert roc["sensitivity"] == pytest.approx(82.5)
        assert roc["specificity"] == pytest.approx(68.2)

    def test_multiple_roc_blocks(self):
        """Two ROC headers should yield two parsed entries."""
        text = (
            "PRIMARY ROC ANALYSIS results for BL\n"
            "  AUC = 0.78\n"
            "PRIMARY ROC ANALYSIS results for W2\n"
            "  AUC = 0.65\n"
        )
        result = parse_stats_predictive(text)
        assert len(result["roc_analyses"]) == 2
        assert result["roc_analyses"][0]["auc"] == pytest.approx(0.78)
        assert result["roc_analyses"][1]["auc"] == pytest.approx(0.65)

    def test_empty_text(self):
        """Empty input yields empty lists for both feature selections and ROC analyses."""
        result = parse_stats_predictive("")
        assert result["feature_selections"] == []
        assert result["roc_analyses"] == []


# ---------------------------------------------------------------------------
# parse_survival
# ---------------------------------------------------------------------------

class TestParseSurvival:
    """Verify extraction of hazard ratios, global LRT, and IPCW weights."""

    def test_hazard_ratios(self):
        """Whitespace-delimited hazard ratio rows: covariate, HR, CI_lo, CI_hi, p-value.

        Two rows are parsed, verifying that all five fields (covariate name,
        hazard ratio, lower/upper 95% CI bounds, p-value) are correctly extracted.
        """
        text = (
            "  mean_adc   1.250   0.980   1.590   0.0680\n"
            "  delta_d    0.750   0.550   1.020   0.0340\n"
        )
        result = parse_survival(text)
        assert len(result["hazard_ratios"]) == 2
        hr0 = result["hazard_ratios"][0]
        assert hr0["covariate"] == "mean_adc"
        assert hr0["hr"] == pytest.approx(1.250)
        assert hr0["ci_lo"] == pytest.approx(0.980)
        assert hr0["ci_hi"] == pytest.approx(1.590)
        assert hr0["p"] == pytest.approx(0.0680)

    def test_global_lrt(self):
        """Global likelihood ratio test line: degrees of freedom, chi-squared statistic, p-value."""
        text = "Global LRT: chi2(2) = 7.82, p = 0.0200\n"
        result = parse_survival(text)
        lrt = result["global_lrt"]
        assert lrt is not None
        assert lrt["df"] == 2
        assert lrt["chi2"] == pytest.approx(7.82)
        assert lrt["p"] == pytest.approx(0.02)

    def test_ipcw_weights(self):
        """IPCW (inverse probability of censoring weights) range is extracted."""
        text = "IPCW weights applied ... range [0.85, 1.42]\n"
        result = parse_survival(text)
        ipcw = result["ipcw"]
        assert ipcw is not None
        assert ipcw["min_weight"] == pytest.approx(0.85)
        assert ipcw["max_weight"] == pytest.approx(1.42)

    def test_empty_text(self):
        """Empty input yields empty hazard ratios and None for LRT/IPCW."""
        result = parse_survival("")
        assert result["hazard_ratios"] == []
        assert result["global_lrt"] is None
        assert result["ipcw"] is None


# ---------------------------------------------------------------------------
# parse_baseline
# ---------------------------------------------------------------------------

class TestParseBaseline:
    """Verify extraction of outlier flags, totals, and baseline exclusions."""

    def test_outlier_flags(self):
        """Per-metric outlier flag lines with LF/LC/CR group breakdowns are parsed.

        Verifies extraction of metric name, total flagged count, and per-group
        counts (LF = local failure, LC = local control, CR = complete response).
        """
        text = (
            "Outlier flag (mean_adc): 3 flagged (LF=2, LC=1, CR=0)\n"
            "Outlier flag (mean_d): 1 flagged (LF=0, LC=0, CR=1)\n"
        )
        result = parse_baseline(text)
        assert len(result["outlier_flags"]) == 2
        f0 = result["outlier_flags"][0]
        assert f0["metric"] == "mean_adc"
        assert f0["n_flagged"] == 3
        assert f0["n_lf"] == 2
        assert f0["n_lc"] == 1
        assert f0["n_cr"] == 0

    def test_total_outliers(self):
        """Total outlier summary line: count removed, total, and percentage."""
        text = "Total outliers removed: 4 / 42 (9.52%)\n"
        result = parse_baseline(text)
        t = result["total_outliers"]
        assert t is not None
        assert t["n_removed"] == 4
        assert t["n_total"] == 42
        assert t["pct"] == pytest.approx(9.52)

    def test_baseline_exclusion_with_lf_rate(self):
        """Baseline exclusion with an LF rate comparison between included and excluded patients."""
        text = (
            "Excluded 6/48 patients due to missing baseline\n"
            "LF rate: included=35.7%, excluded=50.0%\n"
        )
        result = parse_baseline(text)
        be = result["baseline_exclusion"]
        assert be is not None
        assert be["n_excluded"] == 6
        assert be["n_total"] == 48
        assert be["lf_rate_included"] == pytest.approx(35.7)
        assert be["lf_rate_excluded"] == pytest.approx(50.0)

    def test_baseline_exclusion_without_lf_rate(self):
        """Baseline exclusion without a following LF rate line omits rate fields."""
        text = "Excluded 6/48 patients due to missing baseline\n"
        result = parse_baseline(text)
        be = result["baseline_exclusion"]
        assert be is not None
        assert "lf_rate_included" not in be

    def test_empty_text(self):
        """Empty input yields empty outlier flags and None for totals/exclusion."""
        result = parse_baseline("")
        assert result["outlier_flags"] == []
        assert result["total_outliers"] is None
        assert result["baseline_exclusion"] is None


# ---------------------------------------------------------------------------
# _read_log
# ---------------------------------------------------------------------------

class TestReadLog:
    """Verify graceful handling of missing and present log files."""

    def test_reads_existing_file(self, tmp_path: Path):
        """An existing log file's contents are returned as a string."""
        p = tmp_path / "test.log"
        p.write_text("hello\n", encoding="utf-8")
        assert _read_log(p) == "hello\n"

    def test_returns_empty_for_missing(self, tmp_path: Path):
        """A missing log file returns an empty string instead of raising."""
        assert _read_log(tmp_path / "nonexistent.log") == ""


# ---------------------------------------------------------------------------
# parse_all_logs (integration)
# ---------------------------------------------------------------------------

class TestParseAllLogs:
    """Integration test: parse_all_logs reads from DWI subfolders."""

    def test_parses_standard_logs(self, saved_files_with_logs: Path):
        """All four Standard log files are read and parsed into the expected structure.

        Verifies that the orchestrator correctly dispatches each log to its
        dedicated parser and assembles the results under the 'Standard' key.
        """
        results = parse_all_logs(saved_files_with_logs)
        assert "Standard" in results

        std = results["Standard"]
        # Stats comparisons: 2 GLME interaction p-values and 5 excluded patients
        assert len(std["stats_comparisons"]["glme_interactions"]) == 2
        assert std["stats_comparisons"]["glme_excluded"]["n_excluded"] == 5

        # Predictive: 2 elastic net timepoints (BL, W2) and 2 ROC blocks
        assert len(std["stats_predictive"]["feature_selections"]) == 2
        assert len(std["stats_predictive"]["roc_analyses"]) == 2

        # Survival: 2 hazard ratio rows and a global LRT with df=2
        assert len(std["survival"]["hazard_ratios"]) == 2
        assert std["survival"]["global_lrt"]["df"] == 2

        # Baseline: 2 outlier flag entries and 4 total outliers removed
        assert len(std["baseline"]["outlier_flags"]) == 2
        assert std["baseline"]["total_outliers"]["n_removed"] == 4

    def test_skips_missing_dwi_dirs(self, saved_files_dir: Path):
        """DWI directories without log files are skipped gracefully."""
        # saved_files_dir has empty Standard/dnCNN/IVIMnet subdirs
        results = parse_all_logs(saved_files_dir)
        # All dirs exist but have no logs → they still appear with empty results
        for dwi_type in ("Standard", "dnCNN", "IVIMnet"):
            if dwi_type in results:
                assert results[dwi_type]["baseline"]["outlier_flags"] == []

    def test_empty_folder(self, tmp_path: Path):
        """A folder with no DWI subfolders returns an empty dict."""
        empty = tmp_path / "saved_files_empty"
        empty.mkdir()
        results = parse_all_logs(empty)
        assert results == {}
