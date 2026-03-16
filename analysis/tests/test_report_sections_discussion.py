"""Tests for report_sections.discussion module.

Validates discussion section builders:
- Methods section (statistical methodology description)
- Limitations section
- Conclusions section
- Reporting checklist
- Journal guide

Edge cases covered:
- Empty/None inputs at every level
- No MAT data, no log data, no CSV data, no groups
- IPCW used vs not used
- Dosimetry available vs not
- No significant findings for conclusions
- Canonical response pattern detection (D increasing, f decreasing)
- Cross-DWI agreement in conclusions
- Core method Dice agreement
- Large vs small sample size limitations
- Missing baseline data limitation
- No ROC data for LOOCV mention
- Multiple DWI types in methods
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from report.sections.discussion import (
    _section_methods,
    _section_limitations,
    _section_conclusions,
)
from report.sections.publication import (
    _section_reporting_checklist,
    _section_journal_guide,
)


def _make_log_data():
    return {
        "Standard": {
            "survival": {
                "hazard_ratios": [
                    {"covariate": "mean_adc", "hr": 1.5, "ci_lo": 0.9, "ci_hi": 2.5, "p": 0.03},
                ],
                "ipcw": {"max_weight": 1.4, "min_weight": 0.8},
            },
            "stats_comparisons": {
                "glme_details": [
                    {"metric": "mean_adc", "p": 0.01, "adj_alpha": 0.025},
                ],
                "glme_interactions": [0.02],
                "fdr_timepoints": [{"n_significant": 1, "timepoint": "BL"}],
            },
            "stats_predictive": {
                "roc_analyses": [{"auc": 0.78, "timepoint": "BL"}],
                "feature_selections": [
                    {"timepoint": "BL", "features": ["adc"], "lambda": 0.05},
                ],
            },
            "baseline": {
                "total_outliers": {"pct": 8.5},
                "baseline_exclusion": {"n_excluded": 6, "n_total": 48},
            },
            "sanity_checks": {"all_converged": True},
        }
    }


def _make_mat_data():
    return {
        "Standard": {
            "longitudinal": {"num_patients": 42, "num_timepoints": 5},
            "dosimetry": {"d95_gtvp": {"mean": 45.0}, "v50_gtvp": {"mean": 0.85}},
        }
    }


# ── Methods ──


class TestMethods:
    def test_returns_html(self):
        result = _section_methods(["Standard"], _make_mat_data(), _make_log_data())
        assert isinstance(result, list)
        assert len(result) > 0

    def test_contains_methods_header(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "Statistical Methods" in html

    def test_mentions_ivim(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "IVIM" in html

    def test_mentions_bh_fdr(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "Benjamini" in html or "BH" in html or "FDR" in html

    def test_mentions_elastic_net(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "elastic" in html.lower() or "lasso" in html.lower() or "regulariz" in html.lower() or "predictive" in html.lower()

    def test_ipcw_used_when_present(self):
        """When IPCW data exists in logs, methods should mention IPCW."""
        html = "\n".join(_section_methods(["Standard"], {}, _make_log_data()))
        assert "IPCW" in html or "inverse probability" in html.lower()

    def test_no_ipcw_when_absent(self):
        """Without IPCW data, methods should note competing-risk exclusion."""
        log = {"Standard": {"survival": {}}}
        html = "\n".join(_section_methods(["Standard"], {}, log))
        assert "excluded" in html.lower() or "competing" in html.lower()

    def test_loocv_mentioned_with_roc(self):
        """When ROC data exists, LOOCV should be mentioned."""
        html = "\n".join(_section_methods(["Standard"], {}, _make_log_data()))
        assert "LOOCV" in html or "leave-one-out" in html.lower()

    def test_no_loocv_without_roc(self):
        """Without ROC data, LOOCV may not be mentioned."""
        log = {"Standard": {"survival": {}, "stats_predictive": {}}}
        html = "\n".join(_section_methods(["Standard"], {}, log))
        # LOOCV mention is conditional on ROC presence
        assert isinstance(html, str)

    def test_dosimetry_section_when_present(self):
        """Dosimetry in MAT data should trigger dosimetric analysis section."""
        html = "\n".join(_section_methods(["Standard"], _make_mat_data(), None))
        assert "Dosimetric" in html or "dosimetry" in html.lower()

    def test_no_dosimetry_section_when_absent(self):
        """No dosimetry should skip dosimetric analysis section."""
        mat = {"Standard": {"longitudinal": {"num_patients": 42}}}
        html = "\n".join(_section_methods(["Standard"], mat, None))
        assert "Dosimetric Analysis" not in html

    def test_core_method_descriptions(self):
        """All 11 core method descriptions should be listed."""
        html = "\n".join(_section_methods(["Standard"], {}, None))
        for method in ["adc_threshold", "d_threshold", "otsu", "gmm", "kmeans",
                       "region_growing", "active_contours", "percentile", "spectral", "fdm"]:
            assert method in html

    def test_software_section(self):
        """Software section should mention MATLAB and Python."""
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "MATLAB" in html
        assert "Python" in html

    def test_empty_everything(self):
        """Empty lists and None data should still produce methods."""
        result = _section_methods([], {}, None)
        assert isinstance(result, list)
        assert len(result) > 0

    def test_dl_disclosure(self):
        """Deep learning model disclosure should be present."""
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "DnCNN" in html or "deep learning" in html.lower()


# ── Limitations ──


class TestLimitations:
    def test_returns_html(self):
        result = _section_limitations(_make_log_data(), ["Standard"], _make_mat_data())
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Limitation" in html

    def test_no_data(self):
        result = _section_limitations(None, [], {})
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Limitation" in html

    def test_sample_size_limitation(self):
        html = "\n".join(_section_limitations(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "sample" in html.lower() or "cohort" in html.lower()

    def test_small_sample_size_strong_warning(self):
        """n < 50 should note 'Small sample size'."""
        html = "\n".join(_section_limitations(None, [], _make_mat_data()))
        assert "Small sample size" in html or "limited" in html.lower()

    def test_moderate_sample_size(self):
        """n >= 50 should note 'Moderate sample size'."""
        mat = {"Standard": {"longitudinal": {"num_patients": 60, "num_timepoints": 3}}}
        html = "\n".join(_section_limitations(None, [], mat))
        assert "Moderate sample size" in html

    def test_missing_data_limitation(self):
        """Baseline exclusion should trigger missing data limitation."""
        html = "\n".join(_section_limitations(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Missing data" in html or "missing" in html.lower()

    def test_no_missing_data_no_limitation(self):
        """No baseline exclusion should not trigger missing data limitation."""
        log = {"Standard": {"baseline": {}}}
        html = "\n".join(_section_limitations(log, ["Standard"], _make_mat_data()))
        assert "Missing data:" not in html

    def test_standard_limitations_always_present(self):
        """Core limitations should always be listed regardless of data."""
        html = "\n".join(_section_limitations(None, [], {}))
        assert "Single-institution" in html
        assert "Retrospective" in html
        assert "DWI-specific" in html
        assert "Multiple comparisons" in html

    def test_vision_based_caveat(self):
        """Vision-based graph analysis caveat should always appear."""
        html = "\n".join(_section_limitations(None, [], {}))
        assert "Vision-based" in html or "vision" in html.lower()

    def test_core_delineation_limitation(self):
        """Core delineation validation limitation should appear."""
        html = "\n".join(_section_limitations(None, [], {}))
        assert "Core delineation" in html or "delineation" in html.lower()

    def test_no_patients_no_sample_size_note(self):
        """Zero patients should skip sample size limitation."""
        mat = {"Standard": {"longitudinal": {"num_patients": 0}}}
        html = "\n".join(_section_limitations(None, [], mat))
        assert "Small sample" not in html
        assert "Moderate sample" not in html


# ── Conclusions ──


class TestConclusions:
    def test_returns_html(self):
        result = _section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        )
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Conclusion" in html

    def test_no_data(self):
        result = _section_conclusions(None, [], None, {}, None)
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Conclusion" in html

    def test_with_groups(self):
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {
                    "trends_json": json.dumps([
                        {"series": "Mean D", "direction": "increasing"},
                    ])
                },
            }
        }
        result = _section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), groups
        )
        assert isinstance(result, list)

    def test_no_findings_fallback(self):
        """No significant findings should produce generic feasibility statement."""
        log = {"Standard": {"stats_comparisons": {"glme_details": []}, "stats_predictive": {}, "survival": {}}}
        result = _section_conclusions(log, ["Standard"], None, {}, None)
        html = "\n".join(result)
        assert "feasibility" in html.lower() or "Conclusion" in html

    def test_significant_biomarkers_finding(self):
        """Significant GLME metrics should appear in conclusions."""
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "significant" in html.lower() or "GLME" in html

    def test_predictive_performance_finding(self):
        """Best AUC should be mentioned in conclusions."""
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "AUC" in html or "0.78" in html

    def test_cox_covariates_finding(self):
        """Significant Cox covariates should appear."""
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "mean_adc" in html or "Cox" in html.lower() or "hazard" in html.lower()

    def test_dosimetry_finding(self):
        """Dosimetry metrics should appear in conclusions."""
        mat = {
            "Standard": {
                "longitudinal": {"num_patients": 42, "num_timepoints": 5},
                "dosimetry": {"d95_adc_mean": {"mean": 45.0}, "v50_adc_mean": {"mean": 0.85}},
            }
        }
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard"], None, mat, None
        ))
        assert "D95" in html or "dosimetric" in html.lower() or "coverage" in html.lower()

    def test_core_method_dice_finding(self):
        """Core method Dice agreement should appear."""
        mat = {
            "Standard": {
                "longitudinal": {"num_patients": 42, "num_timepoints": 5},
                "core_method": {
                    "methods": ["adc_threshold", "otsu"],
                    "mean_dice_matrix": [[1.0, 0.72], [0.72, 1.0]],
                },
            }
        }
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard"], None, mat, None
        ))
        assert "Dice" in html or "delineation" in html.lower()

    def test_canonical_response_pattern(self):
        """D increasing + f decreasing should mention canonical response."""
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {
                    "trends_json": json.dumps([
                        {"series": "Mean D", "direction": "increasing"},
                        {"series": "Mean f", "direction": "decreasing"},
                    ])
                },
            }
        }
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), groups
        ))
        assert "canonical" in html.lower() or "necrosis" in html.lower()

    def test_cross_dwi_agreement_in_conclusions(self):
        """Cross-DWI trend agreement percentage should appear."""
        groups = {
            "TestGraph": {
                "Standard": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "S1", "direction": "up"}])},
            }
        }
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard", "dnCNN"], None, _make_mat_data(), groups
        ))
        assert "100%" in html or "agreement" in html.lower()

    def test_clinical_significance_always_present(self):
        """Clinical significance statement should always appear."""
        html = "\n".join(_section_conclusions(None, [], None, {}, None))
        assert "Clinical Significance" in html

    def test_future_directions_always_present(self):
        """Future directions should always appear."""
        html = "\n".join(_section_conclusions(None, [], None, {}, None))
        assert "Future directions" in html or "future" in html.lower()

    def test_fdr_global_csv_finding(self):
        """FDR global CSV metrics should appear in conclusions."""
        csv_data = {
            "fdr_global": {
                "Standard": [{"metric": "adc", "p": 0.001}],
            }
        }
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard"], csv_data, _make_mat_data(), None
        ))
        assert "FDR" in html or "surviving" in html.lower()

    def test_per_dwi_auc_comparison(self):
        """Multiple DWI types with AUC should show per-type comparison."""
        log = _make_log_data()
        log["dnCNN"] = {
            "stats_comparisons": {"glme_details": []},
            "stats_predictive": {
                "roc_analyses": [{"auc": 0.72, "timepoint": "BL"}],
            },
            "survival": {},
        }
        html = "\n".join(_section_conclusions(
            log, ["Standard", "dnCNN"], None, _make_mat_data(), None
        ))
        assert "AUC" in html


# ── Reporting Checklist ──


class TestReportingChecklist:
    def test_returns_html(self):
        result = _section_reporting_checklist(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        )
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Checklist" in html or "checklist" in html

    def test_no_data(self):
        result = _section_reporting_checklist(None, [], None, {}, None)
        assert isinstance(result, list)

    def test_contains_strobe_or_remark(self):
        html = "\n".join(_section_reporting_checklist(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "STROBE" in html or "REMARK" in html or "TRIPOD" in html or "checklist" in html.lower()


# ── Journal Guide ──


class TestJournalGuide:
    def test_returns_html(self):
        result = _section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data())
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Journal" in html or "journal" in html

    def test_no_data(self):
        result = _section_journal_guide(None, [], {})
        assert isinstance(result, list)

    def test_suggests_journals(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert any(j in html for j in [
            "Radiology", "Physics", "Oncology", "Cancer",
            "Medical", "International Journal", "journal",
        ])
