"""Tests for report_sections.discussion module.

Validates discussion section builders:
- Methods section (statistical methodology description)
- Limitations section
- Conclusions section
- Reporting checklist
- Journal guide
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
            "dosimetry": {
                "d95_adc_mean": 45.0,
                "v50_adc_mean": 0.85,
                "d95_d_mean": 44.0,
                "v50_d_mean": 0.80,
            },
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

    def test_mentions_wilcoxon(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "Wilcoxon" in html

    def test_mentions_glme(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "GLME" in html or "mixed-effects" in html.lower()

    def test_mentions_cox(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "Cox" in html

    def test_mentions_knn_imputation(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "KNN" in html or "k-nearest" in html.lower()

    def test_mentions_loocv_when_roc_present(self):
        html = "\n".join(_section_methods(["Standard"], {}, _make_log_data()))
        assert "LOOCV" in html or "leave-one-out" in html.lower()

    def test_no_loocv_without_roc(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["roc_analyses"] = []
        html = "\n".join(_section_methods(["Standard"], {}, log))
        assert "LOOCV" not in html

    def test_mentions_ipcw_when_present(self):
        html = "\n".join(_section_methods(["Standard"], {}, _make_log_data()))
        assert "IPCW" in html

    def test_no_ipcw_without_survival_data(self):
        log = _make_log_data()
        del log["Standard"]["survival"]["ipcw"]
        html = "\n".join(_section_methods(["Standard"], {}, log))
        # Should mention competing-risk patients excluded instead
        assert "excluded" in html.lower() or "IPCW" not in html

    def test_dosimetry_section_when_present(self):
        html = "\n".join(_section_methods(["Standard"], _make_mat_data(), None))
        assert "Dosimetric" in html or "dosimetry" in html.lower()

    def test_no_dosimetry_section_when_absent(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "Dosimetric Analysis" not in html

    def test_core_method_descriptions(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "adc_threshold" in html
        assert "otsu" in html
        assert "spectral" in html
        assert "fdm" in html

    def test_deep_learning_disclosure(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "Deep Learning" in html or "DnCNN" in html

    def test_data_quality_assurance(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "sanity check" in html.lower() or "Quality Assurance" in html

    def test_software_section(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "MATLAB" in html
        assert "Python" in html

    def test_firth_mentioned(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "Firth" in html

    def test_collinearity_handling(self):
        html = "\n".join(_section_methods(["Standard"], {}, None))
        assert "collinear" in html.lower() or "correlated" in html.lower()

    def test_empty_inputs(self):
        result = _section_methods([], {}, None)
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Statistical Methods" in html


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

    def test_sample_size_limitation_small(self):
        html = "\n".join(_section_limitations(_make_log_data(), ["Standard"], _make_mat_data()))
        # n=42 < 50, should say "Small sample size"
        assert "Small sample size" in html or "sample" in html.lower()

    def test_sample_size_limitation_moderate(self):
        mat = {"Standard": {"longitudinal": {"num_patients": 60, "num_timepoints": 5}}}
        html = "\n".join(_section_limitations(_make_log_data(), ["Standard"], mat))
        assert "Moderate sample size" in html

    def test_single_institution_limitation(self):
        html = "\n".join(_section_limitations(None, [], {}))
        assert "Single-institution" in html or "single institution" in html.lower()

    def test_retrospective_limitation(self):
        html = "\n".join(_section_limitations(None, [], {}))
        assert "Retrospective" in html or "retrospective" in html.lower()

    def test_missing_data_limitation(self):
        html = "\n".join(_section_limitations(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Missing data" in html or "missing" in html.lower()

    def test_missing_data_shows_counts(self):
        html = "\n".join(_section_limitations(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "6/48" in html

    def test_no_missing_data_limitation_when_no_exclusions(self):
        log = _make_log_data()
        log["Standard"]["baseline"]["baseline_exclusion"]["n_excluded"] = 0
        html = "\n".join(_section_limitations(log, ["Standard"], _make_mat_data()))
        # Should not mention missing data exclusion specifically
        assert "6/48" not in html

    def test_dwi_specific_limitations(self):
        html = "\n".join(_section_limitations(None, [], {}))
        assert "DWI" in html or "IVIM" in html

    def test_tumour_delineation_limitation(self):
        html = "\n".join(_section_limitations(None, [], {}))
        assert "delineation" in html.lower() or "contour" in html.lower()

    def test_multiple_comparisons_limitation(self):
        html = "\n".join(_section_limitations(None, [], {}))
        assert "multiple comparisons" in html.lower() or "FDR" in html

    def test_proportional_hazards_limitation(self):
        html = "\n".join(_section_limitations(None, [], {}))
        assert "Proportional hazards" in html

    def test_knn_imputation_limitation(self):
        html = "\n".join(_section_limitations(None, [], {}))
        assert "KNN" in html or "imputation" in html.lower()

    def test_competing_risk_limitation(self):
        html = "\n".join(_section_limitations(None, [], {}))
        assert "Competing-risk" in html or "competing" in html.lower()

    def test_core_delineation_validation_limitation(self):
        html = "\n".join(_section_limitations(None, [], {}))
        assert "Core delineation" in html or "ground-truth" in html.lower()

    def test_vision_analysis_limitation(self):
        html = "\n".join(_section_limitations(None, [], {}))
        assert "Vision" in html or "vision" in html.lower()

    def test_no_sample_size_when_no_mat_data(self):
        html = "\n".join(_section_limitations(None, [], {}))
        assert "n = " not in html

    def test_limitations_are_list_items(self):
        html = "\n".join(_section_limitations(None, [], {}))
        assert "<li>" in html
        assert "limitation-list" in html


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

    def test_no_data_fallback_message(self):
        result = _section_conclusions(None, [], None, {}, None)
        html = "\n".join(result)
        assert "feasibility" in html.lower() or "Detailed findings" in html

    def test_with_groups(self):
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {
                    "trends_json": json.dumps([
                        {"series": "Mean D", "direction": "increasing"},
                    ])
                },
            },
        }
        result = _section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), groups
        )
        assert isinstance(result, list)

    def test_glme_significant_findings(self):
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "GLME" in html or "interaction" in html.lower()

    def test_fdr_global_findings(self):
        csv_data = {
            "fdr_global": {
                "Standard": [{"metric": "adc"}],
            }
        }
        html = "\n".join(_section_conclusions(
            None, [], csv_data, {}, None
        ))
        assert "FDR" in html or "fdr" in html.lower()

    def test_predictive_performance_acceptable(self):
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        # AUC 0.78 => "acceptable"
        assert "acceptable" in html.lower()

    def test_predictive_performance_excellent(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["roc_analyses"] = [
            {"auc": 0.85, "timepoint": "BL"},
        ]
        html = "\n".join(_section_conclusions(log, ["Standard"], None, _make_mat_data(), None))
        assert "excellent" in html.lower()

    def test_predictive_performance_limited(self):
        log = _make_log_data()
        log["Standard"]["stats_predictive"]["roc_analyses"] = [
            {"auc": 0.55, "timepoint": "BL"},
        ]
        html = "\n".join(_section_conclusions(log, ["Standard"], None, _make_mat_data(), None))
        assert "limited" in html.lower()

    def test_cox_significant_covariates(self):
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "Cox" in html or "mean_adc" in html

    def test_dosimetry_adequate(self):
        html = "\n".join(_section_conclusions(
            None, [], None, _make_mat_data(), None
        ))
        assert "adequate" in html.lower()

    def test_dosimetry_suboptimal(self):
        mat = _make_mat_data()
        mat["Standard"]["dosimetry"]["d95_adc_mean"] = 40.0
        html = "\n".join(_section_conclusions(None, [], None, mat, None))
        assert "sub-optimal" in html.lower() or "suboptimal" in html.lower()

    def test_core_method_agreement(self):
        mat = _make_mat_data()
        mat["Standard"]["core_method"] = {
            "methods": ["adc_threshold", "otsu", "gmm"],
            "mean_dice_matrix": [
                [1.0, 0.65, 0.55],
                [0.65, 1.0, 0.60],
                [0.55, 0.60, 1.0],
            ],
        }
        html = "\n".join(_section_conclusions(None, [], None, mat, None))
        assert "Dice" in html or "agreement" in html.lower()

    def test_cross_dwi_trend_agreement(self):
        groups = {
            "SomeGraph": {
                "Standard": {"trends_json": json.dumps([{"series": "x", "direction": "increasing"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "x", "direction": "increasing"}])},
            },
        }
        html = "\n".join(_section_conclusions(None, ["Standard", "dnCNN"], None, {}, groups))
        assert "agreement" in html.lower() or "100%" in html

    def test_cross_dwi_trend_disagreement(self):
        groups = {
            "SomeGraph": {
                "Standard": {"trends_json": json.dumps([{"series": "x", "direction": "increasing"}])},
                "dnCNN": {"trends_json": json.dumps([{"series": "x", "direction": "decreasing"}])},
            },
        }
        html = "\n".join(_section_conclusions(None, ["Standard", "dnCNN"], None, {}, groups))
        assert "0%" in html or "limited" in html.lower()

    def test_canonical_response_pattern(self):
        groups = {
            "Longitudinal_Mean_Metrics": {
                "Standard": {
                    "trends_json": json.dumps([
                        {"series": "Mean D", "direction": "increasing"},
                        {"series": "Mean f", "direction": "decreasing"},
                    ])
                },
            },
        }
        html = "\n".join(_section_conclusions(None, ["Standard"], None, {}, groups))
        assert "canonical" in html.lower() or "necrosis" in html.lower()

    def test_clinical_significance_statement(self):
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "Clinical Significance" in html

    def test_future_directions(self):
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "Future directions" in html or "future" in html.lower()

    def test_findings_as_ordered_list(self):
        html = "\n".join(_section_conclusions(
            _make_log_data(), ["Standard"], None, _make_mat_data(), None
        ))
        assert "<ol>" in html

    def test_per_dwi_auc_comparison(self):
        log = _make_log_data()
        log["dnCNN"] = {
            "stats_comparisons": {"glme_details": []},
            "stats_predictive": {
                "roc_analyses": [{"auc": 0.72, "timepoint": "BL"}],
            },
            "survival": {},
            "baseline": {},
        }
        html = "\n".join(_section_conclusions(
            log, ["Standard", "dnCNN"], None, _make_mat_data(), None
        ))
        assert "Per-type comparison" in html or "Standard" in html

    def test_invalid_trends_json_handled(self):
        groups = {
            "SomeGraph": {
                "Standard": {"trends_json": "not valid json"},
                "dnCNN": {"trends_json": json.dumps([{"series": "x", "direction": "up"}])},
            },
        }
        # Should not crash
        result = _section_conclusions(None, ["Standard", "dnCNN"], None, {}, groups)
        assert isinstance(result, list)


# ── Reporting Checklist ──


class TestReportingChecklist:
    def test_returns_html(self):
        result = _section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        )
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Checklist" in html or "checklist" in html

    def test_no_data(self):
        result = _section_reporting_checklist(None, [], {}, None, None)
        assert isinstance(result, list)

    def test_contains_strobe_or_remark(self):
        html = "\n".join(_section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        ))
        assert "STROBE" in html or "REMARK" in html or "TRIPOD" in html or "checklist" in html.lower()

    def test_stat_cards_for_counts(self):
        html = "\n".join(_section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        ))
        assert "Addressed" in html
        assert "Partial" in html

    def test_full_data_increases_done_count(self):
        result_full = _section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        )
        result_empty = _section_reporting_checklist(None, [], {}, None, None)
        html_full = "\n".join(result_full)
        html_empty = "\n".join(result_empty)
        # Full data should have more "Addressed" items
        assert html_full.count("checklist-done") >= html_empty.count("checklist-done")

    def test_na_items_shown_when_present(self):
        html = "\n".join(_section_reporting_checklist(None, [], {}, None, None))
        assert "N/A" in html

    def test_checklist_table_structure(self):
        html = "\n".join(_section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        ))
        assert "<thead>" in html
        assert "Item" in html
        assert "Requirement" in html
        assert "Status" in html

    def test_strobe_items_present(self):
        html = "\n".join(_section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        ))
        assert "STROBE 1" in html
        assert "STROBE 12" in html
        assert "STROBE 19" in html

    def test_remark_items_present(self):
        html = "\n".join(_section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, None
        ))
        assert "REMARK 1" in html
        assert "REMARK 5" in html

    def test_with_csv_data(self):
        csv_data = {"fdr_global": {"Standard": [{"metric": "adc"}]}}
        result = _section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), csv_data, None
        )
        assert isinstance(result, list)

    def test_with_graph_rows(self):
        rows = [{"graph_type": "scatter", "dwi_type": "Standard"}]
        result = _section_reporting_checklist(
            _make_log_data(), ["Standard"], _make_mat_data(), None, rows
        )
        assert isinstance(result, list)

    def test_methods_section_always_done(self):
        html = "\n".join(_section_reporting_checklist(None, [], {}, None, None))
        # STROBE 12 "Statistical methods" should always be "done"
        assert "Full methods section auto-generated" in html


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

    def test_survival_adds_acta_oncologica(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Acta Oncologica" in html

    def test_predictive_adds_european_radiology(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "European Radiology" in html

    def test_no_survival_no_acta(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": []},
            "stats_predictive": {"roc_analyses": []},
        }}
        html = "\n".join(_section_journal_guide(log, ["Standard"], _make_mat_data()))
        assert "Acta Oncologica" not in html

    def test_no_predictive_no_european_radiology(self):
        log = {"Standard": {
            "survival": {"hazard_ratios": [{"covariate": "x", "hr": 1.5, "p": 0.03}]},
            "stats_predictive": {"roc_analyses": []},
        }}
        html = "\n".join(_section_journal_guide(log, ["Standard"], _make_mat_data()))
        assert "European Radiology" not in html

    def test_manuscript_checklist(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Manuscript Preparation" in html
        assert "Title page" in html
        assert "Abstract" in html
        assert "IRB" in html

    def test_suggested_keywords(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "diffusion-weighted imaging" in html
        assert "IVIM" in html
        assert "pancreatic cancer" in html

    def test_keywords_include_elastic_net_when_predictive(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "elastic net" in html

    def test_keywords_include_survival_when_present(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "survival analysis" in html

    def test_keywords_include_deep_learning_multi_dwi(self):
        html = "\n".join(_section_journal_guide(
            _make_log_data(), ["Standard", "dnCNN"], _make_mat_data()
        ))
        assert "deep learning" in html.lower()

    def test_no_deep_learning_keyword_single_dwi(self):
        html = "\n".join(_section_journal_guide(
            _make_log_data(), ["Standard"], _make_mat_data()
        ))
        assert "deep learning denoising" not in html

    def test_copy_button_for_keywords(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "copy-btn" in html or "Copy" in html

    def test_journal_table_structure(self):
        html = "\n".join(_section_journal_guide(_make_log_data(), ["Standard"], _make_mat_data()))
        assert "Word Limit" in html
        assert "Article Type" in html
        assert "Scope Match" in html

    def test_base_journals_always_present(self):
        html = "\n".join(_section_journal_guide(None, [], {}))
        assert "Radiotherapy and Oncology" in html
        assert "Medical Physics" in html
        assert "Physics in Medicine" in html

    def test_empty_dwi_types(self):
        result = _section_journal_guide(None, [], {})
        assert isinstance(result, list)
        html = "\n".join(result)
        assert "Journal" in html
