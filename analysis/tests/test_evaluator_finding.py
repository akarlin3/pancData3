"""Tests for the Finding Pydantic model in evaluator.py."""

import sys
import os

import pytest
from pydantic import ValidationError

# Ensure the repo root is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from improvement_loop.evaluator import Finding, should_continue_loop  # type: ignore


# ── Helpers ──────────────────────────────────────────────────────────────────

def _make_finding(**overrides) -> Finding:
    """Build a valid Finding, optionally overriding fields."""
    defaults = {
        "dimension": "correctness",
        "file": "pipeline/core/fit_models.m",
        "function_name": "fit_ivim",
        "description": "IVIM fit uses default TolFun which causes non-convergence",
        "fix": "Set TolFun=1e-4 and add convergence check",
        "importance": 5,
        "branch_name": "improvement/fix-ivim-tolfun",
        "status": "pending",
    }
    defaults.update(overrides)
    return Finding(**defaults)


# ── Valid construction ───────────────────────────────────────────────────────

class TestFindingConstruction:
    def test_valid_finding(self):
        f = _make_finding()
        assert f.dimension == "correctness"
        assert f.file == "pipeline/core/fit_models.m"
        assert f.function_name == "fit_ivim"
        assert f.importance == 5
        assert f.branch_name == "improvement/fix-ivim-tolfun"
        assert f.status == "pending"

    def test_optional_fields_default_none(self):
        f = _make_finding(function_name=None, status=None)
        assert f.function_name is None
        assert f.status is None

    def test_all_valid_dimensions(self):
        dims = [
            "performance", "correctness", "error_handling", "modularity",
            "memory", "code_quality", "test_coverage", "security",
            "cross_platform",
        ]
        for dim in dims:
            f = _make_finding(dimension=dim)
            assert f.dimension == dim

    def test_all_valid_statuses(self):
        for status in ("pending", "implemented", "merged", None):
            f = _make_finding(status=status)
            assert f.status == status

    def test_importance_boundaries(self):
        assert _make_finding(importance=1).importance == 1
        assert _make_finding(importance=10).importance == 10


# ── Invalid dimension ────────────────────────────────────────────────────────

class TestInvalidDimension:
    def test_invalid_dimension_raises(self):
        with pytest.raises(ValidationError):
            _make_finding(dimension="style")

    def test_empty_dimension_raises(self):
        with pytest.raises(ValidationError):
            _make_finding(dimension="")


# ── Invalid importance ───────────────────────────────────────────────────────

class TestInvalidImportance:
    def test_importance_zero_raises(self):
        with pytest.raises(ValidationError, match="importance"):
            _make_finding(importance=0)

    def test_importance_eleven_raises(self):
        with pytest.raises(ValidationError, match="importance"):
            _make_finding(importance=11)

    def test_negative_importance_raises(self):
        with pytest.raises(ValidationError):
            _make_finding(importance=-1)


# ── Invalid branch_name ──────────────────────────────────────────────────────

class TestInvalidBranchName:
    def test_spaces_raises(self):
        with pytest.raises(ValidationError, match="spaces"):
            _make_finding(branch_name="improvement/fix ivim tolfun")

    def test_double_slash_raises(self):
        with pytest.raises(ValidationError, match="slashes"):
            _make_finding(branch_name="improvement/fix//tolfun")

    def test_extra_slash_raises(self):
        with pytest.raises(ValidationError, match="slashes"):
            _make_finding(branch_name="improvement/fix/tolfun")

    def test_missing_prefix_raises(self):
        with pytest.raises(ValidationError, match="improvement/"):
            _make_finding(branch_name="fix-ivim-tolfun")

    def test_empty_slug_raises(self):
        with pytest.raises(ValidationError, match="empty"):
            _make_finding(branch_name="improvement/")

    def test_slug_too_long_raises(self):
        with pytest.raises(ValidationError, match="50"):
            _make_finding(branch_name="improvement/" + "a" * 51)

    def test_tilde_raises(self):
        with pytest.raises(ValidationError, match="invalid"):
            _make_finding(branch_name="improvement/fix~tolfun")

    def test_caret_raises(self):
        with pytest.raises(ValidationError, match="invalid"):
            _make_finding(branch_name="improvement/fix^tolfun")


# ── to_log_dict ──────────────────────────────────────────────────────────────

class TestToLogDict:
    def test_output_matches_legacy_format(self):
        f = _make_finding()
        d = f.to_log_dict()
        assert d == {
            "dimension": "correctness",
            "file": "pipeline/core/fit_models.m",
            "function_name": "fit_ivim",
            "description": "IVIM fit uses default TolFun which causes non-convergence",
            "fix": "Set TolFun=1e-4 and add convergence check",
            "importance": 5,
            "branch_name": "improvement/fix-ivim-tolfun",
            "status": "pending",
        }

    def test_omits_none_function_name(self):
        f = _make_finding(function_name=None, status=None)
        d = f.to_log_dict()
        assert "function_name" not in d
        assert "status" not in d

    def test_includes_status_when_set(self):
        f = _make_finding(status="merged")
        d = f.to_log_dict()
        assert d["status"] == "merged"

    def test_dict_is_plain_dict(self):
        """Ensure to_log_dict returns a plain dict, not a Pydantic object."""
        d = _make_finding().to_log_dict()
        assert type(d) is dict


# ── should_continue_loop ─────────────────────────────────────────────────────

class TestShouldContinueLoop:
    GOOD_SCORES = {
        "specificity": 8, "accuracy": 8, "coverage": 8,
        "prioritization": 8, "domain_appropriateness": 8,
        "overall": 8, "flags": [], "reasoning": "Good."
    }
    LOW_COVERAGE_SCORES = {
        "specificity": 8, "accuracy": 8, "coverage": 4,
        "prioritization": 8, "domain_appropriateness": 8,
        "overall": 7, "flags": [], "reasoning": "Low coverage."
    }

    def test_returns_true_when_high_importance_finding(self):
        findings = [_make_finding(importance=2)]
        assert should_continue_loop(self.GOOD_SCORES, findings) is True

    def test_returns_false_when_all_below_threshold(self):
        findings = [_make_finding(importance=1)]
        assert should_continue_loop(self.GOOD_SCORES, findings) is False

    def test_returns_false_with_empty_findings(self):
        assert should_continue_loop(self.GOOD_SCORES, []) is False

    def test_returns_true_when_coverage_low(self):
        findings = [_make_finding(importance=1)]
        assert should_continue_loop(self.LOW_COVERAGE_SCORES, findings) is True

    def test_importance_exactly_two_continues(self):
        """importance=2 is the threshold — should continue."""
        findings = [_make_finding(importance=2)]
        assert should_continue_loop(self.GOOD_SCORES, findings) is True

    def test_importance_exactly_one_stops(self):
        """importance=1 is below threshold — should stop (with good scores)."""
        findings = [_make_finding(importance=1)]
        assert should_continue_loop(self.GOOD_SCORES, findings) is False

    def test_returns_true_when_flags_present(self):
        flagged_scores = {**self.GOOD_SCORES, "flags": ["LEAKAGE_RISK"]}
        assert should_continue_loop(flagged_scores, []) is True

    def test_evaluation_failed_flag_ignored(self):
        """EVALUATION_FAILED flag alone should not force continuation."""
        failed_scores = {**self.GOOD_SCORES, "flags": ["EVALUATION_FAILED"]}
        assert should_continue_loop(failed_scores, []) is False
