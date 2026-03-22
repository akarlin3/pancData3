"""Tests for improvement_loop.agents.reviewer module."""

import json
import os
import sys

import pytest

# Ensure repo root is on sys.path for absolute imports
_repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from improvement_loop.agents.reviewer import (
    ReviewVerdict,
    _generate_diff,
    _parse_review_verdict,
    review,
)
from improvement_loop.evaluator import Finding
from improvement_loop.loop_config import LoopConfig, reset_config


@pytest.fixture(autouse=True)
def _disable_rag(monkeypatch):
    """Prevent review() from initialising ChromaDB via RAG retriever."""
    monkeypatch.setattr(
        "improvement_loop.agents.reviewer._get_loop_config",
        lambda: LoopConfig(rag_enabled=False),
    )
    yield
    reset_config()


def _make_finding(**overrides) -> Finding:
    """Helper to create a Finding with sensible defaults."""
    defaults = {
        "dimension": "correctness",
        "file": "pipeline/core/fit_models.m",
        "description": "Off-by-one in loop bound at line 42",
        "fix": "Change < to <= in the for-loop condition.",
        "importance": 5,
        "branch_name": "improvement/test-review",
    }
    defaults.update(overrides)
    return Finding(**defaults)


# ---------------------------------------------------------------------------
# review() — dry_run mode
# ---------------------------------------------------------------------------

class TestReviewDryRun:
    """review() in dry_run mode returns approve."""

    def test_returns_approve(self):
        finding = _make_finding()
        verdict = review(finding, "old content", "new content", dry_run=True)
        assert verdict.verdict == "approve"

    def test_dry_run_reasoning(self):
        finding = _make_finding()
        verdict = review(finding, "", "", dry_run=True)
        assert verdict.reasoning == "dry-run"

    def test_no_risk_flags(self):
        finding = _make_finding()
        verdict = review(finding, "", "", dry_run=True)
        assert verdict.risk_flags == []


# ---------------------------------------------------------------------------
# _parse_review_verdict() — valid JSON
# ---------------------------------------------------------------------------

class TestParseReviewVerdictValid:
    """_parse_review_verdict with well-formed JSON input."""

    def test_approve(self):
        raw = json.dumps({
            "verdict": "approve",
            "reasoning": "Change looks correct.",
            "risk_flags": [],
        })
        result = _parse_review_verdict(raw)
        assert result is not None
        assert result.verdict == "approve"
        assert result.reasoning == "Change looks correct."
        assert result.risk_flags == []

    def test_reject_with_flags(self):
        raw = json.dumps({
            "verdict": "reject",
            "reasoning": "Modifies dependencies.",
            "risk_flags": ["DEPS_MODIFIED"],
        })
        result = _parse_review_verdict(raw)
        assert result is not None
        assert result.verdict == "reject"
        assert "DEPS_MODIFIED" in result.risk_flags

    def test_request_changes(self):
        raw = json.dumps({
            "verdict": "request_changes",
            "reasoning": "Missing edge case handling.",
            "risk_flags": [],
        })
        result = _parse_review_verdict(raw)
        assert result is not None
        assert result.verdict == "request_changes"

    def test_markdown_fenced(self):
        inner = json.dumps({
            "verdict": "approve",
            "reasoning": "Looks good.",
            "risk_flags": [],
        })
        raw = f"```json\n{inner}\n```"
        result = _parse_review_verdict(raw)
        assert result is not None
        assert result.verdict == "approve"

    def test_preamble_before_json(self):
        inner = json.dumps({
            "verdict": "approve",
            "reasoning": "Fine.",
            "risk_flags": [],
        })
        raw = f"Here is my review:\n{inner}"
        result = _parse_review_verdict(raw)
        assert result is not None
        assert result.verdict == "approve"


# ---------------------------------------------------------------------------
# _parse_review_verdict() — invalid JSON
# ---------------------------------------------------------------------------

class TestParseReviewVerdictInvalid:
    """_parse_review_verdict with invalid or incomplete input."""

    def test_garbage_returns_none(self):
        assert _parse_review_verdict("not json at all") is None

    def test_empty_string_returns_none(self):
        assert _parse_review_verdict("") is None

    def test_missing_verdict_key(self):
        raw = json.dumps({
            "reasoning": "Something.",
            "risk_flags": [],
        })
        assert _parse_review_verdict(raw) is None

    def test_missing_reasoning_key(self):
        raw = json.dumps({
            "verdict": "approve",
            "risk_flags": [],
        })
        assert _parse_review_verdict(raw) is None

    def test_missing_risk_flags_key(self):
        raw = json.dumps({
            "verdict": "approve",
            "reasoning": "OK.",
        })
        assert _parse_review_verdict(raw) is None

    def test_invalid_verdict_value(self):
        raw = json.dumps({
            "verdict": "maybe",
            "reasoning": "Unsure.",
            "risk_flags": [],
        })
        assert _parse_review_verdict(raw) is None

    def test_risk_flags_not_list(self):
        raw = json.dumps({
            "verdict": "approve",
            "reasoning": "Fine.",
            "risk_flags": "LEAKAGE_RISK",
        })
        assert _parse_review_verdict(raw) is None


# ---------------------------------------------------------------------------
# LEAKAGE_RISK / PHI_RISK override
# ---------------------------------------------------------------------------

class TestCriticalFlagOverride:
    """Critical risk flags force rejection regardless of LLM verdict."""

    def test_leakage_risk_forces_reject(self, monkeypatch):
        """Even if the LLM says 'approve', LEAKAGE_RISK must force reject."""
        api_response = json.dumps({
            "verdict": "approve",
            "reasoning": "Looks fine to me.",
            "risk_flags": ["LEAKAGE_RISK"],
        })
        monkeypatch.setattr(
            "improvement_loop.agents.reviewer.api_call_with_retry",
            lambda kwargs: api_response,
        )
        finding = _make_finding()
        verdict = review(finding, "old", "new")
        assert verdict.verdict == "reject"
        assert "LEAKAGE_RISK" in verdict.risk_flags

    def test_phi_risk_forces_reject(self, monkeypatch):
        """PHI_RISK must also force reject."""
        api_response = json.dumps({
            "verdict": "request_changes",
            "reasoning": "Minor issue.",
            "risk_flags": ["PHI_RISK"],
        })
        monkeypatch.setattr(
            "improvement_loop.agents.reviewer.api_call_with_retry",
            lambda kwargs: api_response,
        )
        finding = _make_finding()
        verdict = review(finding, "old", "new")
        assert verdict.verdict == "reject"
        assert "PHI_RISK" in verdict.risk_flags

    def test_non_critical_flag_does_not_override(self, monkeypatch):
        """DEPS_MODIFIED alone should not override an approve verdict."""
        api_response = json.dumps({
            "verdict": "approve",
            "reasoning": "Change is fine.",
            "risk_flags": ["DEPS_MODIFIED"],
        })
        monkeypatch.setattr(
            "improvement_loop.agents.reviewer.api_call_with_retry",
            lambda kwargs: api_response,
        )
        finding = _make_finding()
        verdict = review(finding, "old", "new")
        assert verdict.verdict == "approve"
        assert "DEPS_MODIFIED" in verdict.risk_flags


# ---------------------------------------------------------------------------
# _generate_diff()
# ---------------------------------------------------------------------------

class TestGenerateDiff:
    """Diff generation produces expected unified diff output."""

    def test_basic_diff(self):
        original = "line1\nline2\nline3\n"
        new = "line1\nline2_modified\nline3\n"
        diff = _generate_diff(original, new, filename="test.m")
        assert "--- a/test.m" in diff
        assert "+++ b/test.m" in diff
        assert "-line2" in diff
        assert "+line2_modified" in diff

    def test_no_diff_for_identical(self):
        content = "same\ncontent\n"
        diff = _generate_diff(content, content, filename="test.m")
        assert diff == ""

    def test_addition_only(self):
        original = "line1\n"
        new = "line1\nline2\n"
        diff = _generate_diff(original, new, filename="added.py")
        assert "+line2" in diff

    def test_deletion_only(self):
        original = "line1\nline2\n"
        new = "line1\n"
        diff = _generate_diff(original, new, filename="removed.py")
        assert "-line2" in diff


# ---------------------------------------------------------------------------
# review() — parse failure fallback
# ---------------------------------------------------------------------------

class TestReviewParseFallback:
    """review() returns request_changes with EVALUATION_FAILED on parse failure."""

    def test_fallback_on_garbage(self, monkeypatch):
        monkeypatch.setattr(
            "improvement_loop.agents.reviewer.api_call_with_retry",
            lambda kwargs: "this is not valid json at all",
        )
        finding = _make_finding()
        verdict = review(finding, "old", "new")
        assert verdict.verdict == "request_changes"
        assert "EVALUATION_FAILED" in verdict.risk_flags
        assert "parse failed" in verdict.reasoning.lower()
