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
    VALID_VERDICTS,
    _parse_review,
    review,
    get_review_system_prompt,
    DEFAULT_REVIEW_PROMPT,
)
from improvement_loop.evaluator import Finding


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
# _parse_review() — valid JSON
# ---------------------------------------------------------------------------

class TestParseReviewValid:
    """_parse_review with well-formed JSON input."""

    def test_approve(self):
        raw = json.dumps({
            "verdict": "APPROVE",
            "summary": "Change looks correct.",
            "issues": [],
            "reasoning": "All good.",
        })
        result = _parse_review(raw)
        assert result is not None
        assert result["verdict"] == "APPROVE"

    def test_reject_with_issues(self):
        raw = json.dumps({
            "verdict": "REJECT",
            "summary": "Modifies dependencies.",
            "issues": ["Touches read-only dir"],
            "reasoning": "Not allowed.",
        })
        result = _parse_review(raw)
        assert result is not None
        assert result["verdict"] == "REJECT"

    def test_request_changes(self):
        raw = json.dumps({
            "verdict": "REQUEST_CHANGES",
            "summary": "Missing edge case.",
            "issues": ["No guard for empty input"],
            "reasoning": "Needs fix.",
        })
        result = _parse_review(raw)
        assert result is not None
        assert result["verdict"] == "REQUEST_CHANGES"

    def test_markdown_fenced(self):
        inner = json.dumps({
            "verdict": "APPROVE",
            "summary": "Looks good.",
            "issues": [],
            "reasoning": "Fine.",
        })
        raw = f"```json\n{inner}\n```"
        result = _parse_review(raw)
        assert result is not None
        assert result["verdict"] == "APPROVE"

    def test_preamble_before_json_returns_none(self):
        """_parse_review does not handle preamble text before JSON."""
        inner = json.dumps({
            "verdict": "APPROVE",
            "summary": "OK.",
            "issues": [],
            "reasoning": "Fine.",
        })
        raw = f"Here is my review:\n{inner}"
        result = _parse_review(raw)
        assert result is None


# ---------------------------------------------------------------------------
# _parse_review() — invalid JSON
# ---------------------------------------------------------------------------

class TestParseReviewInvalid:
    """_parse_review with invalid or incomplete input."""

    def test_garbage_returns_none(self):
        assert _parse_review("not json at all") is None

    def test_empty_string_returns_none(self):
        assert _parse_review("") is None

    def test_invalid_verdict_value(self):
        raw = json.dumps({
            "verdict": "MAYBE",
            "summary": "Unsure.",
            "issues": [],
            "reasoning": "Not sure.",
        })
        assert _parse_review(raw) is None


# ---------------------------------------------------------------------------
# review() — mocked API
# ---------------------------------------------------------------------------

class TestReviewMockedAPI:
    """review() calls the API and returns parsed result."""

    def test_approve_response(self, monkeypatch):
        api_response = json.dumps({
            "verdict": "APPROVE",
            "summary": "Looks fine.",
            "issues": [],
            "reasoning": "Change is correct.",
        })
        monkeypatch.setattr(
            "improvement_loop.agents.reviewer.api_call_with_retry",
            lambda kwargs: api_response,
        )
        finding = _make_finding()
        result = review(finding, "diff text here")
        assert result["verdict"] == "APPROVE"

    def test_reject_response(self, monkeypatch):
        api_response = json.dumps({
            "verdict": "REJECT",
            "summary": "Bad change.",
            "issues": ["Modifies dependencies"],
            "reasoning": "Not allowed.",
        })
        monkeypatch.setattr(
            "improvement_loop.agents.reviewer.api_call_with_retry",
            lambda kwargs: api_response,
        )
        finding = _make_finding()
        result = review(finding, "diff text here")
        assert result["verdict"] == "REJECT"


# ---------------------------------------------------------------------------
# review() — parse failure fallback
# ---------------------------------------------------------------------------

class TestReviewParseFallback:
    """review() returns REQUEST_CHANGES on parse failure."""

    def test_fallback_on_garbage(self, monkeypatch):
        monkeypatch.setattr(
            "improvement_loop.agents.reviewer.api_call_with_retry",
            lambda kwargs: "this is not valid json at all",
        )
        finding = _make_finding()
        result = review(finding, "diff text")
        assert result["verdict"] == "REQUEST_CHANGES"
        assert "parse" in result["summary"].lower() or "parse" in result["reasoning"].lower()


# ---------------------------------------------------------------------------
# get_review_system_prompt
# ---------------------------------------------------------------------------

class TestGetReviewSystemPrompt:
    """Verify system prompt retrieval."""

    def test_returns_string(self):
        prompt = get_review_system_prompt()
        assert isinstance(prompt, str)
        assert len(prompt) > 0

    def test_default_prompt_mentions_reviewer(self):
        assert "reviewer" in DEFAULT_REVIEW_PROMPT.lower()

    def test_valid_verdicts(self):
        assert "APPROVE" in VALID_VERDICTS
        assert "REJECT" in VALID_VERDICTS
        assert "REQUEST_CHANGES" in VALID_VERDICTS
