"""Tests for improvement_loop.agents.implementer module."""

import os
import sys

import pytest

# Ensure repo root is on sys.path for absolute imports
_repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from improvement_loop.agents.implementer import (
    ImplementResult,
    _generate_fix,
    implement,
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
        "branch_name": "improvement/test-fix",
    }
    defaults.update(overrides)
    return Finding(**defaults)


# ---------------------------------------------------------------------------
# ImplementResult dataclass
# ---------------------------------------------------------------------------

class TestImplementResult:
    """Basic sanity checks on the ImplementResult dataclass."""

    def test_success_defaults(self):
        r = ImplementResult(success=True, original_content="a", new_content="b")
        assert r.success is True
        assert r.error == ""

    def test_failure_with_error(self):
        r = ImplementResult(success=False, original_content="", new_content="",
                            error="file not found")
        assert r.success is False
        assert "file not found" in r.error


# ---------------------------------------------------------------------------
# implement() — dry_run mode
# ---------------------------------------------------------------------------

class TestImplementDryRun:
    """implement() in dry_run mode returns success with empty contents."""

    def test_returns_success(self):
        finding = _make_finding()
        result = implement(finding, base_branch="main", dry_run=True)
        assert result.success is True

    def test_empty_contents(self):
        finding = _make_finding()
        result = implement(finding, base_branch="main", dry_run=True)
        assert result.original_content == ""
        assert result.new_content == ""

    def test_no_error(self):
        finding = _make_finding()
        result = implement(finding, base_branch="main", dry_run=True)
        assert result.error == ""


# ---------------------------------------------------------------------------
# implement() — file not found
# ---------------------------------------------------------------------------

class TestImplementFileNotFound:
    """implement() returns failure when the target file doesn't exist."""

    def test_returns_failure(self, monkeypatch):
        # Ensure branch_exists returns False so we reach the file check
        monkeypatch.setattr(
            "improvement_loop.agents.implementer.git_utils.branch_exists",
            lambda name: False,
        )
        finding = _make_finding(file="nonexistent/path/to/file.m")
        result = implement(finding, base_branch="main")
        assert result.success is False
        assert "file not found" in result.error

    def test_empty_contents_on_missing_file(self, monkeypatch):
        monkeypatch.setattr(
            "improvement_loop.agents.implementer.git_utils.branch_exists",
            lambda name: False,
        )
        finding = _make_finding(file="nonexistent/path/to/file.m")
        result = implement(finding, base_branch="main")
        assert result.original_content == ""
        assert result.new_content == ""


# ---------------------------------------------------------------------------
# implement() — branch already exists
# ---------------------------------------------------------------------------

class TestImplementBranchExists:
    """implement() returns failure when the branch already exists."""

    def test_returns_failure(self, monkeypatch):
        monkeypatch.setattr(
            "improvement_loop.agents.implementer.git_utils.branch_exists",
            lambda name: True,
        )
        finding = _make_finding()
        result = implement(finding, base_branch="main")
        assert result.success is False
        assert "branch exists" in result.error


# ---------------------------------------------------------------------------
# _generate_fix() — mocked API
# ---------------------------------------------------------------------------

class TestGenerateFix:
    """_generate_fix() calls the API and returns the response content."""

    def test_returns_api_response(self, monkeypatch):
        expected_content = "% fixed file content\nfunction y = fit_adc(x)\n  y = x + 1;\nend\n"
        monkeypatch.setattr(
            "improvement_loop.agents.implementer.api_call_with_retry",
            lambda kwargs: expected_content,
        )
        finding = _make_finding()
        result = _generate_fix(finding, "% original content")
        assert result == expected_content

    def test_passes_finding_context(self, monkeypatch):
        """Verify the API call includes the finding's file, description, and fix."""
        captured = {}

        def mock_api(kwargs):
            captured.update(kwargs)
            return "fixed"

        monkeypatch.setattr(
            "improvement_loop.agents.implementer.api_call_with_retry",
            mock_api,
        )
        finding = _make_finding(
            file="pipeline/utils/parse_config.m",
            description="Missing default for new_field",
            fix="Add isfield check with default value.",
        )
        _generate_fix(finding, "% original")
        user_content = captured["messages"][0]["content"]
        assert "parse_config.m" in user_content
        assert "Missing default for new_field" in user_content
        assert "Add isfield check" in user_content
        assert "% original" in user_content
