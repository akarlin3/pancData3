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

    def test_happy_path_with_real_file(self, tmp_path):
        """dry_run=True returns new_content reflecting the fix without writing to disk."""
        # Create a real file so the finding references something that exists
        target_file = tmp_path / "pipeline" / "core" / "fit_models.m"
        target_file.parent.mkdir(parents=True, exist_ok=True)
        original = "% original MATLAB content\nfor i = 1:n\n  x(i) = i;\nend\n"
        target_file.write_text(original)

        finding = _make_finding(file=str(target_file))
        result = implement(finding, base_branch="main", dry_run=True)

        assert result.success is True
        assert result.error == ""
        # In dry_run mode the function should succeed without raising
        assert isinstance(result.new_content, str)

    def test_nonexistent_file_dry_run(self):
        """dry_run=True with a non-existent file still succeeds (no disk access needed)."""
        finding = _make_finding(file="totally/fake/nonexistent_file.m")
        result = implement(finding, base_branch="main", dry_run=True)
        # dry_run should short-circuit before checking file existence
        assert result.success is True
        assert result.error == ""

    def test_no_filesystem_side_effects(self, tmp_path):
        """dry_run=True must not create any files, directories, or branches."""
        # Snapshot the tmp_path contents before the call
        before_contents = set(tmp_path.rglob("*"))

        finding = _make_finding(
            file=str(tmp_path / "should_not_be_created.m"),
            branch_name="improvement/dry-run-no-side-effects",
        )
        result = implement(finding, base_branch="main", dry_run=True)

        # Snapshot after
        after_contents = set(tmp_path.rglob("*"))

        assert result.success is True
        assert before_contents == after_contents, (
            "dry_run=True should not create any filesystem artifacts"
        )

    def test_dry_run_does_not_call_git(self, monkeypatch):
        """dry_run=True must not invoke any git operations."""
        git_called = {"called": False}

        def _fail_on_git_call(*args, **kwargs):
            git_called["called"] = True
            raise AssertionError("git should not be called in dry_run mode")

        # Patch all git_utils functions that implement() might use
        monkeypatch.setattr(
            "improvement_loop.agents.implementer.git_utils.branch_exists",
            _fail_on_git_call,
        )

        finding = _make_finding()
        result = implement(finding, base_branch="main", dry_run=True)

        assert result.success is True
        assert not git_called["called"], "git_utils should not be invoked during dry_run"

    def test_dry_run_does_not_call_api(self, monkeypatch):
        """dry_run=True must not invoke the LLM API."""
        def _fail_on_api_call(*args, **kwargs):
            raise AssertionError("API should not be called in dry_run mode")

        monkeypatch.setattr(
            "improvement_loop.agents.implementer.api_call_with_retry",
            _fail_on_api_call,
        )

        finding = _make_finding()
        result = implement(finding, base_branch="main", dry_run=True)

        assert result.success is True


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
        # Disable RAG to avoid ChromaDB initialization during tests
        monkeypatch.setattr(
            "improvement_loop.agents.implementer._get_loop_config",
            lambda: type("Cfg", (), {"rag_enabled": False, "fix_model": "test", "fix_max_tokens": 100})(),
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
        # Disable RAG to avoid ChromaDB initialization during tests
        monkeypatch.setattr(
            "improvement_loop.agents.implementer._get_loop_config",
            lambda: type("Cfg", (), {"rag_enabled": False, "fix_model": "test", "fix_max_tokens": 100})(),
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