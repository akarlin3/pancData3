"""Tests for improvement_loop.agents.implementer module."""

import os
import sys

import pytest

# Ensure repo root is on sys.path for absolute imports
_repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from improvement_loop.agents.implementer import (
    apply_fix,
    get_fix_system_prompt,
    DEFAULT_FIX_PROMPT,
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
# get_fix_system_prompt
# ---------------------------------------------------------------------------

class TestGetFixSystemPrompt:
    """Verify system prompt retrieval."""

    def test_returns_string(self):
        prompt = get_fix_system_prompt()
        assert isinstance(prompt, str)
        assert len(prompt) > 0

    def test_default_prompt_content(self):
        assert "code fixer" in DEFAULT_FIX_PROMPT.lower()


# ---------------------------------------------------------------------------
# apply_fix — file not found
# ---------------------------------------------------------------------------

class TestApplyFixFileNotFound:
    """apply_fix prints a warning and returns when the target file doesn't exist."""

    def test_missing_file_warns(self, tmp_path, capsys):
        finding = _make_finding(file="nonexistent/path/to/file.m")
        apply_fix(finding, repo_root=str(tmp_path))
        captured = capsys.readouterr()
        assert "not found" in captured.out.lower() or "WARNING" in captured.out


# ---------------------------------------------------------------------------
# apply_fix — successful fix (mocked API)
# ---------------------------------------------------------------------------

class TestApplyFixMockedAPI:
    """apply_fix reads, calls API, and writes the file."""

    def test_updates_file(self, tmp_path, monkeypatch):
        # Create a target file
        target_dir = tmp_path / "pipeline" / "core"
        target_dir.mkdir(parents=True)
        target_file = target_dir / "fit_models.m"
        target_file.write_text("% original content\n", encoding="utf-8")

        expected_content = "% fixed content\nfunction y = fit_adc(x)\n  y = x + 1;\nend\n"
        monkeypatch.setattr(
            "improvement_loop.agents.implementer.api_call_with_retry",
            lambda kwargs: expected_content,
        )
        finding = _make_finding()
        apply_fix(finding, repo_root=str(tmp_path))
        assert target_file.read_text(encoding="utf-8") == expected_content

    def test_passes_finding_context(self, tmp_path, monkeypatch):
        """Verify the API call includes the finding's file, description, and fix."""
        target_dir = tmp_path / "pipeline" / "core"
        target_dir.mkdir(parents=True)
        target_file = target_dir / "fit_models.m"
        target_file.write_text("% original\n", encoding="utf-8")

        captured = {}

        def mock_api(kwargs):
            captured.update(kwargs)
            return "% fixed"

        monkeypatch.setattr(
            "improvement_loop.agents.implementer.api_call_with_retry",
            mock_api,
        )
        finding = _make_finding(
            file="pipeline/core/fit_models.m",
            description="Missing default for new_field",
            fix="Add isfield check with default value.",
        )
        apply_fix(finding, repo_root=str(tmp_path))
        user_content = captured["messages"][0]["content"]
        assert "fit_models.m" in user_content
        assert "Missing default for new_field" in user_content
        assert "Add isfield check" in user_content
        assert "% original" in user_content
