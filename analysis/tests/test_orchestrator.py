"""
Smoke tests for orchestrator_v1.py.

Verifies dry-run mode: audit is called, findings are logged,
and no git operations are performed.
"""
import json
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Ensure repo root is importable
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

import orchestrator_v1  # noqa: E402
from evaluator import Finding  # noqa: E402


# ── Fixtures ────────────────────────────────────────────────────────────────

SAMPLE_FINDINGS_JSON = json.dumps([
    {
        "dimension": "correctness",
        "file": "pipeline/core/fit_models.m",
        "function_name": "fit_adc",
        "description": "ADC fit does not check for empty b-value vector",
        "fix": "Add an early-return guard when bvals is empty",
        "importance": 5,
        "branch_name": "improvement/guard-empty-bvals",
        "status": "pending",
    },
    {
        "dimension": "code_quality",
        "file": "pipeline/utils/text_progress_bar.m",
        "function_name": None,
        "description": "Progress bar rebuilds string via concatenation each tick",
        "fix": "Pre-allocate a char buffer and overwrite in-place",
        "importance": 2,
        "branch_name": "improvement/prealloc-progress-bar",
        "status": "pending",
    },
])

SAMPLE_FINDINGS = [
    Finding(**item)
    for item in json.loads(SAMPLE_FINDINGS_JSON)
]


def _make_mock_response(text: str) -> MagicMock:
    """Build a mock Anthropic API response."""
    block = MagicMock()
    block.text = text
    resp = MagicMock()
    resp.content = [block]
    return resp


# ── Tests ───────────────────────────────────────────────────────────────────


class TestOrchestratorDryRun:
    """Orchestrator dry-run mode should audit + log but never touch git."""

    @patch.dict(os.environ, {"ANTHROPIC_API_KEY": "test-key-fake"})
    @patch("orchestrator_v1.loop_tracker")
    @patch("orchestrator_v1.git_utils")
    @patch("orchestrator_v1.anthropic.Anthropic")
    def test_dry_run_calls_log_and_skips_git(
        self, mock_anthropic_cls, mock_git, mock_tracker
    ):
        """
        In dry_run mode:
        - The API is called for the audit (Phase 1)
        - loop_tracker.log_iteration is called with the findings
        - No git operations (create_branch, commit_all, merge_branch) are called
        """
        # Set up the mock Anthropic client
        mock_client = MagicMock()
        mock_anthropic_cls.return_value = mock_client

        # Phase 1 audit returns our sample findings
        mock_client.messages.create.return_value = _make_mock_response(
            SAMPLE_FINDINGS_JSON
        )

        # Phase 0 context
        mock_tracker.get_context_for_next_iteration.return_value = (
            "This is the first iteration. No prior history."
        )

        # Phase 5 log returns an entry with exit_condition_met=True
        mock_tracker.log_iteration.return_value = {
            "iteration": 1,
            "exit_condition_met": True,
            "audit_scores": {"overall": 8.0},
            "findings_count": 2,
            "high_priority_findings": 1,
            "branches_created": [],
            "branches_merged": [],
            "tests_passed": True,
        }

        mock_tracker.print_full_summary = MagicMock()

        # Run
        orchestrator_v1.run_loop(max_iterations=1, dry_run=True, single_iteration=True)

        # Verify Phase 1: API was called for audit
        assert mock_client.messages.create.call_count >= 1

        # Verify Phase 5: log_iteration was called
        mock_tracker.log_iteration.assert_called_once()
        call_kwargs = mock_tracker.log_iteration.call_args
        logged_findings = call_kwargs.kwargs.get("findings") or call_kwargs[1].get("findings")
        # If called positionally
        if logged_findings is None:
            logged_findings = call_kwargs[0][1] if len(call_kwargs[0]) > 1 else None
        assert logged_findings is not None
        assert len(logged_findings) == 2
        assert logged_findings[0].dimension == "correctness"
        assert logged_findings[1].dimension == "code_quality"

        # Verify NO git operations in dry_run
        mock_git.create_branch.assert_not_called()
        mock_git.commit_all.assert_not_called()
        mock_git.merge_branch.assert_not_called()
        mock_git.switch_branch.assert_not_called()

    @patch.dict(os.environ, {"ANTHROPIC_API_KEY": "test-key-fake"})
    @patch("orchestrator_v1.loop_tracker")
    @patch("orchestrator_v1.git_utils")
    @patch("orchestrator_v1.anthropic.Anthropic")
    def test_dry_run_single_iteration_exits(
        self, mock_anthropic_cls, mock_git, mock_tracker
    ):
        """Single-iteration mode exits after one loop even if exit condition is not met."""
        mock_client = MagicMock()
        mock_anthropic_cls.return_value = mock_client

        mock_client.messages.create.return_value = _make_mock_response(
            SAMPLE_FINDINGS_JSON
        )
        mock_tracker.get_context_for_next_iteration.return_value = "First iteration."
        mock_tracker.log_iteration.return_value = {
            "iteration": 1,
            "exit_condition_met": False,  # would normally continue
            "audit_scores": {"overall": 5.0},
            "findings_count": 2,
            "high_priority_findings": 1,
            "branches_created": [],
            "branches_merged": [],
            "tests_passed": True,
        }
        mock_tracker.print_full_summary = MagicMock()

        # single_iteration=True should cause exit after 1 loop
        orchestrator_v1.run_loop(max_iterations=20, dry_run=True, single_iteration=True)

        # Only one call to log_iteration (one iteration)
        assert mock_tracker.log_iteration.call_count == 1


class TestDependenciesGuard:
    """Findings targeting pipeline/dependencies/ must be rejected."""

    def test_guard_raises_on_dependencies_path(self):
        finding = Finding(
            dimension="correctness",
            file="pipeline/dependencies/some_script.m",
            description="Modify dependency",
            fix="Change it",
            importance=3,
            branch_name="improvement/bad-deps-change",
        )
        with pytest.raises(RuntimeError, match="pipeline/dependencies/"):
            orchestrator_v1._guard_dependencies([finding])

    def test_guard_passes_for_normal_paths(self):
        finding = Finding(
            dimension="correctness",
            file="pipeline/core/fit_models.m",
            description="Improve fitting",
            fix="Add guard",
            importance=5,
            branch_name="improvement/fit-guard",
        )
        # Should not raise
        orchestrator_v1._guard_dependencies([finding])


class TestApiKeyCheck:
    """Orchestrator must exit if ANTHROPIC_API_KEY is not set."""

    @patch.dict(os.environ, {}, clear=True)
    def test_missing_api_key_exits(self):
        # Remove the key if present
        os.environ.pop("ANTHROPIC_API_KEY", None)
        with pytest.raises(SystemExit) as exc_info:
            orchestrator_v1._check_api_key()
        assert exc_info.value.code == 1


class TestPhasePlanning:
    """Phase 2 planning groups findings by file overlap."""

    def test_non_overlapping_files_same_group(self):
        f1 = Finding(
            dimension="correctness",
            file="pipeline/core/a.m",
            description="Fix A",
            fix="Fix",
            importance=5,
            branch_name="improvement/fix-a",
        )
        f2 = Finding(
            dimension="code_quality",
            file="pipeline/core/b.m",
            description="Fix B",
            fix="Fix",
            importance=3,
            branch_name="improvement/fix-b",
        )
        groups = orchestrator_v1.phase2_plan([f1, f2])
        # Non-overlapping files can go in same group
        assert len(groups) == 1
        assert len(groups[0]) == 2

    def test_overlapping_files_separate_groups(self):
        f1 = Finding(
            dimension="correctness",
            file="pipeline/core/a.m",
            description="Fix A1",
            fix="Fix",
            importance=5,
            branch_name="improvement/fix-a1",
        )
        f2 = Finding(
            dimension="code_quality",
            file="pipeline/core/a.m",
            description="Fix A2",
            fix="Fix",
            importance=3,
            branch_name="improvement/fix-a2",
        )
        groups = orchestrator_v1.phase2_plan([f1, f2])
        assert len(groups) == 2
