"""Tests for orchestrator_v2: agent pipeline, state dataclasses, dry-run, rejection."""

import os
import sys
import json
import pytest
from unittest.mock import patch, MagicMock

# Add repo root so modules are importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from improvement_loop import loop_tracker
from improvement_loop import orchestrator_v2
from improvement_loop.orchestrator_v2 import (
    FindingState,
    IterationState,
    run_loop,
    _print_agent_summary,
)
from improvement_loop.evaluator import Finding
from improvement_loop.agents.implementer import ImplementResult
from improvement_loop.agents.reviewer import ReviewVerdict


@pytest.fixture(autouse=True)
def _isolate_log(tmp_path, monkeypatch):
    """Point loop_tracker at a temp log file for every test."""
    log_file = str(tmp_path / "test_log.json")
    monkeypatch.setattr(loop_tracker, "LOG_FILE", log_file)


def _make_finding(**overrides) -> Finding:
    defaults = {
        "dimension": "correctness",
        "file": "pipeline/core/fit_models.m",
        "description": "Off-by-one in loop bound",
        "fix": "Change < to <=.",
        "importance": 5,
        "branch_name": "improvement/test-fix",
    }
    defaults.update(overrides)
    return Finding(**defaults)


# ---------------------------------------------------------------------------
# Dataclass construction
# ---------------------------------------------------------------------------

class TestFindingState:
    def test_defaults(self):
        f = _make_finding()
        fs = FindingState(finding=f)
        assert fs.implement_result is None
        assert fs.review_verdict is None
        assert fs.merged is False
        assert fs.error == ""

    def test_with_all_fields(self):
        f = _make_finding()
        ir = ImplementResult(success=True, original_content="a", new_content="b")
        rv = ReviewVerdict(verdict="approve", reasoning="ok", risk_flags=[])
        fs = FindingState(finding=f, implement_result=ir, review_verdict=rv,
                          merged=True, error="")
        assert fs.merged is True
        assert fs.review_verdict.verdict == "approve"


class TestIterationState:
    def test_defaults(self):
        state = IterationState(iteration=1, context="ctx")
        assert state.iteration == 1
        assert state.context == "ctx"
        assert state.findings == []
        assert state.all_tests_passed is True

    def test_with_findings(self):
        f = _make_finding()
        fs = FindingState(finding=f)
        state = IterationState(iteration=2, context="ctx", findings=[fs])
        assert len(state.findings) == 1


# ---------------------------------------------------------------------------
# run_loop dry_run — full pipeline no-op
# ---------------------------------------------------------------------------

class TestRunLoopDryRun:
    """In dry_run mode, audit returns [], all phases are no-ops, exits immediately."""

    def test_completes_without_error(self):
        entries = run_loop(max_iterations=5, dry_run=True)
        assert isinstance(entries, list)
        assert len(entries) >= 1

    def test_exit_condition_met(self):
        entries = run_loop(max_iterations=5, dry_run=True)
        assert entries[0]["exit_condition_met"] is True

    def test_returns_single_entry(self):
        """Dry run should exit on first iteration."""
        entries = run_loop(max_iterations=10, dry_run=True)
        assert len(entries) == 1


# ---------------------------------------------------------------------------
# run_loop — exit on first iteration (mocked log_iteration)
# ---------------------------------------------------------------------------

class TestRunLoopExitCondition:
    def test_stops_on_first_exit(self):
        exit_entry = {
            "iteration": 1,
            "exit_condition_met": True,
            "audit_scores": {"overall": 9.0, "flags": []},
            "findings": [],
            "findings_count": 0,
            "high_priority_findings": 0,
            "branches_created": [],
            "branches_merged": [],
            "tests_passed": True,
        }
        with patch.object(
            loop_tracker, "log_iteration", return_value=exit_entry
        ) as mock_log:
            entries = run_loop(max_iterations=10, dry_run=True)
        assert mock_log.call_count == 1
        assert len(entries) == 1
        assert entries[0]["exit_condition_met"] is True

    def test_stops_after_three_iterations(self):
        continuing_entry = {
            "iteration": 1,
            "exit_condition_met": False,
            "audit_scores": {"overall": 5.0, "flags": []},
            "findings": [],
            "findings_count": 0,
            "high_priority_findings": 0,
            "branches_created": [],
            "branches_merged": [],
            "tests_passed": True,
        }
        call_count = 0

        def side_effect(*args, **kwargs):
            nonlocal call_count
            call_count += 1
            if call_count <= 2:
                return {**continuing_entry, "iteration": call_count}
            return {**continuing_entry, "iteration": 3, "exit_condition_met": True}

        with patch.object(
            loop_tracker, "log_iteration", side_effect=side_effect
        ):
            entries = run_loop(max_iterations=10, dry_run=True)
        assert len(entries) == 3
        assert entries[-1]["exit_condition_met"] is True


# ---------------------------------------------------------------------------
# Rejected finding does not get merged
# ---------------------------------------------------------------------------

class TestRejectedFindingNotMerged:
    """When the reviewer rejects a finding, it must not be merged."""

    def test_rejected_finding_status_pending(self, monkeypatch):
        """Mock the full pipeline: audit returns 1 finding, implement succeeds,
        reviewer rejects → finding.status should be 'pending', not 'merged'."""
        finding = _make_finding(branch_name="improvement/reject-test")

        # Mock audit to return one finding
        monkeypatch.setattr(
            "improvement_loop.orchestrator_v2._audit",
            lambda iteration, context, dry_run: [finding],
        )

        # Mock implement to succeed
        monkeypatch.setattr(
            "improvement_loop.orchestrator_v2._implement",
            lambda f, base_branch, dry_run: ImplementResult(
                success=True, original_content="old", new_content="new"
            ),
        )

        # Mock reviewer to reject
        monkeypatch.setattr(
            "improvement_loop.orchestrator_v2._review",
            lambda f, original_content, new_content, dry_run: ReviewVerdict(
                verdict="reject", reasoning="Bad change", risk_flags=["DEPS_MODIFIED"]
            ),
        )

        # Mock git operations
        monkeypatch.setattr(
            "improvement_loop.orchestrator_v2.git_utils.current_branch",
            lambda: "v2.2-dev",
        )
        monkeypatch.setattr(
            "improvement_loop.orchestrator_v2.git_utils.checkout",
            lambda branch: None,
        )
        monkeypatch.setattr(
            "improvement_loop.orchestrator_v2.git_utils.run_python_tests",
            lambda: True,
        )
        monkeypatch.setattr(
            "improvement_loop.orchestrator_v2.git_utils.merge_branch",
            lambda source, target, delete_after: None,
        )

        # Make it exit after one iteration
        exit_entry = {
            "iteration": 1,
            "exit_condition_met": True,
            "audit_scores": {"overall": 5.0, "flags": []},
            "findings": [],
            "findings_count": 1,
            "high_priority_findings": 0,
            "branches_created": ["improvement/reject-test"],
            "branches_merged": [],
            "tests_passed": True,
        }
        monkeypatch.setattr(
            loop_tracker, "log_iteration",
            lambda **kwargs: exit_entry,
        )

        entries = run_loop(max_iterations=1, dry_run=False)

        # The finding should not have been merged
        assert finding.status == "pending"
        assert len(entries) == 1


# ---------------------------------------------------------------------------
# _print_agent_summary — smoke test
# ---------------------------------------------------------------------------

class TestPrintAgentSummary:
    """_print_agent_summary runs without error on various states."""

    def test_empty_state(self, capsys):
        state = IterationState(iteration=1, context="")
        _print_agent_summary(state)
        captured = capsys.readouterr()
        assert "Audit:" in captured.out
        assert "0 findings" in captured.out

    def test_with_findings(self, capsys):
        f = _make_finding(importance=8)
        ir = ImplementResult(success=True, original_content="a", new_content="b")
        rv = ReviewVerdict(verdict="approve", reasoning="ok", risk_flags=[])
        fs = FindingState(finding=f, implement_result=ir, review_verdict=rv,
                          merged=True)
        state = IterationState(iteration=1, context="", findings=[fs])
        _print_agent_summary(state)
        captured = capsys.readouterr()
        assert "1 findings (1 high-priority)" in captured.out
        assert "1 succeeded" in captured.out
        assert "1 approved" in captured.out
        assert "1 merged" in captured.out

    def test_with_risk_flags(self, capsys):
        f = _make_finding()
        ir = ImplementResult(success=True, original_content="a", new_content="b")
        rv = ReviewVerdict(verdict="reject", reasoning="bad",
                           risk_flags=["LEAKAGE_RISK", "PHI_RISK"])
        fs = FindingState(finding=f, implement_result=ir, review_verdict=rv,
                          error="Rejected")
        state = IterationState(iteration=1, context="", findings=[fs])
        _print_agent_summary(state)
        captured = capsys.readouterr()
        assert "LEAKAGE_RISK" in captured.out
        assert "PHI_RISK" in captured.out
