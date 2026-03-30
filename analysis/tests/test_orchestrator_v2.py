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
    _print_run_summary,
)
from improvement_loop.evaluator import Finding


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
        assert fs.diff == ""
        assert fs.review_verdict == ""
        assert fs.merged is False
        assert fs.error is None

    def test_with_all_fields(self):
        f = _make_finding()
        fs = FindingState(
            finding=f,
            diff="--- a/file\n+++ b/file\n",
            review_verdict="APPROVE",
            review_detail={"verdict": "APPROVE", "reasoning": "ok"},
            merged=True,
            error=None,
        )
        assert fs.merged is True
        assert fs.review_verdict == "APPROVE"


class TestIterationState:
    def test_defaults(self):
        state = IterationState(iteration=1, dry_run=False)
        assert state.iteration == 1
        assert state.dry_run is False
        assert state.finding_states == []
        assert state.all_tests_passed is True

    def test_with_findings(self):
        f = _make_finding()
        fs = FindingState(finding=f)
        state = IterationState(iteration=2, dry_run=True, finding_states=[fs])
        assert len(state.finding_states) == 1


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
# _print_run_summary — smoke test
# ---------------------------------------------------------------------------

class TestPrintRunSummary:
    """_print_run_summary runs without error on various inputs."""

    def test_empty_entries(self, capsys):
        _print_run_summary([])
        captured = capsys.readouterr()
        assert "SUMMARY" in captured.out

    def test_with_entries(self, capsys):
        entries = [{
            "iteration": 1,
            "exit_condition_met": True,
            "audit_scores": {"overall": 8.0, "flags": []},
            "findings": [
                {"dimension": "correctness", "status": "merged"},
                {"dimension": "performance", "status": "pending"},
            ],
            "findings_count": 2,
            "high_priority_findings": 1,
            "branches_created": ["improvement/fix-1"],
            "branches_merged": ["improvement/fix-1"],
            "tests_passed": True,
        }]
        _print_run_summary(entries)
        captured = capsys.readouterr()
        assert "SUMMARY" in captured.out
        assert "1 iteration" in captured.out
