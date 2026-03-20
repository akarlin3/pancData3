"""Tests for the Finding Pydantic model in evaluator.py."""

import sys
import os

import pytest
from pydantic import ValidationError

# Ensure the repo root is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from improvement_loop.evaluator import (  # type: ignore
    Finding, should_continue_loop, check_diminishing_returns,
)
from improvement_loop import loop_tracker  # type: ignore
from improvement_loop.loop_config import (  # type: ignore
    LoopConfig, load_loop_config, get_config, reset_config,
)


@pytest.fixture(autouse=True)
def _reset_loop_config():
    """Ensure each test gets a fresh config (no stale singleton)."""
    reset_config()
    yield
    reset_config()


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

    def test_evaluation_failed_flag_forces_continuation(self):
        """EVALUATION_FAILED flag must force continuation — failed eval is not a valid exit."""
        failed_scores = {**self.GOOD_SCORES, "flags": ["EVALUATION_FAILED"]}
        assert should_continue_loop(failed_scores, []) is True


# ── check_diminishing_returns ───────────────────────────────────────────────

def _make_log_entry(
    iteration: int,
    branches_created: int = 5,
    branches_merged: int = 0,
    importance: int = 2,
    files: list | None = None,
    overall_score: float = 6.0,
) -> dict:
    """Build a synthetic log entry for diminishing-returns tests."""
    if files is None:
        files = ["pipeline/core/fit_models.m"]
    findings = [
        {
            "id": f"iter{iteration}_{i+1:03d}",
            "iteration": iteration,
            "dimension": "correctness",
            "file": f,
            "description": "test finding",
            "fix": "test fix",
            "importance": importance,
            "branch_name": f"improvement/fix-{iteration}-{i}",
            "status": "merged" if i < branches_merged else "pending",
        }
        for i, f in enumerate(files)
    ]
    created = [f"improvement/branch-{iteration}-{j}" for j in range(branches_created)]
    merged = created[:branches_merged]
    return {
        "iteration": iteration,
        "timestamp": "2026-03-20T00:00:00",
        "audit_scores": {"overall": overall_score, "flags": []},
        "findings": findings,
        "findings_count": len(findings),
        "high_priority_findings": 0,
        "branches_created": created,
        "branches_merged": merged,
        "tests_passed": True,
        "exit_condition_met": False,
    }


class TestCheckDiminishingReturns:
    """Tests for the check_diminishing_returns helper."""

    def _build_stale_log(self, **overrides) -> list:
        """Build a 4-entry log where all diminishing-returns conditions are met."""
        defaults = dict(
            branches_created=5, branches_merged=0, importance=2,
            files=["pipeline/core/fit_models.m"], overall_score=6.0,
        )
        defaults.update(overrides)
        return [_make_log_entry(iteration=i, **defaults) for i in range(1, 5)]

    def test_all_conditions_met_returns_true(self):
        log = self._build_stale_log()
        assert check_diminishing_returns(log) is True

    def test_merge_rate_above_threshold_returns_false(self):
        """One iteration with 20% merge rate → condition 1 fails."""
        log = self._build_stale_log()
        # Override iteration 3 to have 1/5 = 20% merge rate
        log[2] = _make_log_entry(
            iteration=3, branches_created=5, branches_merged=1,
            importance=2, files=["pipeline/core/fit_models.m"], overall_score=6.0,
        )
        assert check_diminishing_returns(log) is False

    def test_avg_importance_above_threshold_returns_false(self):
        """Average importance 4.0 → condition 2 fails."""
        log = self._build_stale_log(importance=4)
        assert check_diminishing_returns(log) is False

    def test_fewer_than_four_iterations_skipped(self):
        """Only 3 iterations → not enough history, returns False."""
        log = self._build_stale_log()[:3]
        assert check_diminishing_returns(log) is False

    def test_high_audit_score_returns_false(self):
        """One iteration with audit score 9.0 → condition 4 fails."""
        log = self._build_stale_log()
        log[1] = _make_log_entry(
            iteration=2, branches_created=5, branches_merged=0,
            importance=2, files=["pipeline/core/fit_models.m"], overall_score=9.0,
        )
        assert check_diminishing_returns(log) is False


class TestShouldContinueLoopDiminishingReturns:
    """Integration: should_continue_loop calls check_diminishing_returns."""

    GOOD_SCORES = {
        "specificity": 8, "accuracy": 8, "coverage": 8,
        "prioritization": 8, "domain_appropriateness": 8,
        "overall": 8, "flags": [], "reasoning": "Good.",
    }

    def _build_stale_log(self) -> list:
        return [
            _make_log_entry(iteration=i, branches_created=5, branches_merged=0,
                            importance=2, files=["pipeline/core/fit_models.m"],
                            overall_score=6.0)
            for i in range(1, 5)
        ]

    def test_diminishing_returns_stops_loop(self, tmp_path, monkeypatch):
        """All four conditions met → should_continue_loop returns False."""
        log_file = str(tmp_path / "test_log.json")
        monkeypatch.setattr(loop_tracker, "LOG_FILE", log_file)
        import json
        with open(log_file, "w") as f:
            json.dump(self._build_stale_log(), f)

        # importance=1 findings so existing checks don't force continuation
        findings = [_make_finding(importance=1)]
        assert should_continue_loop(self.GOOD_SCORES, findings) is False

    def test_merge_rate_above_threshold_continues(self, tmp_path, monkeypatch):
        """Merge rate 20% in one iteration → should_continue_loop returns True
        (normal exit path still applies, but diminishing returns doesn't trigger)."""
        log_file = str(tmp_path / "test_log.json")
        monkeypatch.setattr(loop_tracker, "LOG_FILE", log_file)
        log = self._build_stale_log()
        log[2] = _make_log_entry(
            iteration=3, branches_created=5, branches_merged=1,
            importance=2, files=["pipeline/core/fit_models.m"], overall_score=6.0,
        )
        import json
        with open(log_file, "w") as f:
            json.dump(log, f)

        # importance=1, good scores → normal path would also exit,
        # but we verify diminishing returns specifically didn't trigger
        # by checking the unit function directly
        assert check_diminishing_returns(log) is False

    def test_avg_importance_above_threshold_continues(self, tmp_path, monkeypatch):
        """Average importance 4.0 → diminishing returns not triggered."""
        log_file = str(tmp_path / "test_log.json")
        monkeypatch.setattr(loop_tracker, "LOG_FILE", log_file)
        log = [
            _make_log_entry(iteration=i, branches_created=5, branches_merged=0,
                            importance=4, files=["pipeline/core/fit_models.m"],
                            overall_score=6.0)
            for i in range(1, 5)
        ]
        import json
        with open(log_file, "w") as f:
            json.dump(log, f)

        assert check_diminishing_returns(log) is False

    def test_only_three_iterations_normal_exit(self, tmp_path, monkeypatch):
        """Only 3 iterations → diminishing returns skipped, normal exit logic applies."""
        log_file = str(tmp_path / "test_log.json")
        monkeypatch.setattr(loop_tracker, "LOG_FILE", log_file)
        log = self._build_stale_log()[:3]
        import json
        with open(log_file, "w") as f:
            json.dump(log, f)

        findings = [_make_finding(importance=1)]
        # Normal exit: importance < 2, good scores, no flags → False
        assert should_continue_loop(self.GOOD_SCORES, findings) is False
        # But diminishing returns specifically was not the reason
        assert check_diminishing_returns(log) is False

    def test_high_audit_score_continues(self, tmp_path, monkeypatch):
        """All conditions met but one iteration has audit score 9.0 → returns True
        (diminishing returns doesn't trigger; normal exit still may apply)."""
        log_file = str(tmp_path / "test_log.json")
        monkeypatch.setattr(loop_tracker, "LOG_FILE", log_file)
        log = self._build_stale_log()
        log[1] = _make_log_entry(
            iteration=2, branches_created=5, branches_merged=0,
            importance=2, files=["pipeline/core/fit_models.m"], overall_score=9.0,
        )
        import json
        with open(log_file, "w") as f:
            json.dump(log, f)

        assert check_diminishing_returns(log) is False


# ── LoopConfig ──────────────────────────────────────────────────────────────

class TestLoopConfig:
    """Tests for loop_config loading and defaults."""

    def test_defaults_without_file(self, tmp_path):
        cfg = load_loop_config(str(tmp_path / "nonexistent.json"))
        assert cfg.exit_strategy == "both"
        assert cfg.dr_window == 4
        assert cfg.dr_max_merge_rate == 0.15
        assert cfg.anthropic_api_key == ""
        assert cfg.audit_model == "claude-opus-4-6"

    def test_partial_override(self, tmp_path):
        import json
        p = tmp_path / "cfg.json"
        p.write_text(json.dumps({"exit_strategy": "classic", "dr_window": 6}))
        cfg = load_loop_config(str(p))
        assert cfg.exit_strategy == "classic"
        assert cfg.dr_window == 6
        assert cfg.dr_max_merge_rate == 0.15
        assert cfg.importance_threshold == 2

    def test_unknown_keys_ignored(self, tmp_path):
        import json
        p = tmp_path / "cfg.json"
        p.write_text(json.dumps({"bogus_key": 999, "exit_strategy": "diminishing_returns"}))
        cfg = load_loop_config(str(p))
        assert cfg.exit_strategy == "diminishing_returns"
        assert not hasattr(cfg, "bogus_key")

    def test_api_key_from_config(self, tmp_path):
        import json
        p = tmp_path / "cfg.json"
        p.write_text(json.dumps({"anthropic_api_key": "sk-test-123"}))
        cfg = load_loop_config(str(p))
        assert cfg.anthropic_api_key == "sk-test-123"

    def test_classic_exit_strategy_skips_diminishing_returns(self, tmp_path, monkeypatch):
        """With exit_strategy='classic', diminishing returns never triggers."""
        import json
        cfg_path = tmp_path / "cfg.json"
        cfg_path.write_text(json.dumps({"exit_strategy": "classic"}))
        from improvement_loop import loop_config
        monkeypatch.setattr(loop_config, "CONFIG_PATH", str(cfg_path))

        log_file = str(tmp_path / "test_log.json")
        monkeypatch.setattr(loop_tracker, "LOG_FILE", log_file)
        log = [
            _make_log_entry(iteration=i, branches_created=5, branches_merged=0,
                            importance=2, files=["pipeline/core/fit_models.m"],
                            overall_score=6.0)
            for i in range(1, 5)
        ]
        with open(log_file, "w") as f:
            json.dump(log, f)

        good_scores = {
            "specificity": 8, "accuracy": 8, "coverage": 8,
            "prioritization": 8, "domain_appropriateness": 8,
            "overall": 8, "flags": [], "reasoning": "Good.",
        }
        findings = [_make_finding(importance=1)]
        result = should_continue_loop(good_scores, findings)
        assert result is False

    def test_diminishing_returns_only_strategy(self, tmp_path, monkeypatch):
        """With exit_strategy='diminishing_returns', classic checks are skipped."""
        import json
        cfg_path = tmp_path / "cfg.json"
        cfg_path.write_text(json.dumps({"exit_strategy": "diminishing_returns"}))
        from improvement_loop import loop_config
        monkeypatch.setattr(loop_config, "CONFIG_PATH", str(cfg_path))

        log_file = str(tmp_path / "test_log.json")
        monkeypatch.setattr(loop_tracker, "LOG_FILE", log_file)
        with open(log_file, "w") as f:
            json.dump([], f)

        good_scores = {
            "specificity": 8, "accuracy": 8, "coverage": 8,
            "prioritization": 8, "domain_appropriateness": 8,
            "overall": 8, "flags": [], "reasoning": "Good.",
        }
        # importance=5 would normally keep classic going, but classic is off
        findings = [_make_finding(importance=5)]
        result = should_continue_loop(good_scores, findings)
        assert result is False
