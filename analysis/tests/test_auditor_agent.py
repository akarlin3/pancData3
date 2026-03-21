"""Tests for improvement_loop.agents.auditor module."""

import json
import os
import sys

import pytest

# Ensure repo root is on sys.path for absolute imports
_repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from improvement_loop.agents.auditor import audit, collect_source_files, parse_findings


# ---------------------------------------------------------------------------
# parse_findings
# ---------------------------------------------------------------------------

VALID_FINDING_JSON = json.dumps([
    {
        "dimension": "correctness",
        "file": "pipeline/core/fit_models.m",
        "function_name": "fit_adc",
        "description": "Off-by-one in loop bound at line 42",
        "fix": "Change < to <= in the for-loop condition.",
        "importance": 7,
        "branch_name": "improvement/fix-fit-models-obo",
    },
    {
        "dimension": "security",
        "file": "pipeline/utils/escape_shell_arg.m",
        "description": "Missing backtick escaping on Windows PowerShell",
        "fix": "Add backtick to the replacement list.",
        "importance": 5,
        "branch_name": "improvement/escape-backtick",
    },
])


class TestParseFindingsValidJSON:
    """parse_findings with well-formed JSON input."""

    def test_returns_correct_count(self):
        findings = parse_findings(VALID_FINDING_JSON)
        assert len(findings) == 2

    def test_first_finding_fields(self):
        findings = parse_findings(VALID_FINDING_JSON)
        f = findings[0]
        assert f.dimension == "correctness"
        assert f.file == "pipeline/core/fit_models.m"
        assert f.function_name == "fit_adc"
        assert f.importance == 7

    def test_second_finding_fields(self):
        findings = parse_findings(VALID_FINDING_JSON)
        f = findings[1]
        assert f.dimension == "security"
        assert f.branch_name == "improvement/escape-backtick"


class TestParseFindingsTruncatedJSON:
    """parse_findings with truncated (incomplete) JSON input."""

    def test_recovers_complete_findings(self):
        # Truncate mid-way through the second finding
        truncated = VALID_FINDING_JSON[:VALID_FINDING_JSON.rfind('"security"')]
        findings = parse_findings(truncated)
        # Should recover the first finding at minimum
        assert len(findings) >= 1
        assert findings[0].dimension == "correctness"

    def test_returns_empty_on_garbage(self):
        findings = parse_findings("this is not json at all {{{")
        assert findings == []


class TestParseFindingsMarkdownFenced:
    """parse_findings with markdown-fenced JSON."""

    def test_strips_json_fence(self):
        fenced = f"```json\n{VALID_FINDING_JSON}\n```"
        findings = parse_findings(fenced)
        assert len(findings) == 2

    def test_strips_plain_fence(self):
        fenced = f"```\n{VALID_FINDING_JSON}\n```"
        findings = parse_findings(fenced)
        assert len(findings) == 2

    def test_strips_fence_with_preamble(self):
        fenced = f"Here are the findings:\n```json\n{VALID_FINDING_JSON}\n```\nDone."
        findings = parse_findings(fenced)
        assert len(findings) == 2


class TestParseFindingsEmpty:
    """parse_findings with empty or minimal input."""

    def test_empty_string(self):
        assert parse_findings("") == []

    def test_empty_array(self):
        assert parse_findings("[]") == []

    def test_whitespace_only(self):
        assert parse_findings("   \n  ") == []


# ---------------------------------------------------------------------------
# collect_source_files
# ---------------------------------------------------------------------------

class TestCollectSourceFiles:
    """collect_source_files returns context containing known file markers."""

    def test_returns_non_empty(self):
        result = collect_source_files()
        assert len(result) > 0

    def test_contains_fit_models(self):
        result = collect_source_files()
        assert "fit_models" in result

    def test_contains_parse_config(self):
        result = collect_source_files()
        assert "parse_config" in result

    def test_contains_run_analysis(self):
        result = collect_source_files()
        assert "run_analysis" in result


# ---------------------------------------------------------------------------
# audit (dry_run)
# ---------------------------------------------------------------------------

class TestAuditDryRun:
    """audit() in dry_run mode returns an empty list without API calls."""

    def test_returns_empty_list(self):
        result = audit(iteration=1, context="test context", dry_run=True)
        assert result == []

    def test_returns_list_type(self):
        result = audit(iteration=1, context="", dry_run=True)
        assert isinstance(result, list)
