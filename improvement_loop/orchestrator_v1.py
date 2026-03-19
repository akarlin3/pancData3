# orchestrator_v1.py
"""Improvement-loop orchestrator.

Runs repeated audit → fix → evaluate cycles, stopping when the evaluator
signals that no further improvements are warranted.
"""

from . import loop_tracker
from .evaluator import Finding
from typing import List


def _run_audit(iteration: int, context: str, dry_run: bool) -> str:
    """Run or simulate a code audit.  Returns raw audit text."""
    if dry_run:
        return f"[dry-run] audit output for iteration {iteration}"
    # In production this would call Claude to audit the codebase
    raise NotImplementedError("Live audit not yet wired up")


def _parse_findings(audit_output: str, dry_run: bool) -> List[Finding]:
    """Parse or simulate findings from audit output."""
    if dry_run:
        return []
    raise NotImplementedError("Live finding parser not yet wired up")


def _apply_fixes(findings: List[Finding], dry_run: bool) -> bool:
    """Apply fixes and run tests.  Returns True if tests pass."""
    if dry_run:
        return True
    raise NotImplementedError("Live fix application not yet wired up")


def run_loop(max_iterations: int = 10, dry_run: bool = False) -> list:
    """Execute the improvement loop up to *max_iterations* times.

    Stops early when ``loop_tracker.log_iteration`` returns an entry
    whose ``exit_condition_met`` field is ``True``.

    Returns the list of log entries produced during this run.
    """
    entries: list = []

    for i in range(1, max_iterations + 1):
        context = loop_tracker.get_context_for_next_iteration()

        # Phase 1-3: audit
        audit_output = _run_audit(i, context, dry_run)

        # Phase 4: parse findings & apply fixes
        findings = _parse_findings(audit_output, dry_run)
        tests_passed = _apply_fixes(findings, dry_run)

        # Phase 5: log iteration and evaluate exit condition
        entry = loop_tracker.log_iteration(
            audit_output=audit_output,
            findings=findings,
            tests_passed=tests_passed,
            dry_run=dry_run,
        )
        entries.append(entry)

        if entry["exit_condition_met"]:
            break

    return entries


if __name__ == "__main__":
    entries = run_loop(dry_run=True)
    for e in entries:
        print(e)
