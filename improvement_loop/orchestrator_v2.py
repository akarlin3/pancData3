# orchestrator_v2.py
"""Improvement-loop orchestrator v2.

Rewires the improvement loop as a pipeline of discrete agent calls with
shared state: audit → implement → review → test & merge → log.
"""

import argparse
import json
import os
import subprocess
import sys
from dataclasses import dataclass, field
from typing import List, Optional

# Allow direct invocation (python improvement_loop/orchestrator_v2.py)
# by ensuring the repo root is on sys.path for absolute imports.
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from improvement_loop import loop_tracker
from improvement_loop import git_utils
from improvement_loop.evaluator import Finding
from improvement_loop.agents.auditor import audit as _audit
from improvement_loop.agents.implementer import implement as _implement, ImplementResult
from improvement_loop.agents.reviewer import review as _review, ReviewVerdict
from improvement_loop.loop_config import get_config as _get_loop_config

# Safety flags that force rejection regardless of reviewer verdict
_CRITICAL_FLAGS = frozenset({"LEAKAGE_RISK", "PHI_RISK"})


@dataclass
class FindingState:
    """Tracks one finding through the agent pipeline."""
    finding: Finding
    implement_result: Optional[ImplementResult] = None
    review_verdict: Optional[ReviewVerdict] = None
    merged: bool = False
    error: str = ""


@dataclass
class IterationState:
    """Shared state for one iteration of the loop."""
    iteration: int
    context: str
    findings: List[FindingState] = field(default_factory=list)
    all_tests_passed: bool = True


# ---------------------------------------------------------------------------
# Pipeline phases
# ---------------------------------------------------------------------------

def _phase_audit(state: IterationState, dry_run: bool) -> None:
    """Phase 1: Run the auditor agent and populate state.findings."""
    print(f"\n[1/5] Running code audit...")
    raw_findings = _audit(state.iteration, state.context, dry_run)
    state.findings = [FindingState(finding=f) for f in raw_findings]
    print(f"       Found {len(state.findings)} finding(s)")
    for j, fs in enumerate(state.findings, 1):
        f = fs.finding
        print(f"       {j}. [{f.dimension}] {f.description[:80]}"
              f" (importance={f.importance})")


def _phase_implement(state: IterationState, base_branch: str, dry_run: bool) -> None:
    """Phase 2: Implement fixes for each finding."""
    print(f"\n[2/5] Implementing fixes...")
    for fs in state.findings:
        print(f"\n    --- {fs.finding.branch_name} ---")
        print(f"    {fs.finding.dimension}: {fs.finding.description[:80]}")
        try:
            result = _implement(fs.finding, base_branch=base_branch, dry_run=dry_run)
            fs.implement_result = result
            if not result.success:
                fs.error = result.error
                fs.finding.status = "pending"
                if "branch exists" in result.error:
                    print(f"    Branch already exists, skipping")
                else:
                    print(f"    ❌ Implementation failed: {result.error}")
                    state.all_tests_passed = False
            else:
                print(f"    ✅ Implemented")
        except Exception as e:
            fs.error = str(e)
            fs.finding.status = "pending"
            state.all_tests_passed = False
            print(f"    ❌ Error: {e}")
        finally:
            try:
                git_utils.checkout(base_branch)
            except Exception:
                pass


def _phase_review(state: IterationState, base_branch: str, dry_run: bool) -> None:
    """Phase 3: Review implemented changes."""
    print(f"\n[3/5] Reviewing changes...")
    for fs in state.findings:
        if fs.implement_result is None or not fs.implement_result.success:
            continue

        print(f"\n    --- Reviewing {fs.finding.branch_name} ---")
        try:
            verdict = _review(
                fs.finding,
                original_content=fs.implement_result.original_content,
                new_content=fs.implement_result.new_content,
                dry_run=dry_run,
            )
            fs.review_verdict = verdict
            print(f"    Verdict: {verdict.verdict}")
            if verdict.reasoning:
                print(f"    Reasoning: {verdict.reasoning}")

            # Check for critical safety flags
            has_critical = bool(_CRITICAL_FLAGS & set(verdict.risk_flags))

            if verdict.verdict == "reject" or has_critical:
                reason = verdict.reasoning
                if has_critical and verdict.verdict != "reject":
                    reason = (
                        f"Auto-rejected due to critical flags "
                        f"{sorted(_CRITICAL_FLAGS & set(verdict.risk_flags))}. "
                        f"Original: {verdict.reasoning}"
                    )
                fs.error = f"Rejected: {reason}"
                fs.finding.status = "pending"
                state.all_tests_passed = False
                if verdict.risk_flags:
                    print(f"    Risk flags: {verdict.risk_flags}")
                print(f"    ❌ Rejected: {reason}")
                # Delete the rejected branch
                try:
                    git_utils.checkout(base_branch)
                    subprocess.run(
                        ["git", "branch", "-D", fs.finding.branch_name],
                        check=False, capture_output=True, text=True,
                    )
                    print(f"    Deleted branch: {fs.finding.branch_name}")
                except Exception:
                    pass
            elif verdict.verdict == "request_changes":
                fs.error = f"Changes requested: {verdict.reasoning}"
                fs.finding.status = "pending"
                state.all_tests_passed = False
                print(f"    ⚠️  Changes requested (branch preserved)")
            else:
                print(f"    ✅ Approved")

        except Exception as e:
            fs.error = f"Review error: {e}"
            fs.finding.status = "pending"
            state.all_tests_passed = False
            print(f"    ❌ Review error: {e}")


def _phase_test_and_merge(state: IterationState, base_branch: str, dry_run: bool) -> None:
    """Phase 4: Test and merge approved changes."""
    print(f"\n[4/5] Testing and merging...")
    if dry_run:
        return

    cfg = _get_loop_config()

    for fs in state.findings:
        # Only proceed for approved findings
        if fs.review_verdict is None or fs.review_verdict.verdict != "approve":
            continue
        if fs.error:
            continue

        print(f"\n    --- Testing {fs.finding.branch_name} ---")
        try:
            git_utils.checkout(fs.finding.branch_name)

            # Run tests with self-heal retries
            tests_ok = False
            for attempt in range(1, cfg.max_self_heal_attempts + 1):
                print(f"    Running Python tests (attempt {attempt})...")
                py_ok = git_utils.run_python_tests()
                if py_ok:
                    tests_ok = True
                    break
                if attempt < cfg.max_self_heal_attempts:
                    print(f"    ⚠️  Tests failed, retrying ({attempt}/{cfg.max_self_heal_attempts})...")

            if not tests_ok:
                fs.finding.status = "pending"
                fs.error = "tests failed"
                state.all_tests_passed = False
                print(f"    ❌ Tests failed after {cfg.max_self_heal_attempts} attempt(s)")
                git_utils.checkout(base_branch)
                continue

            print(f"    ✅ Tests passed")

            # Merge
            print(f"    ⚙️  Merging {fs.finding.branch_name}...")
            try:
                git_utils.merge_branch(
                    fs.finding.branch_name, target=base_branch,
                    delete_after=True,
                )
                fs.merged = True
                fs.finding.status = "merged"
                print(f"    ✅ Merged")

                # Post-merge sanity check
                if not git_utils.run_syntax_check():
                    print(f"    ❌ Post-merge syntax error")
                    state.all_tests_passed = False
                    continue
                post_ok = git_utils.run_python_tests()
                if post_ok:
                    print(f"    ✅ Post-merge tests passed")
                else:
                    print(f"    ❌ Post-merge tests FAILED")
                    state.all_tests_passed = False

            except Exception as e:
                fs.finding.status = "implemented"
                fs.error = f"merge failed: {e}"
                state.all_tests_passed = False
                print(f"    ❌ Merge failed: {e}")

        except Exception as e:
            fs.finding.status = "pending"
            fs.error = str(e)
            state.all_tests_passed = False
            print(f"    ❌ Error: {e}")
        finally:
            try:
                git_utils.checkout(base_branch)
            except Exception:
                pass


def _phase_log(state: IterationState, dry_run: bool) -> dict:
    """Phase 5: Log the iteration and check exit condition."""
    print(f"\n[5/5] Logging iteration...")

    findings = [fs.finding for fs in state.findings]

    # Build audit_output for the evaluator/tracker
    audit_output = json.dumps(
        [f.to_log_dict() for f in findings], indent=2
    ) if findings else "[]"

    entry = loop_tracker.log_iteration(
        audit_output=audit_output,
        findings=findings,
        tests_passed=state.all_tests_passed,
        dry_run=dry_run,
    )

    _print_agent_summary(state)
    return entry


# ---------------------------------------------------------------------------
# Summary printing
# ---------------------------------------------------------------------------

def _print_agent_summary(state: IterationState) -> None:
    """Print per-agent summary for one iteration."""
    n_findings = len(state.findings)
    n_high = sum(1 for fs in state.findings if fs.finding.importance >= 7)

    n_impl_ok = sum(1 for fs in state.findings
                    if fs.implement_result and fs.implement_result.success)
    n_impl_fail = sum(1 for fs in state.findings
                      if fs.implement_result and not fs.implement_result.success)
    # Also count findings where implement was never attempted (error before)
    n_impl_fail += sum(1 for fs in state.findings if fs.implement_result is None)

    n_approved = sum(1 for fs in state.findings
                     if fs.review_verdict and fs.review_verdict.verdict == "approve")
    n_changes = sum(1 for fs in state.findings
                    if fs.review_verdict and fs.review_verdict.verdict == "request_changes")
    n_rejected = sum(1 for fs in state.findings
                     if fs.review_verdict and fs.review_verdict.verdict == "reject")
    # Count critical-flag overrides as rejections too
    n_rejected += sum(1 for fs in state.findings
                      if fs.review_verdict
                      and fs.review_verdict.verdict != "reject"
                      and bool(_CRITICAL_FLAGS & set(fs.review_verdict.risk_flags)))

    n_merged = sum(1 for fs in state.findings if fs.merged)
    n_merge_fail = n_approved - n_merged

    all_flags: list[str] = []
    for fs in state.findings:
        if fs.review_verdict and fs.review_verdict.risk_flags:
            all_flags.extend(fs.review_verdict.risk_flags)

    print(f"\n  ── Iteration {state.iteration} Agent Summary ──")
    print(f"  Audit:      {n_findings} findings ({n_high} high-priority)")
    print(f"  Implement:  {n_impl_ok} succeeded, {n_impl_fail} failed")
    print(f"  Review:     {n_approved} approved, {n_changes} changes_requested, {n_rejected} rejected")
    print(f"  Merge:      {n_merged} merged, {n_merge_fail} failed")
    if all_flags:
        print(f"  Flags:      {sorted(set(all_flags))}")


def _print_run_summary(entries: list) -> None:
    """Print a concise summary of the full run."""
    n = len(entries)
    print(f"\n{'='*60}")
    print(f"  IMPROVEMENT LOOP v2 SUMMARY — {n} iteration(s)")
    print(f"{'='*60}")

    total_findings = 0
    total_implemented = 0
    total_pending = 0
    by_dimension: dict = {}

    for entry in entries:
        for f in entry.get("findings", []):
            total_findings += 1
            status = f.get("status", "unknown")
            dim = f.get("dimension", "unknown")
            by_dimension.setdefault(dim, {"implemented": 0, "pending": 0})
            if status in ("implemented", "merged"):
                total_implemented += 1
                by_dimension[dim]["implemented"] += 1
            else:
                total_pending += 1
                by_dimension[dim]["pending"] += 1

    for entry in entries:
        it = entry["iteration"]
        n_findings = entry["findings_count"]
        n_merged = len(entry.get("branches_merged", []))
        score = entry.get("audit_scores", {}).get("overall", "?")
        tests = "pass" if entry.get("tests_passed") else "FAIL"
        exit_flag = " [EXIT]" if entry.get("exit_condition_met") else ""
        print(f"  Iter {it}: {n_findings} findings, {n_merged} merged, "
              f"score={score}/10, tests={tests}{exit_flag}")

    print(f"\n  Findings:     {total_findings} total, "
          f"{total_implemented} implemented, {total_pending} pending")

    if by_dimension:
        print(f"\n  By dimension:")
        for dim in sorted(by_dimension):
            counts = by_dimension[dim]
            print(f"    {dim}: {counts['implemented']} implemented, "
                  f"{counts['pending']} pending")

    if entries:
        last = entries[-1]
        if last.get("exit_condition_met"):
            print(f"\n  Status: Converged — all findings below threshold")
        elif last.get("tests_passed") is False:
            print(f"\n  Status: Stopped — test failures remain")
        else:
            print(f"\n  Status: Stopped — max iterations reached")

    print(f"{'='*60}\n")


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

def run_loop(max_iterations: int = 10, dry_run: bool = False) -> list:
    """Execute the improvement loop as a four-agent pipeline.

    Pipeline per iteration: audit → implement → review → test & merge → log.

    Stops early when ``loop_tracker.log_iteration`` returns an entry
    whose ``exit_condition_met`` field is ``True``.

    Returns the list of log entries produced during this run.
    """
    cfg = _get_loop_config()
    mode = "DRY RUN" if dry_run else "LIVE"
    print(f"\n{'='*60}")
    print(f"  Improvement Loop v2 — {mode} (max {max_iterations} iterations)")
    print(f"{'='*60}\n")

    base_branch = git_utils.current_branch() if not dry_run else "main"
    entries: list = []

    # Build / update the RAG index before the first iteration
    if cfg.rag_enabled and not dry_run:
        try:
            from improvement_loop.rag.indexer import (
                build_index, index_improvement_history, REPO_ROOT as _IDX_ROOT,
            )
            print("  Building/updating codebase index...")
            build_index(_IDX_ROOT)
            index_improvement_history()
        except Exception as e:
            print(f"⚠️  RAG index build failed (continuing without): {e}")

    for i in range(1, max_iterations + 1):
        print(f"\n{'─'*60}")
        print(f"  ITERATION {i}/{max_iterations}")
        print(f"{'─'*60}")

        context = loop_tracker.get_context_for_next_iteration()
        state = IterationState(iteration=i, context=context)

        # Agent pipeline
        _phase_audit(state, dry_run)
        _phase_implement(state, base_branch, dry_run)
        _phase_review(state, base_branch, dry_run)
        _phase_test_and_merge(state, base_branch, dry_run)

        # Update RAG index for merged files before logging
        if cfg.rag_enabled and not dry_run:
            merged_files = [
                fs.finding.file for fs in state.findings if fs.merged
            ]
            if merged_files:
                try:
                    from improvement_loop.rag.indexer import update_index_for_files
                    update_index_for_files(merged_files)
                except Exception as e:
                    print(f"⚠️  RAG index update failed (non-fatal): {e}")

        entry = _phase_log(state, dry_run)

        entries.append(entry)

        if entry["exit_condition_met"]:
            print(f"\n  Loop complete after {i} iteration(s).")
            break

    _print_run_summary(entries)
    return entries


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Improvement loop orchestrator v2 (agent pipeline)"
    )
    parser.add_argument("--dry-run", action="store_true",
                        help="Run without API calls or code changes")
    parser.add_argument("--max-iterations", type=int, default=10,
                        help="Maximum number of audit/fix cycles (default: 10)")
    parser.add_argument("--single-iteration", action="store_true",
                        help="Run exactly one iteration (shorthand for --max-iterations 1)")
    args = parser.parse_args()

    max_iter = 1 if args.single_iteration else args.max_iterations
    run_loop(max_iterations=max_iter, dry_run=args.dry_run)
