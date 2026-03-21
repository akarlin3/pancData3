# orchestrator_v1.py
"""Improvement-loop orchestrator.

Runs repeated audit → fix → evaluate cycles, stopping when the evaluator
signals that no further improvements are warranted.
"""

import argparse
import os
import sys

# Allow direct invocation (python improvement_loop/orchestrator_v1.py)
# by ensuring the repo root is on sys.path for absolute imports.
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from improvement_loop import loop_tracker
from improvement_loop import git_utils
from improvement_loop.evaluator import Finding
from improvement_loop.agents.auditor import audit as run_audit
from improvement_loop.agents.implementer import implement as _implement_finding, ImplementResult
from improvement_loop.agents.reviewer import review as _review_change
from typing import List

# Legacy constants — kept as module-level for backward compat but
# actual runtime values come from get_config().
MAX_API_RETRIES = 3
RETRY_BASE_DELAY = 30.0  # seconds — rate limit window is per-minute
MAX_SELF_HEAL_ATTEMPTS = 2  # max retries when a fix causes test failures

# Repo root is one level up from this file's directory
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _apply_fixes(findings: List[Finding], dry_run: bool) -> bool:
    """Apply fixes for each finding on its own branch. Returns True if all tests pass."""
    if dry_run:
        return True

    if not findings:
        return True

    original_branch = git_utils.current_branch()
    all_passed = True

    for finding in findings:
        print(f"\n--- Applying fix: {finding.branch_name} ---")
        print(f"    {finding.dimension}: {finding.description}")

        try:
            # Implement the fix (branch creation, API call, commit, syntax check)
            result = _implement_finding(finding, base_branch=original_branch, dry_run=dry_run)

            if not result.success:
                if "branch exists" in result.error:
                    print(f"    Branch {finding.branch_name} already exists, skipping")
                else:
                    print(f"    ❌ Implementation failed: {result.error}")
                finding.status = "pending"
                if result.error not in ("", "branch exists"):
                    all_passed = False
                continue

            # Run tests
            print("    Running Python tests...")
            py_ok = git_utils.run_python_tests()

            if py_ok:
                print("    ✅ Tests passed on branch")

                # Merge back into the original branch
                print(f"    ⚙️  Attempting merge: {finding.branch_name}")
                try:
                    git_utils.merge_branch(
                        finding.branch_name, target=original_branch,
                        delete_after=True,
                    )
                    finding.status = "merged"
                    print(f"    ✅  Merged: {finding.branch_name}")


                    # Post-merge sanity check
                    print(f"    ⚙️  cwd after merge: {os.getcwd()}")
                    print(f"    ⚙️  Post-merge test run on {original_branch}")
                    print(f"    ⚙️  Post-merge cwd: {os.getcwd()}")
                    if not git_utils.run_syntax_check():
                        print(f"    ❌  Post-merge syntax error — merge may have introduced issues")
                        all_passed = False
                        continue
                    post_ok = git_utils.run_python_tests()
                    print(f"    ⚙️  Post-merge test result: {post_ok}")
                    if post_ok:
                        print(f"    ✅  Post-merge tests passed")
                    else:
                        print(f"    ❌  Post-merge tests FAILED — merge may have introduced issues")
                        all_passed = False
                except Exception as e:
                    print(f"    ❌  Merge failed: {finding.branch_name} — {e}")
                    finding.status = "implemented"
                    all_passed = False
            else:
                finding.status = "pending"
                all_passed = False
                print("    ❌ Tests failed — fix needs review")

        except Exception as e:
            print(f"    ❌ Error applying fix: {e}")
            finding.status = "pending"
            all_passed = False
        finally:
            # Return to original branch for next finding
            try:
                git_utils.checkout(original_branch)
            except Exception:
                pass

    return all_passed


def run_loop(max_iterations: int = 10, dry_run: bool = False) -> list:
    """Execute the improvement loop up to *max_iterations* times.

    Stops early when ``loop_tracker.log_iteration`` returns an entry
    whose ``exit_condition_met`` field is ``True``.

    Returns the list of log entries produced during this run.
    """
    mode = "DRY RUN" if dry_run else "LIVE"
    print(f"\n{'='*60}")
    print(f"  Improvement Loop — {mode} (max {max_iterations} iterations)")
    print(f"{'='*60}\n")

    entries: list = []

    for i in range(1, max_iterations + 1):
        print(f"\n{'─'*60}")
        print(f"  ITERATION {i}/{max_iterations}")
        print(f"{'─'*60}")

        print(f"\n[1/4] Gathering context from prior iterations...")
        context = loop_tracker.get_context_for_next_iteration()

        print(f"[2/4] Running code audit via Claude API...")
        findings = run_audit(i, context, dry_run)
        print(f"       Found {len(findings)} valid finding(s)")
        for j, f in enumerate(findings, 1):
            print(f"       {j}. [{f.dimension}] {f.description[:80]}"
                  f" (importance={f.importance})")

        print(f"[3/4] Applying fixes and running tests...")
        tests_passed = _apply_fixes(findings, dry_run)

        print(f"\n[4/4] Logging iteration and evaluating exit condition...")

        # Build a synthetic audit_output for the tracker/evaluator
        # (they need the raw text for scoring; reconstruct from findings)
        import json
        audit_output = json.dumps(
            [f.to_log_dict() for f in findings], indent=2
        ) if findings else "[]"

        entry = loop_tracker.log_iteration(
            audit_output=audit_output,
            findings=findings,
            tests_passed=tests_passed,
            dry_run=dry_run,
        )
        entries.append(entry)

        if entry["exit_condition_met"]:
            print(f"\n  Loop complete after {i} iteration(s).")
            break

    # Print end-of-run summary
    _print_run_summary(entries)

    return entries


def _print_run_summary(entries: list) -> None:
    """Print a concise summary of this run's iterations and findings."""
    n = len(entries)
    print(f"\n{'='*60}")
    print(f"  IMPROVEMENT LOOP SUMMARY — {n} iteration(s)")
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

    # Per-iteration line
    for entry in entries:
        it = entry["iteration"]
        n_findings = entry["findings_count"]
        n_merged = len(entry.get("branches_merged", []))
        score = entry.get("audit_scores", {}).get("overall", "?")
        tests = "pass" if entry.get("tests_passed") else "FAIL"
        exit_flag = " [EXIT]" if entry.get("exit_condition_met") else ""
        print(f"  Iter {it}: {n_findings} findings, {n_merged} merged, "
              f"score={score}/10, tests={tests}{exit_flag}")

    # Totals
    print(f"\n  Findings:     {total_findings} total, "
          f"{total_implemented} implemented, {total_pending} pending")

    # By dimension
    if by_dimension:
        print(f"\n  By dimension:")
        for dim in sorted(by_dimension):
            counts = by_dimension[dim]
            print(f"    {dim}: {counts['implemented']} implemented, "
                  f"{counts['pending']} pending")

    # Final status
    if entries:
        last = entries[-1]
        if last.get("exit_condition_met"):
            print(f"\n  Status: Converged — all findings below threshold")
        elif last.get("tests_passed") is False:
            print(f"\n  Status: Stopped — test failures remain")
        else:
            print(f"\n  Status: Stopped — max iterations reached")

    print(f"{'='*60}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Improvement loop orchestrator")
    parser.add_argument("--dry-run", action="store_true",
                        help="Run without API calls or code changes")
    parser.add_argument("--max-iterations", type=int, default=10,
                        help="Maximum number of audit/fix cycles (default: 10)")
    parser.add_argument("--single-iteration", action="store_true",
                        help="Run exactly one iteration (shorthand for --max-iterations 1)")
    args = parser.parse_args()

    max_iter = 1 if args.single_iteration else args.max_iterations
    run_loop(max_iterations=max_iter, dry_run=args.dry_run)
