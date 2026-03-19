#!/usr/bin/env python3
"""
orchestrator_v1.py — Programmatic driver for the pancData3 improvement loop.

Replaces the natural-language loop instructions in CLAUDE_WORKFLOWS.md.
Claude is called as a subprocess component for audit and implementation;
this script owns loop control, git operations, logging, and exit decisions.
"""
import argparse
import json
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import List, Optional

import anthropic  # type: ignore

from evaluator import Finding
import loop_tracker
import git_utils

# ── Constants ────────────────────────────────────────────────────────────────

MODEL = "claude-sonnet-4-20250514"
REPO_ROOT = Path(__file__).resolve().parent
CLAUDE_MD_PATH = REPO_ROOT / "CLAUDE.md"
DEPENDENCIES_PREFIX = "pipeline/dependencies/"

# ── Helpers ──────────────────────────────────────────────────────────────────


def _load_claude_md() -> str:
    """Read CLAUDE.md from disk for inclusion in audit prompts."""
    return CLAUDE_MD_PATH.read_text(encoding="utf-8")


def _check_api_key() -> None:
    """Exit immediately if ANTHROPIC_API_KEY is not set."""
    if not os.environ.get("ANTHROPIC_API_KEY"):
        print("❌ ANTHROPIC_API_KEY is not set. Export it before running the orchestrator.")
        sys.exit(1)


def _guard_dependencies(findings: List[Finding]) -> None:
    """Raise RuntimeError if any finding targets pipeline/dependencies/."""
    for f in findings:
        if f.file.startswith(DEPENDENCIES_PREFIX):
            raise RuntimeError(
                f"Finding targets protected path '{f.file}' — "
                f"pipeline/dependencies/ must never be modified."
            )


def _read_file_content(file_path: str) -> Optional[str]:
    """Safely read a file relative to the repo root. Returns None on failure."""
    try:
        full = REPO_ROOT / file_path
        return full.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return None


# ── Phase implementations ────────────────────────────────────────────────────


def phase0_load_context() -> str:
    """Phase 0 — Load context from prior iterations."""
    print("⚙️  Phase 0: Loading iteration context …")
    ctx = loop_tracker.get_context_for_next_iteration()
    print(f"   Context loaded ({len(ctx)} chars)")
    return ctx


def phase1_audit(client: anthropic.Anthropic, claude_md: str, history_context: str) -> List[Finding]:
    """
    Phase 1 — Audit: call Claude with codebase context and history.
    Returns a validated list of Findings.
    """
    print("⚙️  Phase 1: Running codebase audit …")

    system_prompt = (
        "You are an expert code auditor for a MATLAB + Python medical-imaging "
        "research pipeline (pancreatic DWI analysis). Below is the full project "
        "guide (CLAUDE.md) followed by iteration history.\n\n"
        "## CLAUDE.md\n\n"
        f"{claude_md}\n\n"
        "## Iteration history\n\n"
        f"{history_context}\n\n"
        "## Instructions\n\n"
        "Audit the codebase across all 9 dimensions: performance, correctness, "
        "error_handling, modularity, memory, code_quality, test_coverage, "
        "security, cross_platform.\n\n"
        "Return ONLY a JSON array of findings. Each finding must be an object with:\n"
        '  "dimension": one of the 9 dimensions above\n'
        '  "file": relative path from repo root\n'
        '  "function_name": (optional) function name\n'
        '  "description": what the issue is\n'
        '  "fix": how to fix it\n'
        '  "importance": integer 1-10\n'
        '  "branch_name": "improvement/<short-slug>" (no spaces, max 50 char slug)\n\n'
        "Do NOT suggest changes to pipeline/dependencies/.\n"
        "Do NOT repeat already-implemented improvements from the history.\n"
        "Return valid JSON only — no markdown fences, no commentary."
    )

    response = client.messages.create(
        model=MODEL,
        max_tokens=4096,
        system=system_prompt,
        messages=[{"role": "user", "content": "Audit the pancData3 repository now."}],
    )

    raw = response.content[0].text.strip()
    # Strip markdown fences if present
    if raw.startswith("```"):
        raw = raw.split("```")[1]
        if raw.startswith("json"):
            raw = raw[4:]
        raw = raw.strip()

    raw_findings = json.loads(raw)
    findings: List[Finding] = []
    for item in raw_findings:
        item.setdefault("status", "pending")
        findings.append(Finding(**item))

    _guard_dependencies(findings)
    print(f"   ✅ Audit returned {len(findings)} findings")
    return findings


def phase2_plan(findings: List[Finding]) -> List[List[Finding]]:
    """
    Phase 2 — Plan: group findings into work units.
    Findings with no overlapping files run in the same parallel group;
    findings touching the same file are serialized into separate groups.
    """
    print("⚙️  Phase 2: Planning work units …")

    groups: List[List[Finding]] = []
    assigned_files: List[set] = []

    for finding in sorted(findings, key=lambda f: -f.importance):
        placed = False
        for i, group in enumerate(groups):
            if finding.file not in assigned_files[i]:
                group.append(finding)
                assigned_files[i].add(finding.file)
                placed = True
                break
        if not placed:
            groups.append([finding])
            assigned_files.append({finding.file})

    for i, group in enumerate(groups):
        descs = ", ".join(f.branch_name for f in group)
        print(f"   Group {i + 1}: {descs}")

    print(f"   ✅ {len(findings)} findings in {len(groups)} groups")
    return groups


def phase3_implement(
    client: anthropic.Anthropic,
    findings: List[Finding],
    dry_run: bool,
    base_branch: str,
) -> None:
    """
    Phase 3 — Implement each finding.
    In dry_run mode, skip git and API implementation calls.
    """
    print("⚙️  Phase 3: Implementing findings …")

    for finding in findings:
        print(f"   ⚙️  {finding.branch_name}: {finding.description[:72]}")

        if dry_run:
            print(f"      (dry run — skipped)")
            finding.status = "pending"
            continue

        try:
            # Switch back to base branch before creating a new one
            git_utils.switch_branch(base_branch)
            if git_utils.branch_exists(finding.branch_name):
                print(f"      ⚠️  Branch already exists, checking out existing branch: {finding.branch_name}")
                git_utils.checkout(finding.branch_name)
            else:
                git_utils.create_branch(finding.branch_name)

            # Read the relevant file
            file_content = _read_file_content(finding.file) or "(file not found)"

            # Ask Claude for the implementation
            impl_prompt = (
                f"Implement the following improvement.\n\n"
                f"File: {finding.file}\n"
                f"Function: {finding.function_name or '(module-level)'}\n"
                f"Description: {finding.description}\n"
                f"Fix: {finding.fix}\n\n"
                f"Current file content:\n```\n{file_content}\n```\n\n"
                f"Return ONLY the complete updated file content. "
                f"No markdown fences, no commentary."
            )

            response = client.messages.create(
                model=MODEL,
                max_tokens=8192,
                messages=[{"role": "user", "content": impl_prompt}],
            )

            new_content = response.content[0].text.strip()
            # Strip markdown fences if present
            if new_content.startswith("```"):
                lines = new_content.split("\n")
                # Remove first and last fence lines
                if lines[0].startswith("```"):
                    lines = lines[1:]
                if lines and lines[-1].strip() == "```":
                    lines = lines[:-1]
                new_content = "\n".join(lines)

            # Write the changes
            target = REPO_ROOT / finding.file
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(new_content, encoding="utf-8")

            # Commit
            git_utils.commit_all(f"[{finding.dimension}] {finding.description[:72]}")
            finding.status = "implemented"
            print(f"      ✅ Implemented and committed")

        except Exception as e:
            print(f"      ⚠️  Implementation failed: {e}")
            finding.status = "pending"
            # Try to get back to base branch
            try:
                git_utils.switch_branch(base_branch)
            except Exception:
                pass


def phase4_test_and_merge(findings: List[Finding], dry_run: bool, base_branch: str) -> bool:
    """
    Phase 4 — Test each implemented branch and merge if tests pass.
    Returns True if all tests passed across all branches.
    """
    print("⚙️  Phase 4: Testing and merging …")

    if dry_run:
        print("   (dry run — skipped)")
        return True

    all_passed = True
    for finding in findings:
        if finding.status != "implemented":
            continue

        print(f"   ⚙️  Testing {finding.branch_name} …")

        try:
            git_utils.switch_branch(finding.branch_name)
            passed = git_utils.run_python_tests()

            if passed:
                git_utils.switch_branch(base_branch)
                git_utils.merge_branch(finding.branch_name)
                finding.status = "merged"
                print(f"      ✅ Tests passed — merged")
            else:
                finding.status = "failed"
                all_passed = False
                print(f"      ⚠️  Tests failed — branch left unmerged")
                git_utils.switch_branch(base_branch)

        except Exception as e:
            print(f"      ⚠️  Test/merge failed: {e}")
            finding.status = "failed"
            all_passed = False
            try:
                git_utils.switch_branch(base_branch)
            except Exception:
                pass

    return all_passed


def phase5_log(findings: List[Finding], audit_text: str, tests_passed: bool) -> dict:
    """Phase 5 — Log the iteration via loop_tracker."""
    print("⚙️  Phase 5: Logging iteration …")
    entry = loop_tracker.log_iteration(
        audit_output=audit_text,
        findings=findings,
        tests_passed=tests_passed,
    )
    print(f"   ✅ Iteration {entry['iteration']} logged")
    return entry


def phase6_check_exit(entry: dict) -> bool:
    """Phase 6 — Check if exit condition is met. Returns True to exit."""
    if entry["exit_condition_met"]:
        print("✅ Exit condition met — loop complete")
        return True
    print("💡 Exit condition not met — continuing to next iteration")
    return False


# ── Main loop ────────────────────────────────────────────────────────────────


def run_loop(
    max_iterations: int = 20,
    dry_run: bool = False,
    single_iteration: bool = False,
) -> None:
    """
    Run the full improvement loop.

    Args:
        max_iterations: Hard cap on iterations (safety guard).
        dry_run: If True, skip git operations and API implementation calls.
        single_iteration: If True, run exactly one iteration then stop.
    """
    _check_api_key()
    claude_md = _load_claude_md()
    client = anthropic.Anthropic()

    base_branch = git_utils.current_branch() if not dry_run else "main"

    print(f"🚀 Starting improvement loop (max={max_iterations}, "
          f"dry_run={dry_run}, single={single_iteration})")

    for iteration in range(1, max_iterations + 1):
        print(f"\n{'='*60}")
        print(f"🚀 Iteration {iteration}")
        print(f"{'='*60}")

        try:
            # Phase 0
            history_context = phase0_load_context()

            # Phase 1
            findings = phase1_audit(client, claude_md, history_context)

            # Build audit text for logging
            audit_text = json.dumps([f.to_log_dict() for f in findings], indent=2)

            # Phase 2
            groups = phase2_plan(findings)

            # Phase 3
            for group in groups:
                phase3_implement(client, group, dry_run, base_branch)

            # Phase 4
            tests_passed = phase4_test_and_merge(findings, dry_run, base_branch)

            # Phase 5
            entry = phase5_log(findings, audit_text, tests_passed)

            # Phase 6
            should_exit = phase6_check_exit(entry)
            if should_exit or single_iteration:
                break

        except Exception as e:
            print(f"❌ Iteration {iteration} failed: {e}")
            print("⚠️  Continuing to next iteration …")
            continue

    print(f"\n{'='*60}")
    print("🚀 Loop complete — printing summary")
    print(f"{'='*60}")
    loop_tracker.print_full_summary()


# ── CLI ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Programmatic driver for the pancData3 improvement loop."
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=20,
        help="Hard cap on loop iterations (default: 20)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Skip git operations and API implementation calls; only audit and log",
    )
    parser.add_argument(
        "--single-iteration",
        action="store_true",
        help="Run exactly one full iteration regardless of exit condition",
    )
    args = parser.parse_args()

    run_loop(
        max_iterations=args.max_iterations,
        dry_run=args.dry_run,
        single_iteration=args.single_iteration,
    )


if __name__ == "__main__":
    main()
