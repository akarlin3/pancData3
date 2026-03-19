# orchestrator_v1.py
"""Improvement-loop orchestrator.

Runs repeated audit → fix → evaluate cycles, stopping when the evaluator
signals that no further improvements are warranted.
"""

import argparse
import json
import os
import anthropic  # type: ignore

from . import loop_tracker
from . import git_utils
from .evaluator import Finding
from typing import List

# Repo root is one level up from this file's directory
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

client = anthropic.Anthropic()

AUDIT_SYSTEM_PROMPT = """\
You are a senior code auditor for a MATLAB + Python medical-imaging research pipeline \
(pancreatic DWI analysis at Memorial Sloan Kettering).

## Repository overview
- MATLAB pipeline in `pipeline/` — IVIM and ADC diffusion model fitting, survival analysis, \
competing risks, treatment response prediction.
- Python analysis in `analysis/` — post-hoc parsing, cross-DWI comparison, HTML/PDF reports.
- Config-driven via `config.json`; backwards-compatible defaults in `pipeline/utils/parse_config.m`.
- `pipeline/dependencies/` is READ-ONLY — never suggest modifications there.
- Patient data must NEVER be exposed; no PHI in code, logs, or suggestions.

## Your task
Audit the codebase and return a JSON array of findings. Each finding must have:
{
  "dimension": one of ("performance", "correctness", "error_handling", "modularity", \
"memory", "code_quality", "test_coverage", "security", "cross_platform"),
  "file": "path/to/file.m or .py",
  "function_name": "optional — specific function if applicable",
  "description": "what the problem is — cite specific lines",
  "fix": "concrete fix description with code snippets if helpful",
  "importance": integer 1-10 (8-10 = critical bugs/data integrity, 4-7 = moderate, 1-3 = minor),
  "branch_name": "improvement/<slug>" (slug: lowercase, hyphens, no spaces, max 50 chars after prefix)
}

## Rules
- Be specific: cite files, functions, and line numbers.
- Do NOT suggest changes to `pipeline/dependencies/`.
- Do NOT suggest changes that would introduce data leakage (removing patient-stratified CV, \
cross-timepoint imputation, etc.).
- Do NOT inflate importance scores — style nits are 1-2, not 7-8.
- Cover multiple dimensions if possible (performance, correctness, security, tests, etc.).
- Return ONLY the JSON array — no markdown fences, no commentary.
"""


def _collect_source_files() -> str:
    """Collect key source files as context for the audit."""
    key_files = [
        "pipeline/run_dwi_pipeline.m",
        "pipeline/execute_all_workflows.m",
        "pipeline/core/load_dwi_data.m",
        "pipeline/core/sanity_checks.m",
        "pipeline/core/compute_summary_metrics.m",
        "pipeline/core/metrics_baseline.m",
        "pipeline/core/metrics_survival.m",
        "pipeline/core/metrics_stats_predictive.m",
        "pipeline/core/fit_models.m",
        "pipeline/utils/parse_config.m",
        "pipeline/utils/knn_impute_train_test.m",
        "pipeline/utils/extract_tumor_core.m",
        "pipeline/utils/safe_load_mask.m",
        "pipeline/utils/escape_shell_arg.m",
        "pipeline/utils/build_td_panel.m",
        "analysis/run_analysis.py",
        "analysis/shared.py",
        "analysis/parsers/batch_graph_analysis.py",
        "analysis/report/generate_report.py",
    ]

    parts = []
    for rel_path in key_files:
        full_path = os.path.join(REPO_ROOT, rel_path)
        if not os.path.exists(full_path):
            continue
        try:
            with open(full_path, "r", encoding="utf-8", errors="replace") as f:
                content = f.read()
            # Truncate very large files to avoid token limits
            if len(content) > 15000:
                content = content[:15000] + "\n... [truncated]"
            parts.append(f"=== {rel_path} ===\n{content}")
        except OSError:
            continue
    return "\n\n".join(parts)


def _run_audit(iteration: int, context: str, dry_run: bool) -> str:
    """Run or simulate a code audit.  Returns raw audit text."""
    if dry_run:
        return f"[dry-run] audit output for iteration {iteration}"

    source_context = _collect_source_files()
    user_message = (
        f"## Iteration context\n{context}\n\n"
        f"## Source files\n{source_context}\n\n"
        "Return your findings as a JSON array."
    )

    response = client.messages.create(
        model="claude-sonnet-4-20250514",
        max_tokens=4096,
        system=AUDIT_SYSTEM_PROMPT,
        messages=[{"role": "user", "content": user_message}],
    )
    return response.content[0].text


def _parse_findings(audit_output: str, dry_run: bool) -> List[Finding]:
    """Parse or simulate findings from audit output."""
    if dry_run:
        return []

    # Strip markdown fences if present
    text = audit_output.strip()
    if text.startswith("```"):
        text = text.split("```")[1]
        if text.startswith("json"):
            text = text[4:]
        text = text.strip()

    try:
        raw_list = json.loads(text)
    except json.JSONDecodeError as e:
        print(f"WARNING: Could not parse audit findings as JSON: {e}")
        print(f"Raw output (first 500 chars): {audit_output[:500]}")
        return []

    if not isinstance(raw_list, list):
        print(f"WARNING: Expected JSON array, got {type(raw_list).__name__}")
        return []

    findings = []
    for i, raw in enumerate(raw_list):
        try:
            finding = Finding(**raw)
            findings.append(finding)
        except Exception as e:
            print(f"WARNING: Skipping finding {i}: {e}")
    return findings


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
            # Create branch for this fix
            if git_utils.branch_exists(finding.branch_name):
                print(f"    Branch {finding.branch_name} already exists, skipping")
                finding.status = "pending"
                continue

            git_utils.create_branch(finding.branch_name, base=original_branch)

            # Use Claude to generate and apply the fix
            _apply_single_fix(finding)

            # Commit changes
            git_utils.commit_all(
                f"improvement: {finding.dimension} — {finding.description[:60]}"
            )

            # Run tests
            print("    Running Python tests...")
            py_ok = git_utils.run_python_tests()

            if py_ok:
                finding.status = "implemented"
                print("    ✅ Tests passed")
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


def _apply_single_fix(finding: Finding) -> None:
    """Use Claude to generate and apply a code fix for a single finding."""
    file_path = os.path.join(REPO_ROOT, finding.file)
    if not os.path.exists(file_path):
        print(f"    WARNING: File {finding.file} not found, skipping")
        return

    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        original_content = f.read()

    response = client.messages.create(
        model="claude-sonnet-4-20250514",
        max_tokens=8192,
        system=(
            "You are a code fixer. Given the original file and a description of the fix, "
            "return ONLY the complete updated file content. No markdown fences, no commentary, "
            "no explanation — just the raw file content ready to be written to disk."
        ),
        messages=[{
            "role": "user",
            "content": (
                f"File: {finding.file}\n"
                f"Problem: {finding.description}\n"
                f"Fix: {finding.fix}\n\n"
                f"Original file content:\n{original_content}"
            ),
        }],
    )

    new_content = response.content[0].text
    with open(file_path, "w", encoding="utf-8", newline="") as f:
        f.write(new_content)
    print(f"    Updated {finding.file}")


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
    parser = argparse.ArgumentParser(description="Improvement loop orchestrator")
    parser.add_argument("--dry-run", action="store_true",
                        help="Run without API calls or code changes")
    parser.add_argument("--max-iterations", type=int, default=10,
                        help="Maximum number of audit/fix cycles (default: 10)")
    args = parser.parse_args()

    entries = run_loop(max_iterations=args.max_iterations, dry_run=args.dry_run)
    for e in entries:
        print(e)
