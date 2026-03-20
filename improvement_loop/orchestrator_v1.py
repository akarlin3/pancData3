# orchestrator_v1.py
"""Improvement-loop orchestrator.

Runs repeated audit → fix → evaluate cycles, stopping when the evaluator
signals that no further improvements are warranted.
"""

import argparse
import json
import os
import time

import anthropic  # type: ignore

from . import loop_tracker
from . import git_utils
from .evaluator import Finding
from typing import List

from .loop_config import get_config as _get_loop_config

# Legacy constants — kept as module-level for backward compat but
# actual runtime values come from get_config().
MAX_API_RETRIES = 3
RETRY_BASE_DELAY = 30.0  # seconds — rate limit window is per-minute
MAX_SELF_HEAL_ATTEMPTS = 2  # max retries when a fix causes test failures

# Repo root is one level up from this file's directory
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def _get_client() -> anthropic.Anthropic:
    """Return an Anthropic client, using the config API key if set."""
    cfg = _get_loop_config()
    kwargs = {}
    if cfg.anthropic_api_key:
        kwargs["api_key"] = cfg.anthropic_api_key
    return anthropic.Anthropic(**kwargs)

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
  "fix": "concise fix description (1-3 sentences, NO full code blocks — just describe the change)",
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


def _api_call_with_retry(create_kwargs: dict) -> str:
    """Call client.messages.create with rate-limit retry. Returns response text.

    Uses streaming to avoid the SDK's 10-minute timeout guard for large
    max_tokens values.
    """
    cfg = _get_loop_config()
    client = _get_client()
    for attempt in range(1, cfg.max_api_retries + 1):
        try:
            # Use streaming to avoid SDK timeout for large max_tokens
            collected: list[str] = []
            with client.messages.stream(**create_kwargs) as stream:
                for text in stream.text_stream:
                    collected.append(text)
            return "".join(collected)
        except anthropic.RateLimitError:
            delay = cfg.retry_base_delay * attempt
            print(f"    Rate limited (attempt {attempt}/{cfg.max_api_retries}), "
                  f"waiting {delay:.0f}s...")
            time.sleep(delay)
            if attempt == cfg.max_api_retries:
                raise
        except anthropic.APIError:
            raise
    return ""  # unreachable


def _collect_source_files() -> str:
    """Collect key source files as context for the audit."""
    cfg = _get_loop_config()
    # Smaller focused set to stay under rate limits
    key_files = [
        "pipeline/run_dwi_pipeline.m",
        "pipeline/core/fit_models.m",
        "pipeline/core/metrics_survival.m",
        "pipeline/core/metrics_stats_predictive.m",
        "pipeline/utils/parse_config.m",
        "pipeline/utils/knn_impute_train_test.m",
        "pipeline/utils/safe_load_mask.m",
        "pipeline/utils/escape_shell_arg.m",
        "analysis/run_analysis.py",
        "analysis/shared.py",
    ]

    parts = []
    for rel_path in key_files:
        full_path = os.path.join(REPO_ROOT, rel_path)
        if not os.path.exists(full_path):
            continue
        try:
            with open(full_path, "r", encoding="utf-8", errors="replace") as f:
                content = f.read()
            if len(content) > cfg.max_file_chars:
                content = content[:cfg.max_file_chars] + "\n... [truncated]"
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

    cfg = _get_loop_config()
    return _api_call_with_retry({
        "model": cfg.audit_model,
        "max_tokens": cfg.audit_max_tokens,
        "system": AUDIT_SYSTEM_PROMPT,
        "messages": [{"role": "user", "content": user_message}],
    })


def _parse_findings(audit_output: str, dry_run: bool) -> List[Finding]:
    """Parse or simulate findings from audit output."""
    if dry_run:
        return []

    # Strip markdown fences if present (handles preamble before the fence)
    text = audit_output.strip()
    if "```" in text:
        # Extract content between first pair of ``` fences
        parts = text.split("```")
        if len(parts) >= 3:
            # parts[1] is the content between the first pair of fences
            inner = parts[1]
            if inner.startswith("json"):
                inner = inner[4:]
            text = inner.strip()
        elif text.startswith("```"):
            text = parts[1]
            if text.startswith("json"):
                text = text[4:]
            text = text.strip()

    # If text still doesn't look like JSON, try to find a JSON array directly
    if not text.startswith("["):
        bracket_pos = text.find("[")
        if bracket_pos >= 0:
            text = text[bracket_pos:]



    # Truncation guard: check if the JSON array appears complete
    if not text.rstrip().endswith("]"):
        print("⚠️  Audit response appears truncated — consider increasing max_tokens further")
        # Attempt repair: find the last complete finding object by looking
        # for "},\n  {" or "}\n  {" boundaries (the split between findings).
        # Then try progressively shorter substrings.
        import re
        # Find all positions where a finding object boundary occurs
        boundaries = [m.end() - 1 for m in re.finditer(r'\}\s*,\s*\{', text)]
        # Also try the simple last-} approach as fallback
        last_brace = text.rfind("}")
        candidates = sorted(set(boundaries + ([last_brace] if last_brace >= 0 else [])),
                            reverse=True)
        for pos in candidates:
            # Take everything up to and including the } at this boundary
            # For boundary matches, pos points to the comma; use the } before it
            chunk = text[:pos + 1].rstrip().rstrip(",")
            if not chunk.lstrip().startswith("["):
                chunk = "[" + chunk
            chunk = chunk + "]"
            try:
                raw_list = json.loads(chunk)
                if not isinstance(raw_list, list):
                    continue
                print(f"⚠️  Recovered {len(raw_list)} findings from truncated response")
                findings = []
                for i, raw in enumerate(raw_list):
                    try:
                        finding = Finding(**raw)
                        findings.append(finding)
                    except Exception as e:
                        print(f"WARNING: Skipping finding {i}: {e}")
                return findings
            except json.JSONDecodeError:
                continue
        print("WARNING: Could not recover any findings from truncated response")
        print(f"Raw output (last 200 chars): {text[-200:]}")
        return []

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
            print("    Running syntax check...")
            if not git_utils.run_syntax_check():
                print("    ❌ Syntax error detected — skipping tests")
                finding.status = "pending"
                all_passed = False
                continue

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


def _apply_single_fix(finding: Finding) -> None:
    """Use Claude to generate and apply a code fix for a single finding."""
    file_path = os.path.join(REPO_ROOT, finding.file)
    if not os.path.exists(file_path):
        print(f"    WARNING: File {finding.file} not found, skipping")
        return

    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        original_content = f.read()

    cfg = _get_loop_config()
    new_content = _api_call_with_retry({
        "model": cfg.fix_model,
        "max_tokens": cfg.fix_max_tokens,
        "system": (
            "You are a code fixer. Given the original file and a description of the fix, "
            "return ONLY the complete updated file content. No markdown fences, no commentary, "
            "no explanation — just the raw file content ready to be written to disk."
        ),
        "messages": [{
            "role": "user",
            "content": (
                f"File: {finding.file}\n"
                f"Problem: {finding.description}\n"
                f"Fix: {finding.fix}\n\n"
                f"Original file content:\n{original_content}"
            ),
        }],
    })
    with open(file_path, "w", encoding="utf-8", newline="") as f:
        f.write(new_content)
    print(f"    Updated {finding.file}")


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

        print(f"\n[1/5] Gathering context from prior iterations...")
        context = loop_tracker.get_context_for_next_iteration()

        print(f"[2/5] Running code audit via Claude API...")
        audit_output = _run_audit(i, context, dry_run)
        print(f"       Audit response: {len(audit_output)} chars")

        print(f"[3/5] Parsing findings...")
        findings = _parse_findings(audit_output, dry_run)
        print(f"       Found {len(findings)} valid finding(s)")
        for j, f in enumerate(findings, 1):
            print(f"       {j}. [{f.dimension}] {f.description[:80]}"
                  f" (importance={f.importance})")

        print(f"[4/5] Applying fixes and running tests...")
        tests_passed = _apply_fixes(findings, dry_run)

        print(f"\n[5/5] Logging iteration and evaluating exit condition...")

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
