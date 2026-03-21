# implementer.py
"""Implementer agent — generates and applies code fixes for audit findings."""

import os
from dataclasses import dataclass, field
from typing import List

from improvement_loop import git_utils
from improvement_loop.agents._api import api_call_with_retry
from improvement_loop.evaluator import Finding
from improvement_loop.loop_config import get_config as _get_loop_config

# Repo root is two levels up from this file's directory
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

FIX_SYSTEM_PROMPT = (
    "You are a code fixer. Given the original file and a description of the fix, "
    "return ONLY the complete updated file content. No markdown fences, no commentary, "
    "no explanation — just the raw file content ready to be written to disk."
)


@dataclass
class ImplementResult:
    """Result of applying a code fix for a single finding."""
    success: bool
    original_content: str
    new_content: str
    error: str = ""


def _generate_fix(finding: Finding, original_content: str) -> str:
    """Call the Claude API to generate a code fix. Returns new file content."""
    cfg = _get_loop_config()
    return api_call_with_retry({
        "model": cfg.fix_model,
        "max_tokens": cfg.fix_max_tokens,
        "system": FIX_SYSTEM_PROMPT,
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


def implement(finding: Finding, base_branch: str, dry_run: bool = False) -> ImplementResult:
    """Create a branch, generate and apply a code fix, commit. Returns result with file contents."""
    if dry_run:
        return ImplementResult(success=True, original_content="", new_content="")

    # Check if branch already exists
    if git_utils.branch_exists(finding.branch_name):
        return ImplementResult(
            success=False, original_content="", new_content="",
            error=f"branch exists: {finding.branch_name}",
        )

    # Read original file content
    file_path = os.path.join(REPO_ROOT, finding.file)
    if not os.path.exists(file_path):
        return ImplementResult(
            success=False, original_content="", new_content="",
            error=f"file not found: {finding.file}",
        )

    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        original_content = f.read()

    # Create branch
    git_utils.create_branch(finding.branch_name, base=base_branch)

    # Generate fix via API
    new_content = _generate_fix(finding, original_content)

    # Write to disk
    with open(file_path, "w", encoding="utf-8", newline="") as f:
        f.write(new_content)
    print(f"    Updated {finding.file}")

    # Commit
    git_utils.commit_all(
        f"improvement: {finding.dimension} — {finding.description[:60]}"
    )

    # Syntax check
    print("    Running syntax check...")
    if not git_utils.run_syntax_check():
        print("    ❌ Syntax error detected")
        return ImplementResult(
            success=False, original_content=original_content,
            new_content=new_content, error="syntax check failed",
        )

    return ImplementResult(
        success=True, original_content=original_content, new_content=new_content,
    )
