#!/usr/bin/env python3
"""
git_utils.py — Git helper functions for the pancData3 improvement loop.

Provides branch creation, committing, test execution, and merging
as thin wrappers around subprocess git calls.
"""
import os
import subprocess
import sys


def _run(cmd: list[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run a git command and return the result."""
    return subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=check,
        cwd=os.path.dirname(os.path.abspath(__file__)) or ".",
    )


def create_branch(branch_name: str) -> None:
    """Create and switch to a new branch off the current HEAD."""
    _run(["git", "checkout", "-b", branch_name])


def commit_all(message: str) -> None:
    """Stage all changes and commit with the given message."""
    _run(["git", "add", "-A"])
    _run(["git", "commit", "-m", message])


def run_python_tests() -> bool:
    """Run the Python test suite. Returns True if all tests pass."""
    result = _run(
        [sys.executable, "-m", "pytest", "analysis/tests", "-q", "--tb=short"],
        check=False,
    )
    return result.returncode == 0


def merge_branch(branch_name: str, target: str = "HEAD~0") -> None:
    """
    Merge the given branch into the current branch.

    Switches back to the branch we were on before the improvement branch
    was created (typically dev or main), merges, then returns.
    """
    # Remember current branch to return to after merge
    current = _run(["git", "rev-parse", "--abbrev-ref", "HEAD"]).stdout.strip()

    # If we're already on the improvement branch, switch to the parent
    if current == branch_name:
        # Go back to the branch we branched from
        _run(["git", "checkout", "-"])

    _run(["git", "merge", "--no-ff", branch_name, "-m", f"Merge {branch_name}"])


def switch_branch(branch_name: str) -> None:
    """Switch to an existing branch."""
    _run(["git", "checkout", branch_name])


def current_branch() -> str:
    """Return the name of the current branch."""
    return _run(["git", "rev-parse", "--abbrev-ref", "HEAD"]).stdout.strip()


def delete_branch(branch_name: str) -> None:
    """Delete a local branch."""
    _run(["git", "branch", "-d", branch_name])
