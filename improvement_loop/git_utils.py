"""Git utility functions for the pancData3 workflow orchestrator.

Replaces natural-language git instructions with direct subprocess calls.
"""

import logging
import pathlib
import shlex  # noqa: F401 — imported as a guard for future shell-string construction
import shutil
import subprocess
import sys
from typing import List

# Repository root — one level up from this file's directory
REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent

logger = logging.getLogger(__name__)


def _run(args: List[str], *, check: bool = True) -> subprocess.CompletedProcess:
    """Run a git (or other) command, optionally raising on failure."""
    result = subprocess.run(args, check=False, capture_output=True, text=True)
    if check and result.returncode != 0:
        raise RuntimeError(
            f"Command {args!r} failed (exit {result.returncode}):\n{result.stderr}"
        )
    return result


# ---------------------------------------------------------------------------
# Branch helpers
# ---------------------------------------------------------------------------

def current_branch() -> str:
    """Return the name of the current git branch."""
    result = _run(["git", "rev-parse", "--abbrev-ref", "HEAD"])
    return result.stdout.strip()


def branch_exists(branch_name: str) -> bool:
    """Return True if the branch exists locally or on origin."""
    local = _run(
        ["git", "rev-parse", "--verify", branch_name], check=False
    )
    if local.returncode == 0:
        return True
    remote = _run(
        ["git", "rev-parse", "--verify", f"origin/{branch_name}"], check=False
    )
    return remote.returncode == 0


def create_branch(branch_name: str, base: str = "v2.1-dev") -> None:
    """Create and checkout a new branch off *base*.

    Raises RuntimeError if the branch already exists.
    """
    if branch_exists(branch_name):
        raise RuntimeError(f"Branch '{branch_name}' already exists.")
    _run(["git", "checkout", "-b", branch_name, base])


def checkout(branch_name: str) -> None:
    """Checkout an existing branch.

    Raises RuntimeError if the branch does not exist.
    """
    if not branch_exists(branch_name):
        raise RuntimeError(f"Branch '{branch_name}' does not exist.")
    _run(["git", "checkout", branch_name])


# Alias used by orchestrator_v1.py
switch_branch = checkout


def merge_branch(
    source: str, target: str = "v2.1-dev", delete_after: bool = True
) -> None:
    """Merge *source* into *target* using ``--no-ff``.

    Checks out *target* first, merges, then optionally deletes *source*.
    Raises RuntimeError on merge conflict.
    """
    _run(["git", "checkout", target])
    result = _run(
        ["git", "merge", "--no-ff", source, "-m", f"Merge {source} into {target}"],
        check=False,
    )
    if result.returncode != 0:
        # Abort the failed merge so the repo is not left in a dirty state.
        _run(["git", "merge", "--abort"], check=False)
        raise RuntimeError(
            f"Merge conflict merging '{source}' into '{target}':\n{result.stderr}"
        )
    if delete_after:
        _run(["git", "branch", "-d", source])


# ---------------------------------------------------------------------------
# Test runners
# ---------------------------------------------------------------------------

def run_syntax_check() -> bool:
    """Run py_compile on all .py files in analysis/ to catch syntax errors
    before running the full test suite.

    Returns True if all files compile, False on the first syntax error.
    """
    analysis_dir = REPO_ROOT / "analysis"
    py_files = sorted(analysis_dir.rglob("*.py"))

    for py_file in py_files:
        result = subprocess.run(
            [sys.executable, "-m", "py_compile", str(py_file)],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(f"    ❌  Syntax error in {py_file.relative_to(REPO_ROOT)}:")
            print(f"        {result.stderr.strip()}")
            return False

    print(f"    ✅  Syntax check passed ({len(py_files)} files)")
    return True


def run_python_tests(capture_output: bool = False) -> "bool | tuple[bool, str]":
    """Run the Python test suite via pytest.

    If *capture_output* is False (default), streams to stdout and returns bool.
    If *capture_output* is True, returns ``(passed, output_text)`` so callers
    can inspect failure details for self-healing.
    """
    print(f"    [debug] Running tests from: {REPO_ROOT}")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "analysis/tests/", "-q", "--tb=short",
            "--ignore=analysis/tests/test_evaluator_finding.py",
            "--ignore=analysis/tests/test_git_utils.py",
            "--ignore=analysis/tests/test_loop_tracker.py",
            "--ignore=analysis/tests/test_orchestrator.py"],
        cwd=REPO_ROOT,
        check=False,
        capture_output=capture_output,
        text=capture_output,
    )
    print(f"    [debug] pytest exit code: {result.returncode}")
    if capture_output:
        combined = (result.stdout or "") + "\n" + (result.stderr or "")
        return result.returncode == 0, combined
    return result.returncode == 0


def run_matlab_tests() -> bool:
    """Run the MATLAB test suite.

    Returns True if exit code 0, False otherwise.
    If ``matlab`` is not on PATH, logs a warning and returns True
    (non-blocking for environments without MATLAB).
    """
    if shutil.which("matlab") is None:
        logger.warning("matlab not found on PATH — skipping MATLAB tests")
        return True

    result = subprocess.run(
        [
            "matlab",
            "-batch",
            (
                "results = runtests('pipeline/tests'); "
                "if any([results.Failed]), exit(1); else exit(0); end"
            ),
        ],
        cwd=REPO_ROOT,
        check=False,
    )
    return result.returncode == 0


# ---------------------------------------------------------------------------
# Branch slug / staging helpers
# ---------------------------------------------------------------------------

def sanitize_branch_slug(text: str, max_len: int = 50) -> str:
    """Convert arbitrary text to a valid git branch slug.

    * Lowercase
    * Replace spaces and non-alphanumeric characters with hyphens
    * Collapse consecutive hyphens
    * Strip leading/trailing hyphens
    * Truncate to *max_len*
    * Never returns an empty string — falls back to ``"fix"``.
    """
    import re

    slug = text.lower()
    slug = re.sub(r"[^a-z0-9]+", "-", slug)
    slug = re.sub(r"-{2,}", "-", slug)
    slug = slug.strip("-")
    slug = slug[:max_len].rstrip("-")
    return slug or "fix"


def get_staged_files() -> list[str]:
    """Return list of files currently staged for commit."""
    result = _run(["git", "diff", "--cached", "--name-only"])
    lines = result.stdout.strip()
    return lines.splitlines() if lines else []


def commit_all(message: str) -> None:
    """Stage all modified/new tracked files and commit with *message*.

    No-op if nothing to commit.
    """
    _run(["git", "add", "-A"])
    status = _run(["git", "status", "--porcelain"])
    if not status.stdout.strip():
        return  # nothing to commit
    _run(["git", "commit", "-m", message])
