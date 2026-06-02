#!/usr/bin/env python3
"""Subprocess and output helpers for the pancData3 analysis orchestrator.

This module holds the stateless support machinery used by
:mod:`run_analysis`:

- ``ANALYSIS_DIR`` — absolute path to the ``analysis/`` directory.
- ``_handle_windows_path`` — Windows long/UNC path normalisation.
- :class:`TeeWriter` — tees stdout/stderr to a terminal and a log file.
- ``_check_requirements`` / ``_run_tests`` — pre-flight checks.
- ``_mask_key_value`` / ``_mask_api_keys`` — redact API keys from output.
- ``_run_script`` — run a sibling analysis script as a subprocess.

All names defined here are re-exported from :mod:`run_analysis` for backward
compatibility (tests import several of them via ``from run_analysis import``).
"""

from __future__ import annotations

import io
import os
import re
import subprocess
import sys
import time
from pathlib import Path

# Absolute path to the analysis/ directory (where this file lives, alongside
# run_analysis.py).
ANALYSIS_DIR = Path(__file__).resolve().parent


def _handle_windows_path(path: Path) -> Path:
    """Handle Windows-specific path issues including UNC paths and long paths.

    On Windows, applies pathlib.Path.resolve() and handles:
    - Long path prefix (\\\\?\\) for paths longer than 260 characters
    - UNC path detection and appropriate handling
    - Network path normalization

    Parameters
    ----------
    path : Path
        The path to handle

    Returns
    -------
    Path
        The processed path with Windows-specific handling applied
    """
    if os.name != 'nt':
        return path.resolve()

    try:
        # First resolve the path normally
        resolved = path.resolve()
        path_str = str(resolved)

        # Determine path type
        is_unc = path_str.startswith('\\\\')
        path_type = 'local' if not is_unc else 'unc'

        # Handle UNC paths
        if path_type == 'unc':
            if len(path_str) > 260 and not path_str.startswith('\\\\?\\UNC\\'):
                return Path('\\\\?\\UNC\\' + path_str[2:])
        # Handle local paths
        else:
            if len(path_str) > 260 and not path_str.startswith('\\\\?\\'):
                return Path('\\\\?\\' + path_str)

        return resolved

    except (OSError, ValueError) as e:
        print(f"  [WARN] Path handling issue: {e}")
        # Fall back to the original path if resolution fails
        return path


class TeeWriter:
    """Write to both a terminal stream and a log file simultaneously.

    Implements enough of the file-like interface (``write``, ``flush``,
    ``isatty``, ``encoding``, ``fileno``) to be usable as a drop-in
    ``sys.stdout`` / ``sys.stderr`` replacement.
    """

    def __init__(self, terminal: io.TextIOBase, log_file: io.TextIOBase):
        self.terminal = terminal
        self.log_file = log_file
        self.encoding = getattr(terminal, "encoding", "utf-8")

    def write(self, message: str) -> int:
        self.terminal.write(message)
        self.log_file.write(message)
        return len(message)

    def flush(self) -> None:
        self.terminal.flush()
        self.log_file.flush()

    def isatty(self) -> bool:
        """Return the terminal stream's isatty status."""
        return hasattr(self.terminal, "isatty") and self.terminal.isatty()

    def fileno(self) -> int:
        """Delegate fileno to the terminal stream (for subprocess piping).

        Raises OSError when the terminal stream lacks a file descriptor
        (e.g., StringIO during pytest capture), matching the standard
        io.UnsupportedOperation behavior.
        """
        try:
            return self.terminal.fileno()
        except (AttributeError, OSError):
            raise OSError("TeeWriter: underlying stream has no fileno")


def _check_requirements() -> bool:
    """Verify that all packages listed in requirements.txt are installed.

    Returns True if every requirement is satisfied, False otherwise.
    Missing packages are printed to stdout.

    Uses ``packaging.requirements.Requirement`` to properly parse each
    requirement line, correctly handling extras syntax (e.g.,
    ``package[extra]>=1.0``), URL-based requirements (e.g.,
    ``git+https://...``), and complex version specifiers.
    """
    req_file = ANALYSIS_DIR / "requirements.txt"
    if not req_file.exists():
        print("  [WARN] requirements.txt not found — skipping requirement check")
        return True

    missing: list[str] = []
    import importlib.metadata as _meta
    from packaging.requirements import InvalidRequirement, Requirement
    from packaging.utils import canonicalize_name

    # Get all installed packages once (normalize names for reliable lookup)
    installed_packages = {canonicalize_name(dist.metadata["name"]): dist for dist in _meta.distributions()}

    for line in req_file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Skip options lines (e.g., ``--index-url``, ``-f``, etc.)
        if line.startswith("-"):
            continue
        try:
            req = Requirement(line)
            pkg_name = req.name
        except InvalidRequirement:
            # If packaging cannot parse the line (e.g., malformed or
            # very unusual syntax), warn and skip rather than crash.
            print(f"  [WARN] Could not parse requirement line, skipping: {line}")
            continue

        if canonicalize_name(pkg_name) not in installed_packages:
            missing.append(line)

    if missing:
        print("  Missing Python packages:")
        for m in missing:
            print(f"    - {m}")
        print(f"\n  Install with:  pip install -r {req_file}")
        return False

    print("  All requirements satisfied.")
    return True


def _run_tests() -> bool:
    """Run the analysis test suite via pytest.

    Returns True if all tests pass, False otherwise.
    """
    test_dir = ANALYSIS_DIR / "tests"
    if not test_dir.is_dir():
        print("  [WARN] tests/ directory not found — skipping test run")
        return True

    print("  Running analysis test suite (pytest) ...")
    t0 = time.time()
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-v", str(test_dir)],
        cwd=str(ANALYSIS_DIR),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    elapsed = time.time() - t0

    # Show the last portion of output (summary) to keep output manageable.
    if result.stdout:
        # Print full output so failures are visible.
        sys.stdout.write(result.stdout)
    if result.stderr:
        sys.stderr.write(result.stderr)

    if result.returncode == 0:
        print(f"  All tests passed ({elapsed:.1f}s)")
        return True
    else:
        print(f"  Tests FAILED (exit code {result.returncode}, {elapsed:.1f}s)")
        return False


def _mask_key_value(key_value: str) -> str:
    """Mask an API key value aggressively based on its length.

    For keys shorter than 12 characters, the entire value is fully masked
    (no characters revealed) to prevent information leakage on short keys.
    For longer keys, at most the first 2 characters are shown followed by
    ``****`` and a bracketed character count indicating the total key
    length, e.g. ``sk****[39 chars]``.

    This avoids revealing well-known provider prefixes (e.g., ``sk-ant-``,
    ``AIza``) and reduces the information available to an attacker for
    keys of any length.

    Parameters
    ----------
    key_value : str
        The raw API key string to mask.

    Returns
    -------
    str
        The masked key string.
    """
    key_len = len(key_value)
    if key_len < 12:
        # Fully mask short keys — showing any characters is too revealing
        return f"****[{key_len} chars]"
    # For longer keys, show at most the first 2 characters
    return key_value[:2] + f"****[{key_len} chars]"


def _mask_api_keys(text: str) -> str:
    """Redact API key values that may appear in subprocess output.

    Catches common patterns including:
    - ``GEMINI_API_KEY=sk-...`` and similar ENV-style dumps
    - ``gemini_api_key: sk-...`` and ``anthropic_api_key: sk-ant-...``
      (suffix-style key names as found in analysis_config.json)
    - Bearer tokens in HTTP headers
    - Known provider key prefixes (``sk-ant-``, ``AIza``, ``sk-``)
      appearing as standalone values

    Keys shorter than 8 characters are left untouched (unlikely to be real).
    All matched keys are aggressively masked: keys shorter than 12 characters
    are fully redacted, and longer keys show at most 2 leading characters
    plus a length indicator.
    """
    # 1. Mask patterns where the field name ends with _KEY, _SECRET, _TOKEN
    #    (case-insensitive), covering both PREFIX_API_KEY and prefix_api_key
    #    style names.  Matches ``=`` or ``:`` as the separator.
    text = re.sub(
        r'((?:API_KEY|SECRET|TOKEN)\s*[=:]\s*)["\']?([A-Za-z0-9_\-]{8,})["\']?',
        lambda m: m.group(1) + _mask_key_value(m.group(2)),
        text,
        flags=re.IGNORECASE,
    )

    # 2. Mask patterns where the field name ends with _api_key (suffix style),
    #    e.g. ``gemini_api_key: AIza...`` or ``anthropic_api_key = sk-ant-...``
    text = re.sub(
        r'(api_key\s*[=:]\s*)["\']?([A-Za-z0-9_\-]{8,})["\']?',
        lambda m: m.group(1) + _mask_key_value(m.group(2)),
        text,
        flags=re.IGNORECASE,
    )

    # 3. Mask Bearer tokens in HTTP headers
    text = re.sub(
        r'(Bearer\s+)([A-Za-z0-9_\-]{8,})',
        lambda m: m.group(1) + _mask_key_value(m.group(2)),
        text,
        flags=re.IGNORECASE,
    )

    # 4. Mask known provider key prefixes that may appear as standalone values
    #    (e.g. in tracebacks printing repr of a dict value).
    #    - sk-ant-  : Anthropic keys
    #    - AIza     : Google/Gemini keys
    #    - sk-      : OpenAI-style keys (also matches sk-ant- but we run this
    #                 after the more specific pattern above, so duplicates are
    #                 harmless — already masked text won't re-match).
    text = re.sub(
        r'(sk-ant-[A-Za-z0-9_\-]{8,})',
        lambda m: _mask_key_value(m.group(1)),
        text,
    )
    text = re.sub(
        r'(AIza[A-Za-z0-9_\-]{8,})',
        lambda m: _mask_key_value(m.group(1)),
        text,
    )
    text = re.sub(
        r'(sk-[A-Za-z0-9_\-]{8,})',
        lambda m: _mask_key_value(m.group(1)),
        text,
    )

    return text


def _run_script(
    name: str,
    folder: Path,
    log_file: io.TextIOBase | None = None,
    timeout: float | None = None,
) -> bool:
    """Run a sibling Python script as a subprocess.

    Parameters
    ----------
    name : str
        Filename of the script to run (e.g. ``"parse_log_metrics.py"``).
        Must reside in the same directory as this orchestrator.
    folder : Path
        Path to the ``saved_files_*`` output folder, passed as the first
        positional argument to the child script.
    log_file : io.TextIOBase or None
        If provided, subprocess stdout and stderr are captured and written
        to both the terminal and this log file.
    timeout : float or None
        Maximum wall-clock seconds to allow the subprocess to run.  If
        exceeded the process is killed and the step is reported as failed.

    Returns
    -------
    bool
        ``True`` if the script exited with code 0, ``False`` otherwise
        (including when the script file is not found or timed out).
    """
    script = ANALYSIS_DIR / name
    if not script.exists():
        print(f"  [SKIP] {name} not found")
        return False

    print(f"\n  Running {name} ...")
    t0 = time.time()
    # Run the child script using the same Python interpreter as this process.
    # Apply Windows path handling to the folder path
    folder_handled = _handle_windows_path(folder)

    try:
        result = subprocess.run(
            [sys.executable, str(script), str(folder_handled)],
            cwd=str(ANALYSIS_DIR),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        msg = f"  FAILED: {name} (timed out after {elapsed:.1f}s)"
        print(msg)
        if log_file is not None:
            log_file.write(msg + "\n")  # type: ignore
        return False
    elapsed = time.time() - t0

    # Emit captured output to terminal and log file, masking any API keys
    # that may have leaked into stdout/stderr (e.g., verbose HTTP logging).
    for stream_output in (result.stdout, result.stderr):
        if stream_output:
            masked = _mask_api_keys(stream_output)
            sys.stdout.write(masked)
            if log_file is not None:
                log_file.write(masked)  # type: ignore

    if result.returncode == 0:
        print(f"  Done: {name} ({elapsed:.1f}s)")
        return True
    else:
        print(f"  FAILED: {name} (exit code {result.returncode}, {elapsed:.1f}s)")
        return False
