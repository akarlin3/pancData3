"""Unit tests for run_analysis.py — orchestrator utility functions.

Tests the non-main() helper functions:
- _mask_key_value: aggressive API key masking by length
- _mask_api_keys: redaction of API key patterns in free text
- TeeWriter: dual-stream writer for stdout/log file tee
- _handle_windows_path: Windows long-path handling
- _check_requirements: requirements.txt satisfaction check
- _run_script: subprocess script execution with timeout handling

main() is NOT tested here (that is integration-level).
"""

from __future__ import annotations

import io
import os
import subprocess
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Ensure analysis/ is importable (conftest.py already does this, but be safe).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from run_analysis import (
    ANALYSIS_DIR,
    TeeWriter,
    _check_requirements,
    _handle_windows_path,
    _mask_api_keys,
    _mask_key_value,
    _run_script,
)


# ---------------------------------------------------------------------------
# _mask_key_value
# ---------------------------------------------------------------------------

class TestMaskKeyValue:
    """Verify aggressive key masking based on key length."""

    def test_mask_short_key(self):
        """Keys shorter than 12 chars should be fully masked with char count."""
        result = _mask_key_value("abc123")
        assert result == "****[6 chars]"
        # No original characters should be visible.
        assert "abc" not in result

    def test_mask_long_key(self):
        """Keys >= 12 chars show first 2 chars plus masked length indicator."""
        key = "sk-ant-abc123def456ghi789jkl012mno345pqr"
        result = _mask_key_value(key)
        assert result.startswith("sk")
        assert "****" in result
        assert f"[{len(key)} chars]" in result
        # The full key should NOT appear in the output.
        assert key not in result

    def test_mask_empty_key(self):
        """Empty string produces a fully masked indicator with 0 chars."""
        result = _mask_key_value("")
        assert result == "****[0 chars]"

    def test_mask_exactly_12_chars(self):
        """A key of exactly 12 chars should show the first 2 characters."""
        key = "abcdefghijkl"
        result = _mask_key_value(key)
        assert result.startswith("ab")
        assert "****[12 chars]" in result

    def test_mask_11_chars_fully_masked(self):
        """A key of 11 chars (< 12) should be fully masked."""
        result = _mask_key_value("12345678901")
        assert result == "****[11 chars]"


# ---------------------------------------------------------------------------
# _mask_api_keys
# ---------------------------------------------------------------------------

class TestMaskApiKeys:
    """Verify redaction of various API key patterns in text."""

    def test_mask_env_style_key(self):
        """ENV-style API_KEY=value should have the value masked."""
        text = "API_KEY=sk-ant-abc123def456"
        result = _mask_api_keys(text)
        assert "abc123def456" not in result
        assert "API_KEY=" in result
        assert "****" in result

    def test_mask_bearer_token(self):
        """Bearer token in HTTP header should be masked."""
        text = "Authorization: Bearer sk1234567890abcdef"
        result = _mask_api_keys(text)
        assert "sk1234567890abcdef" not in result
        assert "Bearer" in result
        assert "****" in result

    def test_mask_anthropic_key(self):
        """Standalone Anthropic key prefix (sk-ant-...) should be masked."""
        text = "key is sk-ant-abcdefghijklmnop"
        result = _mask_api_keys(text)
        assert "abcdefghijklmnop" not in result
        assert "****" in result

    def test_mask_gemini_key(self):
        """Google/Gemini API key prefix (AIza...) should be masked."""
        text = "key is AIzaSyABCDEFGHIJKLMNOP"
        result = _mask_api_keys(text)
        assert "ABCDEFGHIJKLMNOP" not in result
        assert "****" in result

    def test_no_masking_short_values(self):
        """Values shorter than 8 chars should NOT be masked (not real keys)."""
        text = "API_KEY=short"
        result = _mask_api_keys(text)
        # "short" is only 5 chars — no regex should match it.
        assert "short" in result

    def test_mask_config_style_key(self):
        """Config-style 'api_key: value' should have the value masked."""
        text = "anthropic_api_key: sk-ant-test1234567890"
        result = _mask_api_keys(text)
        assert "test1234567890" not in result
        assert "****" in result

    def test_preserves_unrelated_text(self):
        """Text without API key patterns should pass through unchanged."""
        text = "This is a normal log line with no secrets."
        result = _mask_api_keys(text)
        assert result == text


# ---------------------------------------------------------------------------
# TeeWriter
# ---------------------------------------------------------------------------

class TestTeeWriter:
    """Verify dual-stream writer behavior."""

    def test_tee_write_both_streams(self):
        """write() should send data to both terminal and log file streams."""
        terminal = io.StringIO()
        log_file = io.StringIO()
        tee = TeeWriter(terminal, log_file)
        tee.write("hello world")
        assert terminal.getvalue() == "hello world"
        assert log_file.getvalue() == "hello world"

    def test_tee_write_returns_length(self):
        """write() should return the number of characters written."""
        tee = TeeWriter(io.StringIO(), io.StringIO())
        n = tee.write("test")
        assert n == 4

    def test_tee_flush_both(self):
        """flush() should call flush on both underlying streams."""
        terminal = MagicMock()
        log_file = MagicMock()
        tee = TeeWriter(terminal, log_file)
        tee.flush()
        terminal.flush.assert_called_once()
        log_file.flush.assert_called_once()

    def test_tee_isatty_false(self):
        """isatty() returns False when the terminal is a StringIO (no tty)."""
        tee = TeeWriter(io.StringIO(), io.StringIO())
        assert tee.isatty() is False

    def test_tee_fileno_raises(self):
        """fileno() raises OSError when the terminal lacks a file descriptor."""
        tee = TeeWriter(io.StringIO(), io.StringIO())
        with pytest.raises(OSError, match="no fileno"):
            tee.fileno()

    def test_tee_encoding_attribute(self):
        """TeeWriter should expose an encoding attribute."""
        terminal = io.StringIO()
        tee = TeeWriter(terminal, io.StringIO())
        # TeeWriter uses getattr(terminal, "encoding", "utf-8").
        # StringIO.encoding is None on most Python versions, so the
        # getattr returns None (attribute exists but is None).
        assert hasattr(tee, "encoding")


# ---------------------------------------------------------------------------
# _handle_windows_path
# ---------------------------------------------------------------------------

class TestHandleWindowsPath:
    """Verify Windows-specific and non-Windows path handling."""

    def test_non_windows_returns_resolved(self, tmp_path):
        """On non-Windows platforms, _handle_windows_path returns path.resolve()."""
        # tmp_path is already resolved; create a subdirectory for a clean test.
        sub = tmp_path / "test_dir"
        sub.mkdir()
        with patch("run_analysis.os.name", "posix"):
            result = _handle_windows_path(sub)
        assert result == sub.resolve()
        assert isinstance(result, Path)


# ---------------------------------------------------------------------------
# _check_requirements
# ---------------------------------------------------------------------------

class TestCheckRequirements:
    """Verify requirements satisfaction checking."""

    def test_check_requirements_when_file_missing(self, tmp_path):
        """When requirements.txt does not exist, return True (skip check)."""
        fake_dir = tmp_path / "no_such_dir"
        fake_dir.mkdir()
        fake_req = fake_dir / "requirements.txt"
        # Temporarily point ANALYSIS_DIR at the fake location.
        with patch("run_analysis.ANALYSIS_DIR", fake_dir):
            result = _check_requirements()
        assert result is True

    def test_check_requirements_all_satisfied(self, tmp_path):
        """When all packages are installed, return True."""
        req_file = tmp_path / "requirements.txt"
        req_file.write_text("pytest\ntqdm\n")
        with patch("run_analysis.ANALYSIS_DIR", tmp_path):
            result = _check_requirements()
        assert result is True

    def test_check_requirements_missing_package(self, tmp_path):
        """When a package is missing, return False."""
        req_file = tmp_path / "requirements.txt"
        req_file.write_text("nonexistent_package_xyz_999\n")
        with patch("run_analysis.ANALYSIS_DIR", tmp_path):
            result = _check_requirements()
        assert result is False

    def test_check_requirements_skips_comments(self, tmp_path):
        """Comment lines and blank lines should be ignored."""
        req_file = tmp_path / "requirements.txt"
        req_file.write_text("# This is a comment\n\npytest\n")
        with patch("run_analysis.ANALYSIS_DIR", tmp_path):
            result = _check_requirements()
        assert result is True


# ---------------------------------------------------------------------------
# _run_script
# ---------------------------------------------------------------------------

class TestRunScript:
    """Verify subprocess script execution and error handling."""

    def test_run_nonexistent_script(self, capsys):
        """When the script file does not exist, return False and print SKIP."""
        result = _run_script(
            "totally_nonexistent_script_xyz.py",
            Path("/tmp"),
            log_file=None,
            timeout=10,
        )
        assert result is False
        captured = capsys.readouterr()
        assert "SKIP" in captured.out

    def test_run_script_timeout(self, tmp_path):
        """When subprocess times out, return False."""
        # Create a dummy script that exists (so the file-existence check passes).
        script = tmp_path / "slow_script.py"
        script.write_text("import time; time.sleep(999)")

        with patch("run_analysis.ANALYSIS_DIR", tmp_path), \
             patch("run_analysis.subprocess.run", side_effect=subprocess.TimeoutExpired(cmd="test", timeout=1)):
            result = _run_script(
                "slow_script.py",
                tmp_path,
                log_file=None,
                timeout=1,
            )
        assert result is False

    def test_run_script_success(self, tmp_path):
        """When subprocess succeeds (exit 0), return True."""
        script = tmp_path / "good_script.py"
        script.write_text("print('ok')")

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "ok\n"
        mock_result.stderr = ""

        with patch("run_analysis.ANALYSIS_DIR", tmp_path), \
             patch("run_analysis.subprocess.run", return_value=mock_result):
            result = _run_script(
                "good_script.py",
                tmp_path,
                log_file=None,
                timeout=30,
            )
        assert result is True
