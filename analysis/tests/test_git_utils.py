"""Tests for git_utils module."""

import subprocess
import sys
from unittest import mock

import pytest

# Ensure repo root is importable
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2]))

from improvement_loop import git_utils  # noqa: E402


# ---------------------------------------------------------------------------
# sanitize_branch_slug
# ---------------------------------------------------------------------------

class TestSanitizeBranchSlug:
    def test_empty_string(self):
        assert git_utils.sanitize_branch_slug("") == "fix"

    def test_all_special_chars(self):
        assert git_utils.sanitize_branch_slug("!!!@@@###") == "fix"

    def test_already_valid(self):
        assert git_utils.sanitize_branch_slug("my-feature") == "my-feature"

    def test_too_long(self):
        result = git_utils.sanitize_branch_slug("a" * 100, max_len=50)
        assert len(result) == 50

    def test_leading_trailing_hyphens(self):
        result = git_utils.sanitize_branch_slug("--hello--world--")
        assert result == "hello-world"

    def test_spaces_replaced(self):
        assert git_utils.sanitize_branch_slug("my cool feature") == "my-cool-feature"

    def test_uppercase_lowered(self):
        assert git_utils.sanitize_branch_slug("Fix-BUG") == "fix-bug"

    def test_consecutive_special_collapsed(self):
        assert git_utils.sanitize_branch_slug("a!!!b") == "a-b"

    def test_truncation_strips_trailing_hyphen(self):
        # "abcde-fgh" truncated to 6 -> "abcde-" -> trailing hyphen stripped
        result = git_utils.sanitize_branch_slug("abcde-fgh", max_len=6)
        assert result == "abcde"


# ---------------------------------------------------------------------------
# branch_exists
# ---------------------------------------------------------------------------

class TestBranchExists:
    def test_nonexistent_branch(self):
        assert git_utils.branch_exists("this-branch-does-not-exist-xyz-999") is False


# ---------------------------------------------------------------------------
# current_branch
# ---------------------------------------------------------------------------

class TestCurrentBranch:
    def test_returns_nonempty_string(self):
        branch = git_utils.current_branch()
        assert isinstance(branch, str)
        assert len(branch) > 0


# ---------------------------------------------------------------------------
# run_python_tests
# ---------------------------------------------------------------------------

class TestRunPythonTests:
    def test_returns_true_on_success(self):
        with mock.patch("subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess([], 0)
            assert git_utils.run_python_tests() is True

    def test_returns_false_on_failure(self):
        with mock.patch("subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess([], 1)
            assert git_utils.run_python_tests() is False

    def test_passes_correct_args(self):
        with mock.patch("subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess([], 0)
            git_utils.run_python_tests()
            args = mock_run.call_args[0][0]
            assert any("pytest" in a for a in args)
            assert "analysis/tests/" in args
            assert "-q" in args
            assert "--tb=short" in args


# ---------------------------------------------------------------------------
# create_branch — mocked
# ---------------------------------------------------------------------------

class TestCreateBranch:
    def test_raises_if_branch_exists(self):
        with mock.patch("improvement_loop.git_utils.branch_exists", return_value=True):
            with pytest.raises(RuntimeError, match="already exists"):
                git_utils.create_branch("existing-branch")

    def test_calls_checkout_b(self):
        with mock.patch("improvement_loop.git_utils.branch_exists", return_value=False):
            with mock.patch("improvement_loop.git_utils._run") as mock_run:
                mock_run.return_value = subprocess.CompletedProcess([], 0)
                git_utils.create_branch("new-branch", base="main")
                mock_run.assert_called_once_with(
                    ["git", "checkout", "-b", "new-branch", "main"]
                )


# ---------------------------------------------------------------------------
# checkout — mocked
# ---------------------------------------------------------------------------

class TestCheckout:
    def test_raises_if_branch_missing(self):
        with mock.patch("improvement_loop.git_utils.branch_exists", return_value=False):
            with pytest.raises(RuntimeError, match="does not exist"):
                git_utils.checkout("nonexistent")

    def test_calls_git_checkout(self):
        with mock.patch("improvement_loop.git_utils.branch_exists", return_value=True):
            with mock.patch("improvement_loop.git_utils._run") as mock_run:
                mock_run.return_value = subprocess.CompletedProcess([], 0)
                git_utils.checkout("my-branch")
                mock_run.assert_called_once_with(["git", "checkout", "my-branch"])


# ---------------------------------------------------------------------------
# switch_branch — alias for checkout
# ---------------------------------------------------------------------------

class TestSwitchBranch:
    def test_is_alias_for_checkout(self):
        assert git_utils.switch_branch is git_utils.checkout

    def test_calls_git_checkout(self):
        with mock.patch("improvement_loop.git_utils.branch_exists", return_value=True):
            with mock.patch("improvement_loop.git_utils._run") as mock_run:
                mock_run.return_value = subprocess.CompletedProcess([], 0)
                git_utils.switch_branch("my-branch")
                mock_run.assert_called_once_with(["git", "checkout", "my-branch"])


# ---------------------------------------------------------------------------
# merge_branch — mocked
# ---------------------------------------------------------------------------

class TestMergeBranch:
    def test_raises_on_conflict(self):
        def side_effect(args, *, check=True):
            if args[0:2] == ["git", "merge"]:
                return subprocess.CompletedProcess(
                    args, 1, stdout="", stderr="CONFLICT"
                )
            return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

        with mock.patch("improvement_loop.git_utils._run", side_effect=side_effect):
            with pytest.raises(RuntimeError, match="Merge conflict"):
                git_utils.merge_branch("feat", target="main")

    def test_deletes_source_by_default(self):
        with mock.patch("improvement_loop.git_utils._run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(
                [], 0, stdout="", stderr=""
            )
            git_utils.merge_branch("feat", target="main", delete_after=True)
            calls = [c.args[0] for c in mock_run.call_args_list]
            assert ["git", "branch", "-d", "feat"] in calls

    def test_skips_delete_when_false(self):
        with mock.patch("improvement_loop.git_utils._run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(
                [], 0, stdout="", stderr=""
            )
            git_utils.merge_branch("feat", target="main", delete_after=False)
            calls = [c.args[0] for c in mock_run.call_args_list]
            assert ["git", "branch", "-d", "feat"] not in calls


# ---------------------------------------------------------------------------
# commit_all — mocked
# ---------------------------------------------------------------------------

class TestCommitAll:
    def test_noop_when_nothing_to_commit(self):
        def side_effect(args, *, check=True):
            return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

        with mock.patch("improvement_loop.git_utils._run", side_effect=side_effect) as mock_run:
            git_utils.commit_all("msg")
            calls = [c.args[0] for c in mock_run.call_args_list]
            assert ["git", "commit", "-m", "msg"] not in calls

    def test_commits_when_changes_present(self):
        def side_effect(args, *, check=True):
            if args[:2] == ["git", "status"]:
                return subprocess.CompletedProcess(
                    args, 0, stdout="M file.py\n", stderr=""
                )
            return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

        with mock.patch("improvement_loop.git_utils._run", side_effect=side_effect) as mock_run:
            git_utils.commit_all("my message")
            calls = [c.args[0] for c in mock_run.call_args_list]
            assert ["git", "commit", "-m", "my message"] in calls


# ---------------------------------------------------------------------------
# run_syntax_check
# ---------------------------------------------------------------------------

class TestRunSyntaxCheck:
    def test_returns_false_on_syntax_error(self, tmp_path, capsys):
        bad_file = tmp_path / "bad.py"
        bad_file.write_text("for\n")

        with mock.patch.object(git_utils, "REPO_ROOT", tmp_path):
            # Create the analysis/ subdirectory that run_syntax_check globs
            analysis_dir = tmp_path / "analysis"
            analysis_dir.mkdir()
            bad_analysis_file = analysis_dir / "broken.py"
            bad_analysis_file.write_text("for\n")

            result = git_utils.run_syntax_check()

        assert result is False
        captured = capsys.readouterr()
        assert "broken.py" in captured.out
        assert "Syntax error" in captured.out

    def test_returns_true_when_all_valid(self, tmp_path, capsys):
        analysis_dir = tmp_path / "analysis"
        analysis_dir.mkdir()
        good_file = analysis_dir / "good.py"
        good_file.write_text("x = 1\n")

        with mock.patch.object(git_utils, "REPO_ROOT", tmp_path):
            result = git_utils.run_syntax_check()

        assert result is True
        captured = capsys.readouterr()
        assert "Syntax check passed" in captured.out
        assert "1 files" in captured.out

    def test_returns_true_when_no_files(self, tmp_path, capsys):
        analysis_dir = tmp_path / "analysis"
        analysis_dir.mkdir()

        with mock.patch.object(git_utils, "REPO_ROOT", tmp_path):
            result = git_utils.run_syntax_check()

        assert result is True
        captured = capsys.readouterr()
        assert "0 files" in captured.out
