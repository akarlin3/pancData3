"""Tests for RAG integration into auditor, implementer, reviewer, and orchestrator."""

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from improvement_loop.loop_config import reset_config, LoopConfig


@pytest.fixture(autouse=True)
def _reset_loop_config():
    reset_config()
    yield
    reset_config()


def _make_config(rag_enabled: bool) -> LoopConfig:
    """Return a LoopConfig with rag_enabled set as specified."""
    return LoopConfig(rag_enabled=rag_enabled)


# ---------------------------------------------------------------------------
# Auditor: RAG-enabled uses get_context_for_audit
# ---------------------------------------------------------------------------

class TestAuditorRAGEnabled:
    def test_calls_get_context_for_audit(self, monkeypatch):
        """With rag_enabled=True, _call_audit_api uses retriever context."""
        monkeypatch.setattr(
            "improvement_loop.agents.auditor._get_loop_config",
            lambda: _make_config(rag_enabled=True),
        )

        rag_called = {"called": False}

        def mock_get_context_for_audit(context):
            rag_called["called"] = True
            return "=== rag_context.m ===\nRAG content here"

        monkeypatch.setattr(
            "improvement_loop.rag.retriever.get_context_for_audit",
            mock_get_context_for_audit,
        )

        captured_kwargs = {}

        def mock_api(kwargs):
            captured_kwargs.update(kwargs)
            return "[]"

        monkeypatch.setattr(
            "improvement_loop.agents.auditor.api_call_with_retry",
            mock_api,
        )

        from improvement_loop.agents.auditor import _call_audit_api
        _call_audit_api(iteration=1, context="test context")

        assert rag_called["called"]
        assert "RAG content here" in captured_kwargs["messages"][0]["content"]


class TestAuditorRAGDisabled:
    def test_calls_collect_source_files(self, monkeypatch):
        """With rag_enabled=False, _call_audit_api uses collect_source_files."""
        monkeypatch.setattr(
            "improvement_loop.agents.auditor._get_loop_config",
            lambda: _make_config(rag_enabled=False),
        )

        static_called = {"called": False}

        def mock_collect():
            static_called["called"] = True
            return "=== static_files ===\nstatic content"

        monkeypatch.setattr(
            "improvement_loop.agents.auditor.collect_source_files",
            mock_collect,
        )

        def mock_api(kwargs):
            return "[]"

        monkeypatch.setattr(
            "improvement_loop.agents.auditor.api_call_with_retry",
            mock_api,
        )

        from improvement_loop.agents.auditor import _call_audit_api
        _call_audit_api(iteration=1, context="test context")

        assert static_called["called"]


class TestAuditorRAGFallback:
    def test_falls_back_on_import_error(self, monkeypatch):
        """If RAG import fails, auditor falls back to collect_source_files."""
        monkeypatch.setattr(
            "improvement_loop.agents.auditor._get_loop_config",
            lambda: _make_config(rag_enabled=True),
        )

        # Make the RAG import fail
        import improvement_loop.rag.retriever as retriever_mod
        original_fn = retriever_mod.get_context_for_audit
        monkeypatch.setattr(
            "improvement_loop.rag.retriever.get_context_for_audit",
            lambda ctx: (_ for _ in ()).throw(RuntimeError("chromadb not available")),
        )

        static_called = {"called": False}

        def mock_collect():
            static_called["called"] = True
            return "=== fallback ===\nfallback content"

        monkeypatch.setattr(
            "improvement_loop.agents.auditor.collect_source_files",
            mock_collect,
        )

        def mock_api(kwargs):
            return "[]"

        monkeypatch.setattr(
            "improvement_loop.agents.auditor.api_call_with_retry",
            mock_api,
        )

        from improvement_loop.agents.auditor import _call_audit_api
        _call_audit_api(iteration=1, context="test")

        assert static_called["called"]


# ---------------------------------------------------------------------------
# Implementer: RAG-enabled adds "## Related context"
# ---------------------------------------------------------------------------

class TestImplementerRAGEnabled:
    def test_adds_related_context_section(self, monkeypatch):
        """With rag_enabled=True, _generate_fix appends ## Related context."""
        monkeypatch.setattr(
            "improvement_loop.agents.implementer._get_loop_config",
            lambda: _make_config(rag_enabled=True),
        )

        monkeypatch.setattr(
            "improvement_loop.rag.retriever.get_context_for_fix",
            lambda desc, path: "=== related.m ===\nrelated code",
        )

        captured_kwargs = {}

        def mock_api(kwargs):
            captured_kwargs.update(kwargs)
            return "fixed content"

        monkeypatch.setattr(
            "improvement_loop.agents.implementer.api_call_with_retry",
            mock_api,
        )

        from improvement_loop.evaluator import Finding
        from improvement_loop.agents.implementer import _generate_fix

        finding = Finding(
            dimension="correctness",
            file="pipeline/core/fit_models.m",
            description="Off-by-one",
            fix="Change < to <=",
            importance=5,
            branch_name="improvement/test-rag",
        )
        _generate_fix(finding, "% original code")

        user_content = captured_kwargs["messages"][0]["content"]
        assert "## Related context" in user_content
        assert "related code" in user_content


class TestImplementerRAGDisabled:
    def test_no_related_context_section(self, monkeypatch):
        """With rag_enabled=False, _generate_fix does not add ## Related context."""
        monkeypatch.setattr(
            "improvement_loop.agents.implementer._get_loop_config",
            lambda: _make_config(rag_enabled=False),
        )

        captured_kwargs = {}

        def mock_api(kwargs):
            captured_kwargs.update(kwargs)
            return "fixed content"

        monkeypatch.setattr(
            "improvement_loop.agents.implementer.api_call_with_retry",
            mock_api,
        )

        from improvement_loop.evaluator import Finding
        from improvement_loop.agents.implementer import _generate_fix

        finding = Finding(
            dimension="correctness",
            file="pipeline/core/fit_models.m",
            description="Off-by-one",
            fix="Change < to <=",
            importance=5,
            branch_name="improvement/test-rag",
        )
        _generate_fix(finding, "% original code")

        user_content = captured_kwargs["messages"][0]["content"]
        assert "## Related context" not in user_content


# ---------------------------------------------------------------------------
# Reviewer: RAG-enabled adds "## Codebase context"
# ---------------------------------------------------------------------------

class TestReviewerRAGEnabled:
    def test_adds_codebase_context_section(self, monkeypatch):
        """With rag_enabled=True, review appends ## Codebase context."""
        monkeypatch.setattr(
            "improvement_loop.agents.reviewer._get_loop_config",
            lambda: _make_config(rag_enabled=True),
        )

        monkeypatch.setattr(
            "improvement_loop.rag.retriever.get_context_for_review",
            lambda desc, path: "=== test_fit.m ===\ntest code",
        )

        captured_kwargs = {}

        def mock_api(kwargs):
            captured_kwargs.update(kwargs)
            return json.dumps({
                "verdict": "approve",
                "reasoning": "Looks good.",
                "risk_flags": [],
            })

        monkeypatch.setattr(
            "improvement_loop.agents.reviewer.api_call_with_retry",
            mock_api,
        )

        from improvement_loop.evaluator import Finding
        from improvement_loop.agents.reviewer import review

        finding = Finding(
            dimension="correctness",
            file="pipeline/core/fit_models.m",
            description="Off-by-one",
            fix="Change < to <=",
            importance=5,
            branch_name="improvement/test-rag",
        )
        review(finding, "old", "new")

        user_content = captured_kwargs["messages"][0]["content"]
        assert "## Codebase context" in user_content
        assert "test code" in user_content


# ---------------------------------------------------------------------------
# Orchestrator: calls build_index at loop start
# ---------------------------------------------------------------------------

class TestOrchestratorRAGInit:
    def test_calls_build_index_when_enabled(self, monkeypatch, tmp_path):
        """run_loop calls build_index at start when rag_enabled=True."""
        from improvement_loop import loop_tracker
        from improvement_loop import orchestrator_v2

        monkeypatch.setattr(
            "improvement_loop.orchestrator_v2._get_loop_config",
            lambda: _make_config(rag_enabled=True),
        )

        build_called = {"called": False}

        def mock_build_index(root, **kwargs):
            build_called["called"] = True

            class FakeCol:
                def count(self):
                    return 0
            return FakeCol()

        monkeypatch.setattr(
            "improvement_loop.rag.indexer.build_index",
            mock_build_index,
        )
        monkeypatch.setattr(
            "improvement_loop.rag.indexer.index_improvement_history",
            lambda **kw: None,
        )

        # Make auditor return no findings so we exit immediately
        monkeypatch.setattr(
            "improvement_loop.orchestrator_v2._audit",
            lambda iteration, context, dry_run: [],
        )
        monkeypatch.setattr(
            "improvement_loop.orchestrator_v2.git_utils.current_branch",
            lambda: "main",
        )

        log_file = str(tmp_path / "test_log.json")
        monkeypatch.setattr(loop_tracker, "LOG_FILE", log_file)

        # Make it exit after 1 iteration
        exit_entry = {
            "iteration": 1,
            "exit_condition_met": True,
            "audit_scores": {"overall": 9.0, "flags": []},
            "findings": [],
            "findings_count": 0,
            "high_priority_findings": 0,
            "branches_created": [],
            "branches_merged": [],
            "tests_passed": True,
        }
        monkeypatch.setattr(
            loop_tracker, "log_iteration",
            lambda **kwargs: exit_entry,
        )

        orchestrator_v2.run_loop(max_iterations=1, dry_run=False)

        assert build_called["called"]

    def test_skips_build_index_in_dry_run(self, monkeypatch, tmp_path):
        """run_loop does NOT call build_index in dry_run mode."""
        from improvement_loop import loop_tracker
        from improvement_loop import orchestrator_v2

        monkeypatch.setattr(
            "improvement_loop.orchestrator_v2._get_loop_config",
            lambda: _make_config(rag_enabled=True),
        )

        build_called = {"called": False}

        def mock_build_index(root, **kwargs):
            build_called["called"] = True

        monkeypatch.setattr(
            "improvement_loop.rag.indexer.build_index",
            mock_build_index,
        )

        log_file = str(tmp_path / "test_log.json")
        monkeypatch.setattr(loop_tracker, "LOG_FILE", log_file)

        orchestrator_v2.run_loop(max_iterations=1, dry_run=True)

        assert not build_called["called"]


# ---------------------------------------------------------------------------
# Graceful fallback: RAG exception doesn't crash the loop
# ---------------------------------------------------------------------------

class TestRAGGracefulFallback:
    def test_implementer_continues_on_rag_exception(self, monkeypatch):
        """If retriever raises, implementer proceeds without RAG context."""
        monkeypatch.setattr(
            "improvement_loop.agents.implementer._get_loop_config",
            lambda: _make_config(rag_enabled=True),
        )

        monkeypatch.setattr(
            "improvement_loop.rag.retriever.get_context_for_fix",
            lambda desc, path: (_ for _ in ()).throw(
                RuntimeError("chromadb broken")
            ),
        )

        captured_kwargs = {}

        def mock_api(kwargs):
            captured_kwargs.update(kwargs)
            return "fixed content"

        monkeypatch.setattr(
            "improvement_loop.agents.implementer.api_call_with_retry",
            mock_api,
        )

        from improvement_loop.evaluator import Finding
        from improvement_loop.agents.implementer import _generate_fix

        finding = Finding(
            dimension="correctness",
            file="pipeline/core/fit_models.m",
            description="Off-by-one",
            fix="Change < to <=",
            importance=5,
            branch_name="improvement/test-rag",
        )
        result = _generate_fix(finding, "% original code")

        assert result == "fixed content"
        # Should NOT have ## Related context since RAG failed
        assert "## Related context" not in captured_kwargs["messages"][0]["content"]

    def test_reviewer_continues_on_rag_exception(self, monkeypatch):
        """If retriever raises, reviewer proceeds without codebase context."""
        monkeypatch.setattr(
            "improvement_loop.agents.reviewer._get_loop_config",
            lambda: _make_config(rag_enabled=True),
        )

        monkeypatch.setattr(
            "improvement_loop.rag.retriever.get_context_for_review",
            lambda desc, path: (_ for _ in ()).throw(
                RuntimeError("chromadb broken")
            ),
        )

        captured_kwargs = {}

        def mock_api(kwargs):
            captured_kwargs.update(kwargs)
            return json.dumps({
                "verdict": "approve",
                "reasoning": "OK.",
                "risk_flags": [],
            })

        monkeypatch.setattr(
            "improvement_loop.agents.reviewer.api_call_with_retry",
            mock_api,
        )

        from improvement_loop.evaluator import Finding
        from improvement_loop.agents.reviewer import review

        finding = Finding(
            dimension="correctness",
            file="pipeline/core/fit_models.m",
            description="Off-by-one",
            fix="Change < to <=",
            importance=5,
            branch_name="improvement/test-rag",
        )
        verdict = review(finding, "old", "new")

        assert verdict.verdict == "approve"
        assert "## Codebase context" not in captured_kwargs["messages"][0]["content"]
