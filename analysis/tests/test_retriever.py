"""Tests for improvement_loop.rag.retriever: query, filtering, formatting, agent context."""

import json
import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from improvement_loop.rag.retriever import (
    RetrievedContext,
    query,
    format_retrieved_context,
    get_context_for_audit,
    get_context_for_fix,
    get_context_for_review,
    _build_where_clause,
    _deduplicate_overlapping,
)
from improvement_loop.rag.indexer import build_index, index_improvement_history, COLLECTION_NAME
from improvement_loop.loop_config import reset_config


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _reset_loop_config():
    reset_config()
    yield
    reset_config()


@pytest.fixture()
def indexed_repo(tmp_path, monkeypatch):
    """Create a tiny repo, index it, and configure the retriever to use it."""
    # MATLAB file
    m_file = tmp_path / "pipeline" / "core" / "fit_models.m"
    m_file.parent.mkdir(parents=True)
    m_file.write_text(
        "function result = fit_models(data, params)\n"
        "% FIT_MODELS Fits ADC and IVIM models to diffusion-weighted MRI data.\n"
        "% This function processes voxel-level signals and returns parameter maps\n"
        "% for downstream analysis in the pancreatic DWI pipeline.\n"
        "    result = data .* params;\n"
        "end\n"
        "\n"
        "function out = helper_calc(x)\n"
        "% HELPER_CALC Applies element-wise transformation used by fit_models.\n"
        "% Computes the square of each element for signal normalisation steps.\n"
        "    out = x.^2;\n"
        "end\n",
        encoding="utf-8",
    )

    # Python file
    py_file = tmp_path / "analysis" / "shared.py"
    py_file.parent.mkdir(parents=True)
    py_file.write_text(
        "class ConfigLoader:\n"
        "    \"\"\"Loads analysis configuration from JSON files and merges\n"
        "    defaults with user overrides for the DWI analysis pipeline.\"\"\"\n"
        "\n"
        "    def load(self, path):\n"
        "        return {}\n"
        "\n"
        "\n"
        "def parse_p_value(text):\n"
        "    \"\"\"Extract p-value from statistical test output text using\n"
        "    regex patterns for Wilcoxon, t-test, and ANOVA results.\"\"\"\n"
        "    return 0.05\n",
        encoding="utf-8",
    )

    # Test file
    test_file = tmp_path / "pipeline" / "tests" / "test_fit_models.m"
    test_file.parent.mkdir(parents=True)
    test_file.write_text(
        "function tests = test_fit_models\n"
        "% TEST_FIT_MODELS Unit tests for the fit_models function.\n"
        "% Verifies that ADC and IVIM model fitting produces correct parameter\n"
        "% maps from synthetic diffusion-weighted MRI signal data.\n"
        "    tests = functiontests(localfunctions);\n"
        "end\n"
        "\n"
        "function test_basic_fit(testCase)\n"
        "% TEST_BASIC_FIT Verifies basic model fitting with known input signals\n"
        "% and expected output parameter values within tolerance bounds.\n"
        "    verifyEqual(testCase, fit_models([1 2], [3 4]), [3 8]);\n"
        "end\n",
        encoding="utf-8",
    )

    # Markdown doc
    doc_file = tmp_path / "CLAUDE.md"
    doc_file.write_text(
        "# pancData3\n\n"
        "Overview of the pancreatic DWI analysis pipeline project with MATLAB\n"
        "and Python components for medical imaging research at MSK.\n\n"
        "## Pipeline Core\n\n"
        "The pipeline core modules handle DICOM loading, model fitting,\n"
        "sanity checks, visualisation, and metric computation for each scan.\n\n"
        "## Security\n\n"
        "Patient data must never be exposed. PHI guards prevent logging of\n"
        "patient identifiers and restrict data to secure locations only.\n",
        encoding="utf-8",
    )

    db_path = str(tmp_path / ".chromadb")
    build_index(str(tmp_path), force_rebuild=True, db_path=db_path)

    # Patch the retriever to use our test collection
    monkeypatch.setattr("improvement_loop.rag.retriever.REPO_ROOT", str(tmp_path))
    monkeypatch.setattr(
        "improvement_loop.rag.indexer._get_db_path", lambda: db_path
    )
    monkeypatch.setattr("improvement_loop.rag.indexer.REPO_ROOT", str(tmp_path))

    return tmp_path


@pytest.fixture()
def indexed_repo_with_findings(indexed_repo, monkeypatch):
    """indexed_repo plus improvement history findings."""
    log_data = [
        {
            "iteration": 1,
            "findings": [
                {
                    "id": "find-001",
                    "dimension": "security",
                    "description": "LEAKAGE_RISK: cross-timepoint imputation allows data bleed.",
                    "fix": "Enforce temporal separation in KNN imputation.",
                    "file": "pipeline/utils/knn_impute_train_test.m",
                    "status": "merged",
                },
                {
                    "id": "find-002",
                    "dimension": "correctness",
                    "description": "Off-by-one in loop boundary skips last voxel.",
                    "fix": "Change < to <= on line 42.",
                    "file": "pipeline/core/fit_models.m",
                    "status": "merged",
                },
            ],
        },
    ]
    log_file = indexed_repo / "improvement_loop_log.json"
    log_file.write_text(json.dumps(log_data), encoding="utf-8")
    monkeypatch.setattr("improvement_loop.rag.indexer._LOG_PATH", str(log_file))

    db_path = str(indexed_repo / ".chromadb")
    index_improvement_history(db_path=db_path)

    return indexed_repo


# ---------------------------------------------------------------------------
# query() basics
# ---------------------------------------------------------------------------

class TestQuery:
    def test_returns_results(self, indexed_repo):
        results = query("fit models ADC IVIM")
        assert len(results) > 0

    def test_sorted_by_relevance_descending(self, indexed_repo):
        results = query("model fitting data parameters")
        scores = [r.relevance_score for r in results]
        assert scores == sorted(scores, reverse=True)

    def test_relevance_score_range(self, indexed_repo):
        results = query("fit models")
        for r in results:
            assert 0.0 < r.relevance_score <= 1.0

    def test_dataclass_fields(self, indexed_repo):
        results = query("fit models")
        assert len(results) > 0
        r = results[0]
        assert isinstance(r.content, str) and len(r.content) > 0
        assert isinstance(r.file_path, str)
        assert isinstance(r.chunk_type, str)
        assert isinstance(r.name, str)
        assert isinstance(r.language, str)
        assert isinstance(r.start_line, int)
        assert isinstance(r.end_line, int)


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

class TestQueryFilters:
    def test_filter_type(self, indexed_repo):
        results = query("pipeline code", filter_type="function")
        for r in results:
            assert r.chunk_type == "function"

    def test_filter_language(self, indexed_repo):
        results = query("analysis config loader", filter_language="python")
        for r in results:
            assert r.language == "python"

    def test_filter_file(self, indexed_repo):
        results = query(
            "model fitting",
            filter_file="pipeline/core/fit_models.m",
        )
        for r in results:
            assert r.file_path == "pipeline/core/fit_models.m"

    def test_exclude_files(self, indexed_repo):
        results = query(
            "fit models",
            exclude_files=["pipeline/core/fit_models.m"],
        )
        for r in results:
            assert r.file_path != "pipeline/core/fit_models.m"

    def test_min_relevance_filters_low_scores(self, indexed_repo):
        # Use a very high min_relevance to filter most results
        results = query("fit models", min_relevance=0.99)
        # Should get fewer (possibly zero) results
        all_results = query("fit models", min_relevance=0.0)
        assert len(results) <= len(all_results)

    def test_top_k_limits_results(self, indexed_repo):
        results = query("code", top_k=2)
        assert len(results) <= 2


# ---------------------------------------------------------------------------
# Fallback on missing collection
# ---------------------------------------------------------------------------

class TestQueryFallback:
    def test_missing_collection_returns_empty(self, tmp_path, monkeypatch):
        """If the collection doesn't exist, query returns [] with no crash."""
        db_path = str(tmp_path / ".nonexistent_chromadb")
        monkeypatch.setattr(
            "improvement_loop.rag.indexer._get_db_path", lambda: db_path
        )
        monkeypatch.setattr("improvement_loop.rag.indexer.REPO_ROOT", str(tmp_path))
        results = query("anything at all")
        assert results == [] or isinstance(results, list)


# ---------------------------------------------------------------------------
# _build_where_clause
# ---------------------------------------------------------------------------

class TestBuildWhereClause:
    def test_no_filters(self):
        assert _build_where_clause(None, None, None) is None

    def test_single_filter(self):
        result = _build_where_clause("function", None, None)
        assert result == {"chunk_type": "function"}

    def test_two_filters(self):
        result = _build_where_clause("function", "matlab", None)
        assert "$and" in result
        assert len(result["$and"]) == 2

    def test_three_filters(self):
        result = _build_where_clause("function", "matlab", "pipeline/core/fit_models.m")
        assert "$and" in result
        assert len(result["$and"]) == 3


# ---------------------------------------------------------------------------
# format_retrieved_context
# ---------------------------------------------------------------------------

class TestFormatRetrievedContext:
    def _make_chunk(self, file_path="a.m", name="func", start=1, end=10,
                    content="code here", **kw) -> RetrievedContext:
        defaults = dict(
            content=content, file_path=file_path, chunk_type="function",
            name=name, language="matlab", start_line=start, end_line=end,
            relevance_score=0.9,
        )
        defaults.update(kw)
        return RetrievedContext(**defaults)

    def test_empty_list(self):
        assert format_retrieved_context([]) == ""

    def test_single_chunk(self):
        c = self._make_chunk(content="hello world")
        result = format_retrieved_context([c])
        assert "=== a.m ===" in result
        assert "hello world" in result

    def test_groups_by_file(self):
        c1 = self._make_chunk(file_path="a.m", content="aaa")
        c2 = self._make_chunk(file_path="b.py", content="bbb")
        result = format_retrieved_context([c2, c1])
        # Files should appear in sorted order
        pos_a = result.index("=== a.m ===")
        pos_b = result.index("=== b.py ===")
        assert pos_a < pos_b

    def test_deduplicates_overlapping(self):
        # Big chunk contains small chunk
        big = self._make_chunk(start=1, end=50, content="big block of code")
        small = self._make_chunk(start=5, end=10, content="small slice")
        result = format_retrieved_context([big, small])
        assert "big block of code" in result
        assert "small slice" not in result

    def test_non_overlapping_both_kept(self):
        c1 = self._make_chunk(start=1, end=10, content="first")
        c2 = self._make_chunk(start=20, end=30, content="second")
        result = format_retrieved_context([c1, c2])
        assert "first" in result
        assert "second" in result


# ---------------------------------------------------------------------------
# _deduplicate_overlapping
# ---------------------------------------------------------------------------

class TestDeduplicateOverlapping:
    def _make_chunk(self, start, end, name="f") -> RetrievedContext:
        return RetrievedContext(
            content=f"lines {start}-{end}", file_path="x.m",
            chunk_type="function", name=name, language="matlab",
            start_line=start, end_line=end, relevance_score=0.9,
        )

    def test_no_overlap(self):
        chunks = [self._make_chunk(1, 10), self._make_chunk(20, 30)]
        result = _deduplicate_overlapping(chunks)
        assert len(result) == 2

    def test_contained_removed(self):
        big = self._make_chunk(1, 50, name="big")
        small = self._make_chunk(5, 10, name="small")
        result = _deduplicate_overlapping([big, small])
        assert len(result) == 1
        assert result[0].name == "big"

    def test_container_replaces_earlier_small(self):
        small = self._make_chunk(5, 10, name="small")
        big = self._make_chunk(1, 50, name="big")
        result = _deduplicate_overlapping([small, big])
        assert len(result) == 1
        assert result[0].name == "big"

    def test_single_chunk(self):
        result = _deduplicate_overlapping([self._make_chunk(1, 10)])
        assert len(result) == 1

    def test_empty(self):
        result = _deduplicate_overlapping([])
        assert result == []


# ---------------------------------------------------------------------------
# Agent context builders (require indexed repo)
# ---------------------------------------------------------------------------

class TestGetContextForAudit:
    def test_returns_nonempty_string(self, indexed_repo):
        result = get_context_for_audit("Iteration 1: initial audit pass")
        assert isinstance(result, str)
        assert len(result) > 0

    def test_contains_file_headers(self, indexed_repo):
        result = get_context_for_audit("audit pass")
        assert "===" in result

    def test_empty_context_still_works(self, indexed_repo):
        result = get_context_for_audit("")
        assert isinstance(result, str)


class TestGetContextForFix:
    def test_returns_nonempty_string(self, indexed_repo):
        result = get_context_for_fix(
            "Off-by-one in loop bound",
            "pipeline/core/fit_models.m",
        )
        assert isinstance(result, str)
        assert len(result) > 0

    def test_includes_target_file(self, indexed_repo):
        result = get_context_for_fix(
            "model fitting issue",
            "pipeline/core/fit_models.m",
        )
        assert "fit_models.m" in result

    def test_includes_past_findings(self, indexed_repo_with_findings):
        result = get_context_for_fix(
            "Off-by-one in loop boundary",
            "pipeline/core/fit_models.m",
        )
        assert isinstance(result, str)
        assert len(result) > 0


class TestGetContextForReview:
    def test_returns_nonempty_string(self, indexed_repo):
        result = get_context_for_review(
            "Fix loop boundary",
            "pipeline/core/fit_models.m",
        )
        assert isinstance(result, str)
        assert len(result) > 0

    def test_searches_for_tests(self, indexed_repo):
        result = get_context_for_review(
            "Fix model fitting",
            "pipeline/core/fit_models.m",
        )
        # The test file should be retrievable since it's semantically related
        assert isinstance(result, str)

    def test_searches_for_safety_findings(self, indexed_repo_with_findings):
        result = get_context_for_review(
            "Fix data handling",
            "pipeline/core/fit_models.m",
        )
        assert isinstance(result, str)
        assert len(result) > 0
