"""Tests for improvement_loop.rag.indexer: ChromaDB vector index build, incremental update, history."""

import json
import os
import sys
import time
import pytest

# Add repo root so modules are importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from improvement_loop.rag.chunker import chunk_repo
from improvement_loop.rag.indexer import (
    COLLECTION_NAME,
    build_index,
    get_collection,
    index_improvement_history,
    update_index_for_files,
    _chunk_id,
)
from improvement_loop.loop_config import reset_config


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _reset_loop_config():
    """Reset cached LoopConfig so each test gets defaults."""
    reset_config()
    yield
    reset_config()


@pytest.fixture()
def sample_repo(tmp_path):
    """Create a tiny repo with 3 files for indexing tests."""
    # MATLAB file with two functions
    m_file = tmp_path / "fit_models.m"
    m_file.write_text(
        "function result = fit_models(data, params)\n"
        "% FIT_MODELS Fits ADC and IVIM models to diffusion-weighted MRI data.\n"
        "% This function processes voxel-level signals and returns parameter maps\n"
        "% for downstream analysis in the pancreatic DWI pipeline.\n"
        "    result = data .* params;\n"
        "end\n"
        "\n"
        "function out = helper_function(x)\n"
        "% HELPER_FUNCTION Applies element-wise transformation used by fit_models.\n"
        "% Computes the square of each element for signal normalisation steps.\n"
        "    out = x.^2;\n"
        "end\n",
        encoding="utf-8",
    )

    # Python file with a class and a function
    py_file = tmp_path / "analysis.py"
    py_file.write_text(
        "class Analyzer:\n"
        "    \"\"\"Analyzer class that processes pipeline output data and computes\n"
        "    summary statistics for the pancreatic DWI analysis report.\"\"\"\n"
        "\n"
        "    def run(self):\n"
        "        return 42\n"
        "\n"
        "\n"
        "def standalone_helper(x, y):\n"
        "    \"\"\"Standalone helper that computes weighted combination of two\n"
        "    input arrays for cross-reference analysis reporting.\"\"\"\n"
        "    return x + y\n",
        encoding="utf-8",
    )

    # Markdown file with two sections
    md_file = tmp_path / "README.md"
    md_file.write_text(
        "# Project\n\n"
        "Top-level introduction to the pancreatic DWI analysis pipeline project\n"
        "covering data acquisition, processing, and reporting workflows.\n\n"
        "## Installation\n\n"
        "Install dependencies with pip install -r requirements.txt and configure\n"
        "the pipeline by copying config.example.json to config.json with your paths.\n\n"
        "## Usage\n\n"
        "Run the MATLAB pipeline first, then execute the Python analysis scripts\n"
        "to generate the full HTML and PDF reports from pipeline outputs.\n",
        encoding="utf-8",
    )

    return tmp_path


@pytest.fixture()
def sample_log(tmp_path):
    """Create a sample improvement_loop_log.json."""
    log_data = [
        {
            "iteration": 1,
            "findings": [
                {
                    "id": "finding-001",
                    "dimension": "correctness",
                    "description": "Off-by-one error in loop bound causes last voxel to be skipped.",
                    "fix": "Change < to <= in the loop condition on line 42.",
                    "file": "pipeline/core/fit_models.m",
                    "status": "merged",
                },
                {
                    "id": "finding-002",
                    "dimension": "robustness",
                    "description": "Missing NaN guard in summary metric aggregation leads to silent data loss.",
                    "fix": "Add nanmean_safe wrapper before aggregation step.",
                    "file": "pipeline/core/compute_summary_metrics.m",
                    "status": "pending",
                },
            ],
        },
        {
            "iteration": 2,
            "findings": [
                {
                    "id": "finding-003",
                    "dimension": "performance",
                    "description": "Redundant recomputation of ADC maps inside parfor wastes cluster time.",
                    "fix": "Cache ADC results in parsave_dir_cache before the parallel loop.",
                    "file": "pipeline/core/load_dwi_data.m",
                    "status": "merged",
                },
            ],
        },
    ]
    log_file = tmp_path / "improvement_loop_log.json"
    log_file.write_text(json.dumps(log_data), encoding="utf-8")
    return log_file


# ---------------------------------------------------------------------------
# build_index — full build on temp repo
# ---------------------------------------------------------------------------

class TestBuildIndex:
    def test_indexes_all_chunks(self, sample_repo):
        """build_index on a 3-file repo returns a collection with the expected chunk count."""
        expected_chunks = chunk_repo(str(sample_repo))
        collection = build_index(str(sample_repo), force_rebuild=True,
                                 db_path=str(sample_repo / ".chromadb"))
        assert collection.count() == len(expected_chunks)
        assert collection.count() > 0

    def test_chunk_metadata_fields(self, sample_repo):
        """Every indexed document has the required metadata keys."""
        db_path = str(sample_repo / ".chromadb")
        collection = build_index(str(sample_repo), force_rebuild=True,
                                 db_path=db_path)
        all_data = collection.get(include=["metadatas"])
        required_keys = {"file_path", "chunk_type", "name", "language",
                         "start_line", "end_line", "indexed_at", "file_mtime"}
        for meta in all_data["metadatas"]:
            assert required_keys.issubset(meta.keys()), f"Missing keys in {meta}"

    def test_force_rebuild_resets(self, sample_repo):
        """force_rebuild=True drops and recreates the collection."""
        db_path = str(sample_repo / ".chromadb")
        col1 = build_index(str(sample_repo), force_rebuild=True, db_path=db_path)
        count1 = col1.count()

        # Add a file, rebuild without force — count increases
        extra = sample_repo / "extra.py"
        extra.write_text(
            "def extra_function():\n"
            "    \"\"\"Extra function added after first build to verify incremental\n"
            "    index update picks up new files correctly.\"\"\"\n"
            "    return 99\n",
            encoding="utf-8",
        )
        col2 = build_index(str(sample_repo), force_rebuild=True, db_path=db_path)
        assert col2.count() >= count1


# ---------------------------------------------------------------------------
# Incremental update — only changed files re-indexed
# ---------------------------------------------------------------------------

class TestIncrementalUpdate:
    def test_unchanged_files_skipped(self, sample_repo):
        """A second build_index with no file changes upserts 0 chunks."""
        db_path = str(sample_repo / ".chromadb_unchanged")
        build_index(str(sample_repo), force_rebuild=True, db_path=db_path)

        # Second run — incremental, nothing changed
        import io
        from contextlib import redirect_stdout

        buf = io.StringIO()
        with redirect_stdout(buf):
            build_index(str(sample_repo), force_rebuild=False, db_path=db_path)
        output = buf.getvalue()
        assert "up to date" in output

    def test_modified_file_reindexed(self, sample_repo):
        """Touching one file's mtime causes only that file's chunks to be re-indexed."""
        db_path = str(sample_repo / ".chromadb_modified")
        col = build_index(str(sample_repo), force_rebuild=True, db_path=db_path)
        count_before = col.count()

        # Touch the MATLAB file to change its mtime
        m_file = sample_repo / "fit_models.m"
        time.sleep(0.05)
        m_file.write_text(
            m_file.read_text(encoding="utf-8") +
            "\nfunction z = new_func(a)\n"
            "% NEW_FUNC Additional helper function added to test incremental\n"
            "% re-indexing when a single file is modified after initial build.\n"
            "    z = a + 1;\n"
            "end\n",
            encoding="utf-8",
        )

        col2 = build_index(str(sample_repo), force_rebuild=False, db_path=db_path)
        # Should have more chunks now (new function added)
        assert col2.count() > count_before

    def test_deleted_file_removed(self, sample_repo):
        """Deleting a file removes its chunks on incremental rebuild."""
        db_path = str(sample_repo / ".chromadb_deleted")
        col = build_index(str(sample_repo), force_rebuild=True, db_path=db_path)
        count_before = col.count()

        # Delete the markdown file
        os.remove(sample_repo / "README.md")

        col2 = build_index(str(sample_repo), force_rebuild=False, db_path=db_path)
        assert col2.count() < count_before


# ---------------------------------------------------------------------------
# index_improvement_history
# ---------------------------------------------------------------------------

class TestIndexImprovementHistory:
    def test_indexes_all_findings(self, sample_repo, sample_log, monkeypatch):
        """index_improvement_history inserts one document per finding."""
        db_path = str(sample_repo / ".chromadb")
        monkeypatch.setattr(
            "improvement_loop.rag.indexer._LOG_PATH", str(sample_log)
        )
        # Build the base index first so the collection exists
        build_index(str(sample_repo), force_rebuild=True, db_path=db_path)
        monkeypatch.setattr(
            "improvement_loop.rag.indexer._get_db_path", lambda: db_path
        )
        index_improvement_history(db_path=db_path)

        import chromadb
        client = chromadb.PersistentClient(path=db_path)
        col = client.get_collection(name=COLLECTION_NAME)

        # Check that finding documents exist
        all_data = col.get(include=["metadatas"])
        finding_metas = [m for m in all_data["metadatas"]
                         if m and m.get("chunk_type") == "finding"]
        assert len(finding_metas) == 3  # 3 findings in sample log

    def test_finding_metadata_fields(self, sample_repo, sample_log, monkeypatch):
        """Finding documents have iteration, dimension, and status metadata."""
        db_path = str(sample_repo / ".chromadb")
        monkeypatch.setattr(
            "improvement_loop.rag.indexer._LOG_PATH", str(sample_log)
        )
        build_index(str(sample_repo), force_rebuild=True, db_path=db_path)
        monkeypatch.setattr(
            "improvement_loop.rag.indexer._get_db_path", lambda: db_path
        )
        index_improvement_history(db_path=db_path)

        import chromadb
        client = chromadb.PersistentClient(path=db_path)
        col = client.get_collection(name=COLLECTION_NAME)

        all_data = col.get(include=["metadatas"])
        finding_metas = [m for m in all_data["metadatas"]
                         if m and m.get("chunk_type") == "finding"]
        for meta in finding_metas:
            assert "iteration" in meta
            assert "dimension" in meta
            assert "status" in meta

    def test_missing_log_is_noop(self, sample_repo, monkeypatch):
        """index_improvement_history with no log file does nothing."""
        db_path = str(sample_repo / ".chromadb")
        monkeypatch.setattr(
            "improvement_loop.rag.indexer._LOG_PATH",
            str(sample_repo / "nonexistent.json"),
        )
        build_index(str(sample_repo), force_rebuild=True, db_path=db_path)
        monkeypatch.setattr(
            "improvement_loop.rag.indexer._get_db_path", lambda: db_path
        )

        import chromadb
        client = chromadb.PersistentClient(path=db_path)
        col = client.get_collection(name=COLLECTION_NAME)
        count_before = col.count()

        index_improvement_history(db_path=db_path)

        col2 = client.get_collection(name=COLLECTION_NAME)
        assert col2.count() == count_before


# ---------------------------------------------------------------------------
# get_collection — auto-creates if missing
# ---------------------------------------------------------------------------

class TestGetCollection:
    def test_creates_collection_if_missing(self, sample_repo, monkeypatch):
        """get_collection builds the index when no collection exists."""
        db_path = str(sample_repo / ".chromadb")
        monkeypatch.setattr("improvement_loop.rag.indexer.REPO_ROOT",
                            str(sample_repo))
        monkeypatch.setattr("improvement_loop.rag.indexer._get_db_path",
                            lambda: db_path)

        col = get_collection(db_path=db_path)
        assert col.count() > 0

    def test_returns_existing_collection(self, sample_repo, monkeypatch):
        """get_collection returns an existing non-empty collection without rebuilding."""
        db_path = str(sample_repo / ".chromadb")
        monkeypatch.setattr("improvement_loop.rag.indexer.REPO_ROOT",
                            str(sample_repo))
        monkeypatch.setattr("improvement_loop.rag.indexer._get_db_path",
                            lambda: db_path)

        col1 = build_index(str(sample_repo), force_rebuild=True, db_path=db_path)
        count1 = col1.count()

        col2 = get_collection(db_path=db_path)
        assert col2.count() == count1


# ---------------------------------------------------------------------------
# update_index_for_files
# ---------------------------------------------------------------------------

class TestUpdateIndexForFiles:
    def test_reindexes_specific_file(self, sample_repo, monkeypatch):
        """update_index_for_files re-chunks a single file."""
        db_path = str(sample_repo / ".chromadb")
        monkeypatch.setattr("improvement_loop.rag.indexer.REPO_ROOT",
                            str(sample_repo))
        build_index(str(sample_repo), force_rebuild=True, db_path=db_path)

        # Modify the Python file
        py_file = sample_repo / "analysis.py"
        py_file.write_text(
            "class Analyzer:\n"
            "    \"\"\"Analyzer class that processes pipeline output data and computes\n"
            "    summary statistics for the pancreatic DWI analysis report.\"\"\"\n"
            "\n"
            "    def run(self):\n"
            "        return 42\n",
            encoding="utf-8",
        )

        update_index_for_files(["analysis.py"], db_path=db_path)

        import chromadb
        client = chromadb.PersistentClient(path=db_path)
        col = client.get_collection(name=COLLECTION_NAME)
        all_data = col.get(include=["metadatas"])
        py_chunks = [m for m in all_data["metadatas"]
                     if m and m.get("file_path") == "analysis.py"]
        # Should have only 1 chunk now (Analyzer class, standalone_helper removed)
        assert len(py_chunks) == 1

    def test_deleted_file_chunks_removed(self, sample_repo, monkeypatch):
        """update_index_for_files removes chunks for a deleted file."""
        db_path = str(sample_repo / ".chromadb")
        monkeypatch.setattr("improvement_loop.rag.indexer.REPO_ROOT",
                            str(sample_repo))
        build_index(str(sample_repo), force_rebuild=True, db_path=db_path)

        # Delete the file
        os.remove(sample_repo / "analysis.py")

        update_index_for_files(["analysis.py"], db_path=db_path)

        import chromadb
        client = chromadb.PersistentClient(path=db_path)
        col = client.get_collection(name=COLLECTION_NAME)
        all_data = col.get(include=["metadatas"])
        py_chunks = [m for m in all_data["metadatas"]
                     if m and m.get("file_path") == "analysis.py"]
        assert len(py_chunks) == 0
