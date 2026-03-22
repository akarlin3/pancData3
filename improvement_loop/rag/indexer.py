# indexer.py
"""ChromaDB-based vector index for RAG over the pancData3 codebase.

Indexes semantic code chunks (from chunker.py) and improvement history
(from improvement_loop_log.json) into a persistent ChromaDB collection
for retrieval-augmented context during audits and fixes.
"""

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from typing import List, Optional

import chromadb  # type: ignore

from improvement_loop.rag.chunker import CodeChunk, chunk_file, chunk_repo
from improvement_loop.loop_config import get_config as _get_loop_config

# Repo root is two levels up from this file's directory
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

COLLECTION_NAME = "pancdata3_codebase"
DB_PATH = os.path.join(REPO_ROOT, ".chromadb")

# Path to the improvement loop log
_LOG_PATH = os.path.join(REPO_ROOT, "improvement_loop_log.json")


def _get_db_path() -> str:
    """Return the ChromaDB path, preferring config override if set."""
    cfg = _get_loop_config()
    if cfg.rag_db_path:
        return cfg.rag_db_path
    return DB_PATH


def _get_client(db_path: Optional[str] = None) -> chromadb.ClientAPI:
    """Return a persistent ChromaDB client."""
    path = db_path or _get_db_path()
    os.makedirs(path, exist_ok=True)
    return chromadb.PersistentClient(path=path)


def _chunk_id(chunk: CodeChunk) -> str:
    """Deterministic ID for a code chunk."""
    return f"{chunk.file_path}::{chunk.name}::L{chunk.start_line}"


def _chunk_metadata(chunk: CodeChunk, file_mtime: float) -> dict:
    """Build the metadata dict stored alongside each document."""
    return {
        "file_path": chunk.file_path,
        "chunk_type": chunk.chunk_type,
        "name": chunk.name,
        "language": chunk.language,
        "start_line": chunk.start_line,
        "end_line": chunk.end_line,
        "indexed_at": datetime.now(timezone.utc).isoformat(),
        "file_mtime": file_mtime,
    }


def build_index(repo_root: str, force_rebuild: bool = False,
                db_path: Optional[str] = None) -> "chromadb.Collection":
    """Build or update the ChromaDB collection from the repository.

    If force_rebuild is True, drops and recreates the collection.
    Otherwise, does an incremental update: only re-indexes files whose
    mtime is newer than the stored metadata timestamp.
    """
    client = _get_client(db_path)

    if force_rebuild:
        try:
            client.delete_collection(COLLECTION_NAME)
        except Exception:
            pass

    collection = client.get_or_create_collection(name=COLLECTION_NAME)

    # Chunk the entire repo
    all_chunks = chunk_repo(repo_root)

    # Build a mapping of file_path -> mtime
    file_mtimes: dict[str, float] = {}
    for chunk in all_chunks:
        if chunk.file_path not in file_mtimes:
            abs_path = os.path.join(repo_root, chunk.file_path)
            try:
                file_mtimes[chunk.file_path] = os.path.getmtime(abs_path)
            except OSError:
                file_mtimes[chunk.file_path] = 0.0

    # Get existing IDs and their mtimes for incremental update
    existing_mtimes: dict[str, float] = {}
    if not force_rebuild and collection.count() > 0:
        # Fetch all existing entries in batches
        all_existing = collection.get(include=["metadatas"])
        if all_existing and all_existing["ids"]:
            for doc_id, meta in zip(all_existing["ids"], all_existing["metadatas"]):
                if meta:
                    existing_mtimes[doc_id] = meta.get("file_mtime", 0.0)

    # Determine which chunks need upserting
    ids_to_upsert: list[str] = []
    docs_to_upsert: list[str] = []
    metas_to_upsert: list[dict] = []
    new_ids: set[str] = set()

    for chunk in all_chunks:
        chunk_id = _chunk_id(chunk)
        new_ids.add(chunk_id)
        mtime = file_mtimes.get(chunk.file_path, 0.0)

        # Skip if already indexed with same mtime
        if not force_rebuild and chunk_id in existing_mtimes:
            if existing_mtimes[chunk_id] == mtime:
                continue

        ids_to_upsert.append(chunk_id)
        docs_to_upsert.append(chunk.content)
        metas_to_upsert.append(_chunk_metadata(chunk, mtime))

    # Upsert in batches (ChromaDB has a batch size limit)
    batch_size = 500
    for i in range(0, len(ids_to_upsert), batch_size):
        end = min(i + batch_size, len(ids_to_upsert))
        collection.upsert(
            ids=ids_to_upsert[i:end],
            documents=docs_to_upsert[i:end],
            metadatas=metas_to_upsert[i:end],
        )

    # Remove stale entries (files deleted from repo)
    if not force_rebuild and existing_mtimes:
        stale_ids = [eid for eid in existing_mtimes if eid not in new_ids
                     and not eid.startswith("finding::")]
        if stale_ids:
            for i in range(0, len(stale_ids), batch_size):
                collection.delete(ids=stale_ids[i:i + batch_size])

    upserted = len(ids_to_upsert)
    if upserted > 0:
        print(f"  Indexed {upserted} chunks ({collection.count()} total in collection)")
    else:
        print(f"  Index up to date ({collection.count()} chunks)")

    return collection


def get_collection(db_path: Optional[str] = None) -> "chromadb.Collection":
    """Return the existing collection, or build it if missing."""
    client = _get_client(db_path)
    try:
        col = client.get_collection(name=COLLECTION_NAME)
        if col.count() > 0:
            return col
    except Exception:
        pass
    # Collection doesn't exist or is empty — build it
    return build_index(REPO_ROOT, db_path=db_path)


def update_index_for_files(file_paths: List[str],
                           db_path: Optional[str] = None) -> None:
    """Re-index specific files (e.g., after a fix is applied).

    Parameters
    ----------
    file_paths : list of str
        Paths relative to the repo root.
    """
    client = _get_client(db_path)
    collection = client.get_or_create_collection(name=COLLECTION_NAME)

    for rel_path in file_paths:
        abs_path = os.path.join(REPO_ROOT, rel_path)
        if not os.path.exists(abs_path):
            # File was deleted — remove its chunks
            existing = collection.get(
                where={"file_path": rel_path},
                include=["metadatas"],
            )
            if existing and existing["ids"]:
                collection.delete(ids=existing["ids"])
            continue

        try:
            with open(abs_path, "r", encoding="utf-8", errors="replace") as f:
                content = f.read()
        except OSError:
            continue

        mtime = os.path.getmtime(abs_path)
        chunks = chunk_file(rel_path, content)

        # Remove old chunks for this file
        existing = collection.get(
            where={"file_path": rel_path},
            include=["metadatas"],
        )
        if existing and existing["ids"]:
            collection.delete(ids=existing["ids"])

        # Insert new chunks
        if chunks:
            collection.upsert(
                ids=[_chunk_id(c) for c in chunks],
                documents=[c.content for c in chunks],
                metadatas=[_chunk_metadata(c, mtime) for c in chunks],
            )


def index_improvement_history(db_path: Optional[str] = None) -> None:
    """Index past findings from improvement_loop_log.json as searchable context."""
    if not os.path.exists(_LOG_PATH):
        return

    with open(_LOG_PATH, "r", encoding="utf-8") as f:
        log = json.load(f)

    if not isinstance(log, list):
        return

    client = _get_client(db_path)
    collection = client.get_or_create_collection(name=COLLECTION_NAME)

    ids: list[str] = []
    docs: list[str] = []
    metas: list[dict] = []

    for entry in log:
        for finding in entry.get("findings", []):
            finding_id = finding.get("id", "")
            if not finding_id:
                continue

            doc_id = f"finding::{finding_id}"
            dimension = finding.get("dimension", "unknown")
            description = finding.get("description", "")
            fix = finding.get("fix", "")
            file_path = finding.get("file", "")
            status = finding.get("status", "unknown")
            iteration = finding.get("iteration", entry.get("iteration", 0))

            content = (
                f"[{dimension}] {description}\n"
                f"Fix: {fix}\n"
                f"File: {file_path}\n"
                f"Status: {status}"
            )

            ids.append(doc_id)
            docs.append(content)
            metas.append({
                "file_path": file_path,
                "chunk_type": "finding",
                "name": finding_id,
                "language": "log",
                "start_line": 0,
                "end_line": 0,
                "indexed_at": datetime.now(timezone.utc).isoformat(),
                "file_mtime": 0.0,
                "iteration": iteration,
                "dimension": dimension,
                "status": status,
            })

    if ids:
        batch_size = 500
        for i in range(0, len(ids), batch_size):
            end = min(i + batch_size, len(ids))
            collection.upsert(
                ids=ids[i:end],
                documents=docs[i:end],
                metadatas=metas[i:end],
            )
        print(f"  Indexed {len(ids)} findings from improvement history")


def _print_stats(db_path: Optional[str] = None) -> None:
    """Print collection statistics."""
    client = _get_client(db_path)
    try:
        collection = client.get_collection(name=COLLECTION_NAME)
    except Exception:
        print("No index found. Run with --force-rebuild to create one.")
        return

    total = collection.count()
    print(f"\n  Collection: {COLLECTION_NAME}")
    print(f"  Total chunks: {total}")

    if total == 0:
        return

    # Fetch all metadata for stats
    all_data = collection.get(include=["metadatas"])
    if not all_data or not all_data["metadatas"]:
        return

    by_language: dict[str, int] = {}
    by_type: dict[str, int] = {}
    latest_indexed = ""

    for meta in all_data["metadatas"]:
        if not meta:
            continue
        lang = meta.get("language", "unknown")
        ctype = meta.get("chunk_type", "unknown")
        by_language[lang] = by_language.get(lang, 0) + 1
        by_type[ctype] = by_type.get(ctype, 0) + 1
        ts = meta.get("indexed_at", "")
        if ts > latest_indexed:
            latest_indexed = ts

    print(f"\n  By language:")
    for lang in sorted(by_language):
        print(f"    {lang}: {by_language[lang]}")

    print(f"\n  By type:")
    for ctype in sorted(by_type):
        print(f"    {ctype}: {by_type[ctype]}")

    # Disk size
    path = db_path or _get_db_path()
    if os.path.isdir(path):
        total_bytes = sum(
            os.path.getsize(os.path.join(dp, fn))
            for dp, _, fns in os.walk(path)
            for fn in fns
        )
        if total_bytes < 1024 * 1024:
            print(f"\n  Disk size: {total_bytes / 1024:.1f} KB")
        else:
            print(f"\n  Disk size: {total_bytes / (1024 * 1024):.1f} MB")

    if latest_indexed:
        print(f"  Last indexed: {latest_indexed}")
    print()


if __name__ == "__main__":
    # Allow direct invocation
    _rr = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    if _rr not in sys.path:
        sys.path.insert(0, _rr)

    parser = argparse.ArgumentParser(description="RAG index builder for pancData3")
    parser.add_argument("--force-rebuild", action="store_true",
                        help="Drop and rebuild the entire index")
    parser.add_argument("--stats", action="store_true",
                        help="Print index statistics and exit")
    args = parser.parse_args()

    if args.stats:
        _print_stats()
    else:
        build_index(REPO_ROOT, force_rebuild=args.force_rebuild)
        index_improvement_history()
        if args.force_rebuild:
            _print_stats()
