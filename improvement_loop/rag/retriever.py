# retriever.py
"""RAG retriever — queries the ChromaDB vector store for relevant code context.

Provides targeted query functions for each agent in the improvement loop
(auditor, implementer, reviewer) and a general-purpose query interface.
"""

import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional

from improvement_loop.loop_config import get_config as _get_loop_config

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


@dataclass
class RetrievedContext:
    """A chunk retrieved from the vector store with relevance score."""
    content: str
    file_path: str
    chunk_type: str
    name: str
    language: str
    start_line: int
    end_line: int
    relevance_score: float


# ---------------------------------------------------------------------------
# Core query
# ---------------------------------------------------------------------------

def query(
    query_text: str,
    top_k: Optional[int] = None,
    min_relevance: Optional[float] = None,
    filter_type: Optional[str] = None,
    filter_language: Optional[str] = None,
    filter_file: Optional[str] = None,
    exclude_files: Optional[List[str]] = None,
) -> List[RetrievedContext]:
    """Query the vector store and return relevant chunks.

    Parameters
    ----------
    query_text : str
        Natural-language query.
    top_k : int, optional
        Max results (defaults to config ``rag_top_k``).
    min_relevance : float, optional
        Minimum similarity score (defaults to config ``rag_min_relevance``).
        ChromaDB distances are L2; we convert to a 0-1 relevance score.
    filter_type : str, optional
        Restrict to a specific ``chunk_type`` (e.g. "function", "doc", "finding").
    filter_language : str, optional
        Restrict to a specific language (e.g. "matlab", "python").
    filter_file : str, optional
        Restrict to chunks from a single file path.
    exclude_files : list of str, optional
        File paths to exclude from results.

    Returns
    -------
    list of RetrievedContext
        Sorted by relevance_score descending.
    """
    cfg = _get_loop_config()
    if top_k is None:
        top_k = cfg.rag_top_k
    if min_relevance is None:
        min_relevance = cfg.rag_min_relevance

    try:
        from improvement_loop.rag.indexer import get_collection
        collection = get_collection()
    except Exception as e:
        print(f"⚠️  RAG query failed (collection unavailable): {e}")
        return []

    # Build ChromaDB where clause
    where = _build_where_clause(filter_type, filter_language, filter_file)

    try:
        results = collection.query(
            query_texts=[query_text],
            n_results=top_k,
            include=["documents", "metadatas", "distances"],
            where=where if where else None,
        )
    except Exception as e:
        print(f"⚠️  RAG query failed: {e}")
        return []

    if not results or not results["ids"] or not results["ids"][0]:
        return []

    chunks: List[RetrievedContext] = []
    ids = results["ids"][0]
    docs = results["documents"][0] if results["documents"] else [""] * len(ids)
    metas = results["metadatas"][0] if results["metadatas"] else [{}] * len(ids)
    dists = results["distances"][0] if results["distances"] else [0.0] * len(ids)

    exclude_set = set(exclude_files) if exclude_files else set()

    for doc, meta, dist in zip(docs, metas, dists):
        meta = meta or {}
        file_path = meta.get("file_path", "")

        if file_path in exclude_set:
            continue

        # Convert L2 distance to a 0-1 relevance score.
        # ChromaDB default embedding distances are non-negative; smaller = more similar.
        # We use 1/(1+d) which maps [0,∞) → (0,1].
        relevance = 1.0 / (1.0 + dist)

        if relevance < min_relevance:
            continue

        chunks.append(RetrievedContext(
            content=doc or "",
            file_path=file_path,
            chunk_type=meta.get("chunk_type", ""),
            name=meta.get("name", ""),
            language=meta.get("language", ""),
            start_line=meta.get("start_line", 0),
            end_line=meta.get("end_line", 0),
            relevance_score=relevance,
        ))

    # Sort by relevance descending
    chunks.sort(key=lambda c: c.relevance_score, reverse=True)
    return chunks


def _build_where_clause(
    filter_type: Optional[str],
    filter_language: Optional[str],
    filter_file: Optional[str],
) -> Optional[dict]:
    """Build a ChromaDB ``where`` filter dict from optional constraints."""
    conditions: List[dict] = []

    if filter_type:
        conditions.append({"chunk_type": filter_type})
    if filter_language:
        conditions.append({"language": filter_language})
    if filter_file:
        conditions.append({"file_path": filter_file})

    if not conditions:
        return None
    if len(conditions) == 1:
        return conditions[0]
    return {"$and": conditions}


# ---------------------------------------------------------------------------
# Format helpers
# ---------------------------------------------------------------------------

def format_retrieved_context(chunks: List[RetrievedContext]) -> str:
    """Format retrieved chunks into ``=== file ===`` block format.

    Groups by file, deduplicates overlapping line ranges,
    sorts by file path then start_line.
    """
    if not chunks:
        return ""

    # Group by file
    by_file: Dict[str, List[RetrievedContext]] = {}
    for chunk in chunks:
        by_file.setdefault(chunk.file_path, []).append(chunk)

    parts: List[str] = []
    for file_path in sorted(by_file):
        file_chunks = by_file[file_path]
        # Sort by start_line
        file_chunks.sort(key=lambda c: c.start_line)

        # Deduplicate overlapping ranges
        merged = _deduplicate_overlapping(file_chunks)

        content_block = "\n".join(c.content for c in merged)
        parts.append(f"=== {file_path} ===\n{content_block}")

    return "\n\n".join(parts)


def _deduplicate_overlapping(chunks: List[RetrievedContext]) -> List[RetrievedContext]:
    """Remove chunks whose line range is fully contained within another chunk."""
    if len(chunks) <= 1:
        return chunks

    result: List[RetrievedContext] = []
    for chunk in chunks:
        # Check if this chunk is fully contained by any already-accepted chunk
        contained = False
        for existing in result:
            if (existing.file_path == chunk.file_path
                    and existing.start_line <= chunk.start_line
                    and existing.end_line >= chunk.end_line):
                contained = True
                break
        if not contained:
            # Remove any previously accepted chunks that are now contained by this one
            result = [
                r for r in result
                if not (r.file_path == chunk.file_path
                        and chunk.start_line <= r.start_line
                        and chunk.end_line >= r.end_line)
            ]
            result.append(chunk)

    return result


# ---------------------------------------------------------------------------
# Agent-specific context builders
# ---------------------------------------------------------------------------

def get_context_for_audit(iteration_context: str) -> str:
    """Build an enriched context string for the auditor agent.

    Runs multiple targeted queries to surface diverse code areas:
    1. Code quality / correctness issues
    2. Performance / memory hot paths
    3. Security / PHI-sensitive code
    4. Test coverage gaps
    5. Iteration-context-driven search for files not yet audited

    Deduplicates results across queries and formats as the ``=== file ===``
    block format the auditor system prompt expects.
    """
    queries = [
        ("code quality bugs errors edge cases validation", None),
        ("performance memory optimization parfor parallel", None),
        ("security injection PHI patient data system shell", None),
        ("test coverage gaps missing assertions", None),
    ]

    # Add iteration-context-driven query if non-empty
    if iteration_context and iteration_context.strip():
        queries.append((iteration_context[:500], None))

    all_chunks: List[RetrievedContext] = []
    seen_keys: set = set()

    for query_text, filt in queries:
        results = query(query_text, top_k=8, filter_type=filt)
        for chunk in results:
            key = (chunk.file_path, chunk.name)
            if key not in seen_keys:
                seen_keys.add(key)
                all_chunks.append(chunk)

    # Sort by relevance so the most relevant chunks appear first
    all_chunks.sort(key=lambda c: c.relevance_score, reverse=True)

    # Cap total to avoid exceeding token limits
    cfg = _get_loop_config()
    max_chars = cfg.max_file_chars * 3
    selected: List[RetrievedContext] = []
    total_chars = 0
    for chunk in all_chunks:
        if total_chars + len(chunk.content) > max_chars:
            break
        selected.append(chunk)
        total_chars += len(chunk.content)

    return format_retrieved_context(selected)


def get_context_for_fix(finding_description: str, file_path: str) -> str:
    """Build context for the implementer agent.

    Queries for:
    1. The target file to get full function context
    2. Functions that call or are called by the target function
    3. Similar past findings that were successfully merged
    4. Relevant documentation sections
    """
    all_chunks: List[RetrievedContext] = []
    seen_keys: set = set()

    def _collect(results: List[RetrievedContext]) -> None:
        for chunk in results:
            key = (chunk.file_path, chunk.name)
            if key not in seen_keys:
                seen_keys.add(key)
                all_chunks.append(chunk)

    # 1. Target file's own chunks
    _collect(query(finding_description, top_k=5, filter_file=file_path))

    # 2. Related functions (callers/callees) via semantic search
    basename = os.path.basename(file_path)
    func_query = f"{basename} calls imports uses references"
    _collect(query(func_query, top_k=5, exclude_files=[file_path]))

    # 3. Past merged findings for similar issues
    _collect(query(finding_description, top_k=3, filter_type="finding"))

    # 4. Documentation
    _collect(query(finding_description, top_k=2, filter_language="markdown"))

    # Cap total context
    cfg = _get_loop_config()
    max_chars = cfg.max_file_chars * 2
    selected: List[RetrievedContext] = []
    total_chars = 0
    for chunk in all_chunks:
        if total_chars + len(chunk.content) > max_chars:
            break
        selected.append(chunk)
        total_chars += len(chunk.content)

    return format_retrieved_context(selected)


def get_context_for_review(finding_description: str, file_path: str) -> str:
    """Build context for the reviewer agent.

    Queries for:
    1. The target file's other functions (side-effect awareness)
    2. Test files that cover the target file
    3. Documentation about the module's expected behavior
    4. Past findings flagged as LEAKAGE_RISK or PHI_RISK (pattern matching)
    """
    all_chunks: List[RetrievedContext] = []
    seen_keys: set = set()

    def _collect(results: List[RetrievedContext]) -> None:
        for chunk in results:
            key = (chunk.file_path, chunk.name)
            if key not in seen_keys:
                seen_keys.add(key)
                all_chunks.append(chunk)

    # 1. Other functions in the same file
    _collect(query(
        f"functions in {os.path.basename(file_path)}",
        top_k=5, filter_file=file_path,
    ))

    # 2. Test files for the target module
    module_name = os.path.splitext(os.path.basename(file_path))[0]
    _collect(query(
        f"test {module_name} verify assert",
        top_k=5, filter_type="function",
    ))

    # 3. Documentation about the module
    _collect(query(
        f"{module_name} documentation behavior specification",
        top_k=2, filter_language="markdown",
    ))

    # 4. Past findings with safety flags
    _collect(query(
        "LEAKAGE_RISK PHI_RISK data leakage patient safety",
        top_k=3, filter_type="finding",
    ))

    # Cap total context
    cfg = _get_loop_config()
    max_chars = cfg.max_file_chars * 2
    selected: List[RetrievedContext] = []
    total_chars = 0
    for chunk in all_chunks:
        if total_chars + len(chunk.content) > max_chars:
            break
        selected.append(chunk)
        total_chars += len(chunk.content)

    return format_retrieved_context(selected)
