# chunker.py
"""Semantic code chunker for the pancData3 repository.

Splits source files into semantically meaningful chunks (functions, classes,
sections) for use in RAG-based code retrieval.
"""

import os
import re
from dataclasses import dataclass, field
from typing import List

# Maximum / minimum chunk size in characters
MAX_CHUNK_CHARS = 4000
MIN_CHUNK_CHARS = 50

# File extensions to index, mapped to language names
_EXTENSION_LANG = {
    ".m": "matlab",
    ".py": "python",
    ".md": "markdown",
    ".json": "json",
}

# JSON files worth indexing (relative to repo root)
_INDEXABLE_JSON = {
    "config.example.json",
    "improvement_loop_config.example.json",
    "analysis_config.example.json",
    "package.json",
}

# Directories and patterns to skip entirely.
# Any path component matching one of these names is excluded, plus
# the compound "pipeline/dependencies" path.
_SKIP_DIRS = {
    ".git",
    ".claude",
    "__pycache__",
    "node_modules",
}

# Multi-segment directory paths to skip (matched as substrings of the
# normalized relative path).
_SKIP_DIR_PATHS = {
    "pipeline/dependencies",
}

# Specific files to skip
_SKIP_FILES = {
    "improvement_loop_log.json",
}

# Binary extensions to skip
_BINARY_EXTENSIONS = {
    ".mat", ".nii", ".dcm", ".png", ".jpg", ".jpeg",
    ".pdf", ".gif", ".bmp", ".tif", ".tiff", ".gz",
    ".zip", ".tar", ".exe", ".dll", ".so", ".pyc",
    ".nii.gz",
}


@dataclass
class CodeChunk:
    """One semantically meaningful piece of source code or documentation."""
    file_path: str          # relative to repo root
    chunk_type: str         # "function", "class", "method", "section", "doc", "test", "config"
    name: str               # function/class name, or section heading, or filename
    content: str            # the actual text
    start_line: int         # 1-indexed
    end_line: int
    language: str           # "matlab", "python", "markdown", "json"
    metadata: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Splitting oversized chunks
# ---------------------------------------------------------------------------

def _split_oversized(chunk: CodeChunk) -> List[CodeChunk]:
    """Split a chunk that exceeds MAX_CHUNK_CHARS at blank-line boundaries."""
    if len(chunk.content) <= MAX_CHUNK_CHARS:
        return [chunk]

    lines = chunk.content.splitlines(keepends=True)
    parts: List[CodeChunk] = []
    current_lines: list[str] = []
    current_len = 0
    part_num = 1
    part_start = chunk.start_line

    for i, line in enumerate(lines):
        # If adding this line would exceed the limit and we have content,
        # try to split at a blank line
        if current_len + len(line) > MAX_CHUNK_CHARS and current_lines:
            text = "".join(current_lines)
            if len(text.strip()) >= MIN_CHUNK_CHARS:
                parts.append(CodeChunk(
                    file_path=chunk.file_path,
                    chunk_type=chunk.chunk_type,
                    name=f"{chunk.name}_part{part_num}",
                    content=text,
                    start_line=part_start,
                    end_line=part_start + len(current_lines) - 1,
                    language=chunk.language,
                    metadata=chunk.metadata.copy(),
                ))
                part_num += 1
            part_start = chunk.start_line + i
            current_lines = []
            current_len = 0

        current_lines.append(line)
        current_len += len(line)

    # Remaining content
    if current_lines:
        text = "".join(current_lines)
        if len(text.strip()) >= MIN_CHUNK_CHARS:
            parts.append(CodeChunk(
                file_path=chunk.file_path,
                chunk_type=chunk.chunk_type,
                name=f"{chunk.name}_part{part_num}" if part_num > 1 or len(parts) > 0 else chunk.name,
                content=text,
                start_line=part_start,
                end_line=chunk.end_line,
                language=chunk.language,
                metadata=chunk.metadata.copy(),
            ))

    return parts if parts else [chunk]


# ---------------------------------------------------------------------------
# MATLAB chunking
# ---------------------------------------------------------------------------

# Matches function declarations at the start of a line (ignoring leading whitespace)
_MATLAB_FUNC_RE = re.compile(r"^\s*function\s+", re.MULTILINE)
_MATLAB_FUNC_SIG_RE = re.compile(
    r"^\s*function\s+(?:\[?[^\]]*\]?\s*=\s*)?(\w+)",
)


def _extract_matlab_metadata(content: str) -> dict:
    """Extract metadata from MATLAB function content."""
    meta: dict = {}
    if re.search(r"\bparfor\b", content):
        meta["uses_parfor"] = True
    if re.search(r"\bsystem\s*\(", content):
        meta["calls_system"] = True
    if re.search(r"\beval\s*\(", content):
        meta["calls_eval"] = True
    # Extract first line as signature
    first_line = content.strip().split("\n")[0] if content.strip() else ""
    if first_line:
        meta["signature"] = first_line.strip()
    return meta


def _chunk_matlab(file_path: str, content: str) -> List[CodeChunk]:
    """Split a MATLAB file on function declarations."""
    lines = content.splitlines(keepends=True)
    if not lines:
        return []

    # Find all function declaration line indices
    func_starts: list[int] = []
    for i, line in enumerate(lines):
        if _MATLAB_FUNC_RE.match(line):
            func_starts.append(i)

    # No functions → whole file is one chunk
    if not func_starts:
        text = content
        if len(text.strip()) < MIN_CHUNK_CHARS:
            return []
        chunk = CodeChunk(
            file_path=file_path,
            chunk_type="section",
            name=os.path.basename(file_path),
            content=text,
            start_line=1,
            end_line=len(lines),
            language="matlab",
            metadata=_extract_matlab_metadata(text),
        )
        return _split_oversized(chunk)

    chunks: List[CodeChunk] = []

    # Header: everything before the first function
    if func_starts[0] > 0:
        header_lines = lines[:func_starts[0]]
        header_text = "".join(header_lines)
        if len(header_text.strip()) >= MIN_CHUNK_CHARS:
            chunks.append(CodeChunk(
                file_path=file_path,
                chunk_type="section",
                name="header",
                content=header_text,
                start_line=1,
                end_line=func_starts[0],
                language="matlab",
                metadata={},
            ))

    # Each function
    for idx, start in enumerate(func_starts):
        end = func_starts[idx + 1] if idx + 1 < len(func_starts) else len(lines)
        func_lines = lines[start:end]
        func_text = "".join(func_lines)

        # Extract function name
        sig_match = _MATLAB_FUNC_SIG_RE.match(func_lines[0])
        func_name = sig_match.group(1) if sig_match else f"anonymous_{idx}"

        if len(func_text.strip()) < MIN_CHUNK_CHARS:
            continue

        chunk = CodeChunk(
            file_path=file_path,
            chunk_type="function",
            name=func_name,
            content=func_text,
            start_line=start + 1,
            end_line=end,
            language="matlab",
            metadata=_extract_matlab_metadata(func_text),
        )
        chunks.extend(_split_oversized(chunk))

    return chunks


# ---------------------------------------------------------------------------
# Python chunking
# ---------------------------------------------------------------------------

_PY_TOPLEVEL_RE = re.compile(r"^(class|def)\s+(\w+)", re.MULTILINE)


def _extract_python_metadata(content: str) -> dict:
    """Extract metadata from Python code content."""
    meta: dict = {}
    imports = []
    decorators = []
    for line in content.splitlines():
        stripped = line.strip()
        if stripped.startswith("import ") or stripped.startswith("from "):
            imports.append(stripped)
        elif stripped.startswith("@"):
            decorators.append(stripped)
    if imports:
        meta["imports"] = imports
    if decorators:
        meta["decorators"] = decorators
    return meta


def _chunk_python(file_path: str, content: str) -> List[CodeChunk]:
    """Split a Python file on top-level class and def declarations."""
    lines = content.splitlines(keepends=True)
    if not lines:
        return []

    # Find top-level class/def declarations (zero indentation)
    toplevel_starts: list[tuple[int, str, str]] = []  # (line_idx, kind, name)
    for i, line in enumerate(lines):
        # Must start at column 0 (no indentation)
        if line and not line[0].isspace():
            m = _PY_TOPLEVEL_RE.match(line)
            if m:
                # Include preceding decorator lines
                start = i
                while start > 0 and lines[start - 1].strip().startswith("@"):
                    start -= 1
                toplevel_starts.append((start, m.group(1), m.group(2)))

    if not toplevel_starts:
        text = content
        if len(text.strip()) < MIN_CHUNK_CHARS:
            return []
        chunk = CodeChunk(
            file_path=file_path,
            chunk_type="section",
            name=os.path.basename(file_path),
            content=text,
            start_line=1,
            end_line=len(lines),
            language="python",
            metadata=_extract_python_metadata(text),
        )
        return _split_oversized(chunk)

    chunks: List[CodeChunk] = []

    # Module header: everything before the first top-level declaration
    first_start = toplevel_starts[0][0]
    if first_start > 0:
        header_lines = lines[:first_start]
        header_text = "".join(header_lines)
        if len(header_text.strip()) >= MIN_CHUNK_CHARS:
            chunks.append(CodeChunk(
                file_path=file_path,
                chunk_type="section",
                name="module_header",
                content=header_text,
                start_line=1,
                end_line=first_start,
                language="python",
                metadata=_extract_python_metadata(header_text),
            ))

    # Each top-level class/def
    for idx, (start, kind, name) in enumerate(toplevel_starts):
        end = toplevel_starts[idx + 1][0] if idx + 1 < len(toplevel_starts) else len(lines)
        block_lines = lines[start:end]
        block_text = "".join(block_lines)

        if len(block_text.strip()) < MIN_CHUNK_CHARS:
            continue

        chunk_type = "class" if kind == "class" else "function"

        chunk = CodeChunk(
            file_path=file_path,
            chunk_type=chunk_type,
            name=name,
            content=block_text,
            start_line=start + 1,
            end_line=end,
            language="python",
            metadata=_extract_python_metadata(block_text),
        )
        chunks.extend(_split_oversized(chunk))

    return chunks


# ---------------------------------------------------------------------------
# Markdown chunking
# ---------------------------------------------------------------------------

_MD_H2_RE = re.compile(r"^##\s+(.+)", re.MULTILINE)


def _chunk_markdown(file_path: str, content: str) -> List[CodeChunk]:
    """Split a Markdown file on ## headings."""
    lines = content.splitlines(keepends=True)
    if not lines:
        return []

    # Find ## heading line indices
    h2_starts: list[tuple[int, str]] = []
    for i, line in enumerate(lines):
        m = _MD_H2_RE.match(line)
        if m:
            h2_starts.append((i, m.group(1).strip()))

    if not h2_starts:
        text = content
        if len(text.strip()) < MIN_CHUNK_CHARS:
            return []
        chunk = CodeChunk(
            file_path=file_path,
            chunk_type="doc",
            name=os.path.basename(file_path),
            content=text,
            start_line=1,
            end_line=len(lines),
            language="markdown",
            metadata={},
        )
        return _split_oversized(chunk)

    chunks: List[CodeChunk] = []

    # Preamble before first ##
    if h2_starts[0][0] > 0:
        preamble_lines = lines[:h2_starts[0][0]]
        preamble_text = "".join(preamble_lines)
        if len(preamble_text.strip()) >= MIN_CHUNK_CHARS:
            chunks.append(CodeChunk(
                file_path=file_path,
                chunk_type="doc",
                name="preamble",
                content=preamble_text,
                start_line=1,
                end_line=h2_starts[0][0],
                language="markdown",
                metadata={},
            ))

    for idx, (start, heading) in enumerate(h2_starts):
        end = h2_starts[idx + 1][0] if idx + 1 < len(h2_starts) else len(lines)
        section_lines = lines[start:end]
        section_text = "".join(section_lines)

        if len(section_text.strip()) < MIN_CHUNK_CHARS:
            continue

        chunk = CodeChunk(
            file_path=file_path,
            chunk_type="doc",
            name=heading,
            content=section_text,
            start_line=start + 1,
            end_line=end,
            language="markdown",
            metadata={},
        )
        chunks.extend(_split_oversized(chunk))

    return chunks


# ---------------------------------------------------------------------------
# JSON chunking
# ---------------------------------------------------------------------------

def _chunk_json(file_path: str, content: str) -> List[CodeChunk]:
    """Index a JSON config file as a single chunk."""
    if len(content.strip()) < MIN_CHUNK_CHARS:
        return []
    chunk = CodeChunk(
        file_path=file_path,
        chunk_type="config",
        name=os.path.basename(file_path),
        content=content,
        start_line=1,
        end_line=content.count("\n") + 1,
        language="json",
        metadata={},
    )
    return _split_oversized(chunk)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def chunk_file(file_path: str, content: str) -> List[CodeChunk]:
    """Split a source file into semantic chunks.

    Parameters
    ----------
    file_path : str
        Path relative to repo root (e.g. "pipeline/core/fit_models.m").
    content : str
        The file's text content.

    Returns
    -------
    list of CodeChunk
    """
    ext = os.path.splitext(file_path)[1].lower()

    if ext == ".m":
        return _chunk_matlab(file_path, content)
    elif ext == ".py":
        return _chunk_python(file_path, content)
    elif ext == ".md":
        return _chunk_markdown(file_path, content)
    elif ext == ".json":
        return _chunk_json(file_path, content)
    else:
        return []


def _should_skip(rel_path: str) -> bool:
    """Return True if a file should be skipped during indexing."""
    # Normalize separators
    norm = rel_path.replace("\\", "/")

    # Skip single-segment directory names anywhere in the path
    parts = norm.split("/")
    for part in parts:
        if part in _SKIP_DIRS:
            return True

    # Skip multi-segment directory paths (e.g. "pipeline/dependencies")
    for skip_path in _SKIP_DIR_PATHS:
        if norm.startswith(skip_path + "/") or ("/" + skip_path + "/") in norm:
            return True
    # Also match "dependencies/" as a bare directory component for robustness
    # (catches cases like worktree copies outside the normal tree)
    if "/dependencies/" in norm or norm.startswith("dependencies/"):
        return True

    # Skip specific files
    basename = os.path.basename(norm)
    if basename in _SKIP_FILES:
        return True

    # Skip binary extensions
    _, ext = os.path.splitext(norm)
    if ext.lower() in _BINARY_EXTENSIONS:
        return True

    # JSON files: only index specific ones
    if ext.lower() == ".json":
        if basename not in _INDEXABLE_JSON:
            return True

    return False


def chunk_repo(repo_root: str) -> List[CodeChunk]:
    """Walk the repository and chunk all indexable files.

    Parameters
    ----------
    repo_root : str
        Absolute path to the repository root.

    Returns
    -------
    list of CodeChunk
    """
    all_chunks: List[CodeChunk] = []

    for dirpath, dirnames, filenames in os.walk(repo_root):
        # Prune skipped directories in-place to avoid descending
        dirnames[:] = [
            d for d in dirnames
            if not _should_skip(os.path.relpath(os.path.join(dirpath, d), repo_root))
        ]

        for filename in filenames:
            abs_path = os.path.join(dirpath, filename)
            rel_path = os.path.relpath(abs_path, repo_root).replace("\\", "/")

            if _should_skip(rel_path):
                continue

            ext = os.path.splitext(filename)[1].lower()
            if ext not in _EXTENSION_LANG:
                continue

            try:
                with open(abs_path, "r", encoding="utf-8", errors="replace") as f:
                    content = f.read()
            except OSError:
                continue

            chunks = chunk_file(rel_path, content)
            all_chunks.extend(chunks)

    return all_chunks
