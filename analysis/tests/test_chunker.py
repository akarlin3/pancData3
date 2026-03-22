"""Tests for improvement_loop.rag.chunker module."""

import os
import sys
import textwrap

import pytest

# Ensure repo root is on sys.path for absolute imports
_repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from improvement_loop.rag.chunker import (
    CodeChunk,
    chunk_file,
    chunk_repo,
    MAX_CHUNK_CHARS,
    MIN_CHUNK_CHARS,
)


# ---------------------------------------------------------------------------
# MATLAB chunking
# ---------------------------------------------------------------------------

SAMPLE_MATLAB = textwrap.dedent("""\
    % File header comment
    % This file contains utility functions.

    function result = add_numbers(a, b)
    % ADD_NUMBERS Adds two numbers together.
        result = a + b;
    end

    function result = multiply_numbers(a, b)
    % MULTIPLY_NUMBERS Multiplies two numbers.
        result = a * b;
    end

    function result = subtract_numbers(a, b)
    % SUBTRACT_NUMBERS Subtracts b from a.
        result = a - b;
    end
""")


class TestChunkMATLAB:
    """chunk_file() with MATLAB files."""

    def test_three_functions_plus_header(self):
        chunks = chunk_file("pipeline/utils/math_ops.m", SAMPLE_MATLAB)
        # header + 3 functions = 4 chunks
        assert len(chunks) == 4

    def test_header_chunk(self):
        chunks = chunk_file("pipeline/utils/math_ops.m", SAMPLE_MATLAB)
        header = chunks[0]
        assert header.chunk_type == "section"
        assert header.name == "header"
        assert "File header comment" in header.content

    def test_function_chunks(self):
        chunks = chunk_file("pipeline/utils/math_ops.m", SAMPLE_MATLAB)
        func_chunks = [c for c in chunks if c.chunk_type == "function"]
        assert len(func_chunks) == 3
        names = {c.name for c in func_chunks}
        assert names == {"add_numbers", "multiply_numbers", "subtract_numbers"}

    def test_language_is_matlab(self):
        chunks = chunk_file("pipeline/utils/math_ops.m", SAMPLE_MATLAB)
        assert all(c.language == "matlab" for c in chunks)

    def test_start_line_end_line(self):
        chunks = chunk_file("pipeline/utils/math_ops.m", SAMPLE_MATLAB)
        func_chunks = [c for c in chunks if c.chunk_type == "function"]
        for c in func_chunks:
            assert c.start_line >= 1
            assert c.end_line >= c.start_line

    def test_single_function_file(self):
        content = textwrap.dedent("""\
            function y = square(x)
            % SQUARE Squares each element of the input array x element-wise.
                y = x .^ 2;
            end
        """)
        chunks = chunk_file("pipeline/utils/square.m", content)
        assert len(chunks) == 1
        assert chunks[0].chunk_type == "function"
        assert chunks[0].name == "square"

    def test_no_functions(self):
        content = textwrap.dedent("""\
            % This is a script, not a function file.
            x = 1:10;
            y = x .^ 2;
            plot(x, y);
        """)
        chunks = chunk_file("pipeline/scripts/plot_stuff.m", content)
        assert len(chunks) == 1
        assert chunks[0].chunk_type == "section"

    def test_metadata_parfor(self):
        content = textwrap.dedent("""\
            function process_all(data)
                parfor i = 1:length(data)
                    data(i) = data(i) * 2;
                end
            end
        """)
        chunks = chunk_file("pipeline/core/proc.m", content)
        assert chunks[0].metadata.get("uses_parfor") is True

    def test_metadata_system(self):
        content = textwrap.dedent("""\
            function run_cmd(cmd)
                [status, output] = system(cmd);
            end
        """)
        chunks = chunk_file("pipeline/utils/run.m", content)
        assert chunks[0].metadata.get("calls_system") is True

    def test_metadata_signature(self):
        chunks = chunk_file("pipeline/utils/math_ops.m", SAMPLE_MATLAB)
        func_chunks = [c for c in chunks if c.name == "add_numbers"]
        assert "function result = add_numbers(a, b)" in func_chunks[0].metadata.get("signature", "")


# ---------------------------------------------------------------------------
# Python chunking
# ---------------------------------------------------------------------------

SAMPLE_PYTHON = textwrap.dedent("""\
    \"\"\"Module docstring.\"\"\"

    import os
    import sys

    CONSTANT = 42


    class MyProcessor:
        \"\"\"Processes data.\"\"\"

        def __init__(self, data):
            self.data = data

        def run(self):
            return self.data * 2


    def helper_function(x):
        \"\"\"A standalone helper.\"\"\"
        return x + 1
""")


class TestChunkPython:
    """chunk_file() with Python files."""

    def test_class_and_function(self):
        chunks = chunk_file("analysis/processor.py", SAMPLE_PYTHON)
        types = [c.chunk_type for c in chunks]
        assert "class" in types
        assert "function" in types

    def test_module_header(self):
        chunks = chunk_file("analysis/processor.py", SAMPLE_PYTHON)
        header = [c for c in chunks if c.name == "module_header"]
        assert len(header) == 1
        assert "import os" in header[0].content
        assert "CONSTANT = 42" in header[0].content

    def test_class_contains_methods(self):
        """Methods should NOT be separate chunks — whole class stays together."""
        chunks = chunk_file("analysis/processor.py", SAMPLE_PYTHON)
        class_chunks = [c for c in chunks if c.chunk_type == "class"]
        assert len(class_chunks) == 1
        assert "__init__" in class_chunks[0].content
        assert "def run" in class_chunks[0].content

    def test_language_is_python(self):
        chunks = chunk_file("analysis/processor.py", SAMPLE_PYTHON)
        assert all(c.language == "python" for c in chunks)

    def test_metadata_imports(self):
        chunks = chunk_file("analysis/processor.py", SAMPLE_PYTHON)
        header = [c for c in chunks if c.name == "module_header"][0]
        assert "import os" in header.metadata.get("imports", [])

    def test_decorator_included(self):
        content = textwrap.dedent("""\
            import functools

            @functools.lru_cache
            def cached_fn(x):
                return x * 2
        """)
        chunks = chunk_file("analysis/cached.py", content)
        func = [c for c in chunks if c.chunk_type == "function"]
        assert len(func) == 1
        assert "@functools.lru_cache" in func[0].content
        assert "decorators" in func[0].metadata


# ---------------------------------------------------------------------------
# Markdown chunking
# ---------------------------------------------------------------------------

SAMPLE_MARKDOWN = textwrap.dedent("""\
    # Main Title

    Some preamble text here that is long enough to be indexed as a chunk by the chunker.

    ## Section One

    Content for section one goes here. It needs to be long enough to pass the minimum chunk size.

    ## Section Two

    Content for section two goes here. Also needs to be long enough to be a valid chunk.

    ## Section Three

    Content for section three goes here. Third section with enough content for the chunker.
""")


class TestChunkMarkdown:
    """chunk_file() with Markdown files."""

    def test_three_sections(self):
        chunks = chunk_file("docs/guide.md", SAMPLE_MARKDOWN)
        doc_chunks = [c for c in chunks if c.name != "preamble"]
        assert len(doc_chunks) == 3

    def test_section_names(self):
        chunks = chunk_file("docs/guide.md", SAMPLE_MARKDOWN)
        names = {c.name for c in chunks if c.name != "preamble"}
        assert names == {"Section One", "Section Two", "Section Three"}

    def test_all_doc_type(self):
        chunks = chunk_file("docs/guide.md", SAMPLE_MARKDOWN)
        assert all(c.chunk_type == "doc" for c in chunks)

    def test_language_is_markdown(self):
        chunks = chunk_file("docs/guide.md", SAMPLE_MARKDOWN)
        assert all(c.language == "markdown" for c in chunks)

    def test_preamble(self):
        chunks = chunk_file("docs/guide.md", SAMPLE_MARKDOWN)
        preambles = [c for c in chunks if c.name == "preamble"]
        assert len(preambles) == 1
        assert "Main Title" in preambles[0].content


# ---------------------------------------------------------------------------
# JSON chunking
# ---------------------------------------------------------------------------

class TestChunkJSON:
    """chunk_file() with JSON config files."""

    def test_single_chunk(self):
        content = '{\n  "key": "value",\n  "number": 42,\n  "nested": {"a": 1}\n}'
        chunks = chunk_file("config.example.json", content)
        assert len(chunks) == 1
        assert chunks[0].chunk_type == "config"
        assert chunks[0].language == "json"

    def test_too_small_skipped(self):
        chunks = chunk_file("tiny.json", '{}')
        assert len(chunks) == 0


# ---------------------------------------------------------------------------
# chunk_repo() — skip rules
# ---------------------------------------------------------------------------

class TestChunkRepoSkips:
    """chunk_repo() correctly skips pipeline/dependencies/ and other exclusions."""

    def test_dependencies_skipped(self, tmp_path):
        """Files under pipeline/dependencies/ must not appear in output."""
        # Create a minimal repo structure
        deps_dir = tmp_path / "pipeline" / "dependencies"
        deps_dir.mkdir(parents=True)
        (deps_dir / "third_party.m").write_text(
            "function y = third_party(x)\n% THIRD_PARTY A third-party function that should be skipped.\n    y = x;\nend\n",
            encoding="utf-8",
        )

        # Also create an indexable file (must exceed MIN_CHUNK_CHARS)
        core_dir = tmp_path / "pipeline" / "core"
        core_dir.mkdir(parents=True)
        (core_dir / "my_func.m").write_text(
            "function y = my_func(x)\n% MY_FUNC Adds one to the input value for testing purposes.\n    y = x + 1;\nend\n",
            encoding="utf-8",
        )

        chunks = chunk_repo(str(tmp_path))
        paths = {c.file_path for c in chunks}
        assert any("my_func" in p for p in paths)
        assert not any("dependencies" in p for p in paths)

    def test_pycache_skipped(self, tmp_path):
        pycache = tmp_path / "__pycache__"
        pycache.mkdir()
        (pycache / "module.cpython-312.pyc").write_text("fake", encoding="utf-8")

        (tmp_path / "real.py").write_text(
            "def real_func():\n    return True\n",
            encoding="utf-8",
        )

        chunks = chunk_repo(str(tmp_path))
        paths = {c.file_path for c in chunks}
        assert not any("__pycache__" in p for p in paths)

    def test_binary_files_skipped(self, tmp_path):
        (tmp_path / "data.mat").write_text("binary", encoding="utf-8")
        (tmp_path / "image.png").write_text("binary", encoding="utf-8")
        chunks = chunk_repo(str(tmp_path))
        assert len(chunks) == 0

    def test_improvement_loop_log_skipped(self, tmp_path):
        (tmp_path / "improvement_loop_log.json").write_text(
            '{"iteration": 1}', encoding="utf-8"
        )
        chunks = chunk_repo(str(tmp_path))
        assert not any("improvement_loop_log" in c.file_path for c in chunks)

    def test_indexable_json_included(self, tmp_path):
        content = '{\n  "dataloc": "/data",\n  "skip_tests": false,\n  "key": "value_here"\n}'
        (tmp_path / "config.example.json").write_text(content, encoding="utf-8")
        chunks = chunk_repo(str(tmp_path))
        assert any("config.example.json" in c.file_path for c in chunks)


# ---------------------------------------------------------------------------
# Oversized chunk splitting
# ---------------------------------------------------------------------------

class TestOversizedSplitting:
    """Chunks exceeding MAX_CHUNK_CHARS are split at logical boundaries."""

    def test_large_function_split(self):
        """A function over 4000 chars should be split into multiple parts."""
        # Create a function with many lines totalling > 4000 chars
        body_lines = [f"    x{i} = {i}; % padding to make this line long enough\n"
                      for i in range(120)]
        content = "function result = big_func(x)\n" + "".join(body_lines) + "end\n"
        assert len(content) > MAX_CHUNK_CHARS

        chunks = chunk_file("pipeline/utils/big.m", content)
        assert len(chunks) >= 2
        assert all(len(c.content) <= MAX_CHUNK_CHARS for c in chunks)

    def test_split_names_have_part_suffix(self):
        body_lines = [f"    x{i} = {i}; % padding to make this line long enough\n"
                      for i in range(120)]
        content = "function result = big_func(x)\n" + "".join(body_lines) + "end\n"
        chunks = chunk_file("pipeline/utils/big.m", content)
        part_names = [c.name for c in chunks]
        assert any("_part" in name for name in part_names)

    def test_small_function_not_split(self):
        content = "function y = small(x)\n% SMALL Returns the input unchanged, used for testing.\n    y = x;\nend\n"
        chunks = chunk_file("pipeline/utils/small.m", content)
        assert len(chunks) == 1
        assert "_part" not in chunks[0].name


# ---------------------------------------------------------------------------
# chunk_repo() on the real repo
# ---------------------------------------------------------------------------

class TestChunkRepoReal:
    """chunk_repo() on the actual pancData3 repository."""

    def test_returns_nonempty(self):
        chunks = chunk_repo(_repo_root)
        assert len(chunks) > 0

    def test_no_dependencies_in_output(self):
        chunks = chunk_repo(_repo_root)
        for c in chunks:
            assert "dependencies/" not in c.file_path, (
                f"dependencies file indexed: {c.file_path}"
            )

    def test_contains_fit_models(self):
        chunks = chunk_repo(_repo_root)
        fit_chunks = [c for c in chunks if "fit_models" in c.file_path]
        assert len(fit_chunks) > 0

    def test_contains_python_chunks(self):
        chunks = chunk_repo(_repo_root)
        py_chunks = [c for c in chunks if c.language == "python"]
        assert len(py_chunks) > 0

    def test_contains_markdown_chunks(self):
        chunks = chunk_repo(_repo_root)
        md_chunks = [c for c in chunks if c.language == "markdown"]
        assert len(md_chunks) > 0


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    """Edge cases for chunk_file()."""

    def test_empty_file(self):
        assert chunk_file("empty.py", "") == []
        assert chunk_file("empty.m", "") == []
        assert chunk_file("empty.md", "") == []

    def test_unknown_extension(self):
        assert chunk_file("readme.rst", "some content") == []

    def test_min_chunk_size(self):
        """Content shorter than MIN_CHUNK_CHARS is skipped."""
        tiny = "x = 1;"
        assert len(tiny) < MIN_CHUNK_CHARS
        chunks = chunk_file("tiny.m", tiny)
        assert len(chunks) == 0
