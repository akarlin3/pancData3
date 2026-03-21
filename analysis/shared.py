#!/usr/bin/env python3
"""Shared utilities for the pancData3 analysis suite.

This module provides common functions used by all analysis scripts:

- **Folder resolution**: locating the latest ``saved_files_*`` output directory.
- **DWI type parsing**: extracting Standard / dnCNN / IVIMnet labels from paths.
- **Statistical extraction**: regex-based extraction of p-values and correlation
  coefficients from free-text summaries produced by the vision model.
- **CSV loading**: reading the vision analysis CSV into structured dicts.
- **Grouping**: organising graph rows by normalised graph name for cross-DWI
  comparison.
- **UTF-8 stdout setup**: prevents encoding crashes on Windows consoles when
  printing emoji characters used in pipeline logging.
"""

from __future__ import annotations

import csv
import glob
import io
import json
import os
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path

# ── Analysis configuration ───────────────────────────────────────────────────
# Centralised defaults for all analysis scripts.  These are overridden by
# values in ``analysis_config.json`` (if present) and, for ``dwi_types``,
# optionally by the MATLAB ``config.json`` at the repository root.

_DEFAULTS: dict = {
    "dwi_types": ["Standard", "dnCNN", "IVIMnet"],
    "vision": {
        "gemini_model": "gemini-2.5-flash",
        "claude_model": "claude-opus-4-6",
        "provider": "claude",
        "gemini_api_key": "",
        "anthropic_api_key": "",
        "max_concurrent_requests": 2,
        "max_retries": 4,
        "max_output_tokens": 10000,
        "backoff_base_seconds": 15,
        "request_timeout_seconds": 120,
        "script_timeout_seconds": 1800,
        "fallback_to_local": True,
    },
    "statistics": {
        "p_highly_significant": 0.001,
        "p_significant": 0.01,
        "p_noteworthy": 0.05,
        "correlation_threshold": 0.3,
        "effect_size_medium": 0.5,
        "effect_size_large": 0.8,
    },
    "priority_graphs": [
        "Dose_vs_Diffusion",
        "Longitudinal_Mean_Metrics",
        "Longitudinal_Mean_Metrics_ByOutcome",
        "Longitudinal_Mean_Metrics_LC",
        "Longitudinal_Mean_Metrics_LF",
        "Feature_BoxPlots",
        "Feature_Histograms",
        "core_method_dice_heatmap",
        "core_method_volume_comparison",
    ],
}

# Resolved at first access via :func:`get_config`.
_config_cache: dict | None = None

# Map from provider name to (config key, environment variable name)
_API_KEY_MAP: dict[str, tuple[str, str]] = {
    "gemini": ("gemini_api_key", "GEMINI_API_KEY"),
    "claude": ("anthropic_api_key", "ANTHROPIC_API_KEY"),
    "anthropic": ("anthropic_api_key", "ANTHROPIC_API_KEY"),
}


def get_api_key(provider: str, cfg: dict | None = None) -> str | None:
    """Resolve an API key for *provider*, checking config then environment.

    Resolution order:
    1. ``vision.<provider>_api_key`` in the analysis config.
    2. The corresponding environment variable (``GEMINI_API_KEY`` or
       ``ANTHROPIC_API_KEY``).

    Returns ``None`` if no key is found.
    """
    mapping = _API_KEY_MAP.get(provider.lower())
    if mapping is None:
        return os.environ.get(f"{provider.upper()}_API_KEY")
    config_key, env_var = mapping
    if cfg is None:
        cfg = get_config()
    key = cfg.get("vision", {}).get(config_key, "")
    if key:
        return key
    return os.environ.get(env_var) or None


def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge *override* into *base*, returning a new dict.

    Scalar values in *override* replace those in *base*; nested dicts are
    merged recursively so that partial overrides work correctly (e.g.
    overriding only ``vision.gemini_model`` without touching the other
    vision keys).  Lists are shallow-copied to prevent mutation of the
    originals.

    If *override* contains a ``None`` value for a key whose *base* value
    is a dict, the ``None`` is skipped and the base dict is preserved.
    More generally, if *override* contains a non-dict value for a key
    whose *base* value is a dict, the override is skipped and a warning
    is emitted.  This prevents downstream code from crashing with
    ``AttributeError`` when it does e.g.
    ``cfg.get('vision', {}).get(...)``.
    """
    import copy
    merged = {}
    for key, val in base.items():
        if isinstance(val, dict):
            merged[key] = _deep_merge(val, {})
        elif isinstance(val, list):
            merged[key] = list(val)  # type: ignore
        else:
            merged[key] = val
    for key, val in override.items():
        if key in merged and isinstance(merged[key], dict) and not isinstance(val, dict):
            # Preserve the base dict — overwriting a section-level dict
            # key with a non-dict value (None, list, scalar, etc.) is not
            # allowed because downstream code expects a dict.
            warnings.warn(
                f"[config] Ignoring override for key '{key}': cannot replace "
                f"a dict value with {type(val).__name__} ({val!r}). "
                f"The base dict is preserved.",
                stacklevel=2,
            )
            continue
        if key in merged and isinstance(merged[key], dict) and isinstance(val, dict):
            merged[key] = _deep_merge(merged[key], val)
        elif isinstance(val, list):
            merged[key] = list(val)  # type: ignore
        else:
            merged[key] = val
    return merged


def _normalize_dwi_type(dwi_type: str, canonical_types: list[str]) -> str:
    """Return the canonical-case version of *dwi_type* if it matches an
    existing entry (case-insensitive), otherwise return *dwi_type* unchanged.

    This prevents duplicate entries like ``['Standard', 'standard']`` when
    the MATLAB ``config.json`` uses a different capitalisation than the
    built-in defaults.
    """
    lower_map = {ct.lower(): ct for ct in canonical_types}
    return lower_map.get(dwi_type.lower(), dwi_type)


def _strip_json_comments(text: str) -> str:
    """Remove single-line (``//``) and block (``/* */``) comments from JSON
    text while preserving string literals that may contain ``//``.

    This is a token-aware stripper that tracks whether the current
    position is inside a JSON string literal.  Content inside strings
    (delimited by unescaped double quotes) is never modified, so URLs
    like ``"https://example.com/tool"`` are preserved intact.

    Block comments (``/* ... */``) are removed regardless of position
    outside strings.  Single-line comments (``// ...``) are removed from
    the current position to end-of-line when found outside strings.
    """
    result: list[str] = []
    i = 0
    n = len(text)
    while i < n:
        # Check if we're entering a JSON string literal
        if text[i] == '"':
            # Consume the entire string literal, including escaped characters
            result.append('"')
            i += 1
            while i < n:
                if text[i] == '\\' and i + 1 < n:
                    # Escaped character — consume both the backslash and
                    # the following character unconditionally.
                    result.append(text[i])
                    result.append(text[i + 1])
                    i += 2
                elif text[i] == '"':
                    result.append('"')
                    i += 1
                    break
                else:
                    result.append(text[i])
                    i += 1
        # Block comment outside a string
        elif text[i] == '/' and i + 1 < n and text[i + 1] == '*':
            # Skip until closing */
            end = text.find('*/', i + 2)
            if end == -1:
                # Unterminated block comment — skip rest of text
                break
            i = end + 2
        # Single-line comment outside a string
        elif text[i] == '/' and i + 1 < n and text[i + 1] == '/':
            # Skip until end of line
            while i < n and text[i] != '\n':
                i += 1
            # Keep the newline itself to preserve line structure
        else:
            result.append(text[i])
            i += 1
    return ''.join(result)


def _permissive_json_loads(text: str) -> dict:
    """Attempt to parse JSON text, stripping comments and trailing commas.

    Python's :func:`json.loads` is strict and rejects constructs that
    MATLAB's ``jsondecode`` (and many editors) silently accept, such as:

    - Single-line comments (``// ...``)
    - Block comments (``/* ... */``)
    - Trailing commas before closing ``}`` or ``]``

    This helper strips those constructs and retries parsing when the
    initial strict parse fails.  Importantly, comment stripping is
    *token-aware*: ``//`` sequences inside JSON string values (e.g. URLs)
    are preserved.

    Parameters
    ----------
    text : str
        Raw file content that is expected to be JSON (possibly with
        non-standard extensions).

    Returns
    -------
    dict
        Parsed JSON object.

    Raises
    ------
    json.JSONDecodeError
        If the text cannot be parsed even after sanitisation.
    """
    # First, try strict parsing — fast path for well-formed JSON.
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

    # Use token-aware comment stripping that preserves string contents
    sanitised = _strip_json_comments(text)
    # Strip trailing commas before } or ]
    sanitised = re.sub(r",\s*([}\]])", r"\1", sanitised)

    return json.loads(sanitised)


def load_analysis_config(
    config_path: str | Path | None = None,
    matlab_config_path: str | Path | None = None,
) -> dict:
    """Load the analysis configuration with layered overrides.

    Resolution order (later wins):
    1. Built-in defaults (``_DEFAULTS``).
    2. ``analysis_config.json`` — when no explicit *config_path* is given,
       the function searches two locations:

       a. ``analysis/analysis_config.json`` (next to this module) — checked
          first and takes priority.
       b. ``analysis_config.json`` at the repository root (one level above
          ``analysis/``).

       If an explicit *config_path* is provided, only that path is used.
    3. MATLAB ``config.json`` at the repository root — only ``dwi_type``
       is read (mapped to ``dwi_types`` list for consistency).

    Parameters
    ----------
    config_path : str or Path, optional
        Explicit path to the analysis config JSON.  If ``None``, the
        function searches for ``analysis_config.json`` in the
        ``analysis/`` directory first, then falls back to the repository
        root.
    matlab_config_path : str or Path, optional
        Explicit path to the MATLAB pipeline ``config.json``.  If
        ``None``, the function looks for ``config.json`` at the
        repository root (one level above ``analysis/``).

    Returns
    -------
    dict
        Fully resolved configuration dictionary.
    """
    cfg = dict(_DEFAULTS)
    # Deep-copy nested dicts so mutations don't pollute the template.
    cfg = _deep_merge({}, cfg)

    analysis_dir = Path(__file__).resolve().parent

    # ── Layer 2: analysis_config.json ──
    if config_path is None:
        # Check the analysis/ directory first (next to this file), then
        # fall back to the repository root.  This ensures that a config
        # placed alongside the analysis scripts takes priority.
        local_config = analysis_dir / "analysis_config.json"
        root_config = analysis_dir.parent / "analysis_config.json"
        if local_config.is_file():
            config_path = local_config
        else:
            config_path = root_config
    else:
        config_path = Path(config_path)

    if config_path.is_file():
        try:
            with open(config_path, encoding="utf-8") as f:
                overrides = json.load(f)
            cfg = _deep_merge(cfg, overrides)
        except json.JSONDecodeError as e:
            print(f'[WARN] Invalid JSON in {config_path}: {e}')
        except OSError as e:
            print(f'[WARN] Cannot read {config_path}: {e}')

    # ── Layer 3: MATLAB config.json (dwi_type only) ──
    if matlab_config_path is None:
        matlab_config_path = analysis_dir.parent / "config.json"
    else:
        matlab_config_path = Path(matlab_config_path)

    if matlab_config_path.is_file():
        try:
            with open(matlab_config_path, encoding="utf-8") as f:
                raw_text = f.read()
            try:
                matlab_cfg = _permissive_json_loads(raw_text)
            except json.JSONDecodeError as exc:
                print(
                    f"[WARN] Cannot parse MATLAB config {matlab_config_path}: "
                    f"{exc.__class__.__name__}: {exc}"
                )
                matlab_cfg = None

            if matlab_cfg is not None:
                dwi_type = matlab_cfg.get("dwi_type")
                if dwi_type and isinstance(dwi_type, str):
                    # Normalize to canonical case to avoid duplicates like
                    # ['Standard', 'standard'] when the MATLAB config uses
                    # different capitalisation than the built-in defaults.
                    dwi_type = _normalize_dwi_type(dwi_type, cfg["dwi_types"])
                    # The MATLAB config stores a single active type; we don't
                    # override the full list, but we can ensure it's present.
                    existing_lower = [d.lower() for d in cfg["dwi_types"]]
                    if dwi_type.lower() not in existing_lower:
                        cfg["dwi_types"].append(dwi_type)
        except OSError as exc:
            print(
                f"[WARN] Cannot read MATLAB config {matlab_config_path}: "
                f"{exc.__class__.__name__}: {exc}"
            )

    # ── Layer 4: Environment variable overrides ──
    # These allow run_analysis.py to propagate CLI flags (--gemini-model,
    # --claude-model, --provider, --concurrency) to child scripts that run
    # as separate subprocesses.
    env_model = os.environ.get("PANCDATA3_GEMINI_MODEL")
    if env_model:
        cfg["vision"]["gemini_model"] = env_model
    env_claude_model = os.environ.get("PANCDATA3_CLAUDE_MODEL")
    if env_claude_model:
        cfg["vision"]["claude_model"] = env_claude_model
    env_provider = os.environ.get("PANCDATA3_VISION_PROVIDER")
    if env_provider:
        cfg["vision"]["provider"] = env_provider
    env_concurrency = os.environ.get("PANCDATA3_GEMINI_CONCURRENCY")
    if env_concurrency:
        try:
            cfg["vision"]["max_concurrent_requests"] = int(env_concurrency)
        except ValueError:
            pass

    return cfg


def get_config() -> dict:
    """Return the cached analysis configuration (loading it on first call).

    This is the primary entry point for scripts that need config values.
    The result is cached so that repeated calls do not re-read disk.
    """
    global _config_cache
    if _config_cache is None:
        _config_cache = load_analysis_config()
    return _config_cache


def reset_config_cache() -> None:
    """Clear the cached config (useful in tests)."""
    global _config_cache
    _config_cache = None


# Canonical DWI type names in the order the pipeline processes them.
# Used throughout the suite for iteration and labelling.
# NOTE: This is kept as a module-level constant for backward compatibility;
# scripts that need the configurable list should use ``get_config()["dwi_types"]``.
DWI_TYPES = list(_DEFAULTS["dwi_types"])


def setup_utf8_stdout():
    """Force UTF-8 stdout/stderr on Windows to avoid cp1252 emoji crashes.

    On Python 3.7+ the ``reconfigure`` method is available; for older
    interpreters we fall back to wrapping the underlying byte buffer.
    This is a no-op on non-Windows platforms.
    """
    if sys.platform != "win32":
        return
    # Skip when running under test frameworks — modifying stdout/stderr
    # interferes with capture mechanisms and causes teardown errors.
    if "pytest" in sys.modules or "unittest" in sys.modules:
        return
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
            sys.stderr.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass  # Captured/redirected streams may reject reconfigure
    elif (getattr(sys.stdout, "encoding", "") or "").lower() != "utf-8":
            # Wrap the raw binary buffer with a new TextIOWrapper.
            # Guard: stdout/stderr may lack .buffer when captured (e.g.,
            # pytest StringIO or piped contexts).
            if hasattr(sys.stdout, "buffer"):
                sys.stdout = io.TextIOWrapper(
                    sys.stdout.buffer, encoding="utf-8", errors="replace",
                )
            if hasattr(sys.stderr, "buffer"):
                sys.stderr = io.TextIOWrapper(
                    sys.stderr.buffer, encoding="utf-8", errors="replace",
                )


def find_latest_saved_folder(base_dir: str | None = None) -> Path:
    """Return the most recently modified ``saved_files_*`` directory.

    The timestamp suffix in the folder name (``YYYYMMDD_HHMMSS``) ensures
    that lexicographic sorting picks the newest run when sorted in reverse.

    Parameters
    ----------
    base_dir : str, optional
        Directory to search in. Defaults to the repository root
        (one level above the ``analysis/`` folder).

    Returns
    -------
    Path
        Absolute path to the latest output folder.

    Raises
    ------
    SystemExit
        If no matching directories are found.
    """
    if base_dir is None:
        # Navigate from analysis/ up to the repository root.
        base_dir = str(Path(__file__).resolve().parent.parent)
    pattern = os.path.join(base_dir, "saved_files_*")
    # Reverse-sorted so the newest (lexicographically largest) comes first.
    dirs = sorted(glob.glob(pattern), reverse=True)
    if not dirs:
        sys.exit("ERROR: No saved_files_* folders found in " + base_dir)
    return Path(dirs[0])


def resolve_folder(argv: list[str] | None = None) -> Path:
    """Resolve the output folder from CLI args or auto-detect.

    If a folder path is provided as the second element of *argv*
    (i.e. the first positional CLI argument), it is validated and returned.
    Otherwise, :func:`find_latest_saved_folder` is used for auto-detection.

    Parameters
    ----------
    argv : list[str], optional
        Typically ``sys.argv``.  ``argv[1]`` is treated as the folder path.

    Returns
    -------
    Path
        Validated path to the output folder.

    Usage in scripts::

        folder = resolve_folder(sys.argv)
    """
    if argv and len(argv) > 1:
        p = Path(argv[1])
        if p.is_dir():
            return p
        sys.exit(f"ERROR: Specified folder does not exist: {p}")
    return find_latest_saved_folder()


def parse_dwi_info(file_path: str) -> tuple[str, str]:
    """Extract DWI type and normalised base graph name from a file path.

    The function looks for a ``saved_files_*`` component in the path and
    checks whether the next component is one of the known DWI types.
    DWI-type suffixes (e.g. ``_Standard``) and the ``.png`` extension are
    stripped from the filename to produce a normalised ``base_name`` that
    can be used for cross-DWI comparisons.

    Parameters
    ----------
    file_path : str
        Absolute or relative path to a graph image file.

    Returns
    -------
    dwi_type : str
        One of ``"Standard"``, ``"dnCNN"``, ``"IVIMnet"``, or ``"Root"``
        (when the file is not inside a DWI-type subfolder).
    base_name : str
        Filename without DWI suffix or ``.png`` extension.
    """
    # Normalise path separators for cross-platform compatibility.
    fp = file_path.replace("\\", "/")
    parts = fp.split("/")
    dwi_type = "Root"
    base_name = parts[-1]
    # Walk path components to find the saved_files directory.
    for i, p in enumerate(parts):
        if "saved_files" in p and i + 1 < len(parts):
            # The directory immediately after saved_files is the DWI type.
            if parts[i + 1] in DWI_TYPES:  # type: ignore
                dwi_type = parts[i + 1]  # type: ignore
            break
    # Remove DWI-type suffixes so "Longitudinal_Mean_Metrics_Standard"
    # becomes "Longitudinal_Mean_Metrics" for cross-type matching.
    for t in ["_Standard", "_dnCNN", "_IVIMnet"]:  # type: ignore
        base_name = base_name.replace(t, "")
    return dwi_type, base_name.replace(".png", "")


def safe_text(row: dict, *keys: str) -> str:
    """Safely concatenate multiple CSV row fields into one string.

    Returns a single space-joined string.  Missing or ``None``-valued
    fields are silently replaced with empty strings so that downstream
    regex extraction never crashes on incomplete vision CSV rows.
    """
    parts = []
    for k in keys:
        v = row.get(k, "") or ""
        parts.append(v)
    return " ".join(parts)


def extract_pvalues(text: str) -> list[tuple[float, str]]:
    """Extract all p-values from a text string.

    Searches for patterns like ``p = 0.032``, ``p < 0.001``, or
    ``p-value = 1.2e-4`` (case-insensitive).

    Parameters
    ----------
    text : str
        Free-text string to search (typically a graph summary or
        trend description produced by the vision model).

    Returns
    -------
    list[tuple[float, str]]
        Each tuple contains ``(p_value, context_snippet)`` where the
        context snippet is up to 80 characters before and 40 characters
        after the match to provide surrounding context.
    """
    results: list[tuple[float, str]] = []
    patterns = [
        # Matches: p-value = 0.04, p-value<0.001 (checked first, more specific)
        r"p-value\s*[=<>]\s*([\d.]+(?:e[+-]?\d+)?)",
        # Matches: p = 0.03, p<0.001, p > 1.2e-4
        # (?<![a-zA-Z]) prevents matching "up = 2.5", "group = ...", etc.
        r"(?<![a-zA-Z])p\s*[=<>]\s*([\d.]+(?:e[+-]?\d+)?)",
    ]
    seen_spans: set[tuple[int, int]] = set()
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE):
            # Deduplicate overlapping matches (p-value pattern is a superset).
            span = (m.start(), m.end())
            if any(s[0] <= span[0] < s[1] or s[0] < span[1] <= s[1] for s in seen_spans):
                continue
            seen_spans.add(span)
            try:
                val = float(m.group(1))
                # Skip values that are clearly not p-values (> 1.0 or negative).
                if val > 1.0:
                    continue
                # Extract surrounding context for downstream display.
                start = max(0, m.start() - 80)
                end = min(len(text), m.end() + 40)
                context = text[start:end].strip()  # type: ignore
                results.append((val, context))
            except ValueError:
                pass
    return results


def extract_correlations(text: str) -> list[tuple[float, str]]:
    """Extract correlation coefficients (r, rs, r-squared) from text.

    Searches for patterns like ``r = 0.78``, ``rs = -0.45``, or
    ``r\u00b2 = 0.61`` (where ``\\xb2`` is the superscript-2 character).

    Parameters
    ----------
    text : str
        Free-text string to search.

    Returns
    -------
    list[tuple[float, str]]
        Each tuple contains ``(r_value, context_snippet)`` with up to
        60 characters of surrounding context on either side.
    """
    results: list[tuple[float, str]] = []
    patterns = [
        # Spearman rank correlation: "rs = 0.65" (checked first, more specific)
        r"(?<![a-zA-Z])rs\s*=\s*(-?[\d.]+)",
        # Coefficient of determination: "r² = 0.61"  (\xb2 = superscript 2)
        r"(?<![a-zA-Z])r\xb2\s*=\s*([\d.]+)",
        # Pearson r: "r = 0.78" or "r = -0.45"
        # (?<![a-zA-Z]) prevents matching "parameter = ...", "error = ...", etc.
        r"(?<![a-zA-Z])r\s*=\s*(-?[\d.]+)",
    ]
    seen_spans: set[tuple[int, int]] = set()
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE):
            span = (m.start(), m.end())
            if any(s[0] <= span[0] < s[1] or s[0] < span[1] <= s[1] for s in seen_spans):
                continue
            seen_spans.add(span)
            try:
                val = float(m.group(1))
                # Skip values that are clearly not correlation coefficients.
                if abs(val) > 1.0:
                    continue
                start = max(0, m.start() - 60)
                end = min(len(text), m.end() + 60)
                context = text[start:end].strip()  # type: ignore
                results.append((val, context))
            except ValueError:
                pass
    return results


def load_graph_csv(folder: Path) -> list[dict]:
    """Load ``graph_analysis_results.csv`` from a ``saved_files`` folder.

    This CSV is produced by :mod:`batch_graph_analysis` and contains one
    row per analysed graph image with columns for axes, trends, inflection
    points, and a plain-English summary.

    Parameters
    ----------
    folder : Path
        Path to the ``saved_files_*`` output folder.

    Returns
    -------
    list[dict]
        List of row dicts (one per graph).  Returns an empty list if
        the CSV does not exist (graceful degradation).
    """
    csv_path = folder / "graph_analysis_results.csv"
    if not csv_path.exists():
        return []
    rows: list[dict] = []
    with open(csv_path, encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def group_by_graph_name(rows: list[dict]) -> dict[str, dict[str, dict]]:
    """Group CSV rows by normalised graph name, then by DWI type.

    This enables cross-DWI comparison: for the same graph
    (e.g. ``"Longitudinal_Mean_Metrics"``), all three DWI-type rows are
    gathered under a single key.

    Parameters
    ----------
    rows : list[dict]
        Rows loaded from ``graph_analysis_results.csv``.

    Returns
    -------
    dict[str, dict[str, dict]]
        Mapping of ``base_name -> {dwi_type -> row_dict}``.
    """
    groups: dict[str, dict[str, dict]] = defaultdict(dict)
    for r in rows:
        dwi_type, base_name = parse_dwi_info(r["file_path"])
        groups[base_name][dwi_type] = r  # type: ignore
    return dict(groups)