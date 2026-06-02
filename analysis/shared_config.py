#!/usr/bin/env python3
"""Configuration helpers for the pancData3 analysis suite.

This module holds the stateless configuration machinery used by
:mod:`shared`:

- Built-in defaults (``_DEFAULTS``) and the provider→key map (``_API_KEY_MAP``).
- API-key resolution (``get_api_key``).
- Layered config loading (``load_analysis_config``) and its JSON helpers
  (``_deep_merge``, ``_normalize_dwi_type``, ``_is_string_well_formed``,
  ``_strip_json_comments``, ``_permissive_json_loads``).
- The canonical ``DWI_TYPES`` constant.

The mutable config cache and the cached ``get_config`` / ``reset_config_cache``
accessors live in :mod:`shared` so that in-process cache mutation has a single
owner.  All names defined here are re-exported from :mod:`shared` for backward
compatibility.
"""

from __future__ import annotations

import json
import os
import re
import warnings
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
        "per_image_timeout_seconds": 300,
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
        # Lazy import to avoid an import-time cycle with shared.py, which
        # owns the cached get_config accessor.
        from shared import get_config
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


def _is_string_well_formed(text: str) -> bool:
    """Validate that JSON-like text has balanced quotes and proper escaping.

    This function performs a simplified check to determine if string
    literals in JSON-like text are properly formed before attempting
    comment stripping. It counts quote characters while respecting
    escape sequences to detect unbalanced strings.

    Returns
    -------
    bool
        True if strings appear well-formed, False otherwise.
    """
    in_string = False
    i = 0
    n = len(text)

    while i < n:
        char = text[i]
        if char == '"':
            # Check if this quote is escaped by counting preceding backslashes
            backslash_count = 0
            j = i - 1
            while j >= 0 and text[j] == '\\':
                backslash_count += 1
                j -= 1

            # If even number of backslashes (including 0), quote is not escaped
            if backslash_count % 2 == 0:
                in_string = not in_string
        i += 1

    # String should not still be open at end of text
    return not in_string


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

    .. note::

       This function validates string balance before processing to avoid
       corruption from malformed escape sequences. If validation fails,
       the original text is returned unchanged to prevent incorrect
       comment stripping that could corrupt valid JSON content.
    """
    # Validate string balance before proceeding
    if not _is_string_well_formed(text):
        return text

    result: list[str] = []
    i = 0
    n = len(text)
    while i < n:
        # Check if we're entering a JSON string literal
        if text[i] == '"':
            # Consume the entire string literal, including escaped characters.
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

    The function validates string balance before attempting comment
    stripping to prevent corruption from malformed escape sequences.
    If validation fails, the original text is parsed as-is to preserve
    the original error for debugging.

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

    # Use token-aware comment stripping with validation
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


# Canonical DWI type names in the order the pipeline processes them.
# Used throughout the suite for iteration and labelling.
# NOTE: This is kept as a module-level constant for backward compatibility;
# scripts that need the configurable list should use ``get_config()["dwi_types"]``.
DWI_TYPES = list(_DEFAULTS["dwi_types"])
