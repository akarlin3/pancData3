#!/usr/bin/env python3
"""Vision-API integration for batch graph analysis (Gemini + Claude).

This module holds the async vision-model machinery: lazy SDK imports, the
shared rate-limit coordinator, JSON-repair / response-parsing helpers, and the
per-image ``analyze_image`` (Gemini) and ``analyze_image_claude`` (Claude)
coroutines.  Vision-model tuning constants (model names, concurrency, retry /
backoff / timeout settings) are also defined here, loaded from the centralised
analysis config.

Re-exported by :mod:`parsers.batch_graph_analysis` for backwards compatibility.
"""

from __future__ import annotations

import asyncio
import base64
import json
import re
import sys
from pathlib import Path

from tqdm import tqdm  # type: ignore

# Ensure analysis/ root is on sys.path so 'shared' is importable when this
# module is imported by the standalone-script entry point.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import get_config  # type: ignore

from parsers.batch_graph_schemas import (
    GraphAnalysis,
    SYSTEM_PROMPT,
    media_type_for,
)
from pydantic import ValidationError  # type: ignore

# Lazy-import google.genai and anthropic so that the module can be imported
# (and tested) without either package installed.  The actual imports happen
# inside ``analyze_image`` / ``analyze_image_claude`` and ``main``.
genai = None  # type: ignore[assignment]
types = None  # type: ignore[assignment]
anthropic_mod = None  # type: ignore[assignment]


def _ensure_genai():
    """Import google.genai lazily; raise ImportError if unavailable."""
    global genai, types
    if genai is None:
        from google import genai as _genai  # type: ignore
        from google.genai import types as _types  # type: ignore
        genai = _genai
        types = _types


def _ensure_anthropic():
    """Import anthropic lazily; raise ImportError if unavailable."""
    global anthropic_mod
    if anthropic_mod is None:
        import anthropic as _anthropic  # type: ignore
        anthropic_mod = _anthropic


# ── Vision model configuration ──────────────────────────────────────────────
# Loaded from the centralised analysis config (analysis_config.json /
# shared._DEFAULTS).  Values can be overridden without editing source code.
_vision_cfg = get_config()["vision"]

# Model to use for vision analysis.
GEMINI_MODEL = _vision_cfg["gemini_model"]
CLAUDE_MODEL = _vision_cfg.get("claude_model", "claude-opus-4-6")
VISION_PROVIDER = _vision_cfg.get("provider", "gemini")

# ── Async API call ───────────────────────────────────────────────────────────

# Concurrency limit for API requests.  Kept low (2) to stay well within
# Gemini's rate limits.
SEM_LIMIT = _vision_cfg["max_concurrent_requests"]

# Number of retries on rate-limit (429) errors.  Uses exponential
# backoff: 15s, 30s, 60s, 120s.
MAX_RETRIES = _vision_cfg["max_retries"]

# Base seconds for exponential backoff on rate-limit retries.
BACKOFF_BASE = _vision_cfg["backoff_base_seconds"]

# Maximum output tokens for the vision model response.
MAX_OUTPUT_TOKENS = _vision_cfg["max_output_tokens"]

# Per-request timeout in seconds.  If the Gemini API does not respond
# within this window the request is cancelled and treated as a failure.
REQUEST_TIMEOUT = _vision_cfg["request_timeout_seconds"]

# Per-image timeout in seconds.  Caps the total wall-clock time for a
# single image (including all retries and backoff).  Prevents one
# problematic image (e.g. blank parameter maps) from stalling the
# entire analysis run.
PER_IMAGE_TIMEOUT = _vision_cfg.get("per_image_timeout_seconds", 300)


class _RateLimitCoordinator:
    """Coordinates rate-limit backoff across all workers.

    When any worker hits a 429/rate-limit error it calls :meth:`signal`,
    which sets a shared cooldown timestamp.  Before each API call workers
    call :meth:`wait_if_cooling_down` to honour the cooldown.  This
    prevents the thundering-herd problem where multiple workers slam the
    API simultaneously after individual backoffs expire.
    """

    def __init__(self):
        self._resume_at: float = 0.0  # monotonic timestamp
        self._lock = asyncio.Lock()
        self._consecutive_hits: int = 0

    async def signal(self):
        """Record a rate-limit hit and extend the global cooldown."""
        async with self._lock:
            self._consecutive_hits += 1
            # Exponential cooldown: base * 2^(hits-1), capped at 4× base.
            wait = min(BACKOFF_BASE * (2 ** (self._consecutive_hits - 1)),
                       BACKOFF_BASE * 4)
            now = asyncio.get_event_loop().time()
            new_resume = now + wait
            if new_resume > self._resume_at:
                self._resume_at = new_resume
                print(f"  [RATE-LIMIT] Global cooldown: {wait:.0f}s "
                      f"(consecutive hits: {self._consecutive_hits})", flush=True)

    async def wait_if_cooling_down(self):
        """Sleep until the global cooldown expires (no-op if not active)."""
        now = asyncio.get_event_loop().time()
        if now < self._resume_at:
            delay = self._resume_at - now
            await asyncio.sleep(delay)

    async def record_success(self):
        """Reset the consecutive-hit counter after a successful request."""
        async with self._lock:
            self._consecutive_hits = 0


def _repair_truncated_json(raw: str) -> dict | None:
    """Attempt to repair JSON truncated by token limits.

    Progressively closes open strings, arrays, and objects from the end
    of the raw text.  Returns the parsed dict on success, or ``None``
    if the response is beyond salvage.
    """
    text = raw.rstrip()
    # Try closing open structures (up to a handful of attempts).
    closers = ['"', "]", "}", "]", "}"]
    for i in range(1, len(closers) + 1):
        candidate = text + "".join(closers[:i])  # type: ignore
        try:
            data = json.loads(candidate)
            if isinstance(data, dict):
                return data
        except json.JSONDecodeError:
            continue
    # Brute-force: strip back to last complete key-value, close object.
    # Find the last complete value (ends with , or a closing bracket).
    for trim in range(1, min(300, len(text))):
        stub = text[:-trim].rstrip().rstrip(",")  # type: ignore
        for suffix in ["}", "]}", "\"}", "\"]}", "\"]}",
                        "]}", "\"]}", "\"]}}"]:
            try:
                data = json.loads(stub + suffix)
                if isinstance(data, dict):
                    return data
            except json.JSONDecodeError:
                continue
    return None


def _is_rate_limit_error(exc: Exception) -> bool:
    """Check whether *exc* is a rate-limit / quota error."""
    err_str = str(exc).lower()
    return ("429" in err_str or "rate" in err_str
            or "resource" in err_str or "quota" in err_str)


async def analyze_image(
    client,
    image_path: Path,
    rate_limiter: _RateLimitCoordinator,
    progress_bar: tqdm | None = None,
) -> GraphAnalysis:
    """Send one image to Gemini and parse the structured JSON response.

    The function handles rate-limit retries (coordinated via a shared
    :class:`_RateLimitCoordinator`), markdown fence stripping, JSON
    parsing, and Pydantic validation.  On any non-recoverable error a
    fallback ``GraphAnalysis`` with ``graph_type="unknown"`` or ``"error"``
    is returned so that the batch can continue.

    Parameters
    ----------
    client : genai.Client
        Authenticated Google Gen AI client.
    image_path : Path
        Path to the graph image to analyse.
    rate_limiter : _RateLimitCoordinator
        Shared coordinator that pauses all workers on 429 errors.

    Returns
    -------
    GraphAnalysis
        Validated structured analysis, or a fallback object on failure.
    """
    image_bytes = await asyncio.to_thread(image_path.read_bytes)  # type: ignore
    mime = media_type_for(image_path)
    rel_path = str(image_path)

    if progress_bar is not None:
        progress_bar.set_postfix_str(f"analyzing {image_path.name}", refresh=True)

    # Retry loop with exponential backoff for rate-limit errors.
    response = None
    for attempt in range(MAX_RETRIES + 1):
        # Honour any global cooldown before sending a request.
        await rate_limiter.wait_if_cooling_down()

        try:
            response = await asyncio.wait_for(
                client.aio.models.generate_content(
                    model=GEMINI_MODEL,
                    contents=[
                        types.Part.from_bytes(  # type: ignore
                            data=image_bytes,
                            mime_type=mime,
                        ),
                        f"Analyze this graph image. File: {image_path.name}\n"
                        "Return ONLY the JSON object described in your instructions.",
                    ],
                    config=types.GenerateContentConfig(  # type: ignore
                        system_instruction=SYSTEM_PROMPT,
                        max_output_tokens=MAX_OUTPUT_TOKENS,
                    ),
                ),
                timeout=REQUEST_TIMEOUT,
            )
            await rate_limiter.record_success()
            break  # Success: exit the retry loop.
        except asyncio.TimeoutError:
            if attempt < MAX_RETRIES:
                wait = 2 ** attempt * BACKOFF_BASE
                print(f"  [TIMEOUT] {image_path.name}: retry {attempt+1}/{MAX_RETRIES} in {wait}s", flush=True)
                await asyncio.sleep(wait)
            else:
                raise  # Give up after exhausting all retries.
        except Exception as e:
            if _is_rate_limit_error(e) and attempt < MAX_RETRIES:
                await rate_limiter.signal()
                # Worker-local backoff on top of global cooldown.
                wait = 2 ** attempt * BACKOFF_BASE
                print(f"  [RATE-LIMIT] {image_path.name}: retry {attempt+1}/{MAX_RETRIES} in {wait}s", flush=True)
                await asyncio.sleep(wait)
            elif _is_rate_limit_error(e):
                raise  # Give up after exhausting all retries.
            else:
                raise  # Non-rate-limit errors are re-raised immediately.

    # Guard against the retry loop exhausting without a successful response.
    if response is None:
        raise RuntimeError(f"No response received for {image_path.name} after {MAX_RETRIES + 1} attempts")

    raw_text = response.text.strip()  # type: ignore

    # Strip markdown code fences if the model wraps its JSON response
    # (e.g. "```json\n{...}\n```").
    if raw_text.startswith("```"):
        raw_text = re.sub(r"^```(?:json)?\s*", "", raw_text)
        raw_text = re.sub(r"\s*```$", "", raw_text)

    # ── JSON parsing (with truncation repair) ──
    try:
        data = json.loads(raw_text)
    except json.JSONDecodeError:
        # Attempt to repair truncated JSON from token-limited responses.
        data = _repair_truncated_json(raw_text)
        if data is None:
            if progress_bar is not None:
                progress_bar.update(1)  # type: ignore
            else:
                print(f"  ❌  JSON parse error for {image_path.name}", flush=True)
            return GraphAnalysis(
                file_path=rel_path,  # type: ignore
                graph_type="unknown",  # type: ignore
                summary=f"JSON parse error (unrepairable). Raw response: {raw_text[:200]}",  # type: ignore
            )

    # Inject the file path (not returned by the model).
    data["file_path"] = rel_path  # type: ignore

    # ── Pydantic validation ──
    try:
        analysis = GraphAnalysis(**data)
    except ValidationError as e:
        if progress_bar is not None:
            progress_bar.update(1)
        else:
            print(f"  ⚠️  Validation warning for {image_path.name}: {e}", flush=True)
        # Fallback: preserve whatever fields parsed successfully.
        return GraphAnalysis(
            file_path=rel_path,  # type: ignore
            graph_type=data.get("graph_type", "unknown"),  # type: ignore
            summary=data.get("summary", f"Validation error: {e}"),  # type: ignore
        )

    if progress_bar is not None:
        progress_bar.set_postfix_str(image_path.name, refresh=True)
        progress_bar.update(1)
    return analysis


def _parse_vision_response(
    raw_text: str,
    image_path: Path,
    rel_path: str,
    progress_bar: tqdm | None = None,
) -> GraphAnalysis:
    """Parse a raw JSON response string into a validated GraphAnalysis.

    Shared by both Gemini and Claude response handlers.  Handles markdown
    fence stripping, truncation repair, and Pydantic validation.

    Parameters
    ----------
    raw_text : str
        Raw text response from the vision model.
    image_path : Path
        Path to the source image (for error messages).
    rel_path : str
        Relative file path to embed in the result.
    progress_bar : tqdm or None
        Optional progress bar to update on error.

    Returns
    -------
    GraphAnalysis
        Validated analysis or a fallback on parse/validation failure.
    """
    text = raw_text.strip()

    # Strip markdown code fences if the model wraps its JSON response.
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)

    # ── JSON parsing (with truncation repair) ──
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        data = _repair_truncated_json(text)
        if data is None:
            if progress_bar is not None:
                progress_bar.update(1)
            else:
                print(f"  ❌  JSON parse error for {image_path.name}", flush=True)
            return GraphAnalysis(
                file_path=rel_path,
                graph_type="unknown",
                summary=f"JSON parse error (unrepairable). Raw response: {text[:200]}",
            )

    data["file_path"] = rel_path

    # ── Pydantic validation ──
    try:
        analysis = GraphAnalysis(**data)
    except ValidationError as e:
        if progress_bar is not None:
            progress_bar.update(1)
        else:
            print(f"  ⚠️  Validation warning for {image_path.name}: {e}", flush=True)
        return GraphAnalysis(
            file_path=rel_path,
            graph_type=data.get("graph_type", "unknown"),
            summary=data.get("summary", f"Validation error: {e}"),
        )

    if progress_bar is not None:
        progress_bar.set_postfix_str(image_path.name, refresh=True)
        progress_bar.update(1)
    return analysis


# ── Claude API integration ──────────────────────────────────────────────────


def _is_claude_rate_limit_error(exc: Exception) -> bool:
    """Check whether *exc* is a Claude rate-limit error."""
    err_str = str(exc).lower()
    return ("429" in err_str or "rate_limit" in err_str
            or "rate limit" in err_str or "overloaded" in err_str)


async def analyze_image_claude(
    client,
    image_path: Path,
    rate_limiter: _RateLimitCoordinator,
    progress_bar: tqdm | None = None,
) -> GraphAnalysis:
    """Send one image to Claude and parse the structured JSON response.

    Mirrors :func:`analyze_image` but uses the Anthropic Messages API.
    Uses the same system prompt, Pydantic validation, and retry logic.

    Parameters
    ----------
    client : anthropic.AsyncAnthropic
        Authenticated Anthropic async client.
    image_path : Path
        Path to the graph image to analyse.
    rate_limiter : _RateLimitCoordinator
        Shared coordinator that pauses all workers on rate-limit errors.
    progress_bar : tqdm or None
        Optional progress bar for UI updates.

    Returns
    -------
    GraphAnalysis
        Validated structured analysis, or a fallback object on failure.
    """
    image_bytes = await asyncio.to_thread(image_path.read_bytes)
    mime = media_type_for(image_path)
    rel_path = str(image_path)
    b64_data = base64.standard_b64encode(image_bytes).decode("utf-8")

    if progress_bar is not None:
        progress_bar.set_postfix_str(f"analyzing {image_path.name}", refresh=True)

    response = None
    for attempt in range(MAX_RETRIES + 1):
        await rate_limiter.wait_if_cooling_down()

        try:
            response = await asyncio.wait_for(
                client.messages.create(
                    model=CLAUDE_MODEL,
                    max_tokens=MAX_OUTPUT_TOKENS,
                    system=SYSTEM_PROMPT,
                    messages=[
                        {
                            "role": "user",
                            "content": [
                                {
                                    "type": "image",
                                    "source": {
                                        "type": "base64",
                                        "media_type": mime,
                                        "data": b64_data,
                                    },
                                },
                                {
                                    "type": "text",
                                    "text": (
                                        f"Analyze this graph image. File: {image_path.name}\n"
                                        "Return ONLY the JSON object described in your instructions."
                                    ),
                                },
                            ],
                        }
                    ],
                ),
                timeout=REQUEST_TIMEOUT,
            )
            await rate_limiter.record_success()
            break
        except asyncio.TimeoutError:
            if attempt < MAX_RETRIES:
                wait = 2 ** attempt * BACKOFF_BASE
                print(f"  [TIMEOUT] {image_path.name}: retry {attempt+1}/{MAX_RETRIES} in {wait}s", flush=True)
                await asyncio.sleep(wait)
            else:
                raise
        except Exception as e:
            if _is_claude_rate_limit_error(e) and attempt < MAX_RETRIES:
                await rate_limiter.signal()
                wait = 2 ** attempt * BACKOFF_BASE
                print(f"  [RATE-LIMIT] {image_path.name}: retry {attempt+1}/{MAX_RETRIES} in {wait}s", flush=True)
                await asyncio.sleep(wait)
            elif _is_claude_rate_limit_error(e):
                raise
            else:
                raise

    if response is None:
        raise RuntimeError(f"No response received for {image_path.name} after {MAX_RETRIES + 1} attempts")

    # Extract text from Claude's response content blocks.
    raw_text = ""
    for block in response.content:
        if hasattr(block, "text"):
            raw_text += block.text

    return _parse_vision_response(raw_text, image_path, rel_path, progress_bar=None)
