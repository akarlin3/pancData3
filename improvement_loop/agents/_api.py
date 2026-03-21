# _api.py
"""Shared API retry helper for improvement-loop agents."""

import time

import anthropic  # type: ignore

from improvement_loop.loop_config import get_config as _get_loop_config


def _get_client() -> anthropic.Anthropic:
    """Return an Anthropic client, using the config API key if set."""
    cfg = _get_loop_config()
    kwargs = {}
    if cfg.anthropic_api_key:
        kwargs["api_key"] = cfg.anthropic_api_key
    return anthropic.Anthropic(**kwargs)


def api_call_with_retry(create_kwargs: dict) -> str:
    """Call client.messages.create with rate-limit retry. Returns response text.

    Uses streaming to avoid the SDK's 10-minute timeout guard for large
    max_tokens values.
    """
    cfg = _get_loop_config()
    client = _get_client()
    for attempt in range(1, cfg.max_api_retries + 1):
        try:
            # Use streaming to avoid SDK timeout for large max_tokens
            collected: list[str] = []
            with client.messages.stream(**create_kwargs) as stream:
                for text in stream.text_stream:
                    collected.append(text)
            return "".join(collected)
        except anthropic.RateLimitError:
            delay = cfg.retry_base_delay * attempt
            print(f"    Rate limited (attempt {attempt}/{cfg.max_api_retries}), "
                  f"waiting {delay:.0f}s...")
            time.sleep(delay)
            if attempt == cfg.max_api_retries:
                raise
        except anthropic.APIError:
            raise
    return ""  # unreachable
