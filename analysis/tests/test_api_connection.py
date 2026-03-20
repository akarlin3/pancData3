import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path

# Important: We must skip these tests if we don't have the required API keys!
# These tests verify that the google.genai and anthropic clients can
# successfully connect to their respective APIs and return a response.

# They can be run via: pytest -m "api"

try:
    from google import genai
    from google.genai import types
    HAS_GENAI = True
except ImportError:
    HAS_GENAI = False

try:
    import anthropic
    HAS_ANTHROPIC = True
except ImportError:
    HAS_ANTHROPIC = False


@pytest.mark.api
@pytest.mark.skipif(not HAS_GENAI, reason="google-genai package is not installed.")
def test_gemini_api_connection(tmp_path: Path):
    """
    Test the actual connection to the Gemini API.

    This requires a valid GEMINI_API_KEY environment variable. If the key
    is valid, the API should return a response without raising an exception.
    """
    from shared import get_api_key, get_config
    api_key = get_api_key("gemini")
    if not api_key:
        pytest.skip("Gemini API key not configured. Skipping API test.")

    cfg = get_config()
    model = cfg["vision"]["gemini_model"]

    # We use the synchronous client here since pytest-asyncio is not installed
    client = genai.Client(api_key=api_key)

    # Create a simple tiny 1x1 pixel PNG to send, rather than a full graph
    # This minimizes bandwidth and token usage.
    tiny_png = b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01\x08\x06\x00\x00\x00\x1f\x15\xc4\x89\x00\x00\x00\nIDATx\x9cc\x00\x01\x00\x00\x05\x00\x01\r\n-\xb4\x00\x00\x00\x00IEND\xaeB`\x82'

    try:
        # We send a very simple prompt to the model and ask for a tiny response
        response = client.models.generate_content(
            model=model,
            contents=[
                types.Part.from_bytes(data=tiny_png, mime_type="image/png"),
                "Respond with the exact word 'SUCCESS'. Nothing else.",
            ]
        )

        # If we got a response without an exception, the connection works!
        assert response is not None
        assert response.text is not None

        # We can loosely check if the model understood our simple prompt
        assert "SUCCESS" in response.text.upper()

    except Exception as e:
        err_str = str(e).lower()
        if "429" in err_str or "resource_exhausted" in err_str or "quota" in err_str or "rate" in err_str:
            pytest.skip(f"API rate limit reached (quota exhausted): {e}")
        pytest.fail(f"API connection failed with exception: {e}")


@pytest.mark.api
@pytest.mark.skipif(not HAS_ANTHROPIC, reason="anthropic package is not installed.")
def test_claude_api_connection(tmp_path: Path):
    """
    Test the actual connection to the Claude API.

    This requires a valid ANTHROPIC_API_KEY environment variable. If the key
    is valid, the API should return a response without raising an exception.
    """
    from shared import get_api_key, get_config
    api_key = get_api_key("claude")
    if not api_key:
        pytest.skip("Anthropic API key not configured. Skipping API test.")

    cfg = get_config()
    model = cfg["vision"].get("claude_model", "claude-opus-4-6")

    client = anthropic.Anthropic(api_key=api_key)

    # Create a simple tiny 1x1 pixel PNG
    import base64
    tiny_png = b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01\x08\x06\x00\x00\x00\x1f\x15\xc4\x89\x00\x00\x00\nIDATx\x9cc\x00\x01\x00\x00\x05\x00\x01\r\n-\xb4\x00\x00\x00\x00IEND\xaeB`\x82'
    b64_data = base64.standard_b64encode(tiny_png).decode("utf-8")

    try:
        response = client.messages.create(
            model=model,
            max_tokens=50,
            messages=[
                {
                    "role": "user",
                    "content": [
                        {
                            "type": "image",
                            "source": {
                                "type": "base64",
                                "media_type": "image/png",
                                "data": b64_data,
                            },
                        },
                        {
                            "type": "text",
                            "text": "Respond with the exact word 'SUCCESS'. Nothing else.",
                        },
                    ],
                }
            ],
        )

        assert response is not None
        # Extract text from content blocks
        text = "".join(b.text for b in response.content if hasattr(b, "text"))
        assert "SUCCESS" in text.upper()

    except Exception as e:
        err_str = str(e).lower()
        if "429" in err_str or "rate_limit" in err_str or "overloaded" in err_str:
            pytest.skip(f"API rate limit reached: {e}")
        pytest.fail(f"Claude API connection failed with exception: {e}")
