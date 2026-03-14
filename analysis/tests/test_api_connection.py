import os
import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path

# Important: We must skip this test if we don't have an API key!
# This test verifies that the google.genai client can successfully connect
# to the Gemini API and successfully return a response.

# We will mark it explicitly so it can be run via:
# pytest -m "api"

try:
    from google import genai
    from google.genai import types
    HAS_GENAI = True
except ImportError:
    HAS_GENAI = False

pytestmark = pytest.mark.skipif(
    not HAS_GENAI,
    reason="google-genai package is not installed."
)

@pytest.mark.api
def test_gemini_api_connection(tmp_path: Path):
    """
    Test the actual connection to the Gemini API.
    
    This requires a valid GEMINI_API_KEY environment variable. If the key
    is valid, the API should return a response without raising an exception.
    """
    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        pytest.skip("GEMINI_API_KEY environment variable not set. Skipping API test.")

    # Get the configured model from our shared config
    from shared import get_config
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
