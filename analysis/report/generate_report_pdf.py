#!/usr/bin/env python3
"""Markdown-to-HTML and HTML-to-PDF helpers for the analysis report.

Split out from :mod:`report.generate_report`; the names are re-exported
there for backward compatibility.
"""

from __future__ import annotations

import sys
from pathlib import Path

import markdown  # type: ignore

from report.report_formatters import HTML_TEMPLATE  # type: ignore


def markdown_to_html(md_text: str, title: str) -> str:
    """Convert Markdown text to a styled standalone HTML document.

    Uses the ``python-markdown`` library with table, fenced-code, and
    table-of-contents extensions.  The result is wrapped in the
    :data:`report_formatters.HTML_TEMPLATE` for consistent styling.

    Parameters
    ----------
    md_text : str
        Raw Markdown content.
    title : str
        HTML ``<title>`` for the document.

    Returns
    -------
    str
        Complete HTML document string.
    """
    body = markdown.markdown(
        md_text,
        extensions=["tables", "fenced_code", "toc"],
    )
    return HTML_TEMPLATE.format(title=title, body=body)


def _html_to_pdf_weasyprint(html_content: str, pdf_path: Path) -> bool:
    """Attempt PDF generation via WeasyPrint. Returns True on success."""
    try:
        from weasyprint import HTML as WeasyprintHTML  # type: ignore[import-untyped]
        WeasyprintHTML(string=html_content, base_url=".").write_pdf(str(pdf_path))
        return True
    except ImportError:
        print("  WeasyPrint not installed (pip install weasyprint).")
        return False
    except OSError as e:
        # WeasyPrint is installed but its native GTK/Pango libraries are missing.
        # On Windows this requires installing the GTK3 runtime separately.
        print(f"  WeasyPrint native library error: {e}")
        if sys.platform == "win32":
            print("  On Windows, WeasyPrint requires the GTK3 runtime.")
            print("  Install it from: https://github.com/tschoonj/GTK-for-Windows-Runtime-Environment-Installer")
        return False
    except Exception as e:
        print(f"  WeasyPrint failed: {e}")
        return False


def _html_to_pdf_browser(html_content: str, pdf_path: Path) -> bool:
    """Attempt PDF generation via a headless Chrome/Edge/Chromium subprocess.

    Writes the HTML to a temporary file, invokes a headless browser with
    ``--print-to-pdf``, and cleans up the temp file afterwards.

    Returns True on success.
    """
    import subprocess
    import tempfile

    # Candidate browser executables to try, in order of preference.
    candidates: list[str] = []
    if sys.platform == "win32":
        import os
        # Common Chrome/Edge install paths on Windows.
        win_paths = [
            r"C:\Program Files\Google\Chrome\Application\chrome.exe",
            r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
            os.path.expanduser(r"~\AppData\Local\Google\Chrome\Application\chrome.exe"),
            r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
            r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
        ]
        candidates = [p for p in win_paths if Path(p).exists()]
        # Also try commands on PATH.
        candidates += ["chrome", "msedge", "chromium"]
    elif sys.platform == "darwin":
        candidates = [
            "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
            "/Applications/Chromium.app/Contents/MacOS/Chromium",
            "google-chrome",
            "chromium",
        ]
    else:
        candidates = ["google-chrome", "chromium-browser", "chromium", "google-chrome-stable"]

    if not candidates:
        return False

    # Write HTML to a temporary file so the browser can load it via file://.
    tmp = tempfile.NamedTemporaryFile(
        suffix=".html", delete=False, mode="w", encoding="utf-8"
    )
    try:
        tmp.write(html_content)
        tmp.close()
        tmp_path = Path(tmp.name)

        # file:// URI — use forward slashes even on Windows.
        file_uri = "file:///" + tmp_path.as_posix().lstrip("/")

        for browser in candidates:
            try:
                result = subprocess.run(
                    [
                        browser,
                        "--headless",
                        "--disable-gpu",
                        "--no-sandbox",
                        f"--print-to-pdf={pdf_path}",
                        "--print-to-pdf-no-header",
                        file_uri,
                    ],
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
                if result.returncode == 0 and pdf_path.exists() and pdf_path.stat().st_size > 0:
                    print(f"  Browser PDF engine: {browser}")
                    return True
            except (FileNotFoundError, subprocess.TimeoutExpired):
                continue
            except Exception:
                continue
    finally:
        try:
            Path(tmp.name).unlink(missing_ok=True)
        except Exception:
            pass

    return False


def html_to_pdf(html_content: str, pdf_path: Path) -> bool:
    """Convert an HTML report string to a PDF file.

    Tries WeasyPrint first; if that fails (e.g. missing GTK libraries on
    Windows), falls back to a headless Chrome/Edge/Chromium subprocess.

    Parameters
    ----------
    html_content : str
        Complete HTML document string.
    pdf_path : Path
        Output path for the PDF file.

    Returns
    -------
    bool
        True if PDF generation succeeded, False otherwise.
    """
    # Try WeasyPrint (best quality, but requires GTK native libs on Windows).
    print("  Trying WeasyPrint ...")
    if _html_to_pdf_weasyprint(html_content, pdf_path):
        return True

    # Fallback: headless Chrome/Edge/Chromium (no extra Python packages needed).
    print("  Trying headless browser fallback ...")
    if _html_to_pdf_browser(html_content, pdf_path):
        return True

    print("  PDF generation failed. To enable PDF output, either:")
    print("    1. Install GTK3 runtime for WeasyPrint (Windows):")
    print("       https://github.com/tschoonj/GTK-for-Windows-Runtime-Environment-Installer")
    print("    2. Install Google Chrome or Microsoft Edge.")
    print("  The HTML report is still available.")
    return False
