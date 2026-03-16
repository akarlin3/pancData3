"""Tests for figure gallery section in report generation.

Covers:
- _section_figure_gallery: embedded PNG images, DWI type grouping
"""

from __future__ import annotations

import struct
import zlib
from pathlib import Path

import pytest  # type: ignore

from report.generate_report import (  # type: ignore
    _section_figure_gallery,
    reset_numbering,
)


def _make_tiny_png() -> bytes:
    """Create a minimal valid 1x1 pixel PNG image."""
    raw_data = b'\x00\x00\x00\x00'
    compressed = zlib.compress(raw_data)
    def chunk(ctype: bytes, data: bytes) -> bytes:
        c = ctype + data
        return struct.pack('>I', len(data)) + c + struct.pack('>I', zlib.crc32(c) & 0xffffffff)
    return (b'\x89PNG\r\n\x1a\n' +
            chunk(b'IHDR', struct.pack('>IIBBBBB', 1, 1, 8, 2, 0, 0, 0)) +
            chunk(b'IDAT', compressed) +
            chunk(b'IEND', b''))


class TestFigureGallery:
    """Verify the figure gallery section with embedded images."""

    def test_empty_when_no_images(self, saved_files_dir: Path):
        """Gallery is empty when no images exist."""
        reset_numbering()
        result = _section_figure_gallery(saved_files_dir)
        assert result == []

    def test_embeds_png_images(self, saved_files_dir: Path):
        """Gallery embeds PNG images as base64."""
        reset_numbering()
        img_dir = saved_files_dir / "Standard"
        img_path = img_dir / "test_figure.png"
        img_path.write_bytes(_make_tiny_png())

        result = _section_figure_gallery(saved_files_dir)
        html = "\n".join(result)
        assert "Figure Gallery" in html
        assert "data:image/png;base64," in html
        assert "Test Figure" in html  # filename -> title

    def test_groups_by_dwi_type(self, saved_files_dir: Path):
        """Gallery groups figures by DWI type."""
        reset_numbering()
        png = _make_tiny_png()
        (saved_files_dir / "Standard" / "fig1.png").write_bytes(png)
        (saved_files_dir / "dnCNN" / "fig2.png").write_bytes(png)

        result = _section_figure_gallery(saved_files_dir)
        html = "\n".join(result)
        assert "badge-standard" in html
        assert "badge-dncnn" in html
