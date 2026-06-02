#!/usr/bin/env python3
"""Dual-provider comparison CSV for batch graph analysis.

When the vision analysis runs in ``--provider both`` mode, both Gemini and
Claude analyse every image.  This module compares the two providers'
:class:`GraphAnalysis` results per image and writes a
``graph_analysis_comparison.csv`` noting where they agree and differ.

Re-exported by :mod:`parsers.batch_graph_analysis` for backwards compatibility.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from parsers.batch_graph_schemas import GraphAnalysis


# ── Comparison CSV for dual-provider mode ───────────────────────────────────

COMPARISON_COLUMNS = [
    "file_path",
    "graph_type_gemini",
    "graph_type_claude",
    "graph_type_match",
    "num_trends_gemini",
    "num_trends_claude",
    "num_statistical_tests_gemini",
    "num_statistical_tests_claude",
    "summary_gemini",
    "summary_claude",
    "clinical_relevance_gemini",
    "clinical_relevance_claude",
    "figure_quality_gemini",
    "figure_quality_claude",
    "x_axis_label_gemini",
    "x_axis_label_claude",
    "y_axis_label_gemini",
    "y_axis_label_claude",
    "trend_directions_gemini",
    "trend_directions_claude",
    "trend_directions_match",
    "differences",
]


def _compare_results(
    gemini: GraphAnalysis,
    claude: GraphAnalysis,
) -> dict:
    """Compare Gemini and Claude results for one image, returning a flat dict.

    Parameters
    ----------
    gemini : GraphAnalysis
        Result from the Gemini API.
    claude : GraphAnalysis
        Result from the Claude API.

    Returns
    -------
    dict
        Flat dictionary with comparison columns.
    """
    diffs: list[str] = []

    # Graph type comparison
    type_match = gemini.graph_type == claude.graph_type
    if not type_match:
        diffs.append(f"graph_type: Gemini={gemini.graph_type}, Claude={claude.graph_type}")

    # Trend direction comparison
    gem_dirs = sorted(t.direction for t in gemini.trends)
    cla_dirs = sorted(t.direction for t in claude.trends)
    dirs_match = gem_dirs == cla_dirs
    if not dirs_match:
        diffs.append(f"trend_directions differ")
    if len(gemini.trends) != len(claude.trends):
        diffs.append(f"num_trends: Gemini={len(gemini.trends)}, Claude={len(claude.trends)}")

    # Statistical test count
    if len(gemini.statistical_tests) != len(claude.statistical_tests):
        diffs.append(
            f"num_stat_tests: Gemini={len(gemini.statistical_tests)}, "
            f"Claude={len(claude.statistical_tests)}"
        )

    # Axis labels
    gem_x = (gemini.x_axis.label if gemini.x_axis else "")
    cla_x = (claude.x_axis.label if claude.x_axis else "")
    if gem_x.lower() != cla_x.lower() and gem_x and cla_x:
        diffs.append(f"x_axis_label: Gemini={gem_x!r}, Claude={cla_x!r}")

    gem_y = (gemini.y_axis.label if gemini.y_axis else "")
    cla_y = (claude.y_axis.label if claude.y_axis else "")
    if gem_y.lower() != cla_y.lower() and gem_y and cla_y:
        diffs.append(f"y_axis_label: Gemini={gem_y!r}, Claude={cla_y!r}")

    # Figure quality
    if (gemini.figure_quality or "") != (claude.figure_quality or ""):
        diffs.append(
            f"figure_quality: Gemini={gemini.figure_quality}, Claude={claude.figure_quality}"
        )

    return {
        "file_path": gemini.file_path,
        "graph_type_gemini": gemini.graph_type,
        "graph_type_claude": claude.graph_type,
        "graph_type_match": "yes" if type_match else "no",
        "num_trends_gemini": len(gemini.trends),
        "num_trends_claude": len(claude.trends),
        "num_statistical_tests_gemini": len(gemini.statistical_tests),
        "num_statistical_tests_claude": len(claude.statistical_tests),
        "summary_gemini": gemini.summary[:300] if gemini.summary else "",
        "summary_claude": claude.summary[:300] if claude.summary else "",
        "clinical_relevance_gemini": gemini.clinical_relevance or "",
        "clinical_relevance_claude": claude.clinical_relevance or "",
        "figure_quality_gemini": gemini.figure_quality or "",
        "figure_quality_claude": claude.figure_quality or "",
        "x_axis_label_gemini": gem_x,
        "x_axis_label_claude": cla_x,
        "y_axis_label_gemini": gem_y,
        "y_axis_label_claude": cla_y,
        "trend_directions_gemini": json.dumps(gem_dirs),
        "trend_directions_claude": json.dumps(cla_dirs),
        "trend_directions_match": "yes" if dirs_match else "no",
        "differences": "; ".join(diffs) if diffs else "none",
    }


def _write_comparison_csv(
    folder: Path,
    gemini_results: list[GraphAnalysis],
    claude_results: list[GraphAnalysis],
) -> Path:
    """Write a comparison CSV for dual-provider mode.

    Parameters
    ----------
    folder : Path
        Output folder.
    gemini_results : list[GraphAnalysis]
        Results from Gemini (one per image, ordered by image index).
    claude_results : list[GraphAnalysis]
        Results from Claude (same order).

    Returns
    -------
    Path
        Path to the written comparison CSV.
    """
    out_csv = folder / "graph_analysis_comparison.csv"
    rows = []
    agree_count = 0
    for gem, cla in zip(gemini_results, claude_results):
        row = _compare_results(gem, cla)
        rows.append(row)
        if row["differences"] == "none":
            agree_count += 1

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=COMPARISON_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    total = len(rows)
    pct = (agree_count / total * 100) if total else 0
    print(f"\n\U0001f50d Comparison: {agree_count}/{total} images fully agree ({pct:.1f}%)")
    print(f"\U0001f4c1 Comparison CSV: {out_csv}")

    # Print summary of most common differences
    diff_counts: dict[str, int] = {}
    for row in rows:
        if row["differences"] != "none":
            for d in row["differences"].split("; "):
                key = d.split(":")[0] if ":" in d else d
                diff_counts[key] = diff_counts.get(key, 0) + 1
    if diff_counts:
        print("   Most common differences:")
        for diff_type, count in sorted(diff_counts.items(), key=lambda x: -x[1])[:5]:
            print(f"     {diff_type}: {count} images")

    return out_csv
