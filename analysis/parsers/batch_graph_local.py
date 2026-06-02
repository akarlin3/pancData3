#!/usr/bin/env python3
"""Local filename-based fallback analyser for batch graph analysis.

When the configured vision API (Gemini / Claude) is unavailable, this module
infers graph metadata from the structured filenames produced by the MATLAB
pipeline.  The results are lower-fidelity (no visual analysis) but keep the
downstream CSV pipeline functional for report generation and cross-DWI
comparison.

Re-exported by :mod:`parsers.batch_graph_analysis` for backwards compatibility.
"""

from __future__ import annotations

from pathlib import Path

from parsers.batch_graph_schemas import Axis, GraphAnalysis


# ── Local filename-based fallback analyzer ───────────────────────────────────
# When the Gemini API is unavailable, this function infers graph metadata
# from the structured filenames produced by the MATLAB pipeline.  The
# results are lower-fidelity (no visual analysis) but keep the downstream
# CSV pipeline functional for report generation and cross-DWI comparison.

# Filename keywords → graph_type mapping.  Order matters: first match wins.
_GRAPH_TYPE_KEYWORDS: list[tuple[list[str], str]] = [
    (["heatmap", "dice_heatmap", "hausdorff_heatmap"], "heatmap"),
    (["parameter_map", "param_map", "overlay"], "parameter_map"),
    (["histogram", "hist"], "histogram"),
    (["boxplot", "box_plot", "box"], "box"),
    (["scatter", "correlation"], "scatter"),
    (["violin"], "violin"),
    (["bar", "volume_comparison"], "bar"),
    (["longitudinal", "trajectory", "timeseries", "time_series"], "line"),
    (["kaplan", "survival", "km_curve"], "line"),
    (["roc", "auc"], "line"),
    (["dose_vs", "dvh"], "line"),
]

# Filename keywords → axis label / units hints.
_AXIS_HINTS: dict[str, dict[str, str | None]] = {
    "adc": {"label": "ADC", "units": "mm²/s"},
    "ivim": {"label": "IVIM Parameter", "units": None},
    "d_mean": {"label": "D (mean)", "units": "mm²/s"},
    "f_mean": {"label": "f (mean)", "units": None},
    "dstar": {"label": "D* (mean)", "units": "mm²/s"},
    "dose": {"label": "Dose", "units": "Gy"},
    "time": {"label": "Time", "units": "days"},
    "longitudinal": {"label": "Timepoint", "units": None},
    "survival": {"label": "Time", "units": "days"},
    "kaplan": {"label": "Time", "units": "days"},
    "dice": {"label": "Method", "units": None},
    "hausdorff": {"label": "Method", "units": None},
    "volume": {"label": "Method", "units": None},
}

# Filename keywords → comparison_type hints.
_COMPARISON_HINTS: dict[str, str] = {
    "longitudinal": "longitudinal",
    "trajectory": "longitudinal",
    "byoutcome": "unpaired",
    "lf_vs_lc": "unpaired",
    "dose_vs": "dose-response",
    "dvh": "dose-response",
    "repeat": "paired",
    "baseline": "cross-sectional",
}


def _infer_graph_type(name_lower: str) -> str:
    """Infer graph_type from a lowercased filename stem."""
    for keywords, gtype in _GRAPH_TYPE_KEYWORDS:
        for kw in keywords:
            if kw in name_lower:
                return gtype
    return "unknown"


def _infer_axes(name_lower: str) -> tuple[Axis | None, Axis | None]:
    """Infer x/y axis labels from filename keywords."""
    x_axis = None
    y_axis = None
    for keyword, hints in _AXIS_HINTS.items():
        if keyword in name_lower:
            # For time-related x-axes, set x_axis
            if keyword in ("time", "longitudinal", "survival", "kaplan"):
                x_axis = Axis(
                    label=hints["label"],
                    units=hints.get("units"),
                )
            else:
                # For metric keywords, set as y_axis
                y_axis = Axis(
                    label=hints["label"],
                    units=hints.get("units"),
                )
    return x_axis, y_axis


def _infer_comparison_type(name_lower: str) -> str | None:
    """Infer comparison_type from filename keywords."""
    for keyword, ctype in _COMPARISON_HINTS.items():
        if keyword in name_lower:
            return ctype
    return None


def analyze_image_local(image_path: Path) -> GraphAnalysis:
    """Produce a best-effort ``GraphAnalysis`` from filename heuristics.

    This function does **not** inspect the image pixels.  It uses the
    structured naming conventions of the MATLAB pipeline
    (e.g. ``Longitudinal_Mean_Metrics_Standard.png``,
    ``core_method_dice_heatmap.png``) to infer graph type, axis labels,
    and comparison type.

    Parameters
    ----------
    image_path : Path
        Path to the graph image.

    Returns
    -------
    GraphAnalysis
        A structured analysis with ``graph_type``, axis hints, and a
        summary noting that the result was produced by the local fallback.
    """
    rel_path = str(image_path)
    stem = image_path.stem  # filename without extension
    name_lower = stem.lower()

    graph_type = _infer_graph_type(name_lower)
    x_axis, y_axis = _infer_axes(name_lower)
    comparison_type = _infer_comparison_type(name_lower)

    # Build a human-readable summary from the filename.
    pretty_name = stem.replace("_", " ")
    summary = (
        f"Local fallback analysis of '{pretty_name}' "
        f"(inferred type: {graph_type}). "
        "No visual analysis was performed; metadata was inferred from the "
        "filename. Re-run with GEMINI_API_KEY set for full vision analysis."
    )

    # Detect DWI type from path for clinical context.
    dwi_types = {"Standard", "dnCNN", "IVIMnet"}
    clinical_note = None
    for part in image_path.parts:
        if part in dwi_types:
            clinical_note = f"Graph from {part} DWI processing pipeline."
            break

    return GraphAnalysis(
        file_path=rel_path,
        graph_title=pretty_name,
        graph_type=graph_type,
        x_axis=x_axis,
        y_axis=y_axis,
        summary=summary,
        comparison_type=comparison_type,
        clinical_relevance=clinical_note,
        figure_quality="unknown",
    )
