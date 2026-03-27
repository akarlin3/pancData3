"""Report section: cross-pipeline Dice coefficients at fraction 1."""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


def _section_cross_pipeline_dice(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build HTML section showing cross-pipeline Dice agreement at Fx1.

    Displays a table of mean Dice +/- SD for each of the 11 core methods
    across all 3 pipeline pairs (Std vs DnCNN, Std vs IVIMNet, DnCNN vs
    IVIMNet).  Cells with mean Dice < 0.5 are highlighted in red.

    Parameters
    ----------
    mat_data : dict
        Mapping of DWI type to parsed MAT metrics dict.
    dwi_types : list[str]
        List of DWI types present in the data.

    Returns
    -------
    list[str]
        HTML chunks for the cross-pipeline Dice section.
    """
    h: list[str] = []

    # Check if any DWI type has cross_pipeline_dice data
    has_data = False
    for dwi in DWI_TYPES:
        if dwi in (mat_data or {}) and "cross_pipeline_dice" in mat_data[dwi]:
            cpd = mat_data[dwi]["cross_pipeline_dice"]
            if cpd and cpd.get("methods"):
                has_data = True
                break

    if not has_data:
        h.append(_h2("Cross-Pipeline Dice Coefficients", "cross-pipeline-dice"))
        h.append('<p class="meta">Cross-pipeline Dice data not available. '
                 "Enable <code>run_cross_pipeline_dice</code> in config to generate.</p>")
        return h

    h.append(_h2("Cross-Pipeline Dice Coefficients (Fx1)", "cross-pipeline-dice"))
    h.append(
        "<p>Pairwise Dice coefficients between the 3 DWI processing pipelines "
        "(Standard, DnCNN, IVIMNet) for each tumor core delineation method at "
        "fraction 1. High Dice (&gt;0.7) indicates the core location is robust "
        "to pipeline choice; low Dice (&lt;0.5) suggests the denoising/fitting "
        "method materially affects which voxels are classified as resistant.</p>"
    )

    for dwi in DWI_TYPES:
        if dwi not in (mat_data or {}):
            continue
        cpd = mat_data[dwi].get("cross_pipeline_dice", {})
        if not cpd or not cpd.get("methods"):
            continue

        methods = cpd["methods"]
        pair_labels = cpd.get("pair_labels", ["Std vs DnCNN", "Std vs IVIMNet", "DnCNN vs IVIMNet"])
        mean_dice = cpd.get("mean_dice", [])
        std_dice = cpd.get("std_dice", [])

        h.append(f"<h3>{_dwi_badge(dwi)} &mdash; Cross-Pipeline Agreement</h3>")
        h.append("<table><thead><tr><th>Core Method</th>")
        for lbl in pair_labels:
            h.append(f"<th>{_esc(lbl)}</th>")
        h.append("</tr></thead><tbody>")

        for i, method in enumerate(methods):
            h.append(f"<tr><td><strong>{_esc(method)}</strong></td>")
            for j in range(len(pair_labels)):
                mean_val = _get_val(mean_dice, i, j)
                std_val = _get_val(std_dice, i, j)

                if mean_val is None:
                    h.append("<td>&mdash;</td>")
                else:
                    cell_cls = ' class="sig-strong"' if mean_val < 0.5 else ""
                    if std_val is not None:
                        h.append(f"<td{cell_cls}>{mean_val:.2f} &plusmn; {std_val:.2f}</td>")
                    else:
                        h.append(f"<td{cell_cls}>{mean_val:.2f}</td>")
            h.append("</tr>")

        h.append("</tbody></table>")

    return h


def _get_val(matrix: list | None, row: int, col: int):
    """Safely extract a value from a 2D nested list, returning None on failure."""
    if not matrix:
        return None
    try:
        val = matrix[row][col]
        if val is None:
            return None
        return float(val)
    except (IndexError, TypeError, ValueError):
        return None
