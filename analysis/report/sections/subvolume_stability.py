"""Report section: sub-volume stability over fractions."""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


def _section_subvolume_stability(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build HTML section showing sub-volume stability (Dice vs Fx1 baseline).

    Parameters
    ----------
    mat_data : dict
        Mapping of DWI type to parsed MAT metrics dict.
    dwi_types : list[str]
        List of DWI types present in the data.

    Returns
    -------
    list[str]
        HTML chunks for the sub-volume stability section.
    """
    h: list[str] = []

    has_data = False
    for dwi in DWI_TYPES:
        if dwi in (mat_data or {}) and "subvolume_stability" in mat_data[dwi]:
            ss = mat_data[dwi]["subvolume_stability"]
            if ss and ss.get("method_names"):
                has_data = True
                break

    if not has_data:
        h.append(_h2("Sub-Volume Stability Over Fractions", "subvolume-stability"))
        h.append('<p class="meta">Sub-volume stability data not available. '
                 "Enable <code>run_subvolume_stability</code> in config to generate.</p>")
        return h

    h.append(_h2("Sub-Volume Stability Over Fractions", "subvolume-stability"))
    h.append(
        "<p>Dice coefficient between each fraction's core mask and the Fx1 "
        "baseline mask, averaged across patients. Methods with Dice &gt; 0.7 "
        "at Fx5 show spatially stable sub-volumes suitable for adaptive "
        "planning. Methods with rapidly declining Dice indicate unstable "
        "sub-volume definitions that may not track consistent tumor biology.</p>"
    )

    for dwi in DWI_TYPES:
        if dwi not in (mat_data or {}):
            continue
        ss = mat_data[dwi].get("subvolume_stability", {})
        if not ss or not ss.get("method_names"):
            continue

        methods = ss["method_names"]
        dice = ss.get("dice_vs_baseline", [])
        n_tp = ss.get("n_timepoints", 0)

        if not dice or not n_tp:
            continue

        # Build timepoint labels
        tp_labels = []
        for k in range(n_tp):
            if k < 5:
                tp_labels.append(f"Fx{k + 1}")
            else:
                tp_labels.append("Post")

        h.append(f"<h3>{_dwi_badge(dwi)} — Mean Dice vs Fx1 Baseline</h3>")
        h.append("<table><thead><tr>")
        h.append("<th>Core Method</th>")
        for lbl in tp_labels:
            h.append(f"<th>{_esc(lbl)}</th>")
        h.append("</tr></thead><tbody>")

        for m_idx, method in enumerate(methods):
            h.append(f"<tr><td><strong>{_esc(method)}</strong></td>")
            for k in range(n_tp):
                val = _get_mean_dice(dice, m_idx, k)
                h.append(_dice_cell(val))
            h.append("</tr>")

        h.append("</tbody></table>")

    h.append(
        '<p class="meta"><strong>Interpretation:</strong> Methods with Dice '
        "&ge; 0.7 at Fx5 show spatially stable sub-volumes suitable for "
        "adaptive planning.</p>"
    )

    return h


def _get_mean_dice(dice_3d: list, m_idx: int, k_idx: int) -> float | None:
    """Extract mean Dice from [nMethods x nTp x nPatients] nested list.

    Returns the mean across patients (axis 2) for method m_idx at timepoint k_idx.
    """
    try:
        row = dice_3d[m_idx]
        if not row:
            return None
        tp_data = row[k_idx]
        if isinstance(tp_data, list):
            vals = [v for v in tp_data if v is not None]
            if not vals:
                return None
            return sum(vals) / len(vals)
        if tp_data is None:
            return None
        return float(tp_data)
    except (IndexError, TypeError, ValueError):
        return None


def _dice_cell(val: float | None) -> str:
    """Format a Dice value as a colored table cell."""
    if val is None:
        return "<td>&mdash;</td>"
    if val >= 0.7:
        cls = ' style="background:#d4edda"'  # green
    elif val >= 0.5:
        cls = ' style="background:#fff3cd"'  # yellow
    else:
        cls = ' style="background:#f8d7da"'  # red
    return f"<td{cls}>{val:.3f}</td>"
