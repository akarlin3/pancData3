"""Report section: Fx1 spatial repeatability Dice coefficients."""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


_PARAM_LABELS = [
    ("adc", "ADC"),
    ("d", "D"),
    ("f", "f"),
    ("dstar", "D*"),
]


def _dice_cell(val: float | None) -> str:
    """Format a mean Dice value as a coloured cell."""
    if val is None:
        return "<td>&mdash;</td>"
    if val >= 0.7:
        style = ' style="background:#d4edda"'  # green — reproducible
    elif val >= 0.5:
        style = ' style="background:#fff3cd"'  # yellow — partial
    else:
        style = ' style="background:#f8d7da"'  # red — stochastic
    return f"<td{style}>{val:.2f}</td>"


def _section_repeatability_dice(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build an HTML section showing Fx1 repeat Dice per parameter and DWI type.

    Parameters
    ----------
    mat_data : dict
        Mapping of DWI type name to parsed MAT metrics dict.
    dwi_types : list[str]
        List of DWI types present in the data (unused; kept for signature
        parity with other statistics sections).

    Returns
    -------
    list[str]
        HTML chunks for the repeat Dice section, or empty if no data found.
    """
    if not mat_data:
        return []

    # Find which DWI types have repeatability data with at least one
    # parameter containing a real (non-NaN) sample.
    dwi_with_data: list[str] = []
    for dwi in DWI_TYPES:
        entry = (mat_data or {}).get(dwi, {})
        rpt = entry.get("repeatability_dice") if entry else None
        if not rpt:
            continue
        any_samples = any(
            (info or {}).get("n", 0) > 0 and (info or {}).get("mean") is not None
            for info in rpt.values()
        )
        if any_samples:
            dwi_with_data.append(dwi)

    if not dwi_with_data:
        return []

    h: list[str] = []
    h.append(_h2("Spatial Repeatability (Fx1 Repeat Dice)", "repeat-dice"))
    h.append(
        "<p>Mean Dice coefficient between threshold-defined sub-volumes from "
        "paired Fx1 repeat scans. Parameters with repeat Dice "
        "&ge; 0.7 have spatially reproducible sub-volumes suitable for "
        "adaptive planning. Lower values indicate stochastic threshold "
        "behavior.</p>"
    )

    # Table: rows = parameters, columns = DWI types.
    h.append('<table><thead><tr>')
    h.append("<th>Parameter</th>")
    for dwi in dwi_with_data:
        h.append(f"<th>{_dwi_badge(dwi)}</th>")
    h.append("</tr></thead><tbody>")

    for param_key, param_label in _PARAM_LABELS:
        h.append(f"<tr><td><strong>{_esc(param_label)}</strong></td>")
        for dwi in dwi_with_data:
            rpt = mat_data[dwi].get("repeatability_dice", {})
            info = rpt.get(param_key)
            if not info or info.get("n", 0) == 0 or info.get("mean") is None:
                h.append("<td>&mdash;</td>")
                continue
            mean = info.get("mean")
            std = info.get("std")
            n = info.get("n", 0)
            try:
                mean_f = float(mean)
            except (TypeError, ValueError):
                h.append("<td>&mdash;</td>")
                continue
            cell = _dice_cell(mean_f)
            # Inject SD and N into the cell content.
            if std is None:
                inner = f"{mean_f:.2f} (n={n})"
            else:
                try:
                    std_f = float(std)
                    inner = f"{mean_f:.2f} &plusmn; {std_f:.2f} (n={n})"
                except (TypeError, ValueError):
                    inner = f"{mean_f:.2f} (n={n})"
            # Replace the cell's text with mean ± sd (n=…) while keeping colour.
            cell = cell.replace(f">{mean_f:.2f}<", f">{inner}<")
            h.append(cell)
        h.append("</tr>")

    h.append("</tbody></table>")

    h.append(
        "<p class=\"meta\">Interpretation: Dice &ge; 0.70 = reproducible (green); "
        "0.50&ndash;0.70 = partial reproducibility (yellow); "
        "&lt; 0.50 = stochastic sub-volume definition (red).</p>"
    )

    return h
