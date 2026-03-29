"""Report section: per-method Coefficient of Reproducibility."""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


def _section_per_method_cor(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build HTML section showing per-method CoR from Fx1 repeat scans.

    Parameters
    ----------
    mat_data : dict
        Mapping of DWI type to parsed MAT metrics dict.
    dwi_types : list[str]
        List of DWI types present in the data.

    Returns
    -------
    list[str]
        HTML chunks for the per-method CoR section.
    """
    h: list[str] = []

    has_data = False
    for dwi in DWI_TYPES:
        if dwi in (mat_data or {}) and "per_method_cor" in mat_data[dwi]:
            cr = mat_data[dwi]["per_method_cor"]
            if cr and cr.get("method_names"):
                has_data = True
                break

    if not has_data:
        h.append(_h2("Per-Method Coefficient of Reproducibility", "per-method-cor"))
        h.append('<p class="meta">Per-method CoR data not available. '
                 "Enable <code>run_per_method_cor</code> in config to generate.</p>")
        return h

    h.append(_h2("Per-Method Coefficient of Reproducibility", "per-method-cor"))
    h.append(
        "<p>Coefficient of Reproducibility (CoR) for each core method's "
        "sub-volume fraction, derived from Fx1 repeat scans. CoR = 1.96 "
        "&times; &radic;2 &times; median within-subject CV &times; 100. "
        "Methods with CoR &lt; 15% have reproducible sub-volume definitions "
        "suitable for longitudinal monitoring. Methods with CoR &gt; 30% "
        "should be interpreted with caution.</p>"
    )

    for dwi in DWI_TYPES:
        if dwi not in (mat_data or {}):
            continue
        cr = mat_data[dwi].get("per_method_cor", {})
        if not cr or not cr.get("method_names"):
            continue

        methods = cr["method_names"]
        median_wcv = cr.get("median_wcv", [])
        cor_vals = cr.get("cor", [])
        n_pts = cr.get("n_patients_with_repeats", 0)

        h.append(f"<h3>{_dwi_badge(dwi)} — CoR (n={n_pts} patients with repeats)</h3>")
        h.append("<table><thead><tr>")
        h.append("<th>Core Method</th><th>Median wCV (%)</th><th>CoR (%)</th>")
        h.append("</tr></thead><tbody>")

        for i, method in enumerate(methods):
            wcv = _safe_get(median_wcv, i)
            cor = _safe_get(cor_vals, i)
            h.append(f"<tr><td><strong>{_esc(method)}</strong></td>")

            if wcv is not None:
                h.append(f"<td>{wcv * 100:.1f}%</td>")
            else:
                h.append("<td>&mdash;</td>")

            h.append(_cor_cell(cor))
            h.append("</tr>")

        h.append("</tbody></table>")

    return h


def _safe_get(lst: list | None, idx: int) -> float | None:
    """Safely extract a value from a list."""
    if not lst:
        return None
    try:
        val = lst[idx]
        if val is None:
            return None
        return float(val)
    except (IndexError, TypeError, ValueError):
        return None


def _cor_cell(val: float | None) -> str:
    """Format a CoR value as a colored table cell."""
    if val is None:
        return "<td>&mdash;</td>"
    if val < 15:
        cls = ' style="background:#d4edda"'  # green
    elif val <= 30:
        cls = ' style="background:#fff3cd"'  # yellow
    else:
        cls = ' style="background:#f8d7da"'  # red
    return f"<td{cls}>{val:.1f}%</td>"
