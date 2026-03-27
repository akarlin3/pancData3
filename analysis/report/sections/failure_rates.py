"""Report section: core method failure rate breakdown."""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


def _section_failure_rates(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build HTML section showing core method failure rate breakdown.

    Displays a table of failure rates (Fallback %, Empty %, Insufficient %,
    All-NaN %, Total Fail %) for each of the 11 core methods, sorted by
    total failure rate descending.  Cells are color-coded: green (<10%),
    yellow (10-25%), red (>25%).

    Parameters
    ----------
    mat_data : dict
        Mapping of DWI type to parsed MAT metrics dict.
    dwi_types : list[str]
        List of DWI types present in the data.

    Returns
    -------
    list[str]
        HTML chunks for the failure rates section.
    """
    h: list[str] = []

    # Check if any DWI type has failure_rates data
    has_data = False
    for dwi in DWI_TYPES:
        if dwi in (mat_data or {}) and "failure_rates" in mat_data[dwi]:
            fr = mat_data[dwi]["failure_rates"]
            if fr and fr.get("method_names"):
                has_data = True
                break

    if not has_data:
        h.append(_h2("Core Method Failure Rates", "failure-rates"))
        h.append('<p class="meta">Core method failure rate data not available. '
                 "Enable <code>run_core_failure_rates</code> in config to generate.</p>")
        return h

    h.append(_h2("Core Method Failure Rates", "failure-rates"))
    h.append(
        "<p>Detailed failure mode breakdown for each of the 11 tumor core "
        "delineation methods. Four failure categories are tracked: "
        "<strong>Fallback</strong> (method fell back to ADC threshold), "
        "<strong>Empty</strong> (core mask has zero voxels), "
        "<strong>Insufficient</strong> (core has fewer voxels than <code>min_vox_hist</code>), "
        "and <strong>All-NaN</strong> (input parameter vectors were entirely NaN). "
        "Methods with total failure rate &gt;25% are candidates for pruning.</p>"
    )

    for dwi in DWI_TYPES:
        if dwi not in (mat_data or {}):
            continue
        fr = mat_data[dwi].get("failure_rates", {})
        if not fr or not fr.get("method_names"):
            continue

        methods = fr["method_names"]
        pipeline_names = fr.get("pipeline_names", ["Standard", "DnCNN", "IVIMNet"])
        fallback = fr.get("fallback_rate", [])
        empty = fr.get("empty_rate", [])
        insufficient = fr.get("insufficient_rate", [])
        all_nan = fr.get("all_nan_rate", [])
        any_fail = fr.get("any_failure_rate", [])
        median_vox = fr.get("median_core_voxels", [])

        for p_idx, pipe_name in enumerate(pipeline_names):
            h.append(f"<h3>{_dwi_badge(dwi)} &mdash; {_esc(pipe_name)} Pipeline</h3>")

            # Build rows with total failure for sorting
            rows = []
            for i, method in enumerate(methods):
                fb = _get_rate(fallback, i, p_idx)
                em = _get_rate(empty, i, p_idx)
                ins = _get_rate(insufficient, i, p_idx)
                an = _get_rate(all_nan, i, p_idx)
                af = _get_rate(any_fail, i, p_idx)
                mv = _get_val(median_vox, i, p_idx)
                rows.append((method, fb, em, ins, an, af, mv))

            # Sort by total failure descending
            rows.sort(key=lambda r: r[5] if r[5] is not None else -1, reverse=True)

            h.append("<table><thead><tr>")
            h.append("<th>Core Method</th><th>Fallback</th><th>Empty</th>"
                     "<th>Insufficient</th><th>All-NaN</th><th>Total Fail</th>"
                     "<th>Median Voxels</th>")
            h.append("</tr></thead><tbody>")

            for method, fb, em, ins, an, af, mv in rows:
                h.append(f"<tr><td><strong>{_esc(method)}</strong></td>")
                h.append(_rate_cell(fb))
                h.append(_rate_cell(em))
                h.append(_rate_cell(ins))
                h.append(_rate_cell(an))
                h.append(_rate_cell(af))
                if mv is not None:
                    h.append(f"<td>{mv:.0f}</td>")
                else:
                    h.append("<td>&mdash;</td>")
                h.append("</tr>")

            h.append("</tbody></table>")

    return h


def _get_rate(matrix: list | None, row: int, col: int):
    """Safely extract a rate value from a 2D nested list."""
    if not matrix:
        return None
    try:
        val = matrix[row][col]
        if val is None:
            return None
        return float(val)
    except (IndexError, TypeError, ValueError):
        return None


def _get_val(matrix: list | None, row: int, col: int):
    """Safely extract a scalar value from a 2D nested list."""
    return _get_rate(matrix, row, col)


def _rate_cell(val) -> str:
    """Format a rate value as a colored table cell (percentage)."""
    if val is None:
        return "<td>&mdash;</td>"
    pct = val * 100
    if pct > 25:
        cls = ' class="sig-strong"'
    elif pct > 10:
        cls = ' class="sig-moderate"'
    else:
        cls = ""
    return f"<td{cls}>{pct:.1f}%</td>"
