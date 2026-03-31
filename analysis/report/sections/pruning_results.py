"""Report section: core method pruning results."""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


def _section_pruning_results(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build HTML section showing which core methods were pruned and why.

    Parameters
    ----------
    mat_data : dict
        Mapping of DWI type to parsed MAT metrics dict.
    dwi_types : list[str]
        List of DWI types present in the data.

    Returns
    -------
    list[str]
        HTML chunks for the pruning results section.
    """
    h: list[str] = []

    has_data = False
    for dwi in DWI_TYPES:
        if dwi in (mat_data or {}) and "pruning" in mat_data[dwi]:
            pr = mat_data[dwi]["pruning"]
            if pr and pr.get("active_methods") is not None:
                has_data = True
                break

    if not has_data:
        h.append(_h2("Core Method Pruning", "core-pruning"))
        h.append('<p class="meta">No pruning data available. Pruning is applied when '
                 "<code>max_core_failure_rate</code> &lt; 1.0 or "
                 "<code>excluded_core_methods</code> is set.</p>")
        return h

    h.append(_h2("Core Method Pruning Results", "core-pruning"))
    h.append(
        "<p>Methods excluded from downstream analysis based on failure rate "
        "threshold, manual exclusions, and/or minimum voxel count criteria.</p>"
    )

    # Show min_core_voxels threshold if available from any DWI type
    for dwi in DWI_TYPES:
        pr = (mat_data or {}).get(dwi, {}).get("pruning", {})
        mcv = pr.get("min_core_voxels")
        if mcv is not None:
            h.append(
                f'<p class="meta">Insufficient voxels threshold: '
                f'<code>min_core_voxels = {mcv}</code> '
                f'(minimum for meaningful dosimetry; distinct from '
                f'<code>min_vox_hist</code> used for histogram statistics).</p>'
            )
            break

    for dwi in DWI_TYPES:
        if dwi not in (mat_data or {}):
            continue
        pr = mat_data[dwi].get("pruning", {})
        if not pr or pr.get("active_methods") is None:
            continue

        active = pr.get("active_methods", [])
        pruned = pr.get("pruned_info", [])

        h.append(f"<h3>{_dwi_badge(dwi)} &mdash; Pruning Summary</h3>")

        # Pruned methods table
        if pruned:
            h.append("<h4>Pruned Methods</h4>")
            h.append("<table><thead><tr>"
                     "<th>Method</th><th>Reason</th>"
                     "<th>Failure Rate</th><th>Worst Pipeline</th>"
                     "</tr></thead><tbody>")
            for entry in pruned:
                name = entry.get("name", "")
                reason = entry.get("reason", "").replace("_", " ").title()
                fr = entry.get("failure_rate")
                pipe = entry.get("pipeline", "")
                fr_str = f"{fr * 100:.1f}%" if fr is not None and isinstance(fr, (int, float)) else "&mdash;"
                h.append(f"<tr><td><strong>{_esc(name)}</strong></td>"
                         f"<td>{_esc(reason)}</td>"
                         f"<td>{fr_str}</td>"
                         f"<td>{_esc(pipe)}</td></tr>")
            h.append("</tbody></table>")
        else:
            h.append("<p>No methods were pruned.</p>")

        # Retained-with-warning methods (e.g. adc_threshold kept as fallback)
        rw = pr.get("retained_with_warning", [])
        if rw:
            h.append("<h4>Retained with Warning</h4>")
            h.append('<div class="warn-box">')
            for rw_entry in rw:
                rw_name = rw_entry.get("name", "")
                rw_reason = rw_entry.get("reason", "")
                rw_fr = rw_entry.get("failure_rate")
                fr_str = f" ({rw_fr * 100:.1f}% failure rate)" if rw_fr is not None and isinstance(rw_fr, (int, float)) and rw_fr == rw_fr else ""
                h.append(f"<p><strong>{_esc(rw_name)}</strong>{fr_str}: {_esc(rw_reason)}</p>")
            h.append("</div>")

        # Retained methods
        if active:
            h.append(f"<h4>Retained Methods ({len(active)})</h4>")
            h.append("<p>" + ", ".join(f"<code>{_esc(m)}</code>" for m in active) + "</p>")

    return h
