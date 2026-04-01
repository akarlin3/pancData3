"""Report section: GTV volume confounding check."""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


def _section_gtv_confounding(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build HTML section showing GTV volume confounding analysis.

    Parameters
    ----------
    mat_data : dict
        Mapping of DWI type to parsed MAT metrics dict.
    dwi_types : list[str]
        List of DWI types present in the data.

    Returns
    -------
    list[str]
        HTML chunks for the GTV confounding section.
    """
    h: list[str] = []

    has_data = False
    for dwi in DWI_TYPES:
        if dwi in (mat_data or {}) and "gtv_confounding" in mat_data[dwi]:
            gc = mat_data[dwi]["gtv_confounding"]
            if gc and gc.get("method_results"):
                has_data = True
                break

    if not has_data:
        h.append(_h2("GTV Volume Confounding Check", "gtv-confounding"))
        h.append('<p class="meta">GTV confounding data not available. '
                 "Enable <code>run_gtv_confounding</code> in config to generate.</p>")
        return h

    h.append(_h2("GTV Volume Confounding Check", "gtv-confounding"))
    h.append(
        "<p>Tests whether the association between sub-volume dose coverage (D95) "
        "and local control is confounded by GTV volume. Patients with shrinking "
        "GTVs naturally get better D95 coverage, so an apparent &ldquo;better D95 "
        "&rarr; better LC&rdquo; finding might reflect GTV response rather than "
        "sub-volume dose. The confounding flag is set when adjusting for GTV volume "
        "changes the D95 hazard ratio by &gt;10%.</p>"
    )

    for dwi in DWI_TYPES:
        if dwi not in (mat_data or {}):
            continue
        gc = mat_data[dwi].get("gtv_confounding", {})
        if not gc or not gc.get("method_results"):
            continue

        method_results = gc["method_results"]
        summary = gc.get("summary", "")

        h.append(f"<h3>{_dwi_badge(dwi)} — Confounding Analysis</h3>")
        h.append("<table><thead><tr>")
        h.append("<th>Core Method</th><th>D95-GTV &rho;</th>"
                 "<th>Unadjusted HR</th><th>Adjusted HR</th>"
                 "<th>% Change</th><th>Confounded?</th>")
        h.append("</tr></thead><tbody>")

        for mr in method_results:
            name = mr.get("method_name", "")
            rho = mr.get("d95_gtv_correlation")
            unadj_hr = mr.get("unadjusted_hr")
            adj_hr = mr.get("adjusted_hr")
            flag = mr.get("confounding_flag", False)

            # Compute % change in log(HR)
            pct_change = None
            if unadj_hr is not None and adj_hr is not None and unadj_hr > 0 and adj_hr > 0:
                import math
                log_unadj = abs(math.log(unadj_hr))
                log_adj = abs(math.log(adj_hr))
                if log_unadj > 0:
                    pct_change = abs(log_adj - log_unadj) / log_unadj * 100

            h.append(f"<tr><td><strong>{_esc(name)}</strong></td>")

            if rho is not None:
                h.append(f"<td>{rho:.3f}</td>")
            else:
                h.append("<td>&mdash;</td>")

            if unadj_hr is not None:
                h.append(f"<td>{unadj_hr:.3f}</td>")
            else:
                h.append("<td>&mdash;</td>")

            if adj_hr is not None:
                h.append(f"<td>{adj_hr:.3f}</td>")
            else:
                h.append("<td>&mdash;</td>")

            if pct_change is not None:
                h.append(f"<td>{pct_change:.1f}%</td>")
            else:
                h.append("<td>&mdash;</td>")

            if flag:
                h.append('<td class="sig-strong">Yes</td>')
            else:
                h.append("<td>No</td>")

            h.append("</tr>")

        h.append("</tbody></table>")

        if summary:
            if "confounder" in summary.lower():
                h.append(f'<div class="warn-box">{_esc(summary)}</div>')
            else:
                h.append(f'<p class="meta">{_esc(summary)}</p>')

    return h
