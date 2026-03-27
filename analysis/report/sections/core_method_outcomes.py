"""Report section: core method sub-volume dose coverage vs local control."""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


def _section_core_method_outcomes(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build HTML section showing which core methods' dose coverage predicts LC.

    Parameters
    ----------
    mat_data : dict
        Mapping of DWI type to parsed MAT metrics dict.
    dwi_types : list[str]
        List of DWI types present in the data.

    Returns
    -------
    list[str]
        HTML chunks for the core method outcomes section.
    """
    h: list[str] = []

    has_data = False
    for dwi in DWI_TYPES:
        if dwi in (mat_data or {}) and "core_method_outcomes" in mat_data[dwi]:
            cmo = mat_data[dwi]["core_method_outcomes"]
            if cmo and cmo.get("method_results"):
                has_data = True
                break

    if not has_data:
        h.append(_h2("Core Method Outcome Analysis", "core-method-outcomes"))
        h.append('<p class="meta">Core method outcome data not available. '
                 "Enable <code>run_core_method_outcomes</code> and "
                 "<code>run_all_core_methods</code> in config to generate.</p>")
        return h

    h.append(_h2("Core Method Outcome Analysis", "core-method-outcomes"))
    h.append(
        "<p>Univariable Cox PH analysis testing whether dose coverage (D95, V50) "
        "of each core method's resistant sub-volume predicts local control. "
        "Methods are ranked by their best p-value across dose metrics.</p>"
    )

    for dwi in DWI_TYPES:
        if dwi not in (mat_data or {}):
            continue
        cmo = mat_data[dwi].get("core_method_outcomes", {})
        if not cmo or not cmo.get("method_results"):
            continue

        method_results = cmo["method_results"]
        ranking = cmo.get("ranking", [])

        h.append(f"<h3>{_dwi_badge(dwi)} &mdash; Method Ranking</h3>")

        # Build ranking table
        h.append("<table><thead><tr>"
                 "<th>Rank</th><th>Method</th><th>Best Metric</th>"
                 "<th>HR (95% CI)</th><th>p-value</th><th>Events</th>"
                 "</tr></thead><tbody>")

        # Sort by ranking if available, otherwise use order
        ordered = []
        if ranking:
            for mname in ranking:
                for mr in method_results:
                    if mr.get("method_name") == mname:
                        ordered.append(mr)
                        break
        else:
            ordered = method_results

        has_significant = False
        for rank, mr in enumerate(ordered, 1):
            mname = mr.get("method_name", "")
            univar = mr.get("univariable", [])
            km = mr.get("km", {})
            n_events = mr.get("n_events", 0)

            best_met = km.get("best_metric", "") if km else ""
            hr_str = "&mdash;"
            p_str = "&mdash;"
            p_val = None

            if univar:
                # Find best (lowest p) univariable result
                best = min(univar, key=lambda u: u.get("p_value", 1))
                hr = best.get("hr")
                ci = best.get("hr_ci", [])
                p_val = best.get("p_value")

                if hr is not None:
                    if ci and len(ci) >= 2 and ci[0] is not None and ci[1] is not None:
                        hr_str = f"{hr:.2f} ({ci[0]:.2f}&ndash;{ci[1]:.2f})"
                    else:
                        hr_str = f"{hr:.2f}"
                if p_val is not None:
                    p_str = f"{p_val:.4f}"

            row_cls = ""
            if p_val is not None and p_val < 0.05:
                row_cls = ' class="sig-moderate"'
                has_significant = True

            h.append(f"<tr{row_cls}><td>{rank}</td>"
                     f"<td><strong>{_esc(mname)}</strong></td>"
                     f"<td>{_esc(best_met)}</td>"
                     f"<td>{hr_str}</td>"
                     f"<td>{p_str}</td>"
                     f"<td>{n_events}</td></tr>")

        h.append("</tbody></table>")

        # Interpretation
        if has_significant and ordered:
            best = ordered[0]
            mname = best.get("method_name", "")
            univar = best.get("univariable", [])
            if univar:
                best_uv = min(univar, key=lambda u: u.get("p_value", 1))
                hr = best_uv.get("hr", "?")
                p = best_uv.get("p_value", "?")
                ci = best_uv.get("hr_ci", [])
                ci_str = ""
                if ci and len(ci) >= 2 and ci[0] is not None:
                    ci_str = f", 95% CI [{ci[0]:.2f}, {ci[1]:.2f}]"
                h.append(
                    f"<p><strong>Interpretation:</strong> "
                    f"<code>{_esc(mname)}</code> shows the strongest association "
                    f"between sub-volume dose coverage and local control "
                    f"(HR = {hr:.2f}{ci_str}, p = {p:.4f}). "
                    f"Ensuring adequate dose to the {_esc(mname)}-defined resistant "
                    f"sub-volume may improve local control.</p>"
                )
        elif not has_significant:
            h.append(
                "<p>No core method's dose coverage significantly predicted "
                "local control at the 0.05 level.</p>"
            )

    return h
