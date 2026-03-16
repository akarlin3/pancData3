"""Report sections: supplemental MAT file data (dosimetry, core method comparison)."""

from __future__ import annotations

from shared import (  # type: ignore
    DWI_TYPES,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)
from report.sections._helpers import _scalar_gy  # type: ignore


def _section_mat_data(mat_data) -> list[str]:
    """Build the Supplemental Data (MAT Files) section.

    Displays two sub-sections from parsed MAT-file JSON:
    1. **Dosimetry** -- mean D95 and V50 values for ADC and D sub-volumes.
    2. **Core Method Comparison** -- truncated mean-Dice heatmap matrix
       showing agreement between tumor core delineation methods.

    Parameters
    ----------
    mat_data : dict
        Mapping of DWI type to parsed MAT metrics dict (from
        ``parsed_mat_metrics_{dwi}.json``).

    Returns
    -------
    list[str]
        HTML chunks for the supplemental data section.
    """
    # ── 7.5. Core Method & Dosimetry ──
    h: list[str] = []
    if mat_data:
        h.append(_h2("Supplemental Data (MAT Files)", "supplemental"))

        # Helper formatters defined at function scope (before any branch that uses them).
        def _fmt_gy(val):
            """Format a dose value in Gy — supports dict {mean, std}, plain float, None, NaN."""
            if val is None:
                return "\u2014"
            if isinstance(val, dict):
                mean = val.get("mean")
                std = val.get("std")
                if mean is None or not (isinstance(mean, (int, float)) and mean == mean):
                    return "\u2014"
                if std is not None and isinstance(std, (int, float)) and std == std:
                    return f"{mean:.1f} \u00b1 {std:.1f} Gy"
                return f"{mean:.2f} Gy"
            if isinstance(val, (int, float)) and val == val:
                return f"{val:.2f}"
            return "\u2014"

        def _fmt_pct(val):
            """Format a value as percentage — supports dict {mean, std}, plain float, None, NaN."""
            if val is None:
                return "\u2014"
            if isinstance(val, dict):
                mean = val.get("mean")
                std = val.get("std")
                if mean is None or not (isinstance(mean, (int, float)) and mean == mean):
                    return "\u2014"
                pct_mean = mean * 100 if mean <= 1.0 else mean
                if std is not None and isinstance(std, (int, float)) and std == std:
                    pct_std = std * 100 if std <= 1.0 else std
                    return f"{pct_mean:.1f} \u00b1 {pct_std:.1f}%"
                return f"{pct_mean:.1f}%"
            if isinstance(val, (int, float)) and val == val:
                pct = val * 100 if val <= 1.0 else val
                return f"{pct:.1f}%"
            return "\u2014"

        has_dos = any("dosimetry" in d for d in mat_data.values())
        if has_dos:
            h.append("<h3>Dosimetry (Target Coverage)</h3>")
            h.append("<table><thead><tr><th>DWI</th><th>D95 ADC</th><th>V50 ADC</th>"
                     "<th>D95 D</th><th>V50 D</th></tr></thead><tbody>")
            for dt in DWI_TYPES:
                if dt not in mat_data or "dosimetry" not in mat_data[dt]:
                    continue
                dos = mat_data[dt]["dosimetry"]
                if not dos:
                    continue
                h.append(
                    f"<tr><td>{_dwi_badge(dt)}</td>"
                    f"<td>{_fmt_gy(dos.get('d95_adc_mean'))}</td>"
                    f"<td>{_fmt_pct(dos.get('v50_adc_mean'))}</td>"
                    f"<td>{_fmt_gy(dos.get('d95_d_mean'))}</td>"
                    f"<td>{_fmt_pct(dos.get('v50_d_mean'))}</td></tr>"
                )
            h.append("</tbody></table>")

            # Clinical reference context box (change 6)
            h.append(
                '<div class="info-box">'
                "<strong>Clinical Reference:</strong> Standard pancreatic RT prescribes "
                "50\u201354\u00a0Gy to the GTV (RTOG protocol). "
                "D95 &gt; 45\u00a0Gy indicates adequate tumor coverage; "
                "V50 &gt; 90% is the target threshold for radical treatment intent."
                "</div>"
            )

            # Per-DWI D95 pass/fail against 45 Gy threshold (change 6 continued)
            pass_fail_notes: list[str] = []
            for dt in DWI_TYPES:
                if dt not in mat_data or "dosimetry" not in mat_data[dt]:  # type: ignore
                    continue
                dos = mat_data[dt]["dosimetry"]  # type: ignore
                if not dos:
                    continue
                mean_val = _scalar_gy(dos.get("d95_adc_mean"))
                if mean_val is None:
                    continue
                status = "\u2705 PASS" if mean_val >= 45 else "\u274c FAIL"
                pass_fail_notes.append(
                    f"<li>{_esc(dt)}: D95 ADC = {mean_val:.1f}\u00a0Gy "
                    f"\u2014 {_esc(status)} (threshold: 45\u00a0Gy)</li>"
                )
            if pass_fail_notes:
                h.append(
                    '<div class="info-box"><strong>D95 Coverage Check (\u226545\u00a0Gy):</strong>'
                    "<ul>" + "".join(pass_fail_notes) + "</ul></div>"
                )

            # Pre-existing clinical interpretation notes (preserved, now using _scalar_gy)
            dos_notes: list[str] = []
            for dt in DWI_TYPES:
                if dt not in mat_data or "dosimetry" not in mat_data[dt]:  # type: ignore
                    continue
                dos = mat_data[dt]["dosimetry"]  # type: ignore
                if not dos:
                    continue
                d95_adc = _scalar_gy(dos.get("d95_adc_mean"))
                v50_adc = _scalar_gy(dos.get("v50_adc_mean"))
                if d95_adc is not None:
                    if d95_adc < 45:
                        dos_notes.append(
                            f"{dt}: D95 ADC = {d95_adc:.1f} Gy "
                            "(below 45 Gy \u2014 potential under-dosing of resistant sub-volume)"
                        )
                    elif d95_adc >= 50:
                        dos_notes.append(
                            f"{dt}: D95 ADC = {d95_adc:.1f} Gy "
                            "(adequate coverage of resistant sub-volume)"
                        )
                if v50_adc is not None:
                    v50_pct = v50_adc * 100 if v50_adc <= 1.0 else v50_adc
                    if v50_pct < 90:
                        dos_notes.append(
                            f"{dt}: V50 ADC = {v50_pct:.0f}% "
                            "(partial coverage \u2014 may benefit from dose escalation)"
                        )
            if dos_notes:
                h.append('<div class="info-box"><strong>Dosimetric interpretation:</strong><ul>')
                for note in dos_notes:
                    h.append(f"<li>{_esc(note)}</li>")
                h.append("</ul></div>")

        has_core = any("core_method" in d for d in mat_data.values())
        if has_core:
            h.append("<h3>Core Method Comparison (Mean Dice)</h3>")
            for dt in DWI_TYPES:
                if dt not in mat_data or "core_method" not in mat_data[dt]:
                    continue
                core = mat_data[dt]["core_method"]
                if not core or not core.get("methods"):
                    continue
                h.append(f"<h4>{_dwi_badge(dt)} \u2014 {len(core['methods'])} methods</h4>")
                methods = core["methods"]
                matrix = core["mean_dice_matrix"]
                n = len(methods)

                # Mean Dice matrix
                h.append('<div style="overflow-x:auto">')
                h.append('<table class="table-wide"><thead><tr><th>Method</th>')
                for m in methods:
                    h.append(f"<th>{_esc(m)}</th>")
                h.append("</tr></thead><tbody>")
                for i in range(n):
                    h.append(f"<tr><td><strong>{_esc(methods[i])}</strong></td>")
                    for j in range(n):
                        val = matrix[i][j] if i < len(matrix) and j < len(matrix[i]) else None  # type: ignore
                        if isinstance(val, (int, float)):
                            # Colour-code: high Dice = green-ish
                            if i != j:
                                style = (
                                    f" style=\"color: "
                                    f"{'var(--green)' if val >= 0.7 else ('var(--amber)' if val >= 0.5 else 'var(--red)')};"
                                    f" font-weight:600\""
                                )
                            else:
                                style = ""
                            h.append(f"<td{style}>{val:.2f}</td>")
                        else:
                            h.append("<td>-</td>")
                    h.append("</tr>")
                h.append("</tbody></table>")
                h.append("</div>")

                # Hausdorff matrix (change 4) — shown when hausdorff_matrix key is present
                hausdorff_matrix = core.get("hausdorff_matrix")
                if hausdorff_matrix:
                    h.append("<h5>Mean Hausdorff Distance (mm)</h5>")
                    h.append('<div style="overflow-x:auto">')
                    h.append('<table class="table-wide"><thead><tr><th>Method</th>')
                    for m in methods:
                        h.append(f"<th>{_esc(m)}</th>")
                    h.append("</tr></thead><tbody>")
                    for i in range(n):
                        h.append(f"<tr><td><strong>{_esc(methods[i])}</strong></td>")
                        for j in range(n):
                            hval = (  # type: ignore
                                hausdorff_matrix[i][j]  # type: ignore
                                if i < len(hausdorff_matrix) and j < len(hausdorff_matrix[i])  # type: ignore
                                else None
                            )
                            if isinstance(hval, (int, float)) and i != j:
                                col = ("var(--green)" if hval <= 5 else
                                       ("var(--amber)" if hval <= 15 else "var(--red)"))
                                h.append(
                                    f"<td style=\"color:{col};font-weight:600\">"
                                    f"{hval:.1f}</td>"
                                )
                            elif isinstance(hval, (int, float)):
                                h.append(f"<td>{hval:.1f}</td>")
                            else:
                                h.append("<td>-</td>")
                        h.append("</tr>")
                    h.append("</tbody></table>")
                    h.append("</div>")

                # Summary statistics + recommendation (change 7)
                off_diag: list[float] = []
                method_avg: dict[str, float] = {}
                for i in range(n):
                    row_vals: list[float] = []
                    for j in range(n):
                        if i == j:
                            continue
                        if i < len(matrix) and j < len(matrix[i]):  # type: ignore
                            val = matrix[i][j]  # type: ignore
                            if isinstance(val, (int, float)) and val > 0:  # type: ignore
                                off_diag.append(val)
                                row_vals.append(val)
                    if row_vals:
                        method_avg[methods[i]] = sum(row_vals) / len(row_vals)  # type: ignore

                if off_diag:
                    avg_dice = sum(off_diag) / len(off_diag)
                    min_dice = min(off_diag)
                    max_dice = max(off_diag)
                    h.append(
                        f'<div class="info-box">Mean pairwise Dice: '
                        f'<strong>{avg_dice:.3f}</strong> '
                        f'(range: {min_dice:.3f}\u2013{max_dice:.3f} '
                        f'across {len(off_diag)} pairs)</div>'
                    )

                if method_avg:
                    best_method = max(method_avg, key=lambda k: method_avg[k])
                    interchangeable_pairs: list[str] = []
                    for i in range(n):
                        for j in range(i + 1, n):
                            if i < len(matrix) and j < len(matrix[i]):  # type: ignore
                                pval = matrix[i][j]  # type: ignore
                                if isinstance(pval, (int, float)) and pval >= 0.70:
                                    interchangeable_pairs.append(
                                        f"{_esc(methods[i])} \u2194 {_esc(methods[j])} "  # type: ignore
                                        f"(Dice\u2009=\u2009{pval:.2f})"
                                    )
                    rec_parts = [
                        '<div class="info-box">',
                        "<strong>Recommendation:</strong> Methods with mean Dice "
                        "\u2265\u20090.70 relative to all others can be considered "
                        "interchangeable. ",
                    ]
                    if interchangeable_pairs:
                        rec_parts.append(
                            "Interchangeable pairs: <ul>"
                            + "".join(f"<li>{p}</li>" for p in interchangeable_pairs)
                            + "</ul>"
                        )
                    else:
                        rec_parts.append("No method pairs meet the \u22650.70 threshold. ")
                    rec_parts.append(
                        f"The method with the highest average agreement is "
                        f"<strong>{_esc(best_method)}</strong> "
                        f"(mean Dice\u2009=\u2009{method_avg[best_method]:.3f}) "
                        f"\u2014 recommended as the reference method."
                        "</div>"
                    )
                    h.append("".join(rec_parts))

    return h
