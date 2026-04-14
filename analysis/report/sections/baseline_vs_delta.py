"""Report section: Baseline vs treatment-delta predictive performance.

Compares the prognostic value of baseline (pre-treatment) parameter values
against treatment-induced deltas for each DWI parameter and timepoint.
Winner column identifies the more predictive representation by Cox PH
hazard ratio significance and concordance (C-index).
"""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _sig_class,
)


def _scalar(val) -> float | None:
    if val is None:
        return None
    try:
        f = float(val)
        if f != f:
            return None
        return f
    except (TypeError, ValueError):
        return None


def _fmt_hr_p(hr: float | None, p: float | None) -> str:
    """Render 'HR (p=…)' with significance class on p."""
    hr_s = f"{hr:.2f}" if hr is not None else "&mdash;"
    if p is None:
        return f"{hr_s}"
    try:
        pv = float(p)
    except (TypeError, ValueError):
        return f"{hr_s}"
    cls = _sig_class(pv)
    p_txt = f"{pv:.3g}" if pv >= 0.001 else f"{pv:.1e}"
    p_html = f'<span class="{cls}">p={p_txt}</span>' if cls else f"p={p_txt}"
    return f"{hr_s} ({p_html})"


def _fmt_cindex(c: float | None, is_winner: bool) -> str:
    if c is None:
        return "&mdash;"
    s = f"{c:.3f}"
    return f"<strong>{s}</strong>" if is_winner else s


def _section_baseline_vs_delta(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build HTML comparing baseline vs delta predictors."""
    if not mat_data:
        return []

    dwi_with_data: list[str] = []
    for dwi in DWI_TYPES:
        entry = (mat_data or {}).get(dwi) or {}
        rows = entry.get("baseline_vs_delta")
        if rows:
            dwi_with_data.append(dwi)

    if not dwi_with_data:
        return []

    h: list[str] = []
    h.append(_h2("Baseline vs Delta Predictors", "baseline-vs-delta"))
    h.append(
        "<p>Per-parameter comparison of baseline (pre-treatment) values vs "
        "treatment-induced deltas as predictors of local failure. The "
        "'Winner' column identifies the representation with the larger "
        "concordance index; when a delta wins, measuring change during "
        "treatment adds predictive value beyond baseline alone.</p>"
    )

    for dwi in dwi_with_data:
        rows = mat_data[dwi]["baseline_vs_delta"]

        h.append(f"<h3>{_dwi_badge(dwi)}</h3>")
        h.append("<table><thead><tr>"
                 "<th>Parameter</th><th>Timepoint</th>"
                 "<th>Baseline HR (p)</th><th>Delta HR (p)</th>"
                 "<th>Baseline C</th><th>Delta C</th>"
                 "<th>Winner</th>"
                 "</tr></thead><tbody>")

        for row in rows:
            if not isinstance(row, dict):
                continue
            param = str(row.get("parameter", "") or "")
            tp = str(row.get("timepoint", "") or "")
            b_hr = _scalar(row.get("baseline_hr"))
            b_p = _scalar(row.get("baseline_p"))
            b_c = _scalar(row.get("baseline_cindex"))
            d_hr = _scalar(row.get("delta_hr"))
            d_p = _scalar(row.get("delta_p"))
            d_c = _scalar(row.get("delta_cindex"))
            winner = str(row.get("better_predictor", "") or "")

            baseline_wins = winner.lower().startswith("baseline")
            delta_wins = winner.lower().startswith("delta")

            winner_html = ""
            if winner:
                winner_html = (
                    f'<span style="background:#d4edda;padding:2px 6px;'
                    f'border-radius:3px">{_esc(winner)}</span>'
                )

            h.append("<tr>")
            h.append(f"<td><strong>{_esc(param) or '&mdash;'}</strong></td>")
            h.append(f"<td>{_esc(tp) or '&mdash;'}</td>")
            h.append(f"<td>{_fmt_hr_p(b_hr, b_p)}</td>")
            h.append(f"<td>{_fmt_hr_p(d_hr, d_p)}</td>")
            h.append(f"<td>{_fmt_cindex(b_c, baseline_wins)}</td>")
            h.append(f"<td>{_fmt_cindex(d_c, delta_wins)}</td>")
            h.append(f"<td>{winner_html or '&mdash;'}</td>")
            h.append("</tr>")

        h.append("</tbody></table>")

    h.append(
        "<p><strong>Interpretation:</strong> When the delta C-index exceeds "
        "the baseline C-index, treatment-induced changes provide incremental "
        "predictive value beyond pre-treatment imaging alone. This is "
        "clinically important when deciding whether mid-treatment scans "
        "add information sufficient to justify adaptive therapy decisions.</p>"
    )

    return h
