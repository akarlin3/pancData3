"""Report section: Sub-volume vs whole-GTV dose context.

Compares cohort-level dose coverage (D95, Dmean, V50) between the whole
primary GTV and the diffusion-defined resistant sub-volume.  Highlights
large deficits, which indicate that the sub-volume received substantially
less dose than the rest of the target — a candidate for adaptive dose
escalation.
"""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


def _scalar(val) -> float | None:
    """Return a plain float from either a bare scalar or a {mean,...} dict."""
    if val is None:
        return None
    if isinstance(val, dict):
        val = val.get("mean")
    if isinstance(val, (int, float)) and val == val:
        return float(val)
    return None


def _v50_as_pct(val: float | None) -> float | None:
    """V50 is stored as a fraction in [0, 1] or as a percentage; normalise to %."""
    if val is None:
        return None
    return val * 100.0 if 0.0 <= val <= 1.0 else float(val)


def _section_dose_context(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build HTML comparing whole-GTV vs sub-volume dose coverage.

    Returns an empty list if no DWI type has the required whole-GTV +
    sub-volume dose data.
    """
    if not mat_data:
        return []

    rendered_any = False
    h: list[str] = []
    h.append(_h2("Sub-Volume Dose Context", "dose-context"))
    h.append(
        "<p>Comparison of cohort-level dose coverage between the whole "
        "primary GTV and the diffusion-defined resistant sub-volume. Large "
        "deficits (&gt; 5 Gy or &gt; 10% volume at 50 Gy) indicate "
        "candidates for adaptive dose escalation.</p>"
    )

    for dwi in DWI_TYPES:
        entry = (mat_data or {}).get(dwi) or {}
        dosi = entry.get("dosimetry") or {}
        if not dosi:
            continue

        whole_d95 = _scalar(dosi.get("whole_gtv_d95_mean"))
        whole_dmean = _scalar(dosi.get("whole_gtv_dmean_mean"))
        whole_v50 = _v50_as_pct(_scalar(dosi.get("whole_gtv_v50_mean")))
        sub_d95 = _scalar(dosi.get("d95_adc_mean"))
        sub_v50 = _v50_as_pct(_scalar(dosi.get("v50_adc_mean")))
        # Sub-volume mean dose not typically stored; leave None.
        sub_dmean = _scalar(dosi.get("dmean_adc_mean"))

        if whole_d95 is None and whole_v50 is None and whole_dmean is None:
            continue
        if sub_d95 is None and sub_v50 is None and sub_dmean is None:
            continue

        rendered_any = True
        h.append(f"<h3>{_dwi_badge(dwi)} &mdash; Whole GTV vs Sub-Volume</h3>")
        h.append("<table><thead><tr>")
        h.append("<th>Metric</th><th>Whole GTV</th><th>Sub-Volume</th>"
                 "<th>Deficit</th><th></th></tr></thead><tbody>")

        def _row(label: str, whole: float | None, sub: float | None,
                 threshold: float, fmt: str, unit: str) -> None:
            if whole is None and sub is None:
                return
            deficit = None
            if whole is not None and sub is not None:
                deficit = whole - sub
            flagged = deficit is not None and deficit > threshold
            row_style = ' style="background:#f8d7da"' if flagged else ""
            whole_s = (fmt % whole) if whole is not None else "&mdash;"
            sub_s = (fmt % sub) if sub is not None else "&mdash;"
            deficit_s = (fmt % deficit) if deficit is not None else "&mdash;"
            note = ""
            if flagged:
                note = (
                    f"<span class=\"differ\">{_esc('Clinically significant dose deficit')}</span>"
                )
            h.append(
                f"<tr{row_style}><td><strong>{_esc(label)}</strong></td>"
                f"<td>{whole_s}{_esc(' ' + unit) if whole is not None else ''}</td>"
                f"<td>{sub_s}{_esc(' ' + unit) if sub is not None else ''}</td>"
                f"<td>{deficit_s}{_esc(' ' + unit) if deficit is not None else ''}</td>"
                f"<td>{note}</td></tr>"
            )

        _row("D95", whole_d95, sub_d95, 5.0, "%.1f", "Gy")
        _row("Dmean", whole_dmean, sub_dmean, 5.0, "%.1f", "Gy")
        _row("V50", whole_v50, sub_v50, 10.0, "%.1f", "%")

        h.append("</tbody></table>")

    if not rendered_any:
        return []
    return h
