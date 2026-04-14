"""Report section: Sub-volume sizes stratified by local control / local failure.

Displays the absolute volume (cm^3) and GTV fraction (%) of the diffusion-
defined resistant sub-volume at each treatment timepoint, separately for
patients with local control (LC) and local failure (LF).  A Mann-Whitney
U test compares the LC and LF distributions when scipy is available.
"""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _sig_class,
)


def _fmt_mean_std(mean: float | None, std: float | None, scale: float = 1.0) -> str:
    """Render a 'mean ± std' cell; em-dash when both are missing."""
    if mean is None:
        return "&mdash;"
    m = mean * scale
    if std is None:
        return f"{m:.2f}"
    return f"{m:.2f} &plusmn; {(std * scale):.2f}"


def _fmt_pct(val: float | None) -> str:
    """Render a fraction (0-1 or already 0-100) as a percent string."""
    if val is None:
        return "&mdash;"
    pct = val * 100.0 if 0.0 <= val <= 1.0 else float(val)
    return f"{pct:.1f}%"


def _fmt_p(p: float | None) -> str:
    """Render a p-value cell, classed when significant.

    Returns an em-dash when p is None (e.g. scipy unavailable).
    """
    if p is None:
        return "<td>&mdash;</td>"
    try:
        pv = float(p)
    except (TypeError, ValueError):
        return "<td>&mdash;</td>"
    cls = _sig_class(pv)
    cls_attr = f' class="{cls}"' if cls else ""
    txt = f"{pv:.3g}" if pv >= 0.001 else f"{pv:.1e}"
    return f"<td{cls_attr}>{_esc(txt)}</td>"


def _section_subvolume_sizes(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build HTML section showing sub-volume sizes per timepoint / outcome."""
    if not mat_data:
        return []

    dwi_with_data: list[str] = []
    for dwi in DWI_TYPES:
        entry = (mat_data or {}).get(dwi) or {}
        svs = entry.get("subvolume_sizes")
        if svs and svs.get("timepoints"):
            dwi_with_data.append(dwi)

    if not dwi_with_data:
        return []

    h: list[str] = []
    h.append(_h2("Sub-Volume Sizes", "subvol-sizes"))
    h.append(
        "<p>Diffusion-defined resistant sub-volume at each treatment "
        "timepoint, stratified by local control (LC) vs local failure (LF) "
        "outcome. Mean &plusmn; SD of absolute volume (cm\u00b3) and fraction "
        "of the primary GTV. The p-value is a Mann-Whitney U test comparing "
        "LC vs LF at each timepoint.</p>"
    )

    small_fraction_warn = False

    for dwi in dwi_with_data:
        entry = mat_data[dwi]
        svs = entry["subvolume_sizes"]
        tps = svs.get("timepoints") or []
        lc = svs.get("lc") or {}
        lf = svs.get("lf") or {}
        ps = svs.get("wilcoxon_p") or [None] * len(tps)
        lc_vol = lc.get("mean_vol_cm3") or []
        lc_svol = lc.get("std_vol_cm3") or []
        lf_vol = lf.get("mean_vol_cm3") or []
        lf_svol = lf.get("std_vol_cm3") or []
        lc_frac = lc.get("mean_frac_pct") or []
        lf_frac = lf.get("mean_frac_pct") or []

        h.append(f"<h3>{_dwi_badge(dwi)}</h3>")
        h.append("<table><thead><tr>"
                 "<th>Timepoint</th>"
                 "<th>LC Vol (cm\u00b3)</th>"
                 "<th>LF Vol (cm\u00b3)</th>"
                 "<th>LC Frac (%)</th>"
                 "<th>LF Frac (%)</th>"
                 "<th>p-value</th>"
                 "</tr></thead><tbody>")

        for i, tp in enumerate(tps):
            def _get(lst, idx):
                try:
                    return lst[idx]
                except (IndexError, TypeError):
                    return None

            lcv, lcs = _get(lc_vol, i), _get(lc_svol, i)
            lfv, lfs = _get(lf_vol, i), _get(lf_svol, i)
            lcf, lff = _get(lc_frac, i), _get(lf_frac, i)
            p = _get(ps, i)

            # Normalise fractions to 0-100 scale for the 5% clinical
            # significance check.
            def _pct(val):
                if val is None:
                    return None
                return val * 100.0 if 0.0 <= val <= 1.0 else float(val)

            lcf_pct = _pct(lcf)
            lff_pct = _pct(lff)
            for frac in (lcf_pct, lff_pct):
                if frac is not None and frac < 5.0:
                    small_fraction_warn = True

            h.append("<tr>")
            h.append(f"<td><strong>Fx{_esc(str(tp))}</strong></td>")
            h.append(f"<td>{_fmt_mean_std(lcv, lcs)}</td>")
            h.append(f"<td>{_fmt_mean_std(lfv, lfs)}</td>")
            h.append(f"<td>{_fmt_pct(lcf)}</td>")
            h.append(f"<td>{_fmt_pct(lff)}</td>")
            h.append(_fmt_p(p))
            h.append("</tr>")

        h.append("</tbody></table>")

    if small_fraction_warn:
        h.append(
            '<p class="meta"><strong>Note:</strong> Sub-volumes &lt; 5% of '
            "GTV may lack clinical significance for dose escalation.</p>"
        )

    if all(
        all(p is None for p in (mat_data[d].get("subvolume_sizes", {}).get("wilcoxon_p") or []))
        for d in dwi_with_data
    ):
        h.append(
            '<p class="meta"><em>Mann-Whitney U test not available '
            "(scipy not installed).</em></p>"
        )

    return h
