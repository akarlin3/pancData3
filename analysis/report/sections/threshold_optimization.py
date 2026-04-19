"""Report section: ADC threshold optimisation for sub-volume definition.

Displays median Dice, mean Dice, and volume fraction for each candidate
ADC threshold, highlighting the optimal threshold that maximises spatial
reproducibility between Fx1 repeat scans.
"""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


_DEFAULT_THRESH = 0.0010  # current config default (adc_thresh in config.json)


def _scalar(val) -> float | None:
    """Return a plain float or None."""
    if val is None:
        return None
    try:
        f = float(val)
        if f != f:  # NaN
            return None
        return f
    except (TypeError, ValueError):
        return None


def _section_threshold_optimization(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build HTML section with threshold-sweep results per DWI type."""
    if not mat_data:
        return []

    dwi_with_data: list[str] = []
    for dwi in DWI_TYPES:
        entry = (mat_data or {}).get(dwi) or {}
        opt = entry.get("threshold_optimization")
        if not opt or not opt.get("thresholds"):
            continue
        # Require at least one valid Dice value and at least one patient
        # with repeats; otherwise the sweep produced no usable data.
        n_pts = opt.get("n_patients") or 0
        md = opt.get("median_dice") or []
        has_real_dice = any(v is not None for v in md)
        if n_pts > 0 and has_real_dice:
            dwi_with_data.append(dwi)

    if not dwi_with_data:
        return []

    h: list[str] = []
    h.append(_h2("ADC Threshold Optimisation", "threshold-optimization"))
    h.append(
        "<p>Sweep of candidate ADC thresholds against median Dice between "
        "paired Fx1 repeat scans. The optimal threshold maximises spatial "
        "reproducibility of the resistant sub-volume while keeping the "
        "volume fraction clinically meaningful.</p>"
    )

    for dwi in dwi_with_data:
        entry = mat_data[dwi]
        opt = entry["threshold_optimization"]
        thresholds = opt.get("thresholds") or []
        median_dice = opt.get("median_dice") or []
        mean_dice = opt.get("mean_dice") or []
        std_dice = opt.get("std_dice") or []
        vol_frac = opt.get("median_vol_frac") or []
        opt_thresh = _scalar(opt.get("optimal_thresh"))
        opt_dice = _scalar(opt.get("optimal_dice"))
        opt_vf = _scalar(opt.get("optimal_vol_frac"))
        n_pts = opt.get("n_patients")

        header_suffix = f" (n={n_pts} patients with repeats)" if n_pts else ""
        h.append(f"<h3>{_dwi_badge(dwi)}{_esc(header_suffix)}</h3>")
        h.append("<table><thead><tr>"
                 "<th>Threshold (&times;10<sup>&minus;3</sup>)</th>"
                 "<th>Median Dice</th>"
                 "<th>Mean Dice &plusmn; SD</th>"
                 "<th>Vol Frac (%)</th>"
                 "<th></th>"
                 "</tr></thead><tbody>")

        def _get(lst, idx):
            try:
                return lst[idx]
            except (IndexError, TypeError):
                return None

        for i, t in enumerate(thresholds):
            t_val = _scalar(t)
            md = _scalar(_get(median_dice, i))
            avg = _scalar(_get(mean_dice, i))
            sd = _scalar(_get(std_dice, i))
            vf = _scalar(_get(vol_frac, i))

            is_optimal = (
                opt_thresh is not None
                and t_val is not None
                and abs(t_val - opt_thresh) < 1e-9
            )
            row_style = ""
            marker = ""
            if is_optimal:
                row_style = ' style="background:#d4edda;font-weight:bold"'
                marker = "&larr; Optimal"

            t_str = f"{(t_val * 1000):.2f}" if t_val is not None else "&mdash;"
            md_str = f"{md:.3f}" if md is not None else "&mdash;"
            if avg is not None:
                avg_str = f"{avg:.3f} &plusmn; {sd:.3f}" if sd is not None else f"{avg:.3f}"
            else:
                avg_str = "&mdash;"
            vf_str = f"{(vf * 100):.1f}%" if vf is not None else "&mdash;"

            h.append(
                f"<tr{row_style}><td>{t_str}</td><td>{md_str}</td>"
                f"<td>{avg_str}</td><td>{vf_str}</td><td>{marker}</td></tr>"
            )

        h.append("</tbody></table>")

        # Comparison box: current default vs optimal.
        # Estimate current-default Dice / vol_frac by locating the threshold
        # closest to _DEFAULT_THRESH in the sweep.
        default_dice = None
        default_vf = None
        if thresholds:
            closest_i = None
            closest_diff = None
            for i, t in enumerate(thresholds):
                t_val = _scalar(t)
                if t_val is None:
                    continue
                diff = abs(t_val - _DEFAULT_THRESH)
                if closest_diff is None or diff < closest_diff:
                    closest_diff = diff
                    closest_i = i
            if closest_i is not None:
                default_dice = _scalar(_get(median_dice, closest_i))
                default_vf = _scalar(_get(vol_frac, closest_i))

        if opt_thresh is not None and opt_dice is not None:
            h.append('<div class="summary-box">')
            h.append("<ul>")
            if default_dice is not None:
                dvf_str = f", Volume fraction = {(default_vf * 100):.1f}%" if default_vf is not None else ""
                h.append(
                    f"<li>Current default ({_DEFAULT_THRESH:.4f}): Dice = "
                    f"{default_dice:.3f}{dvf_str}</li>"
                )
            ovf_str = f", Volume fraction = {(opt_vf * 100):.1f}%" if opt_vf is not None else ""
            h.append(
                f"<li>Optimal ({opt_thresh:.4f}): Dice = {opt_dice:.3f}{ovf_str}</li>"
            )
            if default_dice is not None and opt_dice is not None:
                d_dice = opt_dice - default_dice
                if default_dice > 0:
                    pct = 100.0 * d_dice / default_dice
                    h.append(
                        f"<li>Improvement: {d_dice:+.3f} Dice ({pct:+.1f}%)</li>"
                    )
                else:
                    h.append(f"<li>Improvement: {d_dice:+.3f} Dice</li>")
            h.append("</ul>")
            h.append("</div>")

    return h
