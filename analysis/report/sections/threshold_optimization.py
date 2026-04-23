"""Report sections: ADC threshold optimisation (three tactics).

The MATLAB optimize_adc_threshold sweeps 13 candidate ADC thresholds
and reports an optimum from each of three independent tactics:

  Tactic 1 — Reproducibility (max median Fx1 repeat Dice).  Most
             test-retest stable cut.
  Tactic 2 — Volume inflection (knee of vol_frac curve).  Biological
             transition between dense tumour core and surrounding tissue.
  Tactic 3 — Outcome significance (min Wilcoxon p of vol_frac LC vs LF).
             Best at separating treatment outcomes.

Each tactic renders as its own H2 section so the TOC carries one entry
per rationale and the analyst can compare them side by side.
"""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


_DEFAULT_THRESH = 0.0010  # current config default (adc_thresh in config.json)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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


def _get(lst, idx):
    try:
        return lst[idx]
    except (IndexError, TypeError):
        return None


def _dwi_with_data(mat_data: dict, predicate) -> list[str]:
    """Return DWI types whose threshold_optimization entry passes `predicate`."""
    out: list[str] = []
    for dwi in DWI_TYPES:
        entry = (mat_data or {}).get(dwi) or {}
        opt = entry.get("threshold_optimization")
        if not opt or not opt.get("thresholds"):
            continue
        if predicate(opt):
            out.append(dwi)
    return out


def _format_p(p: float | None) -> str:
    if p is None:
        return "&mdash;"
    if p < 0.001:
        return "<0.001"
    if p < 0.01:
        return f"{p:.3f}"
    return f"{p:.2f}"


# ---------------------------------------------------------------------------
# Tactic 1 — Reproducibility (Dice)
# ---------------------------------------------------------------------------

def _section_threshold_dice(mat_data: dict) -> list[str]:
    """Tactic 1: max median Fx1 repeat Dice."""
    def has_dice(opt):
        n_pts = opt.get("n_patients") or 0
        md = opt.get("median_dice") or []
        return n_pts > 0 and any(v is not None for v in md)

    dwi_types = _dwi_with_data(mat_data, has_dice)
    if not dwi_types:
        return []

    h: list[str] = [
        _h2("ADC Threshold Optimisation — Reproducibility (Tactic 1)",
            "threshold-optimization-dice"),
        "<p>Sweep of candidate ADC thresholds against median Dice between "
        "paired Fx1 repeat scans. The optimal threshold by this tactic "
        "maximises spatial reproducibility of the resistant sub-volume — "
        "the cut that is most stable on test-retest.</p>",
    ]

    for dwi in dwi_types:
        opt = mat_data[dwi]["threshold_optimization"]
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
            row_style = ' style="background:#d4edda;font-weight:bold"' if is_optimal else ""
            marker = "&larr; Optimal (Dice)" if is_optimal else ""

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
            h.append('<div class="summary-box"><ul>')
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
            h.append("</ul></div>")

    return h


# ---------------------------------------------------------------------------
# Tactic 2 — Volume inflection
# ---------------------------------------------------------------------------

def _section_threshold_inflection(mat_data: dict) -> list[str]:
    """Tactic 2: knee of median vol_frac vs threshold."""
    def has_inflection(opt):
        return _scalar(opt.get("inflection_thresh")) is not None

    dwi_types = _dwi_with_data(mat_data, has_inflection)
    if not dwi_types:
        return []

    h: list[str] = [
        _h2("ADC Threshold Optimisation — Volume Inflection (Tactic 2)",
            "threshold-optimization-inflection"),
        "<p>The median sub-volume fraction grows monotonically with the "
        "ADC threshold. Its <em>knee</em> &mdash; the point of maximum "
        "concave curvature in the smoothed vol_frac vs threshold curve "
        "&mdash; marks the biological transition from dense tumour core "
        "to surrounding tissue. Beyond this threshold each additional "
        "0.0001 mm&sup2;/s adds sharply less voxels.</p>",
        "<p>Inflection is detected by 3-point smoothing of "
        "<code>median_vol_frac</code> followed by a discrete second "
        "derivative; the index that minimises "
        "<code>d&sup2;V/dt&sup2;</code> (most negative curvature) is "
        "reported as the knee. Endpoints have no defined curvature.</p>",
    ]

    for dwi in dwi_types:
        opt = mat_data[dwi]["threshold_optimization"]
        thresholds = opt.get("thresholds") or []
        vol_frac = opt.get("median_vol_frac") or []
        curvature = opt.get("vol_frac_curvature") or []
        inflect = _scalar(opt.get("inflection_thresh"))
        inflect_curv = _scalar(opt.get("inflection_curvature"))
        n_pts = opt.get("n_patients")

        header_suffix = f" (n={n_pts} patients with repeats)" if n_pts else ""
        h.append(f"<h3>{_dwi_badge(dwi)}{_esc(header_suffix)}</h3>")
        h.append("<table><thead><tr>"
                 "<th>Threshold (&times;10<sup>&minus;3</sup>)</th>"
                 "<th>Vol Frac (%)</th>"
                 "<th>Curvature (d&sup2;V/dt&sup2;)</th>"
                 "<th></th>"
                 "</tr></thead><tbody>")

        for i, t in enumerate(thresholds):
            t_val = _scalar(t)
            vf = _scalar(_get(vol_frac, i))
            cv = _scalar(_get(curvature, i))

            is_knee = (
                inflect is not None
                and t_val is not None
                and abs(t_val - inflect) < 1e-9
            )
            row_style = ' style="background:#e1d5f5;font-weight:bold"' if is_knee else ""
            marker = "&larr; Knee" if is_knee else ""

            t_str = f"{(t_val * 1000):.2f}" if t_val is not None else "&mdash;"
            vf_str = f"{(vf * 100):.1f}%" if vf is not None else "&mdash;"
            cv_str = f"{cv:+.4f}" if cv is not None else "&mdash;"

            h.append(
                f"<tr{row_style}><td>{t_str}</td><td>{vf_str}</td>"
                f"<td>{cv_str}</td><td>{marker}</td></tr>"
            )
        h.append("</tbody></table>")

        if inflect is not None:
            curv_str = f" (curvature = {inflect_curv:+.4f})" if inflect_curv is not None else ""
            h.append('<div class="summary-box"><ul>')
            h.append(
                f"<li>Inflection threshold: <strong>{inflect:.4f}</strong>{curv_str}</li>"
            )
            h.append("<li>Interpretation: voxels with ADC below this value "
                     "are most likely dense tumour core; voxels above "
                     "increasingly include surrounding healthy or "
                     "peri-tumoural tissue.</li>")
            h.append("</ul></div>")

    return h


# ---------------------------------------------------------------------------
# Tactic 3 — Outcome significance
# ---------------------------------------------------------------------------

def _section_threshold_significance(mat_data: dict) -> list[str]:
    """Tactic 3: min Wilcoxon p-value of vol_frac between LC and LF groups."""
    def has_significance(opt):
        return _scalar(opt.get("significance_thresh")) is not None

    dwi_types = _dwi_with_data(mat_data, has_significance)
    if not dwi_types:
        return []

    h: list[str] = [
        _h2("ADC Threshold Optimisation — Outcome Significance (Tactic 3)",
            "threshold-optimization-significance"),
        "<p>Per-threshold Wilcoxon rank-sum test of the per-patient "
        "sub-volume fraction between local-control (LC) and local-failure "
        "(LF) patients. The optimal threshold by this tactic is the one "
        "whose threshold-defined sub-volume best discriminates outcomes "
        "&mdash; i.e., the cut at which the difference in sub-volume "
        "fraction between LC and LF is most statistically significant.</p>",
        "<p>P-values shown for each threshold; thresholds with fewer than "
        "3 patients in either group are skipped (insufficient power for "
        "Wilcoxon). Significance is uncorrected for multiple comparisons "
        "across the 13 thresholds &mdash; treat as exploratory ranking, "
        "not confirmatory inference.</p>",
    ]

    for dwi in dwi_types:
        opt = mat_data[dwi]["threshold_optimization"]
        thresholds = opt.get("thresholds") or []
        pvalues = opt.get("significance_pvalues") or []
        sig_thresh = _scalar(opt.get("significance_thresh"))
        sig_p = _scalar(opt.get("significance_pvalue"))
        n_lc = opt.get("significance_n_lc")
        n_lf = opt.get("significance_n_lf")
        metric = opt.get("significance_metric") or "wilcoxon ranksum on vol_frac"

        header_suffix = f" (test: {metric})"
        h.append(f"<h3>{_dwi_badge(dwi)}{_esc(header_suffix)}</h3>")
        h.append("<table><thead><tr>"
                 "<th>Threshold (&times;10<sup>&minus;3</sup>)</th>"
                 "<th>p-value</th>"
                 "<th></th>"
                 "</tr></thead><tbody>")

        for i, t in enumerate(thresholds):
            t_val = _scalar(t)
            p = _scalar(_get(pvalues, i))

            is_best = (
                sig_thresh is not None
                and t_val is not None
                and abs(t_val - sig_thresh) < 1e-9
            )
            row_style = ' style="background:#fde2c4;font-weight:bold"' if is_best else ""
            marker = "&larr; Most significant" if is_best else ""

            t_str = f"{(t_val * 1000):.2f}" if t_val is not None else "&mdash;"
            p_str = _format_p(p)
            if p is not None and p < 0.05 and not is_best:
                p_str = f"<strong>{p_str}</strong>"

            h.append(
                f"<tr{row_style}><td>{t_str}</td><td>{p_str}</td>"
                f"<td>{marker}</td></tr>"
            )
        h.append("</tbody></table>")

        if sig_thresh is not None and sig_p is not None:
            n_str = ""
            if n_lc is not None and n_lf is not None:
                n_str = f" (n_LC = {n_lc}, n_LF = {n_lf})"
            h.append('<div class="summary-box"><ul>')
            h.append(
                f"<li>Most significant threshold: <strong>{sig_thresh:.4f}</strong>"
                f" with p = {_format_p(sig_p)}{n_str}</li>"
            )
            note = ("<li>Below the conventional 0.05 cutoff &mdash; "
                    "exploratory evidence of outcome separation at this cut.</li>"
                    if sig_p < 0.05
                    else "<li>Above the conventional 0.05 cutoff &mdash; no single "
                         "threshold in the swept range produced a significant "
                         "split of LC vs LF in this cohort.</li>")
            h.append(note)
            h.append("</ul></div>")

    return h


# ---------------------------------------------------------------------------
# Backwards-compatible top-level entry point
# ---------------------------------------------------------------------------

def _section_threshold_optimization(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Render all three tactic sections in priority order.

    `dwi_types` is unused (kept for signature compatibility with the prior
    single-section interface); each tactic decides which DWI types have
    enough data to render itself.
    """
    if not mat_data:
        return []
    h: list[str] = []
    h.extend(_section_threshold_dice(mat_data))
    h.extend(_section_threshold_inflection(mat_data))
    h.extend(_section_threshold_significance(mat_data))
    return h
