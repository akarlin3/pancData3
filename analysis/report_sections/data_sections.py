"""Report sections: data sections group."""

from __future__ import annotations

import json
import re
from pathlib import Path

from shared import (
    DWI_TYPES,
    extract_correlations,
    extract_pvalues,
    parse_dwi_info,
)
from report_formatters import (
    _dwi_badge,
    _esc,
    _figure_caption,
    _h2,
    _sig_class,
    _stat_card,
    _table_caption,
    _trend_tag,
)


def _section_cohort_overview(mat_data, log_data, dwi_types_present) -> list[str]:
    """Build the Cohort Overview section.

    Displays patient counts, timepoint counts, and data quality
    summary across all DWI types from MAT file and log data.

    Parameters
    ----------
    mat_data : dict
        Parsed MAT file metrics.
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.

    Returns
    -------
    list[str]
        HTML chunks for the cohort overview section.
    """
    h: list[str] = []
    has_data = False

    # Check if we have any cohort data to show
    if mat_data:
        for dt in DWI_TYPES:
            if dt in mat_data and "longitudinal" in mat_data[dt]:
                has_data = True
                break

    if log_data:
        for dt in dwi_types_present:
            if dt in log_data:
                bl = log_data[dt].get("baseline", {})
                if bl.get("total_outliers") or bl.get("baseline_exclusion"):
                    has_data = True
                    break

    if not has_data:
        return h

    h.append(_h2("Cohort Overview", "cohort"))

    # Longitudinal data summary from MAT files
    if mat_data:
        cohort_rows = []
        for dt in DWI_TYPES:
            if dt not in mat_data or "longitudinal" not in mat_data[dt]:
                continue
            lon = mat_data[dt]["longitudinal"]
            n_pat = lon.get("num_patients", 0)
            n_tp = lon.get("num_timepoints", 0)
            if n_pat > 0 or n_tp > 0:
                cohort_rows.append((dt, n_pat, n_tp))

        if cohort_rows:
            h.append("<h3>Longitudinal Data Dimensions</h3>")
            h.append("<table><thead><tr><th>DWI Type</th><th>Patients</th>"
                     "<th>Timepoints</th></tr></thead><tbody>")
            for dt, n_pat, n_tp in cohort_rows:
                h.append(f"<tr><td>{_dwi_badge(dt)}</td>"
                         f"<td><strong>{n_pat}</strong></td>"
                         f"<td>{n_tp}</td></tr>")
            h.append("</tbody></table>")

            # Warn about small sample sizes that affect statistical power.
            max_patients = max(r[1] for r in cohort_rows)
            if max_patients < 30:
                h.append(
                    '<div class="warn-box">'
                    f"\u26a0\ufe0f <strong>Small sample size (n\u2009=\u2009{max_patients}):</strong> "
                    "Results should be interpreted with caution. With fewer than 30 patients, "
                    "statistical power is limited and effect size estimates may be imprecise. "
                    "Wide confidence intervals are expected, and non-significant results do not "
                    "rule out clinically meaningful effects. Validation in a larger cohort is "
                    "recommended before clinical application."
                    "</div>"
                )

    # Cross-DWI data quality summary
    if log_data:
        quality_rows = []
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            bl = log_data[dt].get("baseline", {})
            total_out = bl.get("total_outliers")
            baseline_exc = bl.get("baseline_exclusion")
            if total_out or baseline_exc:
                out_str = f"{total_out['n_removed']}/{total_out['n_total']} ({total_out['pct']:.1f}%)" if total_out else "\u2014"
                exc_str = f"{baseline_exc['n_excluded']}/{baseline_exc['n_total']}" if baseline_exc else "\u2014"
                lf_inc = baseline_exc.get("lf_rate_included") if baseline_exc else None
                lf_exc = baseline_exc.get("lf_rate_excluded") if baseline_exc else None
                lf_str = f"{lf_inc:.1f}% / {lf_exc:.1f}%" if lf_inc is not None and lf_exc is not None else "\u2014"
                quality_rows.append((dt, out_str, exc_str, lf_str))

        if quality_rows:
            h.append("<h3>Data Quality Summary</h3>")
            h.append("<table><thead><tr><th>DWI Type</th><th>Outliers Removed</th>"
                     "<th>Baseline Excluded</th><th>LF Rate (incl/excl)</th>"
                     "</tr></thead><tbody>")
            for dt, out_str, exc_str, lf_str in quality_rows:
                h.append(f"<tr><td>{_dwi_badge(dt)}</td>"
                         f"<td>{_esc(out_str)}</td>"
                         f"<td>{_esc(exc_str)}</td>"
                         f"<td>{_esc(lf_str)}</td></tr>")
            h.append("</tbody></table>")

    return h



def _section_patient_flow(log_data, dwi_types_present, mat_data) -> list[str]:
    """Build a CONSORT-style patient flow summary for publication.

    Summarises patient attrition from cohort through each exclusion stage:
    initial cohort, baseline exclusion, outlier removal, and competing-risk
    exclusion.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    mat_data : dict
        Parsed MAT file metrics.

    Returns
    -------
    list[str]
        HTML chunks for the patient flow section.
    """
    h: list[str] = []

    # Gather flow data across DWI types
    flow_rows: list[tuple[str, int, int, int, int, float]] = []
    for dt in dwi_types_present:
        n_initial = 0
        n_baseline_exc = 0
        n_outlier_exc = 0
        n_cr_exc = 0
        outlier_pct = 0.0

        if mat_data and dt in mat_data and "longitudinal" in mat_data[dt]:
            n_initial = mat_data[dt]["longitudinal"].get("num_patients", 0)

        if log_data and dt in log_data:
            bl = log_data[dt].get("baseline", {})
            baseline_exc = bl.get("baseline_exclusion")
            if baseline_exc:
                n_baseline_exc = baseline_exc.get("n_excluded", 0)
                if n_initial == 0:
                    n_initial = baseline_exc.get("n_total", 0)

            total_out = bl.get("total_outliers")
            if total_out:
                n_outlier_exc = total_out.get("n_removed", 0)
                outlier_pct = total_out.get("pct", 0)

            sc = log_data[dt].get("stats_comparisons", {})
            glme_excl = sc.get("glme_excluded")
            if glme_excl:
                n_cr_exc = glme_excl.get("n_excluded", 0)

        if n_initial > 0 or n_baseline_exc > 0:
            flow_rows.append((dt, n_initial, n_baseline_exc, n_outlier_exc,
                              n_cr_exc, outlier_pct))

    if not flow_rows:
        return h

    h.append(_h2("Patient Flow", "patient-flow"))
    h.append(
        '<p class="meta">CONSORT-style summary of patient inclusion and exclusion '
        'at each analysis stage. Numbers may differ across DWI types due to '
        'type-specific data availability and model convergence.</p>'
    )

    h.append("<table>")
    h.append(_table_caption(
        "Patient Flow Summary",
        "Attrition from initial cohort through analysis stages."))
    h.append("<thead><tr>"
             "<th>DWI Type</th>"
             "<th>Initial Cohort</th>"
             "<th>Baseline Excluded</th>"
             "<th>Outliers Removed</th>"
             "<th>CR Excluded (GLME)</th>"
             "<th>Analysed (GLME)</th>"
             "</tr></thead><tbody>")
    for dt, n_init, n_bl, n_out, n_cr, out_pct in flow_rows:
        n_analysed = n_init - n_bl - n_cr
        if n_analysed < 0:
            n_analysed = n_init - n_bl  # fallback if CR data missing
        h.append(
            f"<tr><td>{_dwi_badge(dt)}</td>"
            f"<td><strong>{n_init}</strong></td>"
            f"<td>{n_bl}</td>"
            f"<td>{n_out} ({out_pct:.1f}%)</td>"
            f"<td>{n_cr}</td>"
            f"<td><strong>{max(n_analysed, 0)}</strong></td></tr>"
        )
    h.append("</tbody></table>")

    return h



def _section_data_completeness(log_data, dwi_types_present) -> list[str]:
    """Build the Data Completeness section from sanity check log data.

    Reports convergence issues, outlier detections, dimensional alignment
    results, and excessive NaN warnings extracted from the sanity_checks
    log files.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics (must include ``sanity_checks`` key per DWI type).
    dwi_types_present : list[str]
        DWI types found in this pipeline run.

    Returns
    -------
    list[str]
        HTML chunks for the data completeness section.
    """
    h: list[str] = []
    if not log_data:
        return h

    has_data = False
    for dt in dwi_types_present:
        if dt in log_data and log_data[dt].get("sanity_checks"):
            san = log_data[dt]["sanity_checks"]
            if (san.get("total_convergence", 0) > 0
                    or san.get("all_converged")
                    or san.get("dim_mismatches", 0) > 0
                    or san.get("excessive_nan")):
                has_data = True
                break

    if not has_data:
        return h

    h.append(_h2("Data Completeness", "data-completeness"))
    h.append(
        '<p class="meta">Summary of data quality checks from the pipeline\'s '
        'sanity_checks module, covering voxel-level convergence, outlier '
        'detection, and dose\u2013DWI dimensional alignment.</p>'
    )

    for dt in dwi_types_present:
        if dt not in log_data:
            continue
        san = log_data[dt].get("sanity_checks", {})
        if not san:
            continue

        items_to_show = []

        # Convergence summary
        if san.get("all_converged"):
            items_to_show.append(("Convergence", "All converged", "info-box",
                                  "All voxel-level fit values are finite, "
                                  "non-NaN, and non-negative."))
        elif san.get("total_convergence", 0) > 0:
            items_to_show.append(("Convergence Flags",
                                  str(san["total_convergence"]),
                                  "warn-box",
                                  f"{san['total_convergence']} convergence "
                                  f"flag(s) raised (Inf/NaN/Neg voxels). "
                                  f"Review affected scans for fitting failures."))

        # Dimensional alignment
        if san.get("dim_mismatches", 0) > 0:
            items_to_show.append(("Dim. Mismatches",
                                  str(san["dim_mismatches"]),
                                  "warn-box",
                                  f"{san['dim_mismatches']} dose/DWI pairs have "
                                  f"mismatched voxel counts. Dosimetry metrics may "
                                  f"be unreliable for these scans."))

        if san.get("nan_dose_warnings", 0) > 0:
            items_to_show.append(("NaN Dose Warnings",
                                  str(san["nan_dose_warnings"]),
                                  "warn-box",
                                  f"{san['nan_dose_warnings']} scan(s) have >10% "
                                  f"NaN in-mask dose voxels (partial RT dose coverage)."))

        # Excessive NaN parameters
        if san.get("excessive_nan"):
            for en in san["excessive_nan"]:
                items_to_show.append(("Excessive NaN",
                                      f"{en['parameter']}: {en['pct_nan']:.1f}%",
                                      "warn-box",
                                      f"Parameter {en['parameter']} has "
                                      f"{en['pct_nan']:.1f}% NaN values "
                                      f"(threshold: 50%). Check GTV masks and "
                                      f"model fitting for this parameter."))

        if not items_to_show:
            continue

        h.append(f"<h3>{_dwi_badge(dt)}</h3>")

        # Stat cards for quick overview
        cards = []
        if san.get("all_converged"):
            cards.append(_stat_card("Convergence", "Passed"))
        elif san.get("total_convergence", 0) > 0:
            cards.append(_stat_card("Conv. Flags",
                                    str(san["total_convergence"]),
                                    "Inf/NaN/Neg issues"))
        if san.get("total_outliers", 0) > 0:
            cards.append(_stat_card("Sanity Outliers",
                                    str(san["total_outliers"]),
                                    ">3 IQR from median"))
        dm = san.get("dim_mismatches", 0)
        nw = san.get("nan_dose_warnings", 0)
        if dm == 0 and nw == 0:
            cards.append(_stat_card("Alignment", "Passed"))
        else:
            cards.append(_stat_card("Alignment Issues",
                                    str(dm + nw),
                                    f"{dm} mismatch, {nw} NaN"))
        if cards:
            h.append('<div class="stat-grid">')
            h.extend(cards)
            h.append("</div>")

        # Detailed items
        for label, value, box_cls, desc in items_to_show:
            h.append(f'<div class="{box_cls}">{_esc(desc)}</div>')

    return h



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

        has_dos = any("dosimetry" in d for d in mat_data.values())
        if has_dos:
            h.append("<h3>Dosimetry (Target Coverage)</h3>")
            h.append("<table><thead><tr><th>DWI</th><th>D95 ADC (Gy)</th><th>V50 ADC (%)</th>"
                     "<th>D95 D (Gy)</th><th>V50 D (%)</th></tr></thead><tbody>")
            for dt in DWI_TYPES:
                if dt not in mat_data or "dosimetry" not in mat_data[dt]:
                    continue
                dos = mat_data[dt]["dosimetry"]
                if not dos:
                    continue

                def _fmt_gy(val):
                    """Format a dose value in Gy, handling missing data."""
                    if isinstance(val, (int, float)) and val == val:  # NaN check
                        return f"{val:.2f}"
                    return "N/A"

                def _fmt_pct(val):
                    """Format a fractional value as percentage, handling missing data."""
                    if isinstance(val, (int, float)) and val == val:  # NaN check
                        pct = val * 100 if val <= 1.0 else val
                        return f"{pct:.1f}%"
                    return "N/A"

                h.append(
                    f"<tr><td>{_dwi_badge(dt)}</td>"
                    f"<td>{_fmt_gy(dos.get('d95_adc_mean'))}</td>"
                    f"<td>{_fmt_pct(dos.get('v50_adc_mean'))}</td>"
                    f"<td>{_fmt_gy(dos.get('d95_d_mean'))}</td>"
                    f"<td>{_fmt_pct(dos.get('v50_d_mean'))}</td></tr>"
                )
            h.append("</tbody></table>")

            # Clinical interpretation of dosimetry values
            dos_notes = []
            for dt in DWI_TYPES:
                if dt not in mat_data or "dosimetry" not in mat_data[dt]:
                    continue
                dos = mat_data[dt]["dosimetry"]
                if not dos:
                    continue
                d95_adc = dos.get("d95_adc_mean")
                v50_adc = dos.get("v50_adc_mean")
                if isinstance(d95_adc, (int, float)) and d95_adc == d95_adc:
                    if d95_adc < 45:
                        dos_notes.append(f"{dt}: D95 ADC = {d95_adc:.1f} Gy "
                                         "(below 45 Gy \u2014 potential under-dosing of resistant sub-volume)")
                    elif d95_adc >= 50:
                        dos_notes.append(f"{dt}: D95 ADC = {d95_adc:.1f} Gy "
                                         "(adequate coverage of resistant sub-volume)")
                if isinstance(v50_adc, (int, float)) and v50_adc == v50_adc:
                    v50_pct = v50_adc * 100 if v50_adc <= 1.0 else v50_adc
                    if v50_pct < 90:
                        dos_notes.append(f"{dt}: V50 ADC = {v50_pct:.0f}% "
                                         "(partial coverage \u2014 may benefit from dose escalation)")
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
                h.append('<div style="overflow-x:auto">')
                h.append("<table><thead><tr><th>Method</th>")
                for m in methods:
                    h.append(f"<th>{_esc(m)}</th>")
                h.append("</tr></thead><tbody>")
                for i in range(n):
                    h.append(f"<tr><td><strong>{_esc(methods[i])}</strong></td>")
                    for j in range(n):
                        val = matrix[i][j] if i < len(matrix) and j < len(matrix[i]) else None
                        if isinstance(val, (int, float)):
                            # Colour-code: high Dice = green-ish
                            if i != j:
                                style = f" style=\"color: {'var(--green)' if val >= 0.7 else ('var(--amber)' if val >= 0.5 else 'var(--red)')}; font-weight:600\""
                            else:
                                style = ""
                            h.append(f"<td{style}>{val:.2f}</td>")
                        else:
                            h.append("<td>-</td>")
                    h.append("</tr>")
                h.append("</tbody></table>")
                h.append("</div>")

                # Summary statistics for core method agreement
                off_diag = []
                for i in range(n):
                    for j in range(i + 1, n):
                        if i < len(matrix) and j < len(matrix[i]):
                            val = matrix[i][j]
                            if isinstance(val, (int, float)) and val > 0:
                                off_diag.append(val)
                if off_diag:
                    avg_dice = sum(off_diag) / len(off_diag)
                    min_dice = min(off_diag)
                    max_dice = max(off_diag)
                    h.append(f'<div class="info-box">Mean pairwise Dice: <strong>{avg_dice:.3f}</strong> '
                             f'(range: {min_dice:.3f}\u2013{max_dice:.3f} across {len(off_diag)} pairs)</div>')
    return h



def _section_appendix(rows) -> list[str]:
    """Build the Appendix: All Graphs section.

    Lists every analysed graph in a detailed table with columns for
    DWI type, graph type, title, axis details, trend tags, and a
    collapsible full summary.

    Parameters
    ----------
    rows : list[dict]
        Vision CSV rows (may be empty, in which case the section is skipped).

    Returns
    -------
    list[str]
        HTML chunks for the appendix section.
    """
    # ── 8. Appendix ──
    h: list[str] = []
    if rows:
        h.append(_h2("Appendix: All Graphs", "appendix"))
        h.append(f"<p>{len(rows)} graphs analysed. Each graph is shown with full "
                 "metadata, extracted statistics, and clinical interpretation notes.</p>")

        # Group by graph type for organised presentation
        type_groups: dict[str, list[tuple[int, dict]]] = {}
        for i, r in enumerate(rows, 1):
            gt = r.get("graph_type", "unknown")
            type_groups.setdefault(gt, []).append((i, r))

        for gt in sorted(type_groups.keys(), key=lambda k: -len(type_groups[k])):
            gt_rows = type_groups[gt]
            h.append(f"<h3>{_esc(gt)} ({len(gt_rows)} graphs)</h3>")

            h.append("<table><thead><tr>"
                     "<th>#</th><th>DWI</th><th>Graph</th>"
                     "<th>Title</th><th>Axes</th><th>Trends</th>"
                     "<th>Statistics</th><th>Issues</th><th>Details</th>"
                     "</tr></thead><tbody>")

            for i, r in gt_rows:
                dwi_type, base_name = parse_dwi_info(r["file_path"])
                graph_title = r.get("graph_title", "") or ""

                # Axes summary
                x_lbl = r.get("x_axis_label", "")
                x_unit = r.get("x_axis_units", "")
                x_min = r.get("x_axis_range_min", "")
                x_max = r.get("x_axis_range_max", "")
                y_lbl = r.get("y_axis_label", "")
                y_unit = r.get("y_axis_units", "")
                y_min = r.get("y_axis_range_min", "")
                y_max = r.get("y_axis_range_max", "")
                c_lbl = r.get("color_axis_label", "")
                c_unit = r.get("color_axis_units", "")

                axis_lines = []
                if x_lbl:
                    xr = f" [{x_min}\u2013{x_max}]" if x_min and x_max else ""
                    xu = f" ({x_unit})" if x_unit else ""
                    axis_lines.append(f"X: {x_lbl}{xu}{xr}")
                if y_lbl:
                    yr = f" [{y_min}\u2013{y_max}]" if y_min and y_max else ""
                    yu = f" ({y_unit})" if y_unit else ""
                    axis_lines.append(f"Y: {y_lbl}{yu}{yr}")
                if c_lbl:
                    cu = f" ({c_unit})" if c_unit else ""
                    axis_lines.append(f"C: {c_lbl}{cu}")
                axes_cell = '<br>'.join(_esc(a) for a in axis_lines) if axis_lines else "\u2014"

                # Trends cell with descriptions
                trends_str = r.get("trends_json", "[]") or "[]"
                try:
                    trends_list = json.loads(trends_str)
                except Exception:
                    trends_list = []
                trends_cell = ""
                if isinstance(trends_list, list) and trends_list:
                    trend_parts = []
                    for t in trends_list:
                        if isinstance(t, dict):
                            direction = t.get("direction", "")
                            series = t.get("series", "")
                            desc = t.get("description", "")
                            label = f"{series}: {direction}" if series else direction
                            tag = _trend_tag(label)
                            if desc:
                                tag += f' <span class="axis-info">{_esc(desc[:80])}</span>'
                            trend_parts.append(tag)
                    trends_cell = "<br>".join(trend_parts)

                # Statistics cell: extract p-values and correlations
                all_text = (r.get("summary", "") + " " +
                            r.get("trends_json", "") + " " +
                            r.get("inflection_points_json", ""))
                pvals = extract_pvalues(all_text)
                corrs = extract_correlations(all_text)
                stats_parts = []
                for pval, ctx in pvals:
                    cls = _sig_class(pval)
                    cls_attr = f' class="{cls}"' if cls else ""
                    stats_parts.append(f'<span{cls_attr}>p={pval:.4f}</span>')
                for rval, ctx in corrs:
                    if abs(rval) >= 0.3:
                        stats_parts.append(f'r={rval:+.2f}')
                stats_cell = "<br>".join(stats_parts) if stats_parts else "\u2014"

                # Issues cell
                issues_str = r.get("issues_json", "[]") or "[]"
                try:
                    issues_list = json.loads(issues_str)
                except Exception:
                    issues_list = []
                if isinstance(issues_list, list) and issues_list:
                    issues_cell = "<ul>" + "".join(
                        f"<li>{_esc(iss)}</li>" for iss in issues_list
                    ) + "</ul>"
                else:
                    issues_cell = "\u2014"

                # Details: inflection points + summary in collapsible
                detail_parts = []
                # Inflection points inline
                ips_str = r.get("inflection_points_json", "[]") or "[]"
                try:
                    ips_list = json.loads(ips_str)
                except Exception:
                    ips_list = []
                if isinstance(ips_list, list) and ips_list:
                    ip_items = []
                    for ip in ips_list:
                        if isinstance(ip, dict):
                            x = ip.get("approximate_x", "?")
                            y = ip.get("approximate_y", "")
                            desc_ip = ip.get("description", "")
                            coord = f"x={_esc(str(x))}"
                            if y is not None and y != "":
                                coord += f", y={_esc(str(y))}"
                            ip_items.append(f"({coord}): {_esc(str(desc_ip))}")
                    if ip_items:
                        detail_parts.append(
                            "<strong>Inflection points:</strong><ul>" +
                            "".join(f"<li>{item}</li>" for item in ip_items) +
                            "</ul>"
                        )

                # Summary
                summary = r.get("summary", "") or ""
                if summary:
                    if len(summary) > 150:
                        detail_parts.append(
                            f"<details><summary>Full summary ({len(summary)} chars)</summary>"
                            f'<p class="full-summary">{_esc(summary)}</p></details>'
                        )
                    else:
                        detail_parts.append(f'<p class="full-summary">{_esc(summary)}</p>')

                details_cell = "".join(detail_parts) if detail_parts else "\u2014"

                title_display = _esc(graph_title) if graph_title else "\u2014"
                h.append(
                    f"<tr>"
                    f"<td>{i}</td>"
                    f"<td>{_dwi_badge(dwi_type)}</td>"
                    f"<td>{_esc(base_name)}</td>"
                    f'<td class="axis-info">{title_display}</td>'
                    f'<td class="axis-info">{axes_cell}</td>'
                    f"<td>{trends_cell}</td>"
                    f"<td>{stats_cell}</td>"
                    f"<td>{issues_cell}</td>"
                    f"<td>{details_cell}</td>"
                    f"</tr>"
                )
            h.append("</tbody></table>")
    return h



def _section_figure_gallery(folder) -> list[str]:
    """Build a Figure Gallery section with embedded images and captions.

    Embeds all pipeline-generated figures (PNG/JPG) from the output folder
    as base64 data URIs so the report is fully self-contained. Each figure
    receives an auto-numbered caption with a descriptive title derived
    from the filename.

    Parameters
    ----------
    folder : Path
        Path to the ``saved_files_*`` output folder.

    Returns
    -------
    list[str]
        HTML chunks for the figure gallery section.
    """
    import base64
    from pathlib import Path

    h: list[str] = []
    figures: list[tuple[str, str, Path]] = []  # (dwi_type, clean_name, path)

    # Collect all image files from DWI-type subfolders
    for dt in DWI_TYPES:
        dt_dir = Path(folder) / dt
        if not dt_dir.is_dir():
            continue
        for ext in ("*.png", "*.jpg", "*.jpeg"):
            for img_path in sorted(dt_dir.glob(ext)):
                # Skip empty files
                if img_path.stat().st_size == 0:
                    continue
                # Generate a clean name from the filename
                clean = img_path.stem.replace("_", " ").replace("-", " ").title()
                figures.append((dt, clean, img_path))

    # Also check root folder for cross-DWI comparison figures
    root_folder = Path(folder)
    for ext in ("*.png", "*.jpg", "*.jpeg"):
        for img_path in sorted(root_folder.glob(ext)):
            if img_path.stat().st_size < 1024:
                continue
            clean = img_path.stem.replace("_", " ").replace("-", " ").title()
            figures.append(("Root", clean, img_path))

    if not figures:
        return h

    h.append(_h2("Figure Gallery", "figure-gallery"))
    h.append(
        f'<p class="meta">All pipeline-generated figures are embedded below '
        f"with auto-numbered captions. {len(figures)} figure(s) found across "
        f"the output folder. Figures are grouped by DWI type.</p>"
    )

    current_dt = ""
    for dt, clean_name, img_path in figures:
        if dt != current_dt:
            current_dt = dt
            h.append(f"<h3>{_dwi_badge(dt) if dt != 'Root' else 'Cross-DWI'} Figures</h3>")

        # Read and encode the image as base64
        try:
            img_bytes = img_path.read_bytes()
            mime = "image/png" if img_path.suffix.lower() == ".png" else "image/jpeg"
            b64 = base64.b64encode(img_bytes).decode("ascii")
            caption = _figure_caption(
                f"{clean_name} ({dt})" if dt != "Root" else clean_name,
                f"Generated by the pancData3 pipeline."
            )
            h.append("<figure style=\"margin:1rem 0;page-break-inside:avoid\">")
            h.append(
                f'<img src="data:{mime};base64,{b64}" '
                f'alt="{_esc(clean_name)}" '
                f'style="max-width:100%;height:auto;border:1px solid var(--border);'
                f'border-radius:4px">'
            )
            h.append(caption)
            h.append("</figure>")
        except Exception:
            # Skip files that can't be read
            pass

    return h


