"""Report sections: data sections group (cohort, flow, completeness, MAT data)."""

from __future__ import annotations

from shared import (  # type: ignore
    DWI_TYPES,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _stat_card,
    _table_caption,
)
from report.sections._helpers import _scalar_gy  # type: ignore


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

    # ── Outcome Balance ──
    if log_data:
        h.append("<h3>Outcome Balance</h3>")
        any_outcome = False
        for dt in dwi_types_present:
            if dt not in log_data:
                continue
            bl = log_data[dt].get("baseline", {})
            baseline_exc = bl.get("baseline_exclusion", {})

            # Try to get LF/LC counts from survival data first, then compute from rate
            surv = log_data[dt].get("survival", {})
            n_lf = surv.get("n_lf")
            n_lc = surv.get("n_lc")
            lf_rate = baseline_exc.get("lf_rate_included") if baseline_exc else None

            if n_lf is None and lf_rate is not None:
                # Derive n_total from MAT longitudinal data or baseline_exclusion
                n_total_analyzed = None
                if mat_data and dt in mat_data and "longitudinal" in mat_data[dt]:
                    n_total_analyzed = mat_data[dt]["longitudinal"].get("num_patients")
                if n_total_analyzed is None and baseline_exc:
                    n_total_all = baseline_exc.get("n_total", 0)
                    n_excl = baseline_exc.get("n_excluded", 0)
                    n_total_analyzed = n_total_all - n_excl if n_total_all else None
                if n_total_analyzed:
                    n_lf = round(n_total_analyzed * lf_rate / 100)
                    n_lc = n_total_analyzed - n_lf

            if n_lf is not None and n_lc is not None:
                any_outcome = True
                h.append(f"<p>{_dwi_badge(dt)}</p>")
                h.append('<div class="stat-grid">')
                h.append(_stat_card("Local Failure", str(int(n_lf)), "N"))
                h.append(_stat_card("Local Control", str(int(n_lc)), "N"))
                h.append("</div>")
                # Imbalance warning
                total_oc = int(n_lf) + int(n_lc)
                if total_oc > 0:
                    lf_pct = 100.0 * int(n_lf) / total_oc
                    if lf_pct < 20 or lf_pct > 80:
                        h.append(
                            '<div class="warn-box">'
                            "\u26a0\ufe0f <strong>Imbalanced outcome distribution "
                            f"(LF: {lf_pct:.0f}%):</strong> Imbalanced outcome distribution "
                            "may reduce statistical power. Consider class-weighted analysis."
                            "</div>"
                        )

        if not any_outcome:
            h.append('<p class="meta">Outcome balance data (LF/LC counts) not available '
                     "in parsed log outputs.</p>")

    # ── Patient Attrition Across Timepoints ──
    if mat_data:
        for dt in dwi_types_present:
            if dt not in mat_data:  # type: ignore
                continue
            lon = mat_data[dt].get("longitudinal", {})  # type: ignore
            n_tp = lon.get("num_timepoints")
            if n_tp and n_tp > 1:
                # Check whether per-timepoint counts are available
                tp_counts = lon.get("patients_per_timepoint")
                h.append("<h3>Attrition</h3>")
                if tp_counts and isinstance(tp_counts, list):
                    h.append(f"<p>{_dwi_badge(dt)} \u2014 patients at each timepoint:</p>")
                    h.append("<table><thead><tr><th>Timepoint</th>"
                             "<th>Patients</th></tr></thead><tbody>")
                    for idx, cnt in enumerate(tp_counts, 1):
                        h.append(f"<tr><td>Timepoint {idx}</td>"
                                 f"<td><strong>{cnt}</strong></td></tr>")
                    h.append("</tbody></table>")
                else:
                    h.append('<p class="meta">Per-timepoint attrition data not available '
                             "in parsed outputs.</p>")
                break  # Only show once (first DWI type with multi-timepoint data)

    # ── Cross-DWI Patient Count Comparability Warning ──
    if mat_data and len(dwi_types_present) > 1:
        count_by_type: dict[str, int] = {}
        for dt in dwi_types_present:
            if dt in mat_data and "longitudinal" in mat_data[dt]:
                n = mat_data[dt]["longitudinal"].get("num_patients", 0)
                if n > 0:
                    count_by_type[dt] = n
        if len(count_by_type) >= 2:
            counts = list(count_by_type.values())
            if max(counts) - min(counts) > 2:
                parts = ", ".join(f"{dt}: {n}" for dt, n in count_by_type.items())
                h.append(
                    '<div class="warn-box">'
                    "\u26a0\ufe0f <strong>Patient counts differ across DWI types "
                    f"({parts}).</strong> Ensure identical cohorts for cross-DWI comparison."
                    "</div>"
                )

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

    # Parse-failure warning (change 10)
    has_parse_warnings = False
    if log_data is None:
        has_parse_warnings = True
    elif any(
        log_data.get(dt, {}).get("parse_warnings")
        for dt in dwi_types_present
        if dt in log_data
    ):
        has_parse_warnings = True
    if has_parse_warnings:
        h.append(
            '<div class="warn-box">'
            "\u26a0 <strong>Some flow data could not be parsed from logs.</strong> "
            "Values below may be incomplete."
            "</div>"
        )

    h.append(
        '<p class="meta">CONSORT-style summary of patient inclusion and exclusion '
        'at each analysis stage. Numbers may differ across DWI types due to '
        'type-specific data availability and model convergence.</p>'
    )

    h.append("<table>")
    h.append(_table_caption(
        "Patient Flow Summary",
        "Attrition from initial cohort through analysis stages."))
    # Change 9: clarify CR column label with tooltip-style sub-text
    h.append(
        "<thead><tr>"
        "<th>DWI Type</th>"
        "<th>Initial Cohort</th>"
        "<th>Baseline Excluded</th>"
        "<th>Outliers Removed</th>"
        "<th>Competing-Risk Excluded"
        '<br><span style="font-size:0.75em;font-weight:normal;color:var(--muted)">'
        "(patients who died without local failure, excluded from time-dependent "
        "GLME model)</span></th>"
        "<th>Analysed (GLME)</th>"
        "</tr></thead><tbody>"
    )
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

