"""Report sections: data completeness."""

from __future__ import annotations

from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _stat_card,
    _table_caption,
)


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


def build_registration_quality_section(saved_files: dict) -> str:
    """Build the Registration Quality section.

    Renders registration quality metrics (Jacobian determinant statistics,
    NCC, mutual information) and flags patients with poor registration.

    Parameters
    ----------
    saved_files : dict
        Parsed pipeline output (keyed by DWI type at top level, with
        ``registration_quality`` nested within).

    Returns
    -------
    str
        HTML string for the registration quality section.  Empty string
        if no registration data is present.
    """
    if not saved_files:
        return ""

    h: list[str] = []
    found = False

    for dwi_type, dwi_data in saved_files.items():
        if not isinstance(dwi_data, dict):
            continue
        reg = dwi_data.get("registration_quality")
        if not reg or not isinstance(reg, (dict, list)):
            continue

        # Accept either a list of per-patient dicts or a dict with a
        # "patients" key
        if isinstance(reg, dict):
            patients = reg.get("patients", [])
            summary = reg
        else:
            patients = reg
            summary = {}

        if not patients and not summary:
            continue

        if not found:
            h.append(_h2("Registration Quality", "registration-quality"))
            h.append(
                '<p class="meta">Registration quality metrics assess the '
                "accuracy of deformable image registration between DWI and "
                "dose maps. Poor registration (negative Jacobian determinant "
                "voxels or low NCC) may compromise dosimetry correlations.</p>"
            )
            found = True

        h.append(f"<h3>{_dwi_badge(dwi_type)}</h3>")

        # Summary cards from aggregate stats
        cards: list[str] = []
        mean_ncc = summary.get("mean_ncc")
        mean_nmi = summary.get("mean_nmi")
        n_poor = 0

        if isinstance(patients, list):
            for pt in patients:
                if not isinstance(pt, dict):
                    continue
                pt_ncc = pt.get("ncc")
                jac_neg = pt.get("jacobian_negative_voxels", 0)
                if ((isinstance(pt_ncc, (int, float)) and pt_ncc < 0.7)
                        or (isinstance(jac_neg, (int, float)) and jac_neg > 0)):
                    n_poor += 1

        if isinstance(mean_ncc, (int, float)):
            ncc_cls = "agree" if mean_ncc >= 0.7 else "differ"
            cards.append(_stat_card("Mean NCC", f"{mean_ncc:.3f}",
                                    "good" if ncc_cls == "agree" else "poor"))
        if isinstance(mean_nmi, (int, float)):
            cards.append(_stat_card("Mean NMI", f"{mean_nmi:.3f}"))
        if isinstance(patients, list) and patients:
            cards.append(_stat_card(
                "Poor Registration",
                f"{n_poor}/{len(patients)}",
                "NCC < 0.7 or Jacobian < 0",
            ))
        if cards:
            h.append('<div class="stat-grid">')
            h.extend(cards)
            h.append("</div>")

        # Per-patient table
        if isinstance(patients, list) and patients:
            h.append(_table_caption(
                "Per-Patient Registration Quality",
                f"Registration metrics for {_esc(dwi_type)}.",
            ))
            h.append(
                "<table><thead><tr>"
                "<th>Patient</th>"
                "<th>NCC</th>"
                "<th>NMI</th>"
                "<th>Jacobian &lt;0 Voxels</th>"
                "<th>Jacobian Mean</th>"
                "<th>Status</th>"
                "</tr></thead><tbody>"
            )
            for pt in patients:
                if not isinstance(pt, dict):
                    continue
                pid = pt.get("patient_id", "?")
                pt_ncc = pt.get("ncc")
                pt_nmi = pt.get("nmi")
                jac_neg = pt.get("jacobian_negative_voxels")
                jac_mean = pt.get("jacobian_mean")

                ncc_str = f"{pt_ncc:.3f}" if isinstance(pt_ncc, (int, float)) else "N/A"
                nmi_str = f"{pt_nmi:.3f}" if isinstance(pt_nmi, (int, float)) else "N/A"
                jac_neg_str = str(jac_neg) if isinstance(jac_neg, (int, float)) else "N/A"
                jac_mean_str = f"{jac_mean:.3f}" if isinstance(jac_mean, (int, float)) else "N/A"

                is_poor = (
                    (isinstance(pt_ncc, (int, float)) and pt_ncc < 0.7)
                    or (isinstance(jac_neg, (int, float)) and jac_neg > 0)
                )
                if is_poor:
                    status = '<span class="differ">Poor</span>'
                else:
                    status = '<span class="agree">Good</span>'

                h.append(
                    f"<tr><td><code>{_esc(str(pid))}</code></td>"
                    f"<td>{ncc_str}</td>"
                    f"<td>{nmi_str}</td>"
                    f"<td>{jac_neg_str}</td>"
                    f"<td>{jac_mean_str}</td>"
                    f"<td>{status}</td></tr>"
                )
            h.append("</tbody></table>")

        # Warning for poor registrations
        if n_poor > 0:
            h.append(
                f'<div class="warn-box"><strong>{n_poor} patient(s) with poor '
                f"registration quality.</strong> Dosimetry correlations for "
                f"these patients should be interpreted with caution. Consider "
                f"excluding them from dose-response analyses or performing "
                f"sensitivity analysis with and without these patients.</div>"
            )

    if not found:
        return ""
    return "\n".join(h)
