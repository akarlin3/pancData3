"""Report sections: data completeness."""

from __future__ import annotations

from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _stat_card,
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
