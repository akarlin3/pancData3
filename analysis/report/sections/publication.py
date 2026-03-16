"""Report sections: publication support (reporting checklist, journal guide)."""

from __future__ import annotations

from shared import (  # type: ignore
    DWI_TYPES,
)
from report.report_formatters import (  # type: ignore
    _cite,
    _esc,
    _h2,
    _stat_card,
    _table_caption,
)


def _section_reporting_checklist(
    log_data, dwi_types_present, mat_data, csv_data, rows
) -> list[str]:
    """Build a STROBE/REMARK reporting guideline compliance checklist.

    Each checklist item is dynamically assessed against the available
    data to indicate whether the report addresses it (Done), partially
    addresses it (Partial), or leaves it for the researcher (N/A).

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    mat_data : dict
        Parsed MAT file metrics.
    csv_data : dict or None
        Parsed pipeline CSV exports.
    rows : list[dict]
        Vision CSV rows.

    Returns
    -------
    list[str]
        HTML chunks for the reporting checklist section.
    """
    h: list[str] = []
    h.append(_h2("Reporting Guideline Checklist", "reporting-checklist"))
    h.append(
        f'<p class="meta">Compliance with STROBE{_cite("strobe")} and '
        f'REMARK{_cite("remark")} reporting guidelines, automatically '
        f"assessed where possible. Items marked <strong>Partial</strong> "
        f"require researcher input to complete.</p>"
    )

    # Dynamically assess each checklist item
    has_cohort = bool(mat_data) and any(
        "longitudinal" in mat_data.get(dt, {}) for dt in DWI_TYPES
    )
    has_baseline = bool(log_data) and any(
        log_data.get(dt, {}).get("baseline") for dt in dwi_types_present
    )
    has_survival = bool(log_data) and any(
        log_data.get(dt, {}).get("survival", {}).get("hazard_ratios")
        for dt in dwi_types_present
    )
    has_predictive = bool(log_data) and any(
        log_data.get(dt, {}).get("stats_predictive", {}).get("roc_analyses")
        for dt in dwi_types_present
    )
    has_glme = bool(log_data) and any(
        log_data.get(dt, {}).get("stats_comparisons", {}).get("glme_details")
        for dt in dwi_types_present
    )
    has_fdr = bool(csv_data) and bool(csv_data.get("fdr_global"))
    has_graphs = bool(rows)
    has_dosimetry = bool(mat_data) and any(
        "dosimetry" in mat_data.get(dt, {}) for dt in DWI_TYPES
    )

    # Each item: (Section, Item, Status, Note)
    # Status: "done", "partial", "na"
    items: list[tuple[str, str, str, str]] = [
        ("STROBE 1", "Title: Indicate study design",
         "partial", "Template provided; researcher should confirm study type"),
        ("STROBE 2", "Abstract: Structured summary",
         "done" if has_cohort else "partial",
         "Executive summary auto-generated with structured subsections"),
        ("STROBE 3", "Background: Scientific rationale",
         "partial", "Hypothesis section generated; literature review needed"),
        ("STROBE 4", "Objectives: Specific aims",
         "done", "Objective stated in abstract"),
        ("STROBE 5", "Study design: Key elements",
         "done", "Methods section describes retrospective cohort design"),
        ("STROBE 6", "Setting: Dates and locations",
         "partial", "Institution placeholder provided; dates from pipeline"),
        ("STROBE 7", "Participants: Eligibility criteria",
         "partial" if has_cohort else "na",
         "Patient flow section shows inclusion/exclusion" if has_cohort
         else "No cohort data available"),
        ("STROBE 8", "Variables: Outcomes and exposures",
         "done" if has_glme or has_survival else "partial",
         "DWI biomarkers defined; outcome (LF/LC) described"),
        ("STROBE 12", "Statistical methods",
         "done", "Full methods section auto-generated with citations"),
        ("STROBE 13", "Participants: Numbers at each stage",
         "done" if has_baseline else "partial",
         "Patient flow table with attrition stages" if has_baseline
         else "Baseline data needed for flow diagram"),
        ("STROBE 14", "Descriptive data: Characteristics",
         "partial", "Outlier/quality summary provided; demographics needed"),
        ("STROBE 16", "Main results: Unadjusted and adjusted",
         "done" if has_glme else "partial",
         "GLME interaction tests with FDR adjustment reported" if has_glme
         else "Statistical results needed"),
        ("STROBE 17", "Other analyses: Subgroup, sensitivity",
         "done" if has_survival or has_predictive else "partial",
         "Sensitivity analysis section auto-generated"),
        ("STROBE 18", "Key results: Summary",
         "done", "Conclusions section with numbered findings"),
        ("STROBE 19", "Limitations",
         "done", "Contextual limitations auto-generated from data"),
        ("STROBE 20", "Interpretation",
         "done" if has_cohort else "partial",
         "Clinical significance and future directions provided"),
        ("STROBE 22", "Funding",
         "partial", "Data availability statement provided; funding needed"),
        ("REMARK 1", "Study hypothesis",
         "done", "Data-driven hypothesis section auto-generated"),
        ("REMARK 2", "Patient characteristics",
         "partial" if has_cohort else "na",
         "Cohort overview with sample sizes" if has_cohort
         else "Patient data not available"),
        ("REMARK 5", "Assay methods",
         "done", "DWI acquisition and IVIM fitting described in Methods"),
        ("REMARK 7", "Statistical analysis plan",
         "done", "Full statistical methods with multiple comparison correction"),
        ("REMARK 8", "Effect sizes with CIs",
         "done" if has_survival else "partial",
         "Hazard ratios with 95% CI and forest plots" if has_survival
         else "Effect size data requires survival analysis"),
        ("REMARK 9", "Multivariable analysis",
         "done" if has_predictive else "na",
         "Elastic net with LOOCV reported" if has_predictive
         else "No predictive model data"),
        ("REMARK 10", "Missing data handling",
         "done" if has_baseline else "partial",
         "Exclusions and imputation documented" if has_baseline
         else "Missing data handling not yet documented"),
    ]

    n_done = sum(1 for _, _, s, _ in items if s == "done")
    n_partial = sum(1 for _, _, s, _ in items if s == "partial")
    n_na = sum(1 for _, _, s, _ in items if s == "na")

    h.append('<div class="stat-grid">')
    h.append(_stat_card("Addressed", str(n_done),
                        f"of {len(items)} checklist items"))
    h.append(_stat_card("Partial", str(n_partial), "require researcher input"))
    if n_na > 0:
        h.append(_stat_card("N/A", str(n_na), "not applicable"))
    h.append("</div>")

    h.append('<table class="checklist-table"><thead><tr>'
             "<th>Item</th><th>Requirement</th>"
             "<th>Status</th><th>Notes</th>"
             "</tr></thead><tbody>")
    for section, requirement, status, note in items:
        if status == "done":
            s_html = '<span class="checklist-done">Addressed</span>'
        elif status == "partial":
            s_html = '<span class="checklist-partial">Partial</span>'
        else:
            s_html = '<span class="checklist-na">N/A</span>'
        h.append(
            f"<tr><td>{_esc(section)}</td>"
            f"<td>{_esc(requirement)}</td>"
            f"<td>{s_html}</td>"
            f"<td>{_esc(note)}</td></tr>"
        )
    h.append("</tbody></table>")

    return h



def _section_journal_guide(log_data, dwi_types_present, mat_data) -> list[str]:
    """Build a Journal Submission Guidance section.

    Provides target journal recommendations, word count estimates, and
    formatting checklists to accelerate manuscript preparation.

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
        HTML chunks for the journal submission guidance section.
    """
    h: list[str] = []
    h.append(_h2("Journal Submission Guidance", "journal-guide"))
    h.append(
        '<p class="meta">Recommendations for target journals and formatting '
        "requirements based on the study design and findings.</p>"
    )

    # Assess study characteristics for journal matching
    has_survival = bool(log_data) and any(
        log_data.get(dt, {}).get("survival", {}).get("hazard_ratios")
        for dt in dwi_types_present
    )
    has_predictive = bool(log_data) and any(
        log_data.get(dt, {}).get("stats_predictive", {}).get("roc_analyses")
        for dt in dwi_types_present
    )
    n_patients = 0
    if mat_data:
        for dt in DWI_TYPES:
            if dt in mat_data and "longitudinal" in mat_data[dt]:  # type: ignore
                n = mat_data[dt]["longitudinal"].get("num_patients", 0)  # type: ignore
                if n > n_patients:
                    n_patients = n

    # Target journal recommendations
    h.append("<h3>Suggested Target Journals</h3>")
    journals = [
        ("Radiotherapy and Oncology", "3000\u20134000", "Original Article",
         "Pancreatic RT, DWI biomarkers, treatment response"),
        ("International Journal of Radiation Oncology, Biology, Physics (Red Journal)",
         "3500\u20134500", "Clinical Investigation",
         "RT dose-response, imaging biomarkers"),
        ("Physics in Medicine & Biology", "4000\u20136000", "Paper",
         "DWI methodology, IVIM modelling, denoising comparison"),
        ("Medical Physics", "4000\u20136000", "Research Article",
         "Quantitative imaging, model comparison, dosimetric analysis"),
    ]
    if has_survival:
        journals.append((
            "Acta Oncologica", "3000\u20134000", "Original Article",
            "Survival analysis, prognostic biomarkers"
        ))
    if has_predictive:
        journals.append((
            "European Radiology", "3500\u20134500", "Original Article",
            "Predictive modelling, radiomics, machine learning"
        ))

    h.append("<table>")
    h.append(_table_caption("Suggested Target Journals",
                            "Word limits and scope alignment."))
    h.append("<thead><tr><th>Journal</th><th>Word Limit</th>"
             "<th>Article Type</th><th>Scope Match</th></tr></thead><tbody>")
    for name, words, atype, scope in journals:
        h.append(
            f"<tr><td><strong>{_esc(name)}</strong></td>"
            f"<td>{_esc(words)}</td>"
            f"<td>{_esc(atype)}</td>"
            f"<td>{_esc(scope)}</td></tr>"
        )
    h.append("</tbody></table>")

    # Manuscript structure checklist
    h.append("<h3>Manuscript Preparation Checklist</h3>")
    checklist = [
        ("Title page", "Title, all authors with affiliations, corresponding author, "
         "running head, word count, keywords"),
        ("Abstract", "Structured (Purpose/Methods/Results/Conclusions), "
         "250 words max for most journals"),
        ("Introduction", "Scientific rationale, knowledge gap, study objectives "
         "(typically 400\u2013600 words)"),
        ("Methods", "Auto-generated in this report; review for completeness"),
        ("Results", "Draft available in this report; add demographics table"),
        ("Discussion", "Interpret findings, compare with literature, limitations, "
         "future directions (typically 1000\u20131500 words)"),
        ("References", "BibTeX export available in this report; verify journal "
         "formatting style"),
        ("Figures", "Gallery with captions available; ensure resolution "
         "\u2265300 DPI for print"),
        ("Tables", "Auto-numbered in this report; check journal table formatting"),
        ("Supplementary", "Consider moving detailed tables to supplement"),
        ("Cover letter", "Highlight novelty (multi-strategy DWI, longitudinal "
         "design, cross-DWI comparison)"),
        ("ICMJE forms", "All authors must complete ICMJE disclosure forms"),
        ("IRB statement", "Include protocol number in Methods"),
    ]

    h.append('<table><thead><tr><th>Item</th><th>Notes</th>'
             '</tr></thead><tbody>')
    for item, notes in checklist:
        h.append(f"<tr><td><strong>{_esc(item)}</strong></td>"
                 f"<td>{_esc(notes)}</td></tr>")
    h.append("</tbody></table>")

    # Suggested keywords
    h.append("<h3>Suggested Keywords</h3>")
    keywords = [
        "diffusion-weighted imaging", "IVIM", "pancreatic cancer",
        "radiotherapy", "treatment response", "biomarker",
    ]
    if has_predictive:
        keywords.append("elastic net")
        keywords.append("predictive modelling")
    if has_survival:
        keywords.append("survival analysis")
        keywords.append("competing risks")
    if len(dwi_types_present) >= 2:
        keywords.append("deep learning denoising")
    h.append(
        '<div class="manuscript-sentence" data-copy="'
        + _esc("; ".join(keywords)) + '">'
        + _esc("; ".join(keywords))
        + '<button class="copy-btn" onclick="copyText(this)">Copy</button>'
        + "</div>"
    )

    return h
