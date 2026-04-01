"""Report section: risk model vs dose coverage concordance."""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


def _section_risk_dose_concordance(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build HTML section showing risk-dose concordance analysis.

    Parameters
    ----------
    mat_data : dict
        Mapping of DWI type to parsed MAT metrics dict.
    dwi_types : list[str]
        List of DWI types present in the data.

    Returns
    -------
    list[str]
        HTML chunks for the risk-dose concordance section.
    """
    h: list[str] = []

    has_data = False
    for dwi in DWI_TYPES:
        if dwi in (mat_data or {}) and "risk_dose_concordance" in mat_data[dwi]:
            rc = mat_data[dwi]["risk_dose_concordance"]
            if rc and rc.get("method_results"):
                has_data = True
                break

    if not has_data:
        h.append(_h2("Risk Model vs Dose Coverage Concordance", "risk-dose-concordance"))
        h.append('<p class="meta">Risk-dose concordance data not available. '
                 "Enable <code>run_risk_dose_concordance</code> in config to generate.</p>")
        return h

    h.append(_h2("Risk Model vs Dose Coverage Concordance", "risk-dose-concordance"))
    h.append(
        "<p>Comparison of patient-level risk classifications from the elastic net "
        "LOOCV model against D95-based stratification. Cohen's kappa measures "
        "agreement: &lt;0.2 = poor, 0.2&ndash;0.4 = fair, 0.4&ndash;0.6 = moderate, "
        "0.6&ndash;0.8 = substantial, &gt;0.8 = excellent. "
        "Complementary models (low kappa) identify different at-risk patients; "
        "combined use may improve prediction.</p>"
    )

    for dwi in DWI_TYPES:
        if dwi not in (mat_data or {}):
            continue
        rc = mat_data[dwi].get("risk_dose_concordance", {})
        if not rc or not rc.get("method_results"):
            continue

        method_results = rc["method_results"]
        summary = rc.get("summary", "")

        h.append(f"<h3>{_dwi_badge(dwi)} — Concordance Analysis</h3>")
        h.append("<table><thead><tr>")
        h.append("<th>Core Method</th><th>Dose Metric</th>"
                 "<th>Cohen's &kappa;</th><th>Concordance</th>"
                 "<th>Complementary Patients</th><th>Combined AUC</th>")
        h.append("</tr></thead><tbody>")

        for mr in method_results:
            name = mr.get("method_name", "")
            dose_metric = mr.get("best_dose_metric", "")
            kappa = mr.get("cohen_kappa")
            conc_pct = mr.get("concordance_pct")
            n_comp = mr.get("n_complementary")
            combined_auc = mr.get("combined_auc")

            h.append(f"<tr><td><strong>{_esc(name)}</strong></td>")
            h.append(f"<td><code>{_esc(dose_metric)}</code></td>")

            if kappa is not None:
                kappa_label = _kappa_label(kappa)
                h.append(f"<td>{kappa:.3f} ({kappa_label})</td>")
            else:
                h.append("<td>&mdash;</td>")

            if conc_pct is not None:
                h.append(f"<td>{conc_pct:.1f}%</td>")
            else:
                h.append("<td>&mdash;</td>")

            if n_comp is not None:
                h.append(f"<td>{n_comp}</td>")
            else:
                h.append("<td>&mdash;</td>")

            if combined_auc is not None:
                h.append(f"<td>{combined_auc:.3f}</td>")
            else:
                h.append("<td>&mdash;</td>")

            h.append("</tr>")

        h.append("</tbody></table>")

        if summary:
            h.append(f'<p class="meta"><strong>Summary:</strong> {_esc(summary)}</p>')

    return h


def _kappa_label(kappa: float) -> str:
    """Interpret Cohen's kappa value."""
    if kappa >= 0.8:
        return "excellent"
    if kappa >= 0.6:
        return "substantial"
    if kappa >= 0.4:
        return "moderate"
    if kappa >= 0.2:
        return "fair"
    return "poor"
