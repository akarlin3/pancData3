"""Report section: dose-response ROC analysis."""

from __future__ import annotations

from shared import DWI_TYPES  # type: ignore
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


def _section_dose_response_roc(mat_data: dict, dwi_types: list[str]) -> list[str]:
    """Build HTML section showing dose-response ROC results.

    Parameters
    ----------
    mat_data : dict
        Mapping of DWI type to parsed MAT metrics dict.
    dwi_types : list[str]
        List of DWI types present in the data.

    Returns
    -------
    list[str]
        HTML chunks for the dose-response ROC section.
    """
    h: list[str] = []

    has_data = False
    for dwi in DWI_TYPES:
        if dwi in (mat_data or {}) and "dose_response_roc" in mat_data[dwi]:
            rr = mat_data[dwi]["dose_response_roc"]
            if rr and rr.get("method_results"):
                has_data = True
                break

    if not has_data:
        h.append(_h2("Dose-Response ROC Analysis", "dose-response-roc"))
        h.append('<p class="meta">Dose-response ROC data not available. '
                 "Enable <code>run_dose_response_roc</code> in config to generate.</p>")
        return h

    h.append(_h2("Dose-Response ROC Analysis", "dose-response-roc"))
    h.append(
        "<p>ROC analysis on sub-volume dose metrics (D95, V50) to find the "
        "optimal dose cutoff separating local control from local failure. "
        "The Youden index identifies the threshold that maximises sensitivity + "
        "specificity. AUC 95% CI estimated via 1000 bootstrap resamples.</p>"
    )

    for dwi in DWI_TYPES:
        if dwi not in (mat_data or {}):
            continue
        rr = mat_data[dwi].get("dose_response_roc", {})
        if not rr or not rr.get("method_results"):
            continue

        method_results = rr["method_results"]

        h.append(f"<h3>{_dwi_badge(dwi)} — ROC Summary</h3>")
        h.append("<table><thead><tr>")
        h.append("<th>Core Method</th><th>Best Metric</th><th>AUC</th>"
                 "<th>AUC 95% CI</th><th>Optimal Threshold</th>"
                 "<th>Sensitivity</th><th>Specificity</th>")
        h.append("</tr></thead><tbody>")

        for mr in method_results:
            name = mr.get("method_name", "")
            best_metric = mr.get("best_metric", "")
            best_auc = mr.get("best_auc")

            # Find the metrics entry for the best metric
            auc_ci_str = "&mdash;"
            thresh_str = "&mdash;"
            sens_str = "&mdash;"
            spec_str = "&mdash;"

            metrics_list = mr.get("metrics", [])
            for met in (metrics_list if isinstance(metrics_list, list) else [metrics_list]):
                if met and met.get("metric_name") == best_metric:
                    ci = met.get("auc_ci", [])
                    if isinstance(ci, list) and len(ci) >= 2 and ci[0] is not None:
                        auc_ci_str = f"[{ci[0]:.3f}, {ci[1]:.3f}]"
                    thresh = met.get("optimal_threshold")
                    if thresh is not None:
                        thresh_str = f"{thresh:.2f}"
                    sens = met.get("sensitivity")
                    if sens is not None:
                        sens_str = f"{sens * 100:.1f}%"
                    spec = met.get("specificity")
                    if spec is not None:
                        spec_str = f"{spec * 100:.1f}%"
                    break

            auc_cls = ""
            if best_auc is not None and best_auc >= 0.7:
                auc_cls = ' style="background:#d4edda"'

            h.append(f"<tr><td><strong>{_esc(name)}</strong></td>")
            h.append(f"<td><code>{_esc(best_metric)}</code></td>")
            if best_auc is not None:
                h.append(f"<td{auc_cls}>{best_auc:.3f}</td>")
            else:
                h.append("<td>&mdash;</td>")
            h.append(f"<td>{auc_ci_str}</td>")
            h.append(f"<td>{thresh_str}</td>")
            h.append(f"<td>{sens_str}</td>")
            h.append(f"<td>{spec_str}</td>")
            h.append("</tr>")

        h.append("</tbody></table>")

        # Clinical guidance
        for mr in method_results:
            best_auc = mr.get("best_auc")
            if best_auc is not None and best_auc >= 0.7:
                metrics_list = mr.get("metrics", [])
                for met in (metrics_list if isinstance(metrics_list, list) else [metrics_list]):
                    if met and met.get("metric_name") == mr.get("best_metric"):
                        thresh = met.get("optimal_threshold")
                        sens = met.get("sensitivity")
                        spec = met.get("specificity")
                        if thresh is not None and sens is not None and spec is not None:
                            h.append(
                                f'<p class="meta"><strong>Clinical guidance:</strong> '
                                f'Based on ROC analysis, a D95 threshold of {thresh:.1f} Gy '
                                f'to the {_esc(mr["method_name"])}-defined resistant sub-volume '
                                f'achieves {sens * 100:.0f}% sensitivity and {spec * 100:.0f}% '
                                f'specificity for predicting local control.</p>'
                            )
                        break

    return h
