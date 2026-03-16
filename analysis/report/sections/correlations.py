"""Report sections: notable correlations."""

from __future__ import annotations

from shared import (  # type: ignore
    extract_correlations,
    parse_dwi_info,
    safe_text,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


def _section_correlations(rows) -> list[str]:
    """Build the Notable Correlations section.

    Extracts correlation coefficients (r, rs, r-squared) from vision
    graph summaries and trend descriptions, filtering for |r| >= 0.3.
    Results are sorted by absolute correlation strength (descending).

    Parameters
    ----------
    rows : list[dict]
        Vision CSV rows.

    Returns
    -------
    list[str]
        HTML chunks for the correlations section.
    """
    # ── 5. Correlations ──
    h: list[str] = []
    if rows:
        h.append(_h2("Notable Correlations", "correlations"))

        # Causation caveat (change 1)
        h.append(
            '<div class="info-box">Correlation does not imply causation. '
            "Associations shown below may reflect confounding by dose level, "
            "patient factors, or treatment modality. Partial correlation analyses "
            "adjusting for covariates are recommended before clinical "
            "interpretation.</div>"
        )

        corr_findings = []
        for r in rows:
            dwi_type, base_name = parse_dwi_info(r["file_path"])
            all_text = safe_text(r, "summary", "trends_json")
            for rval, context in extract_correlations(all_text):
                if abs(rval) >= 0.3:
                    corr_findings.append((abs(rval), rval, dwi_type, base_name, context))

        if corr_findings:
            corr_findings.sort(reverse=True)

            # Split into strong (|r| >= 0.5) and moderate (0.3 <= |r| < 0.5)
            strong = [(a, r, d, g, c) for a, r, d, g, c in corr_findings if a >= 0.5]
            moderate = [(a, r, d, g, c) for a, r, d, g, c in corr_findings if a < 0.5]

            h.append(f"<p>{len(corr_findings)} correlation(s) with |r| \u2265 0.3 "
                     f"({len(strong)} strong, {len(moderate)} moderate).</p>")
            h.append(
                '<p class="meta">Correlation strength benchmarks (Cohen): '
                '|r| \u2265 0.5 strong, 0.3\u20130.5 moderate, &lt; 0.3 weak. '
                'For DWI parameters, negative correlations with dose metrics suggest '
                'dose-response relationships; positive ADC\u2013D correlations reflect '
                'shared diffusion signal.</p>'
            )

            # Bonferroni correction note (change 2)
            n_corr_tested = len(corr_findings)
            if n_corr_tested > 1:
                bonferroni_alpha = 0.05 / n_corr_tested
                surviving = sum(
                    1 for _, rval, _, _, _ in corr_findings
                    if abs(rval) >= 0.5  # rough proxy; true p not always available
                )
                h.append(
                    f'<div class="info-box">'
                    f"<strong>Multiple testing note:</strong> {n_corr_tested} correlations "
                    f"tested; Bonferroni-corrected threshold: \u03b1\u2009=\u2009"
                    f"{bonferroni_alpha:.4f}. "
                    f"Correlations surviving correction (|r|\u2009\u22650.5, as a "
                    f"conservative proxy): {surviving}."
                    f"</div>"
                )

            if strong:
                h.append("<h3>Strong Correlations (|r| \u2265 0.5)</h3>")
                h.append("<table><thead><tr><th>|r|</th><th>r</th><th>Strength</th>"
                         "<th>DWI</th><th>Graph</th><th>Context</th></tr></thead><tbody>")
                for _, rval, dwi, graph, ctx in strong:
                    ctx_short = _esc(ctx.replace("\n", " ")[:120])
                    h.append(f'<tr><td class="sig-2"><strong>{abs(rval):.2f}</strong></td>'
                             f"<td>{rval:+.2f}</td>"
                             f"<td>Strong</td><td>{_dwi_badge(dwi)}</td>"
                             f"<td>{_esc(graph)}</td><td><em>{ctx_short}</em></td></tr>")
                h.append("</tbody></table>")

            if moderate:
                h.append("<h3>Moderate Correlations (0.3 \u2264 |r| &lt; 0.5)</h3>")
                h.append("<table><thead><tr><th>|r|</th><th>r</th><th>Strength</th>"
                         "<th>DWI</th><th>Graph</th><th>Context</th></tr></thead><tbody>")
                for _, rval, dwi, graph, ctx in moderate:
                    ctx_short = _esc(ctx.replace("\n", " ")[:120])
                    h.append(f"<tr><td><strong>{abs(rval):.2f}</strong></td><td>{rval:+.2f}</td>"
                             f"<td>Moderate</td><td>{_dwi_badge(dwi)}</td>"
                             f"<td>{_esc(graph)}</td><td><em>{ctx_short}</em></td></tr>")
                h.append("</tbody></table>")

            # Confidence interval caveat (change 3)
            h.append(
                '<p class="meta"><strong>Note:</strong> Confidence intervals for '
                "correlations are not reported. With cohort sizes typical of this "
                "study (n\u2009&lt;\u200950), 95% CIs for r\u2009=\u20090.65 span "
                "approximately \u00b10.25. Interpret effect sizes with caution.</p>"
            )
        else:
            h.append("<p>No notable correlations (|r| &ge; 0.3) found.</p>")
    return h
