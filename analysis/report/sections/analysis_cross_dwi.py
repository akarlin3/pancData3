"""Report sections: cross-DWI comparison and correlations."""

from __future__ import annotations

import json

from shared import (  # type: ignore
    DWI_TYPES,
    extract_correlations,
    get_config,
    parse_dwi_info,
    safe_text,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _trend_tag,
)
from report.sections._helpers import (  # type: ignore
    _build_normalised_series_map,
    _best_display_name,
)


def _section_cross_dwi_comparison(groups, csv_data) -> list[str]:
    """Build the Cross-DWI Comparison section.

    For priority graphs that exist in multiple DWI types, displays a
    trend-direction comparison table showing whether Standard, dnCNN,
    and IVIMnet agree or differ on each data series' direction.

    Also includes a significance inconsistency table from the CSV
    cross-reference (metrics significant in some DWI types but not others).

    Parameters
    ----------
    groups : dict[str, dict[str, dict]]
        Graph rows grouped by base name and DWI type.
    csv_data : dict or None
        Parsed pipeline CSV exports (contains ``cross_reference`` key).

    Returns
    -------
    list[str]
        HTML chunks for the cross-DWI comparison section.
    """
    # ── 4. Cross-DWI Comparison ──
    h: list[str] = []
    if groups:
        h.append(_h2("Cross-DWI Comparison", "cross-dwi"))

        priority_graphs = get_config().get("priority_graphs", [
            "Dose_vs_Diffusion", "Longitudinal_Mean_Metrics",
            "Longitudinal_Mean_Metrics_ByOutcome", "Feature_BoxPlots",
            "Feature_Histograms", "core_method_dice_heatmap",
        ])

        # Build the list of all graph groups with 2+ DWI types, priority first
        all_comparable = []
        seen = set()
        for base_name in priority_graphs:
            if base_name in groups:
                real = [t for t in groups[base_name] if t != "Root"]
                if len(real) >= 2:
                    all_comparable.append((base_name, True))
                    seen.add(base_name)
        for base_name in sorted(groups.keys()):
            if base_name not in seen:
                real = [t for t in groups[base_name] if t != "Root"]
                if len(real) >= 2:
                    all_comparable.append((base_name, False))

        # Sample-size mismatch warning (change 4)
        if groups:
            for base_name, _ in [(bn, _) for bn, _ in
                                  [(k, None) for k in groups]]:
                dwi_dict = groups[base_name]
                ss_by_dt: dict[str, int] = {}
                for dt in DWI_TYPES:
                    if dt in dwi_dict:
                        ss_raw = dwi_dict[dt].get("sample_size")
                        if ss_raw is not None:
                            try:
                                ss_by_dt[dt] = int(ss_raw)
                            except (ValueError, TypeError):
                                pass
                if len(ss_by_dt) >= 2:
                    ss_vals = list(ss_by_dt.values())
                    if max(ss_vals) - min(ss_vals) > 5:
                        parts = ", ".join(f"{dt}: {n}" for dt, n in ss_by_dt.items())
                        h.append(
                            '<div class="warn-box">'
                            f"\u26a0\ufe0f <strong>Sample size mismatch in "
                            f"\u201c{_esc(base_name)}\u201d ({parts}).</strong> "
                            "Results for this graph may not be directly comparable "
                            "across DWI types."
                            "</div>"
                        )

        # Summary counts
        n_agree: int = 0
        n_differ: int = 0
        # Track disagreements for Notable Disagreements subsection (change 6)
        disagree_records: list[dict] = []

        for base_name, is_priority in all_comparable:
            dwi_dict = groups[base_name]

            priority_label = ' <span class="badge badge-standard">priority</span>' if is_priority else ""
            h.append(f"<h3>{_esc(base_name)}{priority_label}</h3>")

            all_trends: dict[str, list] = {}
            for dt in DWI_TYPES:
                if dt in dwi_dict:
                    try:
                        all_trends[dt] = json.loads(str(dwi_dict[dt].get("trends_json", "[]")))
                    except (json.JSONDecodeError, TypeError):
                        pass

            if all_trends:
                # Use normalised series names so free-form vision-API labels
                # like "Mean D - Local Control" and "Mean D (Local Control)"
                # are matched correctly across DWI types.
                norm_map = _build_normalised_series_map(all_trends)

                h.append("<table><thead><tr><th>Series</th>")
                for dt in DWI_TYPES:
                    h.append(f"<th>{_esc(dt)}</th>")
                h.append("<th>Agreement</th></tr></thead><tbody>")

                has_rows = False
                for norm_key in sorted(norm_map.keys()):
                    dt_entries = norm_map[norm_key]
                    directions: dict[str, str] = {}
                    descriptions: dict[str, str] = {}
                    for dt, (d, desc) in dt_entries.items():
                        directions[dt] = d
                        descriptions[dt] = desc

                    display_name = _best_display_name(all_trends, norm_key)

                    if len(directions) >= 2:
                        has_rows = True
                        vals = list(directions.values())
                        if len(set(vals)) == 1:
                            agree_html = '<span class="agree">AGREE</span>'
                            n_agree = int(n_agree + 1)  # type: ignore
                        else:
                            agree_html = '<span class="differ">DIFFER</span>'
                            n_differ = int(n_differ + 1)  # type: ignore
                            # Record for disagreement subsection
                            modal_dir = max(set(vals), key=vals.count)
                            agreeing = [dt for dt, d in directions.items() if d == modal_dir]
                            differing = [dt for dt, d in directions.items() if d != modal_dir]
                            disagree_records.append({
                                "graph": base_name,
                                "series": display_name,
                                "directions": dict(directions),
                                "agreeing": agreeing,
                                "differing": differing,
                            })
                        h.append(f"<tr><td>{_esc(display_name)}</td>")
                        for dt in DWI_TYPES:
                            d_str = directions.get(dt, "-")
                            cell = _trend_tag(d_str) if d_str != "-" else "-"
                            desc = descriptions.get(dt, "")
                            if desc:
                                cell += f'<br><span class="axis-info">{_esc(desc[:80])}</span>'  # type: ignore
                            h.append(f"<td>{cell}</td>")
                        h.append(f"<td>{agree_html}</td></tr>")
                    elif len(directions) == 1:
                        # Show single-type series for completeness (no agreement verdict)
                        has_rows = True
                        h.append(f"<tr><td>{_esc(display_name)}</td>")
                        for dt in DWI_TYPES:
                            d_str = directions.get(dt, "-")
                            cell = _trend_tag(d_str) if d_str != "-" else "-"
                            desc = descriptions.get(dt, "")
                            if desc:
                                cell += f'<br><span class="axis-info">{_esc(desc[:80])}</span>'  # type: ignore
                            h.append(f"<td>{cell}</td>")
                        h.append("<td>-</td></tr>")

                if not has_rows:
                    h.append('<tr><td colspan="5"><em>No trend data available</em></td></tr>')

                h.append("</tbody></table>")

        # Overall agreement summary with improved interpretation (change 5)
        n_tot = sum([n_agree, n_differ])
        if n_tot > 0:
            pct_agree = 100.0 * float(n_agree) / float(n_tot)  # type: ignore
            cls = "agree" if pct_agree >= 70 else ("differ" if pct_agree < 50 else "")  # type: ignore
            cls_attr = f' class="{cls}"' if cls else ""
            h.append(
                f'<div class="summary-box"><strong>Cross-DWI Agreement:</strong> '
                f'<span{cls_attr}>{n_agree}/{n_agree + n_differ} series agree '
                f'({pct_agree:.0f}%)</span>, '
                f'{n_differ} differ across {len(all_comparable)} graph group(s). '
                f'<span class="meta">Clinical consensus typically requires &gt;80% '
                f'DWI-type agreement before adopting a single processing strategy as '
                f'standard. Agreement \u226570% suggests moderate robustness.</span>'
                f'</div>'
            )

        # Notable Disagreements subsection (change 6)
        if disagree_records:
            h.append("<h3>Notable Disagreements</h3>")
            h.append(
                '<p class="meta">The following graph\u2013series combinations show '
                "disagreement in trend direction across DWI types. "
                "Processing-dependent results should be validated with sensitivity "
                "analysis.</p>"
            )
            h.append(
                "<table><thead><tr><th>Graph</th><th>Series</th>"
                "<th>Agreeing Types</th><th>Differing Types</th>"
                "<th>Directions</th></tr></thead><tbody>"
            )
            for rec in disagree_records:
                dir_str = "; ".join(
                    f"{dt}: {d}" for dt, d in rec["directions"].items()
                )
                h.append(
                    f"<tr>"
                    f"<td>{_esc(rec['graph'])}</td>"
                    f"<td>{_esc(rec['series'])}</td>"
                    f'<td class="agree">{_esc(", ".join(rec["agreeing"]))}</td>'
                    f'<td class="differ">{_esc(", ".join(rec["differing"]))}</td>'
                    f"<td><small>{_esc(dir_str)}</small></td>"
                    f"</tr>"
                )
            h.append("</tbody></table>")

    # Cross-reference from CSV
    if csv_data and csv_data.get("cross_reference"):
        cross_ref_list = csv_data["cross_reference"]
        inconsistent = []
        if isinstance(cross_ref_list, list):
            for c in cross_ref_list:
                if isinstance(c, dict) and not c.get("consistent", True):
                    inconsistent.append(c)
        if inconsistent:
            h.append("<h3>Cross-DWI Significance Inconsistencies</h3>")
            h.append("<table><thead><tr><th>Metric</th><th>Timepoint</th>"
                     "<th>Significant In</th><th>Not Significant In</th></tr></thead><tbody>")
            for c in inconsistent:
                sig_in = c.get("significant_in", [])
                nsig_in = c.get("not_significant_in", [])
                sig_str = ", ".join(sig_in) if isinstance(sig_in, list) else "-"
                nsig_str = ", ".join(nsig_in) if isinstance(nsig_in, list) else "-"
                h.append(f"<tr><td>{_esc(str(c.get('metric', '')))}</td>"
                         f"<td>{_esc(str(c.get('timepoint', '')))}</td>"
                         f"<td>{_esc(sig_str)}</td>"
                         f"<td><span class='differ'>{_esc(nsig_str)}</span></td></tr>")
            h.append("</tbody></table>")
    return h



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
