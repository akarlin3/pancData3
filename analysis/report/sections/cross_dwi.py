"""Report sections: cross-DWI comparison and feature overlap."""

from __future__ import annotations

import json

from shared import (  # type: ignore
    DWI_TYPES,
    get_config,
    parse_dwi_info,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _stat_card,
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
                    except Exception:
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



def _section_feature_overlap(log_data, dwi_types_present) -> list[str]:
    """Build the Cross-DWI Feature Overlap Analysis section.

    Compares elastic-net selected features across DWI types at each
    timepoint, highlighting features that are consistently selected
    (high confidence) vs type-specific (potentially noise).

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.

    Returns
    -------
    list[str]
        HTML chunks for the feature overlap section.
    """
    h: list[str] = []
    if not log_data or len(dwi_types_present) < 2:
        return h

    # Gather all feature selections keyed by timepoint and DWI type.
    tp_features: dict[str, dict[str, list[str]]] = {}
    for dt in dwi_types_present:
        if dt not in log_data:
            continue
        fs_list = log_data[dt].get("stats_predictive", {}).get("feature_selections", [])
        for fs in fs_list:
            tp = fs.get("timepoint", "?")
            tp_features.setdefault(tp, {})[dt] = fs.get("features", [])

    # Only show if we have at least one timepoint with 2+ DWI types.
    multi_tp = {tp: dts for tp, dts in tp_features.items() if len(dts) >= 2}
    if not multi_tp:
        return h

    h.append(_h2("Cross-DWI Feature Overlap", "feature-overlap"))
    h.append(
        '<p class="meta">Features consistently selected by elastic net across '
        'multiple DWI types at the same timepoint are more likely to reflect true '
        'biological signal rather than processing-specific artefacts.</p>'
    )

    total_shared: int = 0
    total_unique: int = 0

    for tp in sorted(multi_tp.keys()):
        dts = multi_tp[tp]
        all_features: set[str] = set()
        for feats in dts.values():
            all_features.update(feats)

        if not all_features:
            continue

        h.append(f"<h3>Timepoint: <code>{_esc(tp)}</code></h3>")

        # Classify each feature by how many DWI types selected it.
        feature_counts: dict[str, list[str]] = {}
        for feat in sorted(all_features):
            selected_in = [dt for dt, feats in dts.items() if feat in feats]
            feature_counts[feat] = selected_in

        shared = {f: dts_list for f, dts_list in feature_counts.items()
                  if len(dts_list) >= 2}
        unique = {f: dts_list for f, dts_list in feature_counts.items()
                  if len(dts_list) == 1}

        total_shared = int(total_shared + len(shared))  # type: ignore
        total_unique = int(total_unique + len(unique))  # type: ignore

        # Stat cards
        h.append('<div class="stat-grid">')
        h.append(_stat_card("Shared Features", str(len(shared)),
                            f"selected by \u22652 DWI types"))
        h.append(_stat_card("Type-Specific", str(len(unique)),
                            "selected by 1 DWI type only"))
        h.append(_stat_card("Total Unique", str(len(all_features)),
                            "across all DWI types"))
        h.append("</div>")

        # Table
        h.append("<table><thead><tr><th>Feature</th>")
        for dt in dwi_types_present:
            if dt in dts:
                h.append(f"<th>{_esc(dt)}</th>")
        h.append("<th>Consensus</th></tr></thead><tbody>")
        # Shared features first (sorted by count desc), then unique
        for feat in sorted(shared.keys()):
            dts_list = shared[feat]
            h.append(f"<tr><td><code>{_esc(feat)}</code></td>")
            for dt in dwi_types_present:
                if dt in dts:
                    if dt in dts_list:
                        h.append('<td class="agree">\u2713</td>')
                    else:
                        h.append("<td>\u2717</td>")
            h.append(f'<td class="agree"><strong>Shared ({len(dts_list)})</strong></td></tr>')
        for feat in sorted(unique.keys()):
            dts_list = unique[feat]
            h.append(f"<tr><td><code>{_esc(feat)}</code></td>")
            for dt in dwi_types_present:
                if dt in dts:
                    if dt in dts_list:
                        h.append(f"<td>{_dwi_badge(dt)}</td>")
                    else:
                        h.append("<td>\u2014</td>")
            h.append('<td>Type-specific</td></tr>')
        h.append("</tbody></table>")

    tot_st = sum([total_shared, total_unique])
    if tot_st > 0:
        pct_shared = 100.0 * float(total_shared) / float(tot_st)  # type: ignore
        cls = "agree" if pct_shared >= 50 else ("differ" if pct_shared < 25 else "")
        cls_attr = f' class="{cls}"' if cls else ""
        h.append(
            f'<div class="summary-box"><strong>Feature Overlap Summary:</strong> '
            f'<span{cls_attr}>{total_shared}/{total_shared + total_unique} '
            f'features ({pct_shared:.0f}%) are shared</span> across DWI types. '
            f'Shared features are more robust candidates for clinical biomarker '
            f'development.</div>'
        )

    # Feature importance context note (change 8)
    h.append(
        '<div class="info-box">'
        "<strong>Feature importance context:</strong> Features listed are those "
        "selected by elastic net regularization. Selection frequency across "
        "timepoints/DWI-types indicates robustness. Features stable across "
        "\u226575% of timepoints are candidate biomarkers for prospective "
        "validation."
        "</div>"
    )

    # ── Feature Stability Across Timepoints ──
    # Identify features that appear consistently across multiple timepoints
    # within the same DWI type (temporal stability = stronger biomarker).
    tp_stability_data = False
    for dt in dwi_types_present:
        if dt not in log_data:
            continue
        fs_list = log_data[dt].get("stats_predictive", {}).get("feature_selections", [])
        if len(fs_list) < 2:
            continue

        # Count how many timepoints each feature is selected at
        feat_tp_count: dict[str, int] = {}
        for fs in fs_list:
            for feat in fs.get("features", []):
                feat_tp_count[feat] = feat_tp_count.get(feat, 0) + 1

        n_timepoints = len(fs_list)
        stable_feats = {f: c for f, c in feat_tp_count.items() if c >= 2}
        if not stable_feats:
            continue

        if not tp_stability_data:
            h.append("<h3>Feature Stability Across Timepoints</h3>")
            h.append(
                '<p class="meta">Features selected at multiple timepoints within '
                'the same DWI type demonstrate temporal stability, suggesting they '
                'capture consistent biological signal rather than timepoint-specific noise.</p>'
            )
            tp_stability_data = True

        h.append(f"<p>{_dwi_badge(dt)} ({n_timepoints} timepoints):</p>")
        h.append("<table><thead><tr><th>Feature</th><th>Timepoints Selected</th>"
                 "<th>Stability</th></tr></thead><tbody>")
        for feat, count in sorted(stable_feats.items(), key=lambda x: -x[1]):
            pct = 100 * count / n_timepoints
            if pct >= 75:
                stab_cls = "agree"
                stab_lbl = "High"
            elif pct >= 50:
                stab_cls = ""
                stab_lbl = "Moderate"
            else:
                stab_cls = ""
                stab_lbl = "Low"
            s_attr = f' class="{stab_cls}"' if stab_cls else ""
            h.append(
                f"<tr><td><code>{_esc(feat)}</code></td>"
                f"<td>{count}/{n_timepoints} ({pct:.0f}%)</td>"
                f"<td{s_attr}><strong>{_esc(stab_lbl)}</strong></td></tr>"
            )
        h.append("</tbody></table>")

    # Potential duplicate feature detection (change 9)
    all_seen_feats: set[str] = set()
    for dt in dwi_types_present:
        if dt not in log_data:
            continue
        for fs in log_data[dt].get("stats_predictive", {}).get("feature_selections", []):
            all_seen_feats.update(fs.get("features", []))

    if all_seen_feats:
        def _feat_root(name: str) -> str:
            """Normalise a feature name to a canonical root for duplicate detection."""
            return name.lower().replace("_", "").replace("-", "").replace(" ", "")

        feat_list = sorted(all_seen_feats)
        duplicate_pairs: list[str] = []
        roots: dict[str, str] = {}
        for feat in feat_list:
            root = _feat_root(feat)
            if root in roots and roots[root] != feat:
                duplicate_pairs.append(f"{_esc(roots[root])} / {_esc(feat)}")
            else:
                roots[root] = feat

        if duplicate_pairs:
            h.append(
                '<div class="warn-box">'
                "\u26a0\ufe0f <strong>Potential duplicate features detected:</strong> "
                "<ul>"
                + "".join(f"<li><code>{pair}</code></li>" for pair in duplicate_pairs)
                + "</ul>"
                "Verify these are not the same metric exported under different names, "
                "which would inflate apparent feature stability."
                "</div>"
            )

    return h
