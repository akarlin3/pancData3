"""Report section: cross-DWI feature overlap analysis."""

from __future__ import annotations

from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _stat_card,
    _table_caption,
)


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


def build_texture_section(saved_files: dict) -> str:
    """Build the Texture Features summary section.

    Renders a table of the top discriminating texture features (GLCM,
    GLRLM, first-order, shape) ranked by AUC from parsed pipeline output.

    Parameters
    ----------
    saved_files : dict
        Parsed pipeline output (keyed by DWI type at top level, with
        ``texture_features`` nested within).

    Returns
    -------
    str
        HTML string for the texture features section.  Empty string if
        no texture data is present.
    """
    if not saved_files:
        return ""

    h: list[str] = []
    found = False

    for dwi_type, dwi_data in saved_files.items():
        if not isinstance(dwi_data, dict):
            continue
        tex = dwi_data.get("texture_features")
        if not tex or not isinstance(tex, (dict, list)):
            continue

        # Accept either a list of feature dicts or a dict with a "features" key
        if isinstance(tex, dict):
            features = tex.get("features", [])
        else:
            features = tex

        if not features or not isinstance(features, list):
            continue

        if not found:
            h.append(_h2("Texture Features", "texture-features"))
            h.append(
                '<p class="meta">Texture features extracted from DWI parameter '
                "maps (GLCM, GLRLM, first-order statistics, and shape descriptors) "
                "ranked by discriminative ability (AUC). Higher AUC indicates "
                "stronger separation between outcome groups.</p>"
            )
            found = True

        h.append(f"<h3>{_dwi_badge(dwi_type)}</h3>")

        # Sort by AUC descending, take top features
        scored = [
            f for f in features
            if isinstance(f, dict) and isinstance(f.get("auc"), (int, float))
        ]
        scored.sort(key=lambda f: f["auc"], reverse=True)

        # Summary cards
        if scored:
            categories: dict[str, int] = {}
            for f in features:
                if isinstance(f, dict):
                    cat = f.get("category", "unknown")
                    categories[cat] = categories.get(cat, 0) + 1

            cards = [
                _stat_card("Total Features", str(len(features)),
                           "extracted from parameter maps"),
                _stat_card("Top AUC", f"{scored[0]['auc']:.3f}",
                           _esc(str(scored[0].get("name", "?")))),
            ]
            if categories:
                cat_summary = ", ".join(
                    f"{v} {k}" for k, v in sorted(
                        categories.items(), key=lambda x: -x[1]
                    )
                )
                cards.append(_stat_card("Categories", cat_summary))
            h.append('<div class="stat-grid">')
            h.extend(cards)
            h.append("</div>")

        # Table of top features (up to 20)
        top_n = scored[:20]
        if top_n:
            h.append(_table_caption(
                "Top Texture Features by AUC",
                f"Top {len(top_n)} discriminating texture features for "
                f"{_esc(dwi_type)}.",
            ))
            h.append(
                "<table><thead><tr>"
                "<th>Feature</th><th>Category</th>"
                "<th>AUC</th><th>p-value</th>"
                "</tr></thead><tbody>"
            )
            for f in top_n:
                name = f.get("name", "?")
                cat = f.get("category", "?")
                auc = f["auc"]
                p = f.get("p")
                auc_cls = "agree" if auc >= 0.7 else ""
                auc_attr = f' class="{auc_cls}"' if auc_cls else ""
                p_str = f"{p:.4f}" if isinstance(p, (int, float)) else "N/A"
                h.append(
                    f"<tr><td><code>{_esc(str(name))}</code></td>"
                    f"<td>{_esc(str(cat))}</td>"
                    f"<td{auc_attr}><strong>{auc:.3f}</strong></td>"
                    f"<td>{p_str}</td></tr>"
                )
            h.append("</tbody></table>")

    if not found:
        return ""
    return "\n".join(h)
