"""Report section: treatment response / longitudinal trends."""

from __future__ import annotations

import json

from shared import (  # type: ignore
    DWI_TYPES,
    extract_correlations,
    extract_pvalues,
    safe_text,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
    _sig_class,
    _sig_tag,
    _trend_tag,
)


def _section_treatment_response(groups) -> list[str]:
    """Build the Treatment Response / Longitudinal Trends section.

    For each graph with "Longitudinal" in its name, displays per-DWI-type
    trend tags, axis information, inflection points (in collapsible
    ``<details>`` blocks), and full summaries.

    Parameters
    ----------
    groups : dict[str, dict[str, dict]]
        Graph rows grouped by base name and DWI type.

    Returns
    -------
    list[str]
        HTML chunks for the treatment response section.
    """
    # ── 6. Treatment Response ──
    h: list[str] = []
    if groups:
        has_longitudinal = any("Longitudinal" in bn for bn in groups)
        if has_longitudinal:
            h.append(_h2("Treatment Response \u2014 Longitudinal Trends", "treatment"))
            for base_name in sorted(groups.keys()):
                if "Longitudinal" not in base_name:
                    continue
                dwi_dict = groups[base_name]
                real = [t for t in dwi_dict if t != "Root"]
                if not real:
                    continue

                h.append(f"<h3>{_esc(base_name)}</h3>")
                for dt in DWI_TYPES:
                    if dt not in dwi_dict:
                        continue
                    r = dwi_dict[dt]
                    summary = r.get("summary", "")

                    # Axis info with range analysis
                    x_lbl = r.get("x_axis_label", "")
                    x_unit = r.get("x_axis_units", "")
                    x_min = r.get("x_axis_range_min", "")
                    x_max = r.get("x_axis_range_max", "")
                    y_lbl = r.get("y_axis_label", "")
                    y_unit = r.get("y_axis_units", "")
                    y_min = r.get("y_axis_range_min", "")
                    y_max = r.get("y_axis_range_max", "")
                    axis_parts = []
                    if x_lbl:
                        x_range = ""
                        if x_min and x_max:
                            x_range = f" [{x_min}\u2013{x_max}]"
                        axis_parts.append(f"X: {x_lbl}" + (f" ({x_unit})" if x_unit else "") + x_range)
                    if y_lbl:
                        y_range = ""
                        if y_min and y_max:
                            y_range = f" [{y_min}\u2013{y_max}]"
                        axis_parts.append(f"Y: {y_lbl}" + (f" ({y_unit})" if y_unit else "") + y_range)

                    h.append(f"<p>{_dwi_badge(dt)}")
                    if axis_parts:
                        h.append(f' <span class="axis-info">[{" | ".join(axis_parts)}]</span>')

                    # Extract any p-values mentioned in this specific graph
                    all_text = safe_text(r, "summary", "trends_json", "inflection_points_json")
                    graph_pvals = extract_pvalues(all_text)
                    graph_corrs = extract_correlations(all_text)
                    if graph_pvals or graph_corrs:
                        stat_parts = []
                        for pv, _ in graph_pvals:
                            cls = _sig_class(pv)
                            cls_attr = f' class="{cls}"' if cls else ""
                            stat_parts.append(f'<span{cls_attr}>p={pv:.4f}{_sig_tag(pv)}</span>')
                        for rv, _ in graph_corrs:
                            if abs(rv) >= 0.3:
                                stat_parts.append(f"r={rv:+.2f}")
                        if stat_parts:
                            h.append(f' <span class="meta">[{", ".join(stat_parts)}]</span>')
                    h.append("</p>")

                    # Trends
                    trends_str = r.get("trends_json", "[]") if isinstance(r, dict) else "[]"
                    try:
                        trends = json.loads(str(trends_str))
                    except Exception:
                        trends = []
                    if isinstance(trends, list) and trends:
                        for t in trends:
                            if isinstance(t, dict):
                                series = t.get("series", "")
                                direction = t.get("direction", "")
                                desc = t.get("description", "")
                                h.append(f"<p style='margin-left:1rem'>")
                                if series:
                                    h.append(f"<code>{_esc(series)}</code> ")
                                h.append(_trend_tag(direction))
                                if desc:
                                    h.append(f" <em>{_esc(desc)}</em>")
                                h.append("</p>")

                    # Inflection points
                    ips_str = r.get("inflection_points_json", "[]") if isinstance(r, dict) else "[]"
                    try:
                        ips = json.loads(str(ips_str))
                    except Exception:
                        ips = []
                    if isinstance(ips, list) and ips:
                        h.append("<details open><summary>Inflection points</summary><ul>")
                        for ip in ips:
                            if isinstance(ip, dict):
                                x = ip.get("approximate_x", "?")
                                y = ip.get("approximate_y", "")
                                desc_ip = ip.get("description", "")
                                coord = f"x={_esc(str(x))}"
                                if y is not None and y != "":
                                    coord += f", y={_esc(str(y))}"
                                h.append(f"<li>({coord}): {_esc(str(desc_ip))}</li>")
                        h.append("</ul></details>")

                    # Full summary in collapsible
                    if summary:
                        if len(summary) > 200:
                            h.append(f"<details open><summary>Full summary</summary>"
                                     f'<p class="full-summary">{_esc(summary)}</p></details>')
                        else:
                            h.append(f'<p class="full-summary" style="margin-left:1rem">{_esc(summary)}</p>')

                # Cross-DWI consensus for this longitudinal graph
                if len(real) >= 2:
                    all_dirs: dict[str, list[str]] = {}
                    for dt in real:
                        r = dwi_dict[dt]
                        t_str = r.get("trends_json", "[]") if isinstance(r, dict) else "[]"
                        try:
                            t_list = json.loads(str(t_str))
                        except Exception:
                            t_list = []
                        if isinstance(t_list, list):
                            for t in t_list:
                                if isinstance(t, dict):
                                    s = t.get("series") or "overall"
                                    d = str(t.get("direction", ""))
                                    all_dirs.setdefault(s, []).append(d)
                    consensus_parts = []
                    for s, dirs in all_dirs.items():
                        if len(dirs) >= 2:
                            if len(set(dirs)) == 1:
                                consensus_parts.append(f"{s}: all {dirs[0]}")
                            else:
                                consensus_parts.append(f"{s}: mixed ({', '.join(dirs)})")
                    if consensus_parts:
                        cons_cls = "agree" if all("all" in p for p in consensus_parts) else "differ"
                        h.append(
                            f'<div class="info-box"><strong>Cross-DWI consensus:</strong> '
                            f'<span class="{cons_cls}">{"; ".join(consensus_parts)}</span></div>'
                        )

        # Non-longitudinal inflection points
        non_long_ips = []
        for base_name in sorted(groups.keys()):
            if "Longitudinal" in base_name:
                continue
            dwi_dict = groups[base_name]
            for dt in DWI_TYPES:
                if dt not in dwi_dict:
                    continue
                r = dwi_dict[dt]
                ips_str = r.get("inflection_points_json", "[]") if isinstance(r, dict) else "[]"
                try:
                    ips = json.loads(str(ips_str))
                except Exception:
                    ips = []
                if isinstance(ips, list):
                    for ip in ips:
                        if isinstance(ip, dict) and ip.get("description"):
                            non_long_ips.append((base_name, dt, ip))

        if non_long_ips:
            h.append("<h3>Non-Longitudinal Inflection Points</h3>")
            h.append('<p class="meta">Inflection points detected in non-longitudinal graphs '
                     'may indicate dose thresholds, feature boundaries, or regime transitions.</p>')
            h.append("<table><thead><tr><th>Graph</th><th>DWI</th>"
                     "<th>Location</th><th>Description</th></tr></thead><tbody>")
            for gname, dt, ip in non_long_ips:
                x = ip.get("approximate_x", "?")
                y = ip.get("approximate_y", "")
                coord = f"x={_esc(str(x))}"
                if y is not None and y != "":
                    coord += f", y={_esc(str(y))}"
                desc_ip = _esc(str(ip.get("description", "")))
                h.append(f"<tr><td>{_esc(gname)}</td>"
                         f"<td>{_dwi_badge(dt)}</td>"
                         f"<td><code>{coord}</code></td>"
                         f"<td>{desc_ip}</td></tr>")
            h.append("</tbody></table>")
    return h
