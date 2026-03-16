"""Report sections: figure gallery and graph appendix."""

from __future__ import annotations

import json
from pathlib import Path

from shared import (  # type: ignore
    DWI_TYPES,
    extract_correlations,
    extract_pvalues,
    parse_dwi_info,
    safe_text,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _figure_caption,
    _h2,
    _sig_class,
    _trend_tag,
)


def _section_appendix(rows) -> list[str]:
    """Build the Appendix: All Graphs section.

    Provides a compact index of every analysed graph with full metadata,
    extracted statistics, and clinical interpretation notes.  Detailed
    analysis alongside each image is presented in the Figure Gallery
    section; this appendix serves as a quick-reference table grouped by
    graph type.

    Parameters
    ----------
    rows : list[dict]
        Vision CSV rows (may be empty, in which case the section is skipped).

    Returns
    -------
    list[str]
        HTML chunks for the appendix section.
    """
    # ── 8. Appendix ──
    h: list[str] = []
    if rows:
        h.append(_h2("Appendix: All Graphs", "appendix"))
        h.append(f"<p>{len(rows)} graphs analysed. Each graph is shown with full "
                 "metadata, extracted statistics, and clinical interpretation notes. "
                 "See the <a href=\"#figure-gallery\">Figure Gallery</a> for images "
                 "with inline analysis.</p>")

        # Group by graph type for organised presentation
        type_groups: dict[str, list[tuple[int, dict]]] = {}
        for i, r in enumerate(rows, 1):
            gt = r.get("graph_type", "unknown")
            type_groups.setdefault(gt, []).append((i, r))

        for gt in sorted(type_groups.keys(), key=lambda k: -len(type_groups[k])):
            gt_rows = type_groups[gt]
            h.append(f"<h3>{_esc(gt)} ({len(gt_rows)} graphs)</h3>")

            for i, r in gt_rows:
                dwi_type, base_name = parse_dwi_info(r["file_path"])
                graph_title = r.get("graph_title", "") or ""

                # ── Build card header ──
                title_str = f" \u2014 {_esc(graph_title)}" if graph_title else ""
                h.append(f'<div class="graph-card">')
                h.append(f'<div class="graph-card-header">'
                         f'<span>#{i}</span> {_dwi_badge(dwi_type)} '
                         f'<span>{_esc(base_name)}{title_str}</span></div>')
                h.append('<dl class="graph-card-grid">')
                h.extend(_build_graph_analysis_html(r))
                h.append('</dl></div>')
    return h



def _build_graph_analysis_html(r: dict) -> list[str]:
    """Build the analysis details HTML for a single vision CSV row.

    Renders trends, statistics, inflection points, issues, and summary
    in a ``<dl>`` card format identical to the appendix layout.

    Parameters
    ----------
    r : dict
        A single row from the vision CSV.

    Returns
    -------
    list[str]
        HTML chunks for the analysis card body (inside a ``<dl>``).
    """
    h: list[str] = []

    # ── Axes ──
    x_lbl = r.get("x_axis_label", "")
    y_lbl = r.get("y_axis_label", "")
    c_lbl = r.get("color_axis_label", "")
    if x_lbl or y_lbl or c_lbl:
        axis_parts = []
        if x_lbl:
            x_unit = r.get("x_axis_units", "")
            xu = f" ({x_unit})" if x_unit else ""
            x_min = r.get("x_axis_range_min", "")
            x_max = r.get("x_axis_range_max", "")
            xr = f" [{x_min}\u2013{x_max}]" if x_min and x_max else ""
            axis_parts.append(f"X: {_esc(x_lbl)}{_esc(xu)}{_esc(xr)}")
        if y_lbl:
            y_unit = r.get("y_axis_units", "")
            yu = f" ({y_unit})" if y_unit else ""
            y_min = r.get("y_axis_range_min", "")
            y_max = r.get("y_axis_range_max", "")
            yr = f" [{y_min}\u2013{y_max}]" if y_min and y_max else ""
            axis_parts.append(f"Y: {_esc(y_lbl)}{_esc(yu)}{_esc(yr)}")
        if c_lbl:
            c_unit = r.get("color_axis_units", "")
            cu = f" ({c_unit})" if c_unit else ""
            axis_parts.append(f"Color: {_esc(c_lbl)}{_esc(cu)}")
        h.append(f'<dt>Axes</dt><dd>{" &bull; ".join(axis_parts)}</dd>')

    # ── Trends ──
    trends_str = r.get("trends_json", "[]") or "[]"
    try:
        trends_list = json.loads(trends_str)
    except Exception:
        trends_list = []
    if isinstance(trends_list, list) and trends_list:
        trend_parts = []
        for t in trends_list:
            if isinstance(t, dict):
                direction = t.get("direction", "")
                series = t.get("series", "")
                desc = t.get("description", "")
                label = f"{series}: {direction}" if series else direction
                tag = _trend_tag(label)
                if desc:
                    tag += f' <span class="axis-info">{_esc(desc[:100])}</span>'
                trend_parts.append(tag)
        if trend_parts:
            h.append(f'<dt>Trends</dt><dd>{"<br>".join(trend_parts)}</dd>')

    # ── Statistics ──
    all_text = safe_text(r, "summary", "trends_json", "inflection_points_json")
    pvals = extract_pvalues(all_text)
    corrs = extract_correlations(all_text)
    stats_parts = []
    for pval, ctx in pvals:
        cls = _sig_class(pval)
        cls_attr = f' class="{cls}"' if cls else ""
        stats_parts.append(f'<span{cls_attr}>p={pval:.4f}</span>'
                           f' <span class="axis-info">{_esc(ctx[:60])}</span>')
    for rval, ctx in corrs:
        if abs(rval) >= 0.3:
            stats_parts.append(f'r={rval:+.2f}'
                               f' <span class="axis-info">{_esc(ctx[:60])}</span>')
    if stats_parts:
        h.append(f'<dt>Statistics</dt><dd>{"<br>".join(stats_parts)}</dd>')

    # ── Issues ──
    issues_str = r.get("issues_json", "[]") or "[]"
    try:
        issues_list = json.loads(issues_str)
    except Exception:
        issues_list = []
    if isinstance(issues_list, list) and issues_list:
        items = "".join(f"<li>{_esc(iss)}</li>" for iss in issues_list)
        h.append(f'<dt>Issues</dt><dd><ul style="margin:0;padding-left:1.2em">{items}</ul></dd>')

    # ── Inflection points ──
    ips_str = r.get("inflection_points_json", "[]") or "[]"
    try:
        ips_list = json.loads(ips_str)
    except Exception:
        ips_list = []
    if isinstance(ips_list, list) and ips_list:
        ip_items = []
        for ip in ips_list:
            if isinstance(ip, dict):
                x = ip.get("approximate_x", "?")
                y = ip.get("approximate_y", "")
                desc_ip = ip.get("description", "")
                coord = f"x={_esc(str(x))}"
                if y is not None and y != "":
                    coord += f", y={_esc(str(y))}"
                ip_items.append(f"({coord}) {_esc(str(desc_ip))}")
        if ip_items:
            h.append(f'<dt>Inflection pts</dt><dd>{"<br>".join(ip_items)}</dd>')

    # ── Summary ──
    summary = r.get("summary", "") or ""
    if summary and not summary.startswith("JSON parse error"):
        h.append('<dt class="graph-card-full">Summary</dt>')
        if len(summary) > 200:
            h.append(
                f'<dd class="graph-card-full"><details open>'
                f'<summary>{_esc(summary[:150])}\u2026</summary>'
                f'<p>{_esc(summary)}</p></details></dd>'
            )
        else:
            h.append(f'<dd class="graph-card-full">{_esc(summary)}</dd>')

    return h


def _section_figure_gallery(folder, rows=None) -> list[str]:
    """Build a Figure Gallery section with embedded images and analysis.

    Embeds all pipeline-generated figures (PNG/JPG) from the output folder
    as base64 data URIs so the report is fully self-contained.  When vision
    CSV rows are provided, each figure is paired with its extracted analysis
    (trends, statistics, inflection points, and summary) so the reader sees
    each graph alongside its interpretation.

    Parameters
    ----------
    folder : Path
        Path to the ``saved_files_*`` output folder.
    rows : list[dict] or None
        Vision CSV rows (from ``graph_analysis_results.csv``).  When
        provided, each image is matched to its analysis row and the
        analysis is rendered directly below the image.

    Returns
    -------
    list[str]
        HTML chunks for the figure gallery section.
    """
    import base64

    h: list[str] = []
    figures: list[tuple[str, str, Path]] = []  # (dwi_type, clean_name, path)

    # Collect all image files from DWI-type subfolders
    for dt in DWI_TYPES:
        dt_dir = Path(folder) / dt
        if not dt_dir.is_dir():
            continue
        for ext in ("*.png", "*.jpg", "*.jpeg"):
            for img_path in sorted(dt_dir.glob(ext)):
                # Skip empty files
                if img_path.stat().st_size == 0:
                    continue
                # Generate a clean name from the filename
                clean = img_path.stem.replace("_", " ").replace("-", " ").title()
                figures.append((dt, clean, img_path))

    # Also check root folder for cross-DWI comparison figures
    root_folder = Path(folder)
    for ext in ("*.png", "*.jpg", "*.jpeg"):
        for img_path in sorted(root_folder.glob(ext)):
            if img_path.stat().st_size < 1024:
                continue
            clean = img_path.stem.replace("_", " ").replace("-", " ").title()
            figures.append(("Root", clean, img_path))

    if not figures:
        return h

    # Build a lookup from image filename to vision CSV row for matching.
    row_lookup: dict[str, dict] = {}
    if rows:
        for r in rows:
            fp = r.get("file_path", "")
            # Normalise to just the filename for matching.
            fname = Path(fp.replace("\\", "/")).name.lower()
            row_lookup[fname] = r

    h.append(_h2("Figure Gallery", "figure-gallery"))
    h.append(
        f'<p class="meta">All pipeline-generated figures are embedded below '
        f"with auto-numbered captions and analysis. {len(figures)} figure(s) "
        f"found across the output folder. Figures are grouped by DWI type.</p>"
    )

    current_dt = ""
    for dt, clean_name, img_path in figures:
        if dt != current_dt:
            current_dt = dt
            h.append(f"<h3>{_dwi_badge(dt) if dt != 'Root' else 'Cross-DWI'} Figures</h3>")

        # Read and encode the image as base64
        try:
            img_bytes = img_path.read_bytes()
            mime = "image/png" if img_path.suffix.lower() == ".png" else "image/jpeg"
            b64 = base64.b64encode(img_bytes).decode("ascii")
            caption = _figure_caption(
                f"{clean_name} ({dt})" if dt != "Root" else clean_name,
                f"Generated by the pancData3 pipeline."
            )
            h.append("<figure style=\"margin:1rem 0;page-break-inside:avoid\">")
            h.append(
                f'<img src="data:{mime};base64,{b64}" '
                f'alt="{_esc(clean_name)}" '
                f'style="max-width:100%;height:auto;border:1px solid var(--border);'
                f'border-radius:4px">'
            )
            h.append(caption)
            h.append("</figure>")

            # Match this image to its vision CSV analysis row
            match_key = img_path.name.lower()
            matched_row = row_lookup.get(match_key)
            if matched_row:
                graph_title = matched_row.get("graph_title", "") or ""
                graph_type = matched_row.get("graph_type", "") or ""
                type_str = f" ({_esc(graph_type)})" if graph_type else ""
                title_str = f" \u2014 {_esc(graph_title)}" if graph_title else ""
                h.append(f'<div class="graph-card">')
                h.append(
                    f'<div class="graph-card-header">'
                    f'{_dwi_badge(dt) if dt != "Root" else ""} '
                    f'<span>{_esc(clean_name)}{title_str}{type_str}</span></div>'
                )
                h.append('<dl class="graph-card-grid">')
                h.extend(_build_graph_analysis_html(matched_row))
                h.append('</dl></div>')
        except Exception:
            # Skip files that can't be read
            pass

    # Render analysis for vision CSV rows that had no matching image on disk.
    if rows:
        matched_fnames = {img_path.name.lower() for _, _, img_path in figures}
        unmatched = [r for r in rows
                     if Path(r.get("file_path", "").replace("\\", "/")).name.lower()
                     not in matched_fnames]
        if unmatched:
            h.append("<h3>Additional Analysed Graphs (no image on disk)</h3>")
            for i, r in enumerate(unmatched, 1):
                dwi_type, base_name = parse_dwi_info(r["file_path"])
                graph_title = r.get("graph_title", "") or ""
                title_str = f" \u2014 {_esc(graph_title)}" if graph_title else ""
                h.append(f'<div class="graph-card">')
                h.append(
                    f'<div class="graph-card-header">'
                    f'{_dwi_badge(dwi_type)} '
                    f'<span>{_esc(base_name)}{title_str}</span></div>'
                )
                h.append('<dl class="graph-card-grid">')
                h.extend(_build_graph_analysis_html(r))
                h.append('</dl></div>')

    return h
