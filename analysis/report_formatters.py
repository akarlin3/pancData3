"""Formatting utilities and constants for the HTML analysis report.

This module is imported by :mod:`generate_report` and :mod:`report_sections`
to keep presentation logic separate from data loading and section assembly.

Contents:

- **HTML escaping** (:func:`_esc`) and **section headings** (:func:`_h2`,
  :func:`_section`).
- **Significance helpers** (:func:`_sig_tag`, :func:`_sig_class`) --
  map p-values to visual indicators (asterisks and CSS classes).
- **CSS stylesheet** (:data:`CSS`) -- embedded in the ``<style>`` tag of
  the HTML report for self-contained styling.
- **DWI badge** (:func:`_dwi_badge`) -- colour-coded ``<span>`` labels
  for Standard / dnCNN / IVIMnet.
- **Trend tag** (:func:`_trend_tag`) -- directional arrow badges for
  increasing / decreasing / stable / non-monotonic trends.
- **Navigation bar** (:func:`_nav_bar`, :data:`NAV_SECTIONS`) -- sticky
  top-of-page anchor links.
- **Stat card** (:func:`_stat_card`) -- summary metric display cards.
- **Consensus helper** (:func:`_get_consensus`) -- vote-based trend
  consensus across DWI types.
- **HTML template** (:data:`HTML_TEMPLATE`) -- standalone HTML wrapper
  for Markdown-to-HTML conversion.
"""

from __future__ import annotations

import html as _html


def _esc(text: str) -> str:
    """HTML-escape a string to prevent XSS and rendering issues.

    Parameters
    ----------
    text : str
        Raw text to escape.

    Returns
    -------
    str
        Escaped HTML string safe for embedding in ``<td>``, ``<p>``, etc.
    """
    return _html.escape(str(text))


def _section(title: str, level: int = 2) -> str:
    """Return a Markdown-style section heading string.

    This is a utility / test helper used by the legacy Markdown report
    path (now superseded by HTML generation).

    Parameters
    ----------
    title : str
        Section title text.
    level : int
        Heading level (number of ``#`` characters).

    Returns
    -------
    str
        Markdown heading string.
    """
    return f"\n{'#' * level} {title}\n"


def _sig_tag(p: float) -> str:
    """Return asterisk significance markers for a p-value.

    Follows biomedical convention:
    - ``***`` for p < 0.001
    - ``**``  for p < 0.01
    - ``*``   for p < 0.05
    - ``""``  for p >= 0.05

    Parameters
    ----------
    p : float
        P-value.

    Returns
    -------
    str
        Significance marker string.
    """
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def _sig_class(p: float) -> str:
    """Return a CSS class name for the significance level.

    Classes ``sig-1``, ``sig-2``, ``sig-3`` map to amber, red, and bold-red
    styling defined in :data:`CSS`.

    Parameters
    ----------
    p : float
        P-value.

    Returns
    -------
    str
        CSS class name, or empty string if not significant.
    """
    if p < 0.001:
        return "sig-3"
    if p < 0.01:
        return "sig-2"
    if p < 0.05:
        return "sig-1"
    return ""


# ── CSS stylesheet ────────────────────────────────────────────────────────────
# Embedded directly in the HTML report's <style> tag so the report is a
# single self-contained file with no external dependencies.  Uses CSS custom
# properties (variables) for theming consistency.
CSS = """\
:root {
    --bg: #ffffff; --fg: #1a1a2e; --muted: #64748b;
    --accent: #2563eb; --accent-light: #dbeafe;
    --border: #e2e8f0; --row-alt: #f8fafc;
    --green: #16a34a; --red: #dc2626; --amber: #d97706;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    color: var(--fg); background: var(--bg);
    line-height: 1.6; max-width: 1100px; margin: 0 auto; padding: 2rem 1.5rem;
}
h1 { font-size: 1.8rem; border-bottom: 3px solid var(--accent); padding-bottom: 0.5rem; margin-bottom: 1.5rem; }
h2 { font-size: 1.35rem; color: var(--accent); margin: 2rem 0 0.75rem; border-bottom: 1px solid var(--border); padding-bottom: 0.3rem; }
h3 { font-size: 1.1rem; margin: 1.25rem 0 0.5rem; }
h4 { font-size: 0.98rem; margin: 1rem 0 0.4rem; color: var(--muted); }
p, ul { margin-bottom: 0.75rem; }
ul { padding-left: 1.5rem; }
.meta { color: var(--muted); font-size: 0.9rem; margin-bottom: 1.5rem; }
.meta span { display: inline-block; margin-right: 1.5rem; }
code { background: var(--row-alt); border: 1px solid var(--border); border-radius: 3px; padding: 0.1em 0.35em; font-size: 0.88em; }
table { width: 100%; border-collapse: collapse; margin: 0.75rem 0 1.25rem; font-size: 0.92rem; }
th { background: var(--accent); color: #fff; text-align: left; padding: 0.5rem 0.75rem; font-weight: 600; }
td { padding: 0.4rem 0.75rem; border-bottom: 1px solid var(--border); vertical-align: top; }
tr:nth-child(even) td { background: var(--row-alt); }
tr:hover td { background: var(--accent-light); }
.sig-1 { color: var(--amber); font-weight: 600; }
.sig-2 { color: var(--red); font-weight: 600; }
.sig-3 { color: var(--red); font-weight: 700; }
.agree { color: var(--green); font-weight: 600; }
.differ { color: var(--red); font-weight: 700; }
.badge { display: inline-block; padding: 0.15em 0.5em; border-radius: 4px; font-size: 0.82rem; font-weight: 600; }
.badge-standard { background: #dbeafe; color: #1e40af; }
.badge-dncnn { background: #dcfce7; color: #166534; }
.badge-ivimnet { background: #fef3c7; color: #92400e; }
.badge-root { background: #f1f5f9; color: #475569; }
.summary-box { background: var(--accent-light); border-left: 4px solid var(--accent); padding: 1rem 1.25rem; border-radius: 0 6px 6px 0; margin: 1rem 0; }
.warn-box { background: #fef3c7; border-left: 4px solid var(--amber); padding: 0.75rem 1rem; border-radius: 0 6px 6px 0; margin: 0.75rem 0; font-size: 0.9rem; }
.info-box { background: #f0fdf4; border-left: 4px solid var(--green); padding: 0.75rem 1rem; border-radius: 0 6px 6px 0; margin: 0.75rem 0; font-size: 0.9rem; }
footer { margin-top: 3rem; padding-top: 1rem; border-top: 1px solid var(--border); color: var(--muted); font-size: 0.85rem; }
nav.toc {
    position: sticky; top: 0; background: var(--bg);
    border-bottom: 2px solid var(--border); padding: 0.5rem 0 0.4rem;
    margin-bottom: 1.5rem; z-index: 100; overflow-x: auto;
    white-space: nowrap; font-size: 0.83rem;
}
nav.toc a { color: var(--accent); text-decoration: none; margin-right: 1.1rem; }
nav.toc a:hover { text-decoration: underline; }
details { margin: 0.3rem 0; }
details > summary {
    cursor: pointer; color: var(--accent); font-size: 0.88rem;
    padding: 0.2rem 0.1rem; list-style: none; user-select: none;
}
details > summary::before { content: "\\25B6\\00A0"; font-size: 0.7em; }
details[open] > summary::before { content: "\\25BC\\00A0"; }
details[open] > summary { margin-bottom: 0.4rem; }
.axis-info { font-size: 0.82rem; color: var(--muted); line-height: 1.4; }
.trend-tag {
    display: inline-block; padding: 0.1em 0.45em; border-radius: 3px;
    font-size: 0.8rem; margin: 0.1em 0.1em; border: 1px solid var(--border);
}
.trend-incr { background: #dcfce7; border-color: #86efac; color: #166534; }
.trend-decr { background: #fee2e2; border-color: #fca5a5; color: #991b1b; }
.trend-flat { background: var(--row-alt); }
.trend-nm   { background: #f5f3ff; border-color: #c4b5fd; color: #5b21b6; }
.full-summary { font-size: 0.88rem; color: #374151; margin-top: 0.25rem; }
.stat-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(220px, 1fr)); gap: 0.75rem; margin: 0.75rem 0 1rem; }
.stat-card { background: var(--row-alt); border: 1px solid var(--border); border-radius: 6px; padding: 0.75rem 1rem; }
.stat-card .label { font-size: 0.78rem; color: var(--muted); text-transform: uppercase; letter-spacing: 0.04em; }
.stat-card .value { font-size: 1.35rem; font-weight: 700; color: var(--accent); }
.stat-card .sub { font-size: 0.8rem; color: var(--muted); }
.abstract-box { background: #f8fafc; border: 1px solid var(--border); border-radius: 8px; padding: 1.25rem 1.5rem; margin: 1rem 0; }
.abstract-box h4 { color: var(--accent); margin: 0.75rem 0 0.3rem; font-size: 0.95rem; text-transform: uppercase; letter-spacing: 0.03em; }
.abstract-box h4:first-child { margin-top: 0; }
.abstract-box p { margin-bottom: 0.5rem; font-size: 0.93rem; }
.methods-box { background: #fafbfc; border: 1px solid var(--border); border-radius: 6px; padding: 1rem 1.25rem; margin: 0.75rem 0; font-size: 0.92rem; line-height: 1.65; }
.methods-box p { margin-bottom: 0.6rem; }
.methods-box strong { color: var(--fg); }
.forest-row { display: flex; align-items: center; height: 1.2rem; position: relative; }
.forest-bar { height: 3px; background: var(--accent); position: absolute; }
.forest-point { width: 8px; height: 8px; border-radius: 50%; position: absolute; transform: translate(-50%, -50%); top: 50%; }
.forest-point-sig { background: var(--red); }
.forest-point-ns { background: var(--muted); }
.forest-ref { position: absolute; width: 1px; height: 100%; background: #999; top: 0; }
.effect-sm { color: var(--muted); }
.effect-md { color: var(--amber); font-weight: 600; }
.effect-lg { color: var(--red); font-weight: 700; }
.ci-text { font-size: 0.82rem; color: var(--muted); }
.limitation-list li { margin-bottom: 0.4rem; line-height: 1.55; }
.conclusion-box { background: #f0f7ff; border-left: 4px solid var(--accent); padding: 1rem 1.25rem; border-radius: 0 6px 6px 0; margin: 0.75rem 0; }
.conclusion-box li { margin-bottom: 0.35rem; }
.diag-box { background: #fffbeb; border: 1px solid #fde68a; border-radius: 6px; padding: 0.75rem 1rem; margin: 0.5rem 0; font-size: 0.9rem; }
"""


def _dwi_badge(dwi_type: str) -> str:
    """Return a colour-coded HTML badge ``<span>`` for a DWI type.

    Each DWI type gets a distinct background colour:
    - Standard: blue
    - dnCNN: green
    - IVIMnet: amber
    - Other/Root: grey

    Parameters
    ----------
    dwi_type : str
        DWI type name.

    Returns
    -------
    str
        HTML ``<span class="badge ...">`` string.
    """
    cls = {
        "Standard": "badge-standard", "dnCNN": "badge-dncnn",
        "IVIMnet": "badge-ivimnet",
    }.get(dwi_type, "badge-root")
    return f'<span class="badge {cls}">{_esc(dwi_type)}</span>'


def _trend_tag(direction: str) -> str:
    """Return a directional arrow badge ``<span>`` for a trend direction.

    Keyword matching on the direction string determines the arrow and
    colour class:
    - Increasing/up/higher/rising: green with up-arrow
    - Decreasing/down/lower/falling/drop: red with down-arrow
    - Flat/stable/constant: grey with right-arrow
    - Anything else (non-monotonic, U-shaped): purple, no arrow

    Parameters
    ----------
    direction : str
        Trend direction text from the vision model.

    Returns
    -------
    str
        HTML ``<span class="trend-tag ...">`` string.
    """
    d = direction.lower()
    if "increas" in d or "up" in d or "higher" in d or "rising" in d:
        cls = "trend-incr"
        arrow = "\u2191\u00a0"  # Up arrow + non-breaking space
    elif "decreas" in d or "down" in d or "lower" in d or "falling" in d or "drop" in d:
        cls = "trend-decr"
        arrow = "\u2193\u00a0"  # Down arrow
    elif "flat" in d or "stable" in d or "constant" in d:
        cls = "trend-flat"
        arrow = "\u2192\u00a0"  # Right arrow
    else:
        cls = "trend-nm"  # Non-monotonic
        arrow = ""
    return f'<span class="trend-tag {cls}">{arrow}{_esc(direction)}</span>'


# ── Navigation sections ────────────────────────────────────────────────────────
# Defines the order and labels for the sticky top-of-page navigation bar.
# Each tuple is (anchor_id, display_label).

NAV_SECTIONS = [
    ("exec-summary", "Abstract"),
    ("methods", "Methods"),
    ("cohort", "Cohort"),
    ("data-quality", "Data Quality"),
    ("hypothesis", "Hypothesis"),
    ("graph-overview", "Graphs"),
    ("stats-by-type", "By Type"),
    ("significance", "Statistics"),
    ("effect-sizes", "Effect Sizes"),
    ("mult-comp", "Corrections"),
    ("cross-dwi", "Cross-DWI"),
    ("fdr-global", "FDR Global"),
    ("correlations", "Correlations"),
    ("treatment", "Treatment"),
    ("predictive", "Predictive"),
    ("model-diag", "Diagnostics"),
    ("supplemental", "Supplemental"),
    ("limitations", "Limitations"),
    ("conclusions", "Conclusions"),
    ("appendix", "All Graphs"),
]


def _nav_bar() -> str:
    """Build the sticky navigation bar HTML from :data:`NAV_SECTIONS`.

    Returns
    -------
    str
        HTML ``<nav class="toc">`` element with anchor links.
    """
    links = "".join(
        f'<a href="#{anchor}">{_esc(label)}</a>'
        for anchor, label in NAV_SECTIONS
    )
    return f'<nav class="toc">{links}</nav>'


def _h2(text: str, anchor: str) -> str:
    """Return an ``<h2>`` element with an ``id`` attribute for anchor linking.

    Parameters
    ----------
    text : str
        Heading text.
    anchor : str
        HTML ``id`` attribute (used by navigation bar links).

    Returns
    -------
    str
        HTML ``<h2>`` string.
    """
    return f'<h2 id="{anchor}">{_esc(text)}</h2>'


def _stat_card(label: str, value: str, sub: str = "") -> str:
    """Return a summary statistic card HTML block.

    Cards are displayed in a CSS grid and show a label (small caps),
    a large value, and an optional subtitle.

    Parameters
    ----------
    label : str
        Card label (e.g. "Best AUC").
    value : str
        Main display value (e.g. "0.843").
    sub : str, optional
        Subtitle text below the value (e.g. "across all DWI types").

    Returns
    -------
    str
        HTML ``<div class="stat-card">`` block.
    """
    sub_html = f'<div class="sub">{_esc(sub)}</div>' if sub else ""
    return (
        f'<div class="stat-card">'
        f'<div class="label">{_esc(label)}</div>'
        f'<div class="value">{_esc(value)}</div>'
        f'{sub_html}</div>'
    )


def _forest_plot_cell(hr: float, ci_lo: float, ci_hi: float, p: float,
                      log_min: float = -1.5, log_max: float = 1.5) -> str:
    """Return an inline HTML forest plot cell for a single hazard ratio.

    Renders a horizontal bar representing the 95% CI with a point
    estimate dot, plus a vertical reference line at HR=1.0.

    Parameters
    ----------
    hr : float
        Hazard ratio point estimate.
    ci_lo, ci_hi : float
        Lower and upper bounds of the 95% confidence interval.
    p : float
        P-value for significance-based colouring.
    log_min, log_max : float
        Log-scale axis bounds for positioning.

    Returns
    -------
    str
        HTML ``<div>`` representing the forest plot bar.
    """
    import math
    width = 180  # pixels

    def _pos(val: float) -> float:
        lv = math.log(max(val, 0.01))
        return max(0, min(width, (lv - log_min) / (log_max - log_min) * width))

    ref_x = _pos(1.0)
    hr_x = _pos(hr)
    lo_x = _pos(ci_lo)
    hi_x = _pos(ci_hi)
    pt_cls = "forest-point-sig" if p < 0.05 else "forest-point-ns"

    return (
        f'<div class="forest-row" style="width:{width}px">'
        f'<div class="forest-ref" style="left:{ref_x:.0f}px"></div>'
        f'<div class="forest-bar" style="left:{lo_x:.0f}px;width:{max(hi_x - lo_x, 1):.0f}px"></div>'
        f'<div class="forest-point {pt_cls}" style="left:{hr_x:.0f}px"></div>'
        f'</div>'
    )


def _effect_size_class(d: float) -> str:
    """Return a CSS class for an effect size magnitude.

    Uses Cohen's conventions: small (0.2), medium (0.5), large (0.8).

    Parameters
    ----------
    d : float
        Absolute effect size (Cohen's d or similar).

    Returns
    -------
    str
        CSS class name.
    """
    d = abs(d)
    if d >= 0.8:
        return "effect-lg"
    if d >= 0.5:
        return "effect-md"
    return "effect-sm"


def _effect_size_label(d: float) -> str:
    """Return a human-readable label for an effect size magnitude.

    Parameters
    ----------
    d : float
        Absolute effect size.

    Returns
    -------
    str
        One of ``"Large"``, ``"Medium"``, or ``"Small"``.
    """
    d = abs(d)
    if d >= 0.8:
        return "Large"
    if d >= 0.5:
        return "Medium"
    return "Small"


def _get_consensus(trend_list: list[str]) -> str:
    """Determine the consensus trend direction from a list of direction strings.

    Uses simple keyword voting: counts how many entries contain increasing-
    related keywords vs decreasing-related keywords.

    Parameters
    ----------
    trend_list : list[str]
        Direction strings (e.g. ``["increasing", "decreasing", "rising"]``).

    Returns
    -------
    str
        One of ``"increasing"``, ``"decreasing"``, ``"stable"``, or
        ``"unknown"`` (if the list is empty).
    """
    if not trend_list: return "unknown"
    increasers = sum(1 for x in trend_list if "increas" in x or "higher" in x or "up" in x)
    decreasers = sum(1 for x in trend_list if "decreas" in x or "lower" in x or "down" in x)
    if increasers > decreasers: return "increasing"
    if decreasers > increasers: return "decreasing"
    return "stable"


# ── HTML template for Markdown-to-HTML conversion ────────────────────────────
# Used by generate_report.markdown_to_html() to wrap rendered Markdown in a
# styled HTML document.  The ``{title}`` and ``{body}`` placeholders are
# filled via str.format().  Curly braces in CSS are doubled to escape them.
HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
  :root {{
    --bg: #ffffff;
    --fg: #1a1a2e;
    --accent: #0f3460;
    --accent-light: #e8eef6;
    --border: #d0d7de;
    --table-stripe: #f6f8fa;
    --sig-strong: #d32f2f;
    --sig-moderate: #e65100;
    --sig-mild: #f9a825;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
                 "Helvetica Neue", Arial, sans-serif;
    line-height: 1.6;
    color: var(--fg);
    background: var(--bg);
    max-width: 1100px;
    margin: 0 auto;
    padding: 2rem 1.5rem;
  }}
  h1 {{
    font-size: 1.8rem;
    color: var(--accent);
    border-bottom: 3px solid var(--accent);
    padding-bottom: 0.5rem;
    margin-bottom: 1rem;
  }}
  h2 {{
    font-size: 1.4rem;
    color: var(--accent);
    border-bottom: 1px solid var(--border);
    padding-bottom: 0.3rem;
    margin-top: 2rem;
    margin-bottom: 0.8rem;
  }}
  h3 {{
    font-size: 1.15rem;
    margin-top: 1.5rem;
    margin-bottom: 0.5rem;
  }}
  p {{ margin-bottom: 0.6rem; }}
  code {{
    background: var(--accent-light);
    padding: 0.15em 0.4em;
    border-radius: 3px;
    font-size: 0.9em;
  }}
  table {{
    border-collapse: collapse;
    width: 100%;
    margin: 0.8rem 0 1.2rem;
    font-size: 0.92rem;
  }}
  th, td {{
    border: 1px solid var(--border);
    padding: 0.45rem 0.7rem;
    text-align: left;
  }}
  th {{
    background: var(--accent);
    color: #fff;
    font-weight: 600;
  }}
  tr:nth-child(even) {{ background: var(--table-stripe); }}
  tr:hover {{ background: var(--accent-light); }}
  ul, ol {{ margin: 0.5rem 0 0.5rem 1.5rem; }}
  li {{ margin-bottom: 0.25rem; }}
  strong {{ color: var(--accent); }}
  hr {{
    border: none;
    border-top: 1px solid var(--border);
    margin: 2rem 0;
  }}
  em {{ color: #555; }}
  .footer {{
    margin-top: 2rem;
    padding-top: 1rem;
    border-top: 1px solid var(--border);
    font-size: 0.85rem;
    color: #666;
  }}
</style>
</head>
<body>
{body}
</body>
</html>
"""
