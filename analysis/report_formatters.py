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

from shared import get_config


# ── Table / Figure numbering ─────────────────────────────────────────────────
# A simple counter-based numbering system for publication-quality
# cross-referencing.  Each call to ``_table_caption`` or ``_figure_caption``
# increments the counter and returns an HTML ``<caption>`` or ``<figcaption>``.

class _NumberingContext:
    """Thread-local numbering state for tables and figures.

    Call :meth:`reset` at the start of each report generation to begin
    numbering from 1.  Also tracks titles for generating a List of Tables
    and List of Figures at the end of the report.
    """

    def __init__(self):
        self.table_num = 0
        self.figure_num = 0
        self.table_titles: list[tuple[int, str]] = []
        self.figure_titles: list[tuple[int, str]] = []

    def reset(self):
        """Reset all counters (call at the start of report generation)."""
        self.table_num = 0
        self.figure_num = 0
        self.table_titles = []
        self.figure_titles = []

    def next_table(self, title: str = "") -> int:
        """Increment and return the next table number."""
        self.table_num += 1
        if title:
            self.table_titles.append((self.table_num, title))
        return self.table_num

    def next_figure(self, title: str = "") -> int:
        """Increment and return the next figure number."""
        self.figure_num += 1
        if title:
            self.figure_titles.append((self.figure_num, title))
        return self.figure_num


_numbering = _NumberingContext()


def reset_numbering() -> None:
    """Reset table/figure counters for a new report."""
    _numbering.reset()


def get_numbering() -> _NumberingContext:
    """Return the global numbering context (for generating table/figure indices)."""
    return _numbering


def _table_caption(title: str, description: str = "") -> str:
    """Return an HTML ``<caption>`` with an auto-incremented table number.

    Parameters
    ----------
    title : str
        Short caption title (e.g. "Hazard Ratio Effect Sizes").
    description : str, optional
        Additional descriptive text displayed after the title.

    Returns
    -------
    str
        HTML ``<caption>`` element.
    """
    num = _numbering.next_table(title)
    desc_html = f" {_html.escape(description)}" if description else ""
    return (
        f'<caption style="caption-side:top;text-align:left;font-size:0.9rem;'
        f'color:#374151;padding:0.4rem 0;font-weight:600">'
        f'Table {num}. {_html.escape(title)}'
        f'<span style="font-weight:400">{desc_html}</span></caption>'
    )


def _figure_caption(title: str, description: str = "") -> str:
    """Return an HTML ``<figcaption>`` with an auto-incremented figure number.

    Parameters
    ----------
    title : str
        Short caption title (e.g. "Longitudinal ADC Trends").
    description : str, optional
        Additional descriptive text displayed after the title.

    Returns
    -------
    str
        HTML ``<figcaption>`` element.
    """
    num = _numbering.next_figure(title)
    desc_html = f" {_html.escape(description)}" if description else ""
    return (
        f'<figcaption style="text-align:left;font-size:0.9rem;'
        f'color:#374151;padding:0.4rem 0;font-weight:600">'
        f'Figure {num}. {_html.escape(title)}'
        f'<span style="font-weight:400">{desc_html}</span></figcaption>'
    )


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

    Thresholds are read from the centralised analysis config
    (``statistics.p_highly_significant``, ``p_significant``,
    ``p_noteworthy``).

    Parameters
    ----------
    p : float
        P-value.

    Returns
    -------
    str
        Significance marker string.
    """
    stats = get_config()["statistics"]
    if p < stats["p_highly_significant"]:
        return "***"
    if p < stats["p_significant"]:
        return "**"
    if p < stats["p_noteworthy"]:
        return "*"
    return ""


def _sig_class(p: float) -> str:
    """Return a CSS class name for the significance level.

    Classes ``sig-1``, ``sig-2``, ``sig-3`` map to amber, red, and bold-red
    styling defined in :data:`CSS`.  Thresholds are read from the
    centralised analysis config.

    Parameters
    ----------
    p : float
        P-value.

    Returns
    -------
    str
        CSS class name, or empty string if not significant.
    """
    stats = get_config()["statistics"]
    if p < stats["p_highly_significant"]:
        return "sig-3"
    if p < stats["p_significant"]:
        return "sig-2"
    if p < stats["p_noteworthy"]:
        return "sig-1"
    return ""


from report_constants import (  # noqa: F401
    CSS,
    HTML_TEMPLATE,
    REFERENCES,
    REPORT_JS,
    _REF_INDEX,
)



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
    ("manuscript-findings", "Manuscript"),
    ("methods", "Methods"),
    ("cohort", "Cohort"),
    ("patient-flow", "Flow"),
    ("data-quality", "Data Quality"),
    ("data-completeness", "Completeness"),
    ("hypothesis", "Hypothesis"),
    ("graph-overview", "Graphs"),
    ("graph-issues", "Issues"),
    ("stats-by-type", "By Type"),
    ("significance", "Statistics"),
    ("effect-sizes", "Effect Sizes"),
    ("mult-comp", "Corrections"),
    ("cross-dwi", "Cross-DWI"),
    ("fdr-global", "FDR Global"),
    ("correlations", "Correlations"),
    ("treatment", "Treatment"),
    ("predictive", "Predictive"),
    ("feature-overlap", "Features"),
    ("model-diag", "Diagnostics"),
    ("power", "Power"),
    ("sensitivity", "Sensitivity"),
    ("supplemental", "Supplemental"),
    ("limitations", "Limitations"),
    ("conclusions", "Conclusions"),
    ("reporting-checklist", "Checklist"),
    ("table-index", "Tables"),
    ("figure-index", "Figures"),
    ("results-draft", "Results Draft"),
    ("figure-gallery", "Gallery"),
    ("journal-guide", "Journal"),
    ("data-availability", "Data"),
    ("references", "References"),
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

    Thresholds are read from the centralised analysis config
    (``statistics.effect_size_medium``, ``effect_size_large``).

    Parameters
    ----------
    d : float
        Absolute effect size (Cohen's d or similar).

    Returns
    -------
    str
        CSS class name.
    """
    stats = get_config()["statistics"]
    d = abs(d)
    if d >= stats["effect_size_large"]:
        return "effect-lg"
    if d >= stats["effect_size_medium"]:
        return "effect-md"
    return "effect-sm"


def _effect_size_label(d: float) -> str:
    """Return a human-readable label for an effect size magnitude.

    Thresholds are read from the centralised analysis config.

    Parameters
    ----------
    d : float
        Absolute effect size.

    Returns
    -------
    str
        One of ``"Large"``, ``"Medium"``, or ``"Small"``.
    """
    stats = get_config()["statistics"]
    d = abs(d)
    if d >= stats["effect_size_large"]:
        return "Large"
    if d >= stats["effect_size_medium"]:
        return "Medium"
    return "Small"



def _cite(*keys: str) -> str:
    """Return superscript citation(s) for the given reference key(s).

    Parameters
    ----------
    *keys : str
        One or more reference keys (e.g. ``"wilcoxon"``, ``"bh_fdr"``).

    Returns
    -------
    str
        HTML superscript citation, e.g. ``<sup class="cite">[1,2]</sup>``.
    """
    nums = sorted(_REF_INDEX[k] for k in keys if k in _REF_INDEX)
    if not nums:
        return ""
    return f'<sup class="cite">[{",".join(str(n) for n in nums)}]</sup>'


def _copy_button(target_id: str = "", onclick: str = "") -> str:
    """Return an HTML copy-to-clipboard button.

    Parameters
    ----------
    target_id : str, optional
        Element ID to copy from (uses ``copyById``).
    onclick : str, optional
        Custom onclick handler (overrides target_id).

    Returns
    -------
    str
        HTML ``<button>`` element.
    """
    if onclick:
        return f'<button class="copy-btn" onclick="{onclick}">Copy</button>'
    if target_id:
        return f'<button class="copy-btn" onclick="copyById(\'{target_id}\')">Copy</button>'
    return '<button class="copy-btn" onclick="copyText(this)">Copy</button>'


def _manuscript_sentence(text: str) -> str:
    """Return a copyable manuscript-ready sentence block.

    Parameters
    ----------
    text : str
        Publication-ready sentence text.

    Returns
    -------
    str
        HTML block with the sentence and a copy button.
    """
    return (
        f'<div class="manuscript-sentence" data-copy="{_html.escape(text)}">'
        f'{_html.escape(text)}'
        f'<button class="copy-btn" onclick="copyText(this)">Copy</button>'
        f'</div>'
    )


def _references_section() -> list[str]:
    """Build the References section HTML with BibTeX export.

    Returns
    -------
    list[str]
        HTML chunks for the references section.
    """
    h: list[str] = []
    h.append(_h2("References", "references"))
    h.append('<ol class="ref-list">')
    for ref in REFERENCES:
        h.append(f"<li>{_html.escape(ref['text'])}</li>")
    h.append("</ol>")

    # BibTeX export block
    all_bibtex = "\n\n".join(
        ref.get("bibtex", "") for ref in REFERENCES if ref.get("bibtex")
    )
    if all_bibtex:
        escaped_bibtex = _html.escape(all_bibtex)
        h.append('<details><summary><strong>Export BibTeX</strong></summary>')
        h.append(f'<div id="bibtex-block" data-copy="{escaped_bibtex}">')
        h.append(f'<button class="copy-btn" onclick="copyById(\'bibtex-block\')">Copy All BibTeX</button>')
        h.append(f'<pre style="background:#f8fafc;border:1px solid #e2e8f0;'
                 f'border-radius:6px;padding:1rem;font-size:0.82rem;'
                 f'overflow-x:auto;margin-top:0.5rem;white-space:pre-wrap">'
                 f'{escaped_bibtex}</pre>')
        h.append('</div></details>')

    return h


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
    # Keywords must match _trend_tag() to ensure consistent classification.
    increasers = sum(1 for x in trend_list
                     if "increas" in x or "higher" in x or "up" in x or "rising" in x)
    decreasers = sum(1 for x in trend_list
                     if "decreas" in x or "lower" in x or "down" in x
                     or "falling" in x or "drop" in x)
    if increasers > decreasers: return "increasing"
    if decreasers > increasers: return "decreasing"
    return "stable"


