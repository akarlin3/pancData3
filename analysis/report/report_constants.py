"""Large constants extracted from report_formatters.py.

Contains the CSS stylesheet, JavaScript, publication references,
and HTML template used by the analysis report.
"""

from __future__ import annotations

# ── CSS stylesheet ────────────────────────────────────────────────────────────
# Embedded directly in the HTML report's <style> tag so the report is a
# single self-contained file with no external dependencies.  Uses CSS custom
# properties (variables) for theming consistency.
CSS = """\
:root {
    --bg: #ffffff; --fg: #1a1a2e; --muted: #5a6a7e;
    --accent: #1e3569; --accent-light: #e8eff7;
    --border: #cdd5e0; --row-alt: #f5f7fa;
    --green: #166534; --red: #991b1b; --amber: #92400e;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    font-size: 1rem; color: var(--fg); background: var(--bg);
    line-height: 1.7; max-width: 1100px; margin: 0 auto; padding: 2rem 1.5rem;
}
h1, h2, h3, h4 {
    font-family: Georgia, 'Times New Roman', serif;
    font-weight: 700;
}
h1 { font-size: 1.75rem; color: var(--fg); border-bottom: 2px solid var(--accent); padding-bottom: 0.6rem; margin-bottom: 1.5rem; }
h2 { font-size: 1.35rem; color: var(--accent); margin: 2.5rem 0 1rem; border-bottom: 1px solid var(--border); padding-bottom: 0.4rem; scroll-margin-top: 3.5rem; }
h3 { font-size: 1.1rem; color: #1e293b; margin: 1.5rem 0 0.65rem; }
h4 { font-size: 1.0rem; font-weight: 600; margin: 1.25rem 0 0.5rem; color: #374151; }
p, ul { margin-bottom: 0.95rem; }
ul { padding-left: 1.5rem; }
.meta { color: var(--muted); font-size: 0.9rem; margin-bottom: 1.5rem; }
.meta span { display: inline-block; margin-right: 1.5rem; }
code { background: var(--row-alt); border: 1px solid var(--border); border-radius: 3px; padding: 0.1em 0.35em; font-size: 0.88em; }
table { width: 100%; border-collapse: collapse; margin: 0.75rem 0 1.5rem; font-size: 0.95rem; }
th { background: var(--accent); color: #fff; text-align: left; padding: 0.6rem 0.9rem; font-weight: 600; font-family: Georgia, 'Times New Roman', serif; letter-spacing: 0.01em; }
td { padding: 0.6rem 0.9rem; border-bottom: 1px solid var(--border); vertical-align: top; }
tr:nth-child(even) td { background: var(--row-alt); }
tr:hover td { background: var(--accent-light); }
.sig-1 { color: var(--amber); font-weight: 600; }
.sig-2 { color: var(--red); font-weight: 600; }
.sig-3 { color: var(--red); font-weight: 700; }
.agree { color: var(--green); font-weight: 600; }
.differ { color: var(--red); font-weight: 700; }
.badge { display: inline-block; padding: 0.15em 0.55em; border-radius: 2px; font-size: 0.8rem; font-weight: 600; border: 1px solid; }
.badge-standard { background: #eef3fb; color: #1e3569; border-color: #93acd4; }
.badge-dncnn { background: #f0fdf4; color: #14532d; border-color: #86efac; }
.badge-ivimnet { background: #fefce8; color: #713f12; border-color: #fcd34d; }
.badge-root { background: #f8f9fb; color: #374151; border-color: #9ca3af; }
.summary-box { background: #f0f4fb; border-left: 4px solid var(--accent); padding: 1rem 1.25rem; border-radius: 0 4px 4px 0; margin: 1rem 0; }
.warn-box { background: #fdf8f0; border-left: 4px solid #b45309; padding: 0.75rem 1rem; border-radius: 0 4px 4px 0; margin: 0.75rem 0; font-size: 0.9rem; }
.info-box { background: #f0f7f2; border-left: 4px solid #166534; padding: 0.75rem 1rem; border-radius: 0 4px 4px 0; margin: 0.75rem 0; font-size: 0.9rem; }
footer { margin-top: 3rem; padding-top: 1rem; border-top: 1px solid var(--border); color: var(--muted); font-size: 0.85rem; }
nav.toc {
    position: sticky; top: 0; background: var(--bg);
    border-bottom: 2px solid var(--border); padding: 0.35rem 0;
    margin-bottom: 1.75rem; z-index: 100; font-size: 0.82rem;
}
.toc-row { display: flex; flex-wrap: wrap; align-items: baseline; padding: 0.1rem 0; }
.toc-group-label { color: var(--muted); font-size: 0.72rem; font-weight: 700; text-transform: uppercase; letter-spacing: 0.05em; white-space: nowrap; margin-right: 0.5rem; }
nav.toc a { color: var(--accent); text-decoration: none; margin-right: 0.9rem; white-space: nowrap; }
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
.stat-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(220px, 1fr)); gap: 1rem; margin: 1rem 0 1.5rem; }
.stat-card { background: var(--row-alt); border: 1px solid var(--border); border-radius: 6px; padding: 1rem 1.25rem; }
.stat-card .label { font-size: 0.78rem; color: var(--muted); text-transform: uppercase; letter-spacing: 0.04em; }
.stat-card .value { font-size: 1.35rem; font-weight: 700; color: var(--accent); }
.stat-card .sub { font-size: 0.8rem; color: var(--muted); }
.abstract-box { background: #f7f9fc; border: 1px solid var(--border); border-radius: 4px; padding: 1.25rem 1.75rem; margin: 1rem 0; }
.abstract-box h4 { color: var(--fg); margin: 0.9rem 0 0.35rem; font-size: 0.9rem; text-transform: uppercase; letter-spacing: 0.06em; font-family: Georgia, 'Times New Roman', serif; }
.abstract-box h4:first-child { margin-top: 0; }
.abstract-box p { margin-bottom: 0.55rem; font-size: 0.94rem; }
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
.conclusion-box { background: #f0f4fb; border-left: 4px solid var(--accent); padding: 1rem 1.25rem; border-radius: 0 4px 4px 0; margin: 0.75rem 0; }
.conclusion-box li { margin-bottom: 0.35rem; }
.diag-box { background: #fffbeb; border: 1px solid #fde68a; border-radius: 6px; padding: 0.75rem 1rem; margin: 0.5rem 0; font-size: 0.9rem; }
.pub-meta { background: #f8fafc; border: 1px solid var(--border); border-radius: 8px; padding: 1.25rem 1.5rem; margin: 1rem 0; }
.pub-meta h4 { color: var(--accent); margin: 0.5rem 0 0.2rem; font-size: 0.95rem; }
.pub-meta p { margin-bottom: 0.4rem; font-size: 0.92rem; color: var(--muted); }
.pub-meta .placeholder { border-bottom: 1px dashed var(--border); display: inline-block; min-width: 200px; color: var(--muted); font-style: italic; }
.ref-list { counter-reset: ref-counter; list-style: none; padding-left: 0; }
.ref-list li { counter-increment: ref-counter; padding: 0.4rem 0 0.4rem 2.5rem; text-indent: -2.5rem; font-size: 0.9rem; line-height: 1.6; border-bottom: 1px solid #edf0f5; }
.ref-list li::before { content: counter(ref-counter) ". "; font-weight: 700; color: var(--fg); }
.cite { color: var(--accent); font-size: 0.8rem; vertical-align: super; font-weight: 600; cursor: help; }
caption { caption-side: top; text-align: left; font-size: 0.9rem; color: #374151; padding: 0.4rem 0; font-weight: 600; }
.copy-btn {
    display: inline-block; padding: 0.25em 0.6em; border: 1px solid var(--border);
    border-radius: 4px; background: var(--row-alt); color: var(--accent);
    font-size: 0.78rem; cursor: pointer; margin-left: 0.5rem; vertical-align: middle;
    transition: background 0.15s, border-color 0.15s;
}
.copy-btn:hover { background: var(--accent-light); border-color: var(--accent); }
.copy-btn.copied { background: #dcfce7; border-color: var(--green); color: var(--green); }
.manuscript-sentence { background: #f5f8fd; border: 1px solid #c3d0e8; border-radius: 3px; padding: 0.55rem 0.85rem; margin: 0.4rem 0; font-size: 0.9rem; line-height: 1.6; position: relative; }
.manuscript-sentence .copy-btn { position: absolute; top: 0.35rem; right: 0.5rem; }
.checklist-table td:first-child { font-weight: 600; white-space: nowrap; }
.table-wide { font-size: 0.88rem; table-layout: fixed; }
.table-wide td, .table-wide th { padding: 0.4rem 0.6rem; overflow-wrap: break-word; word-break: break-word; }
.table-wide .axis-info { font-size: 0.8rem; }
.graph-card { background: #fafbfc; border: 1px solid var(--border); border-radius: 6px; padding: 0.75rem 1rem; margin: 0.5rem 0; font-size: 0.88rem; line-height: 1.55; }
.graph-card-header { display: flex; align-items: center; gap: 0.5rem; margin-bottom: 0.4rem; font-weight: 600; }
.graph-card-header .badge { flex-shrink: 0; }
.graph-card-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 0.25rem 1.5rem; margin-top: 0.4rem; }
.graph-card-grid dt { font-weight: 600; color: var(--muted); font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.03em; margin-top: 0.3rem; }
.graph-card-grid dd { margin: 0; }
.graph-card-full { grid-column: 1 / -1; }
.table-compact td, .table-compact th { padding: 0.3rem 0.5rem; font-size: 0.85rem; }
.checklist-done { color: var(--green); }
.checklist-partial { color: var(--amber); }
.checklist-na { color: var(--muted); font-style: italic; }
.toc-list { column-count: 2; column-gap: 2rem; font-size: 0.88rem; margin: 0.5rem 0 1rem; }
.toc-list li { break-inside: avoid; margin-bottom: 0.2rem; }
.toc-list a { color: var(--accent); text-decoration: none; }
.toc-list a:hover { text-decoration: underline; }
.print-toc { background: #f8fafc; border: 1px solid var(--border); border-radius: 6px; padding: 1.5rem; margin: 1.5rem 0; column-count: 2; column-gap: 2rem; }
.print-toc h2 { column-span: all; margin-top: 0; border: none; padding: 0; }
.print-toc-group { break-inside: avoid; margin-bottom: 1rem; }
.print-toc-group-label { display: block; color: var(--accent); font-size: 0.78rem; font-weight: 700; text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 0.3rem; }
.print-toc-list { list-style: none; padding: 0; margin: 0; }
.print-toc-list li { margin-bottom: 0.15rem; font-size: 0.9rem; }
.print-toc-list a { color: var(--fg); text-decoration: none; }
.print-toc-list a:hover { color: var(--accent); text-decoration: underline; }
.word-count { display: inline-block; font-size: 0.75rem; color: var(--muted); margin-left: 0.5rem; font-weight: 400; }
.cover-page { display: none; }
.part-break { display: none; }
@page { size: A4; margin: 1.5cm 1.2cm 1.8cm; }
@page { @bottom-left { content: "Generated by pancData3 analysis pipeline"; font-size: 7pt; color: #aaa; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; } }
@page { @bottom-center { content: "Page " counter(page) " of " counter(pages); font-size: 8pt; color: #888; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; } }
@page :first { @bottom-left { content: none; } @bottom-center { content: none; } }
@media print {
    body { max-width: none; padding: 0; font-size: 10.5pt; line-height: 1.55; }
    nav.toc { display: none; }
    .cover-page { display: flex; flex-direction: column; justify-content: center; align-items: center; height: 100vh; page-break-after: always; break-after: page; text-align: center; }
    .cover-page h1 { font-size: 22pt; border: none; margin-bottom: 1.5cm; color: var(--accent); }
    .cover-page .cover-meta { font-size: 10pt; color: #555; line-height: 2; margin-top: 1cm; }
    .cover-page .cover-rule { width: 8cm; border-top: 2pt solid var(--accent); margin: 1cm auto; }
    .part-break { display: block; page-break-before: always; break-before: page; page-break-after: avoid; break-after: avoid; margin: 0 0 1.5em; padding: 0.6em 0; border-bottom: 2pt solid var(--accent); color: var(--accent); font-family: Georgia, 'Times New Roman', serif; font-size: 14pt; font-weight: 700; letter-spacing: 0.04em; }
    h1 { font-size: 15pt; }
    h2 { font-size: 12.5pt; page-break-after: avoid; margin-top: 1.5em; }
    h3 { font-size: 11pt; page-break-after: avoid; }
    table { font-size: 9pt; table-layout: auto; width: 100%; overflow-wrap: break-word; word-break: break-word; }
    td, th { overflow-wrap: break-word; word-break: break-word; padding: 0.3rem 0.45rem; }
    th { background: #333 !important; color: #fff !important; -webkit-print-color-adjust: exact; print-color-adjust: exact; }
    tr:nth-child(even) td { background: #f5f5f5 !important; -webkit-print-color-adjust: exact; print-color-adjust: exact; }
    thead { display: table-header-group; }
    tr { page-break-inside: avoid; }
    .stat-grid { grid-template-columns: repeat(4, 1fr); }
    .stat-card { border: 1px solid #ccc; }
    .summary-box, .warn-box, .info-box, .conclusion-box, .diag-box, .methods-box, .abstract-box, .pub-meta {
        -webkit-print-color-adjust: exact; print-color-adjust: exact;
    }
    .forest-row { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
    .badge { border: 1px solid #999; -webkit-print-color-adjust: exact; print-color-adjust: exact; }
    a { color: var(--accent); text-decoration: none; }
    footer { font-size: 8pt; }
    details { display: block; }
    details[open] > div, details[open] > p, details[open] > pre { display: block; }
    details > summary { list-style: none; }
    details > summary::before { content: ""; }
    .copy-btn { display: none; }
    .word-count { display: none; }
    .print-toc { column-count: 2; page-break-after: always; break-after: page; border: 1px solid #ccc; }
    img { max-width: 100%; height: auto; }
    .table-wide { font-size: 8pt; }
    .table-wide td, .table-wide th { padding: 0.2rem 0.35rem; }
    .graph-card { break-inside: avoid; border: 1px solid #ccc; padding: 0.5rem 0.75rem; margin: 0.35rem 0; font-size: 9pt; }
    .graph-card-grid { display: block; }
    .graph-card-grid dt { display: inline; font-weight: bold; }
    .graph-card-grid dd { display: inline; margin: 0; }
    .graph-card-grid dd::after { content: ""; display: block; }
}
"""

# ── JavaScript for interactive features (copy-to-clipboard) ──────────────────
REPORT_JS = """\
function copyText(el) {
    var target = el.closest('[data-copy]') || el.parentElement;
    var text = target.getAttribute('data-copy') || target.innerText;
    // Strip the button text from the copy content
    text = text.replace(/Copy(?: to clipboard)?\\s*$/m, '').trim();
    navigator.clipboard.writeText(text).then(function() {
        el.textContent = 'Copied!';
        el.classList.add('copied');
        setTimeout(function() { el.textContent = 'Copy'; el.classList.remove('copied'); }, 1500);
    });
}
function copyById(id) {
    var el = document.getElementById(id);
    if (!el) return;
    var text = el.getAttribute('data-copy') || el.innerText;
    navigator.clipboard.writeText(text).then(function() {
        var btn = el.querySelector('.copy-btn');
        if (btn) { btn.textContent = 'Copied!'; btn.classList.add('copied');
            setTimeout(function() { btn.textContent = 'Copy'; btn.classList.remove('copied'); }, 1500);
        }
    });
}
"""

# ── Publication references ────────────────────────────────────────────────────
# Numbered references used in the Methods, Diagnostics, and Effect Size sections
# to support a journal-ready report.

REFERENCES: list[dict[str, str]] = [
    {"key": "wilcoxon",
     "text": "Mann HB, Whitney DR. On a test of whether one of two random "
             "variables is stochastically larger than the other. Ann Math "
             "Stat. 1947;18(1):50\u201360.",
     "bibtex": "@article{mann1947test,\n"
               "  author  = {Mann, H. B. and Whitney, D. R.},\n"
               "  title   = {On a test of whether one of two random variables is "
               "stochastically larger than the other},\n"
               "  journal = {Annals of Mathematical Statistics},\n"
               "  year    = {1947},\n"
               "  volume  = {18},\n"
               "  number  = {1},\n"
               "  pages   = {50--60}\n"
               "}"},
    {"key": "bh_fdr",
     "text": "Benjamini Y, Hochberg Y. Controlling the false discovery rate: "
             "a practical and powerful approach to multiple testing. J R Stat "
             "Soc Series B. 1995;57(1):289\u2013300.",
     "bibtex": "@article{benjamini1995controlling,\n"
               "  author  = {Benjamini, Yoav and Hochberg, Yosef},\n"
               "  title   = {Controlling the false discovery rate: a practical and "
               "powerful approach to multiple testing},\n"
               "  journal = {Journal of the Royal Statistical Society: Series B},\n"
               "  year    = {1995},\n"
               "  volume  = {57},\n"
               "  number  = {1},\n"
               "  pages   = {289--300}\n"
               "}"},
    {"key": "cox_ph",
     "text": "Cox DR. Regression models and life-tables. J R Stat Soc Series B. "
             "1972;34(2):187\u2013220.",
     "bibtex": "@article{cox1972regression,\n"
               "  author  = {Cox, David R.},\n"
               "  title   = {Regression models and life-tables},\n"
               "  journal = {Journal of the Royal Statistical Society: Series B},\n"
               "  year    = {1972},\n"
               "  volume  = {34},\n"
               "  number  = {2},\n"
               "  pages   = {187--220}\n"
               "}"},
    {"key": "elastic_net",
     "text": "Zou H, Hastie T. Regularization and variable selection via the "
             "elastic net. J R Stat Soc Series B. 2005;67(2):301\u2013320.",
     "bibtex": "@article{zou2005regularization,\n"
               "  author  = {Zou, Hui and Hastie, Trevor},\n"
               "  title   = {Regularization and variable selection via the elastic net},\n"
               "  journal = {Journal of the Royal Statistical Society: Series B},\n"
               "  year    = {2005},\n"
               "  volume  = {67},\n"
               "  number  = {2},\n"
               "  pages   = {301--320}\n"
               "}"},
    {"key": "hosmer_lemeshow",
     "text": "Hosmer DW, Lemeshow S. Applied Logistic Regression. 2nd ed. "
             "New York: Wiley; 2000.",
     "bibtex": "@book{hosmer2000applied,\n"
               "  author    = {Hosmer, David W. and Lemeshow, Stanley},\n"
               "  title     = {Applied Logistic Regression},\n"
               "  edition   = {2nd},\n"
               "  publisher = {Wiley},\n"
               "  address   = {New York},\n"
               "  year      = {2000}\n"
               "}"},
    {"key": "ipcw",
     "text": "Robins JM, Finkelstein DM. Correcting for noncompliance and "
             "dependent censoring in an AIDS clinical trial with inverse "
             "probability of censoring weighted (IPCW) log-rank tests. "
             "Biometrics. 2000;56(3):779\u2013788.",
     "bibtex": "@article{robins2000correcting,\n"
               "  author  = {Robins, James M. and Finkelstein, Dianne M.},\n"
               "  title   = {Correcting for noncompliance and dependent censoring in an "
               "AIDS clinical trial with inverse probability of censoring weighted "
               "(IPCW) log-rank tests},\n"
               "  journal = {Biometrics},\n"
               "  year    = {2000},\n"
               "  volume  = {56},\n"
               "  number  = {3},\n"
               "  pages   = {779--788}\n"
               "}"},
    {"key": "firth",
     "text": "Firth D. Bias reduction of maximum likelihood estimates. "
             "Biometrika. 1993;80(1):27\u201338.",
     "bibtex": "@article{firth1993bias,\n"
               "  author  = {Firth, David},\n"
               "  title   = {Bias reduction of maximum likelihood estimates},\n"
               "  journal = {Biometrika},\n"
               "  year    = {1993},\n"
               "  volume  = {80},\n"
               "  number  = {1},\n"
               "  pages   = {27--38}\n"
               "}"},
    {"key": "ivim",
     "text": "Le Bihan D, Breton E, Lallemand D, Aubin ML, Vignaud J, "
             "Laval-Jeantet M. Separation of diffusion and perfusion in "
             "intravoxel incoherent motion MR imaging. Radiology. "
             "1988;168(2):497\u2013505.",
     "bibtex": "@article{lebihan1988separation,\n"
               "  author  = {Le Bihan, Denis and Breton, Eric and Lallemand, "
               "Denis and Aubin, Marie-Louise and Vignaud, Jean and "
               "Laval-Jeantet, Michel},\n"
               "  title   = {Separation of diffusion and perfusion in intravoxel "
               "incoherent motion MR imaging},\n"
               "  journal = {Radiology},\n"
               "  year    = {1988},\n"
               "  volume  = {168},\n"
               "  number  = {2},\n"
               "  pages   = {497--505}\n"
               "}"},
    {"key": "cohen_d",
     "text": "Cohen J. Statistical Power Analysis for the Behavioral Sciences. "
             "2nd ed. Hillsdale, NJ: Lawrence Erlbaum Associates; 1988.",
     "bibtex": "@book{cohen1988statistical,\n"
               "  author    = {Cohen, Jacob},\n"
               "  title     = {Statistical Power Analysis for the Behavioral Sciences},\n"
               "  edition   = {2nd},\n"
               "  publisher = {Lawrence Erlbaum Associates},\n"
               "  address   = {Hillsdale, NJ},\n"
               "  year      = {1988}\n"
               "}"},
    {"key": "dice",
     "text": "Dice LR. Measures of the amount of ecologic association between "
             "species. Ecology. 1945;26(3):297\u2013302.",
     "bibtex": "@article{dice1945measures,\n"
               "  author  = {Dice, Lee R.},\n"
               "  title   = {Measures of the amount of ecologic association between "
               "species},\n"
               "  journal = {Ecology},\n"
               "  year    = {1945},\n"
               "  volume  = {26},\n"
               "  number  = {3},\n"
               "  pages   = {297--302}\n"
               "}"},
    {"key": "dncnn",
     "text": "Zhang K, Zuo W, Chen Y, Meng D, Zhang L. Beyond a Gaussian "
             "denoiser: residual learning of deep CNN for image denoising. "
             "IEEE Trans Image Process. 2017;26(7):3142\u20133155.",
     "bibtex": "@article{zhang2017beyond,\n"
               "  author  = {Zhang, Kai and Zuo, Wangmeng and Chen, Yunjin and "
               "Meng, Deyu and Zhang, Lei},\n"
               "  title   = {Beyond a Gaussian denoiser: residual learning of deep "
               "CNN for image denoising},\n"
               "  journal = {IEEE Transactions on Image Processing},\n"
               "  year    = {2017},\n"
               "  volume  = {26},\n"
               "  number  = {7},\n"
               "  pages   = {3142--3155}\n"
               "}"},
    {"key": "strobe",
     "text": "von Elm E, Altman DG, Egger M, Pocock SJ, G\u00f8tzsche PC, "
             "Vandenbroucke JP. The Strengthening the Reporting of "
             "Observational Studies in Epidemiology (STROBE) statement: "
             "guidelines for reporting observational studies. Ann Intern Med. "
             "2007;147(8):573\u2013577.",
     "bibtex": "@article{vonelm2007strobe,\n"
               "  author  = {von Elm, Erik and Altman, Douglas G. and Egger, Matthias "
               "and Pocock, Stuart J. and G{\\o}tzsche, Peter C. and "
               "Vandenbroucke, Jan P.},\n"
               "  title   = {The {STROBE} statement: guidelines for reporting "
               "observational studies},\n"
               "  journal = {Annals of Internal Medicine},\n"
               "  year    = {2007},\n"
               "  volume  = {147},\n"
               "  number  = {8},\n"
               "  pages   = {573--577}\n"
               "}"},
    {"key": "remark",
     "text": "McShane LM, Altman DG, Sauerbrei W, Taube SE, Gion M, Clark GM. "
             "REporting recommendations for tumour MARKer prognostic studies "
             "(REMARK). Br J Cancer. 2005;93(4):387\u2013391.",
     "bibtex": "@article{mcshane2005remark,\n"
               "  author  = {McShane, Lisa M. and Altman, Douglas G. and "
               "Sauerbrei, Willi and Taube, Sheila E. and Gion, Massimo and "
               "Clark, Gary M.},\n"
               "  title   = {{REMARK}: {REporting} recommendations for tumour "
               "{MARKer} prognostic studies},\n"
               "  journal = {British Journal of Cancer},\n"
               "  year    = {2005},\n"
               "  volume  = {93},\n"
               "  number  = {4},\n"
               "  pages   = {387--391}\n"
               "}"},
]

_REF_INDEX = {r["key"]: i + 1 for i, r in enumerate(REFERENCES)}


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
    --accent: #1e3569;
    --accent-light: #e8eff7;
    --border: #cdd5e0;
    --table-stripe: #f5f7fa;
    --sig-strong: #991b1b;
    --sig-moderate: #92400e;
    --sig-mild: #713f12;
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
  h1, h2, h3 {{
    font-family: Georgia, 'Times New Roman', serif;
    font-weight: 700;
  }}
  h1 {{
    font-size: 1.75rem;
    color: var(--fg);
    border-bottom: 2px solid var(--accent);
    padding-bottom: 0.5rem;
    margin-bottom: 1rem;
  }}
  h2 {{
    font-size: 1.35rem;
    color: var(--accent);
    border-bottom: 1px solid var(--border);
    padding-bottom: 0.35rem;
    margin-top: 2.5rem;
    margin-bottom: 0.9rem;
  }}
  h3 {{
    font-size: 1.1rem;
    color: #1e293b;
    margin-top: 1.5rem;
    margin-bottom: 0.55rem;
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
