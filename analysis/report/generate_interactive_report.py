#!/usr/bin/env python3
"""Generate an interactive HTML analysis report with client-side filtering.

Produces a standalone HTML file with embedded JavaScript that allows the
reader to:

- Filter by patient, DWI type, or core method via a sidebar
- Toggle between visualisations using tabbed panels
- Drill into individual scan results with expandable patient cards
- Compare metrics interactively with Chart.js bar/line charts
- Sort any table column by clicking its header
- Search across all content with a live search box

The report embeds all data as a JSON blob in a ``<script>`` tag and uses
Chart.js (loaded from CDN) for charting.  No server is required.

Usage:
    python generate_interactive_report.py [saved_files_path]
"""

from __future__ import annotations

import html as _html
import json
import sys
from datetime import datetime
from pathlib import Path

from tqdm import tqdm  # type: ignore

# Ensure analysis/ root is on sys.path so sibling packages are importable.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import (  # type: ignore
    DWI_TYPES,
    group_by_graph_name,
    load_graph_csv,
    resolve_folder,
    setup_utf8_stdout,
)
from parsers.parse_csv_results import parse_all_csvs  # type: ignore
from parsers.parse_log_metrics import parse_all_logs  # type: ignore
from report.interactive_constants import INTERACTIVE_CSS, INTERACTIVE_JS  # type: ignore

# Re-export helpers and section builders for backward compatibility (tests and
# other modules import these names from generate_interactive_report).
from report.generate_interactive_report_helpers import (  # noqa: F401  # type: ignore
    _dwi_badge,
    _esc,
    _extract_core_methods,
    _extract_patients,
    _sig_class,
    _trend_tag,
)
from report.generate_interactive_report_sections import (  # noqa: F401  # type: ignore
    _build_core_comparison,
    _build_dosimetry_section,
    _build_graph_explorer,
    _build_overview_section,
    _build_patient_explorer,
    _build_significance_table,
    _build_visualizations_section,
)

setup_utf8_stdout()




# ── Main report assembly ─────────────────────────────────────────────────────

def generate_interactive_report(folder: Path) -> str:
    """Build the full interactive HTML report from all available data sources.

    Parameters
    ----------
    folder : Path
        Path to the ``saved_files_*`` output folder.

    Returns
    -------
    str
        Complete HTML document as a string.
    """
    # Extract timestamp from folder name.
    raw_ts = folder.name.replace("saved_files_", "")
    try:
        _dt = datetime.strptime(raw_ts, "%Y%m%d_%H%M%S")
        timestamp = _dt.strftime(f"%B {_dt.day}, %Y at %H:%M:%S")
    except ValueError:
        timestamp = raw_ts
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # ── Load data sources ──
    report_bar = tqdm(total=10, desc="Building interactive report", unit="step",
                      bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}] {postfix}")

    report_bar.set_postfix_str("loading CSVs", refresh=True)
    try:
        csv_data = parse_all_csvs(folder)
    except Exception as e:
        print(f"  ⚠️  CSV parsing failed ({type(e).__name__}): {e}")
        csv_data = None
    report_bar.update(1)

    report_bar.set_postfix_str("loading logs", refresh=True)
    try:
        log_data = parse_all_logs(folder)
    except Exception as e:
        print(f"  ⚠️  Log parsing failed ({type(e).__name__}): {e}")
        log_data = None
    report_bar.update(1)

    report_bar.set_postfix_str("loading MAT data", refresh=True)
    mat_data: dict = {}
    for dt in DWI_TYPES:
        mat_path = folder / f"parsed_mat_metrics_{dt}.json"
        if mat_path.exists():
            try:
                with open(mat_path, "r", encoding="utf-8") as f:
                    mat_data[dt] = json.load(f)
            except (json.JSONDecodeError, OSError) as e:
                print(f"  ⚠️  Failed to load {mat_path.name} ({type(e).__name__}): {e}")
    report_bar.update(1)

    report_bar.set_postfix_str("loading graphs", refresh=True)
    rows = load_graph_csv(folder)
    groups = group_by_graph_name(rows) if rows else {}
    report_bar.update(1)

    # Detect DWI types present
    dwi_types_present: list[str] = []
    if log_data:
        for d in DWI_TYPES:
            if d in log_data:
                dwi_types_present.append(d)
    if not dwi_types_present:
        dwi_types_present = [d for d in DWI_TYPES if (folder / d).is_dir()]

    # Extract structured data
    report_bar.set_postfix_str("extracting patients", refresh=True)
    patients = _extract_patients(log_data, mat_data)
    core_methods = _extract_core_methods(mat_data)
    report_bar.update(1)

    # Prepare JSON data blob for JavaScript
    report_data = {
        "patients": patients,
        "dwi_types": dwi_types_present,
        "core_methods": core_methods,
        "graph_rows": [
            {k: v for k, v in r.items()
             if isinstance(v, (str, int, float, bool, type(None)))}
            for r in (rows or [])
        ],
    }

    # ── Assemble HTML ──
    h: list[str] = []
    h.append("<!DOCTYPE html>")
    h.append('<html lang="en">')
    h.append("<head>")
    h.append('<meta charset="utf-8">')
    h.append('<meta name="viewport" content="width=device-width, initial-scale=1">')
    h.append(f"<title>Interactive Report \u2014 {_esc(timestamp)}</title>")
    h.append(f"<style>{INTERACTIVE_CSS}</style>")
    # Chart.js from CDN
    h.append('<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.7/dist/chart.umd.min.js"></script>')
    # Inject data as JSON
    h.append("<script>")
    h.append("window.REPORT_DATA = ")
    h.append(json.dumps(report_data, default=str, ensure_ascii=False))
    h.append(";")
    h.append("</script>")
    h.append(f"<script>{INTERACTIVE_JS}</script>")
    h.append("</head>")
    h.append("<body>")

    # ── Sidebar ──
    h.append('<aside class="sidebar">')
    h.append("<h2>Filters</h2>")

    # DWI type filter
    h.append("<h3>DWI Type</h3>")
    h.append('<div id="dwi-filter-checks"></div>')

    # Patient filter
    h.append("<h3>Patient</h3>")
    h.append('<select id="patient-select" multiple size="5">')
    h.append('<option value="">All patients</option>')
    h.append("</select>")

    # Core method filter
    h.append('<div id="core-method-section">')
    h.append("<h3>Core Method</h3>")
    h.append('<div id="core-method-checks"></div>')
    h.append("</div>")

    # Navigation links
    h.append('<div class="nav-links">')
    h.append("<h3>Sections</h3>")
    nav_items = [
        ("overview", "Overview"),
        ("visualizations", "Visualizations"),
        ("patients", "Patient Explorer"),
        ("significance", "Significance"),
        ("graphs", "Graph Analysis"),
        ("core-methods", "Core Methods"),
        ("dosimetry", "Dosimetry"),
    ]
    for anchor, label in nav_items:
        h.append(f'<a href="#{anchor}">{label}</a>')
    h.append("</div>")

    h.append("</aside>")

    # ── Main content ──
    h.append('<main class="main-content">')
    h.append(f"<h1>Interactive Analysis Report \u2014 {_esc(timestamp)}</h1>")
    h.append(f'<div class="meta">Generated: {now} | DWI types: '
             f'{", ".join(dwi_types_present) or "None detected"}</div>')

    # Overview section
    report_bar.set_postfix_str("overview", refresh=True)
    h.extend(_build_overview_section(log_data, dwi_types_present, rows or [], csv_data, mat_data, patients))
    report_bar.update(1)

    # Visualizations section
    report_bar.set_postfix_str("visualizations", refresh=True)
    h.extend(_build_visualizations_section(dwi_types_present))
    report_bar.update(1)

    # Patient explorer
    report_bar.set_postfix_str("patient explorer", refresh=True)
    h.extend(_build_patient_explorer(patients, dwi_types_present))
    report_bar.update(1)

    # Significance tables
    report_bar.set_postfix_str("significance", refresh=True)
    h.extend(_build_significance_table(csv_data, log_data, dwi_types_present))

    # Graph analysis
    h.extend(_build_graph_explorer(rows or [], groups))

    # Core method comparison
    h.extend(_build_core_comparison(mat_data))

    # Dosimetry
    h.extend(_build_dosimetry_section(mat_data, dwi_types_present))
    report_bar.update(1)

    # Footer
    h.append(f"<footer>Interactive report generated by pancData3 on {now}</footer>")
    h.append("</main>")
    h.append("</body>")
    h.append("</html>")

    report_bar.set_postfix_str("complete", refresh=True)
    report_bar.close()

    return "\n".join(h)


def main():
    """CLI entry point: generate the interactive HTML report."""
    folder = resolve_folder(sys.argv)
    html_report = generate_interactive_report(folder)
    out_path = folder / "interactive_report.html"
    out_path.write_text(html_report, encoding="utf-8")
    print(f"Interactive report written to: {out_path}")
    print(f"  Length: {len(html_report)} characters")


if __name__ == "__main__":
    main()
