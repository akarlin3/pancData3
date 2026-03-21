"""Constants for the interactive HTML report.

Contains the CSS stylesheet, JavaScript application code, and HTML
skeleton used by :mod:`generate_interactive_report`.  All assets are
embedded directly so the report is a single self-contained file.
"""

from __future__ import annotations

# ── CSS stylesheet for the interactive report ────────────────────────────────
INTERACTIVE_CSS = """\
:root {
    --bg: #ffffff; --fg: #1a1a2e; --muted: #5a6a7e;
    --accent: #1e3569; --accent-light: #e8eff7;
    --border: #cdd5e0; --row-alt: #f5f7fa;
    --green: #166534; --red: #991b1b; --amber: #92400e;
    --sidebar-bg: #f8fafc; --sidebar-w: 260px;
    --card-shadow: 0 1px 3px rgba(0,0,0,0.08);
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    font-size: 0.95rem; color: var(--fg); background: var(--bg);
    line-height: 1.6; display: flex; min-height: 100vh;
}
/* ── Sidebar ── */
.sidebar {
    width: var(--sidebar-w); background: var(--sidebar-bg);
    border-right: 1px solid var(--border); padding: 1.25rem 1rem;
    position: fixed; top: 0; left: 0; bottom: 0; overflow-y: auto;
    z-index: 100; flex-shrink: 0;
}
.sidebar h2 {
    font-family: Georgia, 'Times New Roman', serif;
    font-size: 1.1rem; color: var(--accent); margin-bottom: 1rem;
    padding-bottom: 0.5rem; border-bottom: 2px solid var(--accent);
}
.sidebar h3 {
    font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.05em;
    color: var(--muted); margin: 1rem 0 0.4rem; font-weight: 700;
}
.sidebar label {
    display: flex; align-items: center; gap: 0.4rem;
    padding: 0.2rem 0; font-size: 0.88rem; cursor: pointer;
}
.sidebar label:hover { color: var(--accent); }
.sidebar input[type="checkbox"] { accent-color: var(--accent); }
.sidebar select {
    width: 100%; padding: 0.35rem 0.5rem; border: 1px solid var(--border);
    border-radius: 4px; font-size: 0.88rem; background: #fff;
    margin-top: 0.25rem;
}
.sidebar .filter-count {
    font-size: 0.78rem; color: var(--muted); margin-top: 0.3rem;
}
.sidebar .nav-links { margin-top: 1rem; }
.sidebar .nav-links a {
    display: block; padding: 0.3rem 0.5rem; margin: 0.1rem 0;
    color: var(--accent); text-decoration: none; font-size: 0.85rem;
    border-radius: 3px; transition: background 0.15s;
}
.sidebar .nav-links a:hover, .sidebar .nav-links a.active {
    background: var(--accent-light);
}

/* ── Main content ── */
.main-content {
    margin-left: var(--sidebar-w); flex: 1; padding: 1.5rem 2rem;
    max-width: 1100px;
}
h1, h2, h3, h4 { font-family: Georgia, 'Times New Roman', serif; font-weight: 700; }
h1 { font-size: 1.6rem; color: var(--fg); border-bottom: 2px solid var(--accent);
     padding-bottom: 0.5rem; margin-bottom: 1.25rem; }
h2 { font-size: 1.25rem; color: var(--accent); margin: 2rem 0 0.75rem;
     border-bottom: 1px solid var(--border); padding-bottom: 0.35rem;
     scroll-margin-top: 1rem; }
h3 { font-size: 1.05rem; color: #1e293b; margin: 1.25rem 0 0.5rem; }
h4 { font-size: 0.95rem; font-weight: 600; color: #374151; margin: 1rem 0 0.4rem; }
p, ul { margin-bottom: 0.75rem; }
ul { padding-left: 1.5rem; }
code { background: var(--row-alt); border: 1px solid var(--border);
       border-radius: 3px; padding: 0.1em 0.3em; font-size: 0.88em; }

/* ── Badges ── */
.badge { display: inline-block; padding: 0.12em 0.5em; border-radius: 2px;
         font-size: 0.78rem; font-weight: 600; border: 1px solid; }
.badge-standard { background: #eef3fb; color: #1e3569; border-color: #93acd4; }
.badge-dncnn { background: #f0fdf4; color: #14532d; border-color: #86efac; }
.badge-ivimnet { background: #fefce8; color: #713f12; border-color: #fcd34d; }

/* ── Stat cards ── */
.stat-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
             gap: 0.75rem; margin: 1rem 0; }
.stat-card { background: var(--row-alt); border: 1px solid var(--border);
             border-radius: 6px; padding: 0.85rem 1rem; box-shadow: var(--card-shadow); }
.stat-card .label { font-size: 0.72rem; color: var(--muted); text-transform: uppercase;
                    letter-spacing: 0.04em; }
.stat-card .value { font-size: 1.2rem; font-weight: 700; color: var(--accent); }
.stat-card .sub { font-size: 0.78rem; color: var(--muted); }

/* ── Tables ── */
table { width: 100%; border-collapse: collapse; margin: 0.5rem 0 1rem; font-size: 0.9rem; }
th { background: var(--accent); color: #fff; text-align: left; padding: 0.5rem 0.75rem;
     font-weight: 600; font-family: Georgia, serif; letter-spacing: 0.01em;
     cursor: pointer; user-select: none; position: relative; }
th:hover { background: #2a4a8a; }
th .sort-arrow { font-size: 0.7em; margin-left: 0.3em; opacity: 0.5; }
th.sorted .sort-arrow { opacity: 1; }
td { padding: 0.5rem 0.75rem; border-bottom: 1px solid var(--border); vertical-align: top; }
tr:nth-child(even) td { background: var(--row-alt); }
tr:hover td { background: var(--accent-light); }
tr.hidden { display: none; }
.sig-1 { color: var(--amber); font-weight: 600; }
.sig-2 { color: var(--red); font-weight: 600; }
.sig-3 { color: var(--red); font-weight: 700; }

/* ── Tabs ── */
.tab-bar { display: flex; gap: 0; border-bottom: 2px solid var(--border); margin: 1rem 0 0; }
.tab-btn {
    padding: 0.5rem 1.25rem; border: none; background: none;
    font-size: 0.88rem; color: var(--muted); cursor: pointer;
    border-bottom: 2px solid transparent; margin-bottom: -2px;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    transition: color 0.15s, border-color 0.15s;
}
.tab-btn:hover { color: var(--fg); }
.tab-btn.active { color: var(--accent); border-bottom-color: var(--accent); font-weight: 600; }
.tab-panel { display: none; padding: 1rem 0; }
.tab-panel.active { display: block; }

/* ── Charts ── */
.chart-container { position: relative; width: 100%; max-width: 800px;
                   margin: 0.5rem auto; }
.chart-container canvas { width: 100% !important; }
.chart-row { display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; margin: 1rem 0; }
.chart-row .chart-container { max-width: none; }

/* ── Patient detail cards ── */
.patient-card {
    background: #fafbfc; border: 1px solid var(--border); border-radius: 6px;
    padding: 1rem 1.25rem; margin: 0.75rem 0; box-shadow: var(--card-shadow);
}
.patient-card-header {
    display: flex; align-items: center; justify-content: space-between;
    margin-bottom: 0.75rem; cursor: pointer;
}
.patient-card-header h3 { margin: 0; font-size: 1rem; }
.patient-card-body { display: none; }
.patient-card.expanded .patient-card-body { display: block; }
.patient-card .toggle-icon { font-size: 0.8rem; color: var(--muted);
                             transition: transform 0.2s; }
.patient-card.expanded .toggle-icon { transform: rotate(90deg); }
.metric-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
               gap: 0.5rem; margin: 0.5rem 0; }
.metric-item { padding: 0.4rem 0.6rem; background: #fff; border: 1px solid var(--border);
               border-radius: 4px; }
.metric-item .m-label { font-size: 0.72rem; color: var(--muted); text-transform: uppercase; }
.metric-item .m-value { font-size: 0.95rem; font-weight: 600; color: var(--accent); }

/* ── Comparison mode ── */
.compare-panel { display: grid; grid-template-columns: 1fr 1fr; gap: 1.5rem; margin: 1rem 0; }
.compare-col { background: #fafbfc; border: 1px solid var(--border); border-radius: 6px;
               padding: 1rem; }
.compare-col h4 { margin-top: 0; }

/* ── Search / filter bar ── */
.search-bar { display: flex; gap: 0.5rem; margin: 0.75rem 0; align-items: center; }
.search-bar input {
    flex: 1; padding: 0.4rem 0.75rem; border: 1px solid var(--border);
    border-radius: 4px; font-size: 0.88rem;
}
.search-bar input:focus { outline: none; border-color: var(--accent);
                          box-shadow: 0 0 0 2px var(--accent-light); }

/* ── Trend tags ── */
.trend-tag { display: inline-block; padding: 0.1em 0.4em; border-radius: 3px;
             font-size: 0.78rem; border: 1px solid var(--border); margin: 0.1em; }
.trend-incr { background: #dcfce7; border-color: #86efac; color: #166534; }
.trend-decr { background: #fee2e2; border-color: #fca5a5; color: #991b1b; }
.trend-flat { background: var(--row-alt); }
.trend-nm { background: #f5f3ff; border-color: #c4b5fd; color: #5b21b6; }

/* ── Utility ── */
.meta { color: var(--muted); font-size: 0.85rem; margin-bottom: 1rem; }
.hidden { display: none !important; }
.summary-box { background: #f0f4fb; border-left: 4px solid var(--accent);
               padding: 0.75rem 1rem; border-radius: 0 4px 4px 0; margin: 0.75rem 0; }
.warn-box { background: #fdf8f0; border-left: 4px solid #b45309;
            padding: 0.6rem 0.85rem; border-radius: 0 4px 4px 0; margin: 0.5rem 0;
            font-size: 0.88rem; }
footer { margin-top: 2rem; padding-top: 0.75rem; border-top: 1px solid var(--border);
         color: var(--muted); font-size: 0.82rem; }

/* ── Responsive ── */
@media (max-width: 900px) {
    .sidebar { width: 100%; position: static; border-right: none;
               border-bottom: 1px solid var(--border); }
    .main-content { margin-left: 0; padding: 1rem; }
    .chart-row { grid-template-columns: 1fr; }
    .compare-panel { grid-template-columns: 1fr; }
}
"""

# ── JavaScript application code ──────────────────────────────────────────────
# Provides filtering, sorting, charting, and drill-down interactivity.
# Chart.js is loaded from a CDN <script> tag in the HTML head; the code
# below assumes ``Chart`` is available globally.
INTERACTIVE_JS = """\
(function() {
"use strict";

// ── Data injected by Python (available as window.REPORT_DATA) ──
var DATA = window.REPORT_DATA || {};
var patients = DATA.patients || [];
var dwiTypes = DATA.dwi_types || [];
var coreMethods = DATA.core_methods || [];
var graphRows = DATA.graph_rows || [];
var logData = DATA.log_data || {};
var csvData = DATA.csv_data || {};
var matData = DATA.mat_data || {};

// ── State ──
var state = {
    selectedDwi: dwiTypes.slice(),         // all selected by default
    selectedPatients: [],                   // empty = all
    searchTerm: "",
    activeTab: {},                          // tab-group -> active tab id
    sortCol: null,
    sortAsc: true,
    expandedPatients: new Set(),
};

// ── Utility ──
function $(sel, ctx) { return (ctx || document).querySelector(sel); }
function $$(sel, ctx) { return Array.from((ctx || document).querySelectorAll(sel)); }
function esc(s) {
    var d = document.createElement("div"); d.textContent = s; return d.innerHTML;
}

// ── Sidebar: DWI type checkboxes ──
function initDwiFilters() {
    var container = $("#dwi-filter-checks");
    if (!container) return;
    dwiTypes.forEach(function(dt) {
        var id = "dwi-" + dt;
        var label = document.createElement("label");
        label.innerHTML = '<input type="checkbox" id="' + id + '" checked value="' + esc(dt) + '"> ' + esc(dt);
        container.appendChild(label);
    });
    container.addEventListener("change", function() {
        state.selectedDwi = $$('#dwi-filter-checks input:checked').map(function(cb) { return cb.value; });
        applyFilters();
    });
}

// ── Sidebar: Patient selector ──
function initPatientFilter() {
    var sel = $("#patient-select");
    if (!sel) return;
    patients.forEach(function(p) {
        var opt = document.createElement("option");
        opt.value = p.id;
        opt.textContent = p.id;
        sel.appendChild(opt);
    });
    sel.addEventListener("change", function() {
        var opts = Array.from(sel.selectedOptions).map(function(o) { return o.value; });
        state.selectedPatients = opts.filter(function(v) { return v !== ""; });
        applyFilters();
    });
}

// ── Sidebar: Core method filter ──
function initCoreMethodFilter() {
    var container = $("#core-method-checks");
    if (!container || coreMethods.length === 0) {
        var section = $("#core-method-section");
        if (section) section.classList.add("hidden");
        return;
    }
    coreMethods.forEach(function(cm) {
        var id = "cm-" + cm.replace(/[^a-zA-Z0-9]/g, "_");
        var label = document.createElement("label");
        label.innerHTML = '<input type="checkbox" id="' + id + '" checked value="' + esc(cm) + '"> ' + esc(cm);
        container.appendChild(label);
    });
}

// ── Sidebar: Search ──
function initSearch() {
    var input = $("#search-input");
    if (!input) return;
    input.addEventListener("input", function() {
        state.searchTerm = input.value.toLowerCase().trim();
        applyFilters();
    });
}

// ── Tab switching ──
function initTabs() {
    $$(".tab-bar").forEach(function(bar) {
        var group = bar.dataset.tabGroup;
        var btns = $$(".tab-btn", bar);
        var panels = $$('.tab-panel[data-tab-group="' + group + '"]');
        if (btns.length > 0) {
            btns[0].classList.add("active");
            state.activeTab[group] = btns[0].dataset.tab;
        }
        if (panels.length > 0) panels[0].classList.add("active");
        btns.forEach(function(btn) {
            btn.addEventListener("click", function() {
                btns.forEach(function(b) { b.classList.remove("active"); });
                panels.forEach(function(p) { p.classList.remove("active"); });
                btn.classList.add("active");
                var target = btn.dataset.tab;
                state.activeTab[group] = target;
                var panel = document.querySelector('.tab-panel[data-tab-group="' + group + '"][data-tab="' + target + '"]');
                if (panel) panel.classList.add("active");
            });
        });
    });
}

// ── Sortable tables ──
function initSortableTables() {
    $$("table.sortable").forEach(function(table) {
        var headers = $$("th", table);
        headers.forEach(function(th, colIdx) {
            th.innerHTML += ' <span class="sort-arrow">\\u25B2</span>';
            th.addEventListener("click", function() {
                var tbody = $("tbody", table);
                if (!tbody) return;
                var rows = $$("tr", tbody);
                var asc = th.classList.contains("sorted") && th.dataset.sortDir === "asc";
                var newDir = asc ? "desc" : "asc";
                headers.forEach(function(h) { h.classList.remove("sorted"); h.dataset.sortDir = ""; });
                th.classList.add("sorted");
                th.dataset.sortDir = newDir;
                $(".sort-arrow", th).textContent = newDir === "asc" ? "\\u25B2" : "\\u25BC";
                rows.sort(function(a, b) {
                    var av = (a.cells[colIdx] || {}).textContent || "";
                    var bv = (b.cells[colIdx] || {}).textContent || "";
                    var an = parseFloat(av), bn = parseFloat(bv);
                    if (!isNaN(an) && !isNaN(bn)) return newDir === "asc" ? an - bn : bn - an;
                    return newDir === "asc" ? av.localeCompare(bv) : bv.localeCompare(av);
                });
                rows.forEach(function(r) { tbody.appendChild(r); });
            });
        });
    });
}

// ── Patient cards toggle ──
function initPatientCards() {
    $$(".patient-card-header").forEach(function(header) {
        header.addEventListener("click", function() {
            var card = header.parentElement;
            card.classList.toggle("expanded");
            var pid = card.dataset.patientId;
            if (card.classList.contains("expanded")) state.expandedPatients.add(pid);
            else state.expandedPatients.delete(pid);
        });
    });
}

// ── Apply all filters ──
function applyFilters() {
    // Filter table rows by DWI type
    $$("tr[data-dwi]").forEach(function(row) {
        var dwi = row.dataset.dwi;
        var show = state.selectedDwi.indexOf(dwi) !== -1;
        if (show && state.searchTerm) {
            show = row.textContent.toLowerCase().indexOf(state.searchTerm) !== -1;
        }
        row.classList.toggle("hidden", !show);
    });
    // Filter patient cards
    $$(".patient-card").forEach(function(card) {
        var pid = card.dataset.patientId;
        var dwiShow = !card.dataset.dwi || state.selectedDwi.indexOf(card.dataset.dwi) !== -1;
        var patientShow = state.selectedPatients.length === 0 || state.selectedPatients.indexOf(pid) !== -1;
        var searchShow = !state.searchTerm || card.textContent.toLowerCase().indexOf(state.searchTerm) !== -1;
        card.classList.toggle("hidden", !(dwiShow && patientShow && searchShow));
    });
    // Update filter count
    updateFilterCount();
    // Rebuild charts if needed
    rebuildCharts();
}

function updateFilterCount() {
    var el = $("#filter-count");
    if (!el) return;
    var visible = $$(".patient-card:not(.hidden)").length;
    var total = $$(".patient-card").length;
    el.textContent = visible + " / " + total + " patients shown";
}

// ── Chart rendering ──
var chartInstances = {};

function destroyChart(id) {
    if (chartInstances[id]) {
        chartInstances[id].destroy();
        delete chartInstances[id];
    }
}

function rebuildCharts() {
    // Longitudinal overview chart
    buildLongitudinalChart();
    // Metric distribution chart
    buildDistributionChart();
    // DWI comparison chart
    buildDwiComparisonChart();
}

function getChartColors(n) {
    var palette = [
        "rgba(30, 53, 105, 0.8)", "rgba(22, 101, 52, 0.8)",
        "rgba(146, 64, 14, 0.8)", "rgba(153, 27, 27, 0.8)",
        "rgba(91, 33, 182, 0.8)", "rgba(3, 105, 161, 0.8)",
    ];
    var out = [];
    for (var i = 0; i < n; i++) out.push(palette[i % palette.length]);
    return out;
}

function buildLongitudinalChart() {
    var canvas = document.getElementById("chart-longitudinal");
    if (!canvas || typeof Chart === "undefined") return;
    destroyChart("longitudinal");

    // Build datasets from patient data
    var datasets = [];
    var labels = [];
    var labelSet = new Set();

    // Collect all timepoints
    patients.forEach(function(p) {
        if (state.selectedPatients.length > 0 && state.selectedPatients.indexOf(p.id) === -1) return;
        (p.timepoints || []).forEach(function(tp) {
            if (!labelSet.has(tp.label)) { labelSet.add(tp.label); labels.push(tp.label); }
        });
    });
    labels.sort(function(a, b) { return parseFloat(a) - parseFloat(b); });

    var colors = getChartColors(state.selectedDwi.length);
    state.selectedDwi.forEach(function(dwi, dIdx) {
        var values = labels.map(function() { return []; });
        patients.forEach(function(p) {
            if (state.selectedPatients.length > 0 && state.selectedPatients.indexOf(p.id) === -1) return;
            (p.timepoints || []).forEach(function(tp) {
                if (tp.dwi_type !== dwi) return;
                var idx = labels.indexOf(tp.label);
                if (idx >= 0 && tp.adc_mean != null) values[idx].push(tp.adc_mean);
            });
        });
        var means = values.map(function(arr) {
            if (arr.length === 0) return null;
            return arr.reduce(function(s, v) { return s + v; }, 0) / arr.length;
        });
        datasets.push({
            label: dwi + " (mean ADC)",
            data: means,
            borderColor: colors[dIdx],
            backgroundColor: colors[dIdx].replace("0.8", "0.15"),
            tension: 0.3,
            spanGaps: true,
        });
    });

    chartInstances["longitudinal"] = new Chart(canvas, {
        type: "line",
        data: { labels: labels, datasets: datasets },
        options: {
            responsive: true,
            plugins: { title: { display: true, text: "Longitudinal ADC by DWI Type" },
                       legend: { position: "bottom" } },
            scales: { x: { title: { display: true, text: "Timepoint" } },
                      y: { title: { display: true, text: "Mean ADC (mm\\u00b2/s)" } } },
        },
    });
}

function buildDistributionChart() {
    var canvas = document.getElementById("chart-distribution");
    if (!canvas || typeof Chart === "undefined") return;
    destroyChart("distribution");

    var allValues = {};
    state.selectedDwi.forEach(function(dwi) { allValues[dwi] = []; });
    patients.forEach(function(p) {
        if (state.selectedPatients.length > 0 && state.selectedPatients.indexOf(p.id) === -1) return;
        (p.timepoints || []).forEach(function(tp) {
            if (state.selectedDwi.indexOf(tp.dwi_type) === -1) return;
            if (tp.adc_mean != null) allValues[tp.dwi_type].push(tp.adc_mean);
        });
    });

    var colors = getChartColors(state.selectedDwi.length);
    var datasets = state.selectedDwi.map(function(dwi, i) {
        return {
            label: dwi,
            data: allValues[dwi],
            backgroundColor: colors[i].replace("0.8", "0.5"),
            borderColor: colors[i],
            borderWidth: 1,
        };
    });

    chartInstances["distribution"] = new Chart(canvas, {
        type: "bar",
        data: {
            labels: state.selectedDwi,
            datasets: [{
                label: "Mean ADC values",
                data: state.selectedDwi.map(function(dwi) {
                    var vals = allValues[dwi];
                    if (vals.length === 0) return 0;
                    return vals.reduce(function(s,v){return s+v;},0) / vals.length;
                }),
                backgroundColor: colors,
            }],
        },
        options: {
            responsive: true,
            plugins: { title: { display: true, text: "Mean ADC by DWI Type" },
                       legend: { display: false } },
            scales: { y: { title: { display: true, text: "Mean ADC (mm\\u00b2/s)" }, beginAtZero: false } },
        },
    });
}

function buildDwiComparisonChart() {
    var canvas = document.getElementById("chart-dwi-compare");
    if (!canvas || typeof Chart === "undefined") return;
    destroyChart("dwi-compare");

    var metrics = ["adc_mean", "adc_median", "d_mean", "f_mean"];
    var metricLabels = ["ADC Mean", "ADC Median", "D Mean", "f Mean"];
    var colors = getChartColors(state.selectedDwi.length);

    var datasets = state.selectedDwi.map(function(dwi, dIdx) {
        var vals = metrics.map(function(metric) {
            var arr = [];
            patients.forEach(function(p) {
                if (state.selectedPatients.length > 0 && state.selectedPatients.indexOf(p.id) === -1) return;
                (p.timepoints || []).forEach(function(tp) {
                    if (tp.dwi_type !== dwi) return;
                    if (tp[metric] != null) arr.push(tp[metric]);
                });
            });
            if (arr.length === 0) return null;
            return arr.reduce(function(s,v){return s+v;},0) / arr.length;
        });
        return { label: dwi, data: vals, backgroundColor: colors[dIdx] };
    });

    chartInstances["dwi-compare"] = new Chart(canvas, {
        type: "bar",
        data: { labels: metricLabels, datasets: datasets },
        options: {
            responsive: true,
            plugins: { title: { display: true, text: "Metric Comparison Across DWI Types" },
                       legend: { position: "bottom" } },
            scales: { y: { title: { display: true, text: "Value" }, beginAtZero: false } },
        },
    });
}

// ── Sidebar navigation highlighting ──
function initNavHighlight() {
    var links = $$(".sidebar .nav-links a");
    if (links.length === 0) return;
    var sections = links.map(function(a) {
        return document.getElementById(a.getAttribute("href").slice(1));
    }).filter(Boolean);
    window.addEventListener("scroll", function() {
        var scrollY = window.scrollY + 80;
        var current = null;
        sections.forEach(function(sec) { if (sec.offsetTop <= scrollY) current = sec; });
        links.forEach(function(a) {
            a.classList.toggle("active", current && a.getAttribute("href") === "#" + current.id);
        });
    });
}

// ── Initialise on DOM ready ──
document.addEventListener("DOMContentLoaded", function() {
    initDwiFilters();
    initPatientFilter();
    initCoreMethodFilter();
    initSearch();
    initTabs();
    initSortableTables();
    initPatientCards();
    initNavHighlight();
    rebuildCharts();
    updateFilterCount();
});

})();
"""
