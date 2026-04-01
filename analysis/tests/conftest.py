"""Shared fixtures for the pancData3 analysis test suite.

Provides reusable temporary directory structures, synthetic log files,
and mock CSV data that mirror the real pipeline output layout expected
by all analysis scripts.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import os
import struct
import sys
import zlib
from pathlib import Path

import pytest  # type: ignore

# Ensure the analysis/ directory is on sys.path so `import shared` works.
ANALYSIS_DIR = Path(__file__).resolve().parent.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))


# ---------------------------------------------------------------------------
# Dependency pre-flight check
# ---------------------------------------------------------------------------
# Verify that required third-party packages are installed before test
# collection begins.  This gives a single, clear error message instead of
# cryptic ImportError tracebacks scattered across individual test files.

_REQUIRED_PACKAGES = ["tqdm", "pydantic", "markdown"]
_missing = [pkg for pkg in _REQUIRED_PACKAGES if importlib.util.find_spec(pkg) is None]
if _missing:
    _requirements_file = ANALYSIS_DIR / "requirements.txt"
    pytest.exit(
        f"Missing required packages: {', '.join(_missing)}. "
        f"Install them with:  pip install -r {_requirements_file}",
        returncode=1,
    )


# ---------------------------------------------------------------------------
# Synthetic data constants
# ---------------------------------------------------------------------------

SAMPLE_GRAPH_CSV_ROWS = [
    # Row 1: Standard box plot with a significant p-value (0.003) and a
    # notable correlation (r=0.65).  The "increasing" trend direction
    # deliberately differs from the dnCNN row below so cross-DWI comparison
    # tests can verify AGREE/DIFFER logic.
    {
        "file_path": "saved_files_20260301_120000/Standard/Feature_BoxPlots_Standard.png",
        "graph_title": "Feature BoxPlots",
        "graph_type": "box",
        "x_axis_label": "Feature",
        "x_axis_units": "",
        "x_axis_range_min": "",
        "x_axis_range_max": "",
        "x_axis_scale_type": "categorical",
        "x_axis_tick_count": "4",
        "x_axis_tick_labels_json": json.dumps(["ADC", "D", "f", "D*"]),
        "y_axis_label": "Value",
        "y_axis_units": "mm²/s",
        "y_axis_range_min": "0",
        "y_axis_range_max": "0.003",
        "y_axis_scale_type": "linear",
        "y_axis_tick_count": "",
        "y_axis_tick_labels_json": "[]",
        "color_axis_label": "",
        "color_axis_units": "",
        "color_axis_range_min": "",
        "color_axis_range_max": "",
        "color_axis_scale_type": "",
        "color_axis_tick_count": "",
        "color_axis_tick_labels_json": "[]",
        "num_trends": "1",
        "trends_json": json.dumps([
            {"series": "ADC", "direction": "increasing", "description": "ADC rises over time",
             "magnitude": "~15% increase", "statistical_significance": "p=0.003",
             "confidence_band": None, "start_value": 0.001, "end_value": 0.0015}
        ]),
        "num_inflection_points": "0",
        "inflection_points_json": "[]",
        "num_statistical_tests": "1",
        "statistical_tests_json": json.dumps([
            {"test_name": "Wilcoxon", "statistic_value": None, "p_value": 0.003,
             "comparison_groups": "LF vs LC"}
        ]),
        "num_outliers": "0",
        "outliers_json": "[]",
        "num_reference_lines": "0",
        "reference_lines_json": "[]",
        "num_issues": "1",
        "issues_json": json.dumps(["Y-axis label partially obscured by legend"]),
        "summary": "Box plot showing p = 0.003 for ADC comparison. r = 0.65 correlation.",
        "sample_size": "n=42",
        "data_series_count": "2",
        "error_bars": "IQR",
        "num_annotations": "1",
        "annotations_json": json.dumps(["**"]),
        "clinical_relevance": "Higher ADC in LF patients suggests less restricted diffusion associated with treatment failure.",
        "data_density": "moderate",
        "spatial_pattern": "",
        "num_legend_items": "2",
        "legend_items_json": json.dumps(["LF", "LC"]),
        "subpanel_count": "",
        "comparison_type": "unpaired",
        "figure_quality": "high",
    },
    # Row 2: dnCNN box plot — same graph name as row 1 but with a
    # non-significant p-value (0.12) and opposite trend direction
    # ("decreasing" vs "increasing").  Used to test cross-DWI
    # inconsistency detection and grouping by graph name.
    {
        "file_path": "saved_files_20260301_120000/dnCNN/Feature_BoxPlots_dnCNN.png",
        "graph_title": "Feature BoxPlots",
        "graph_type": "box",
        "x_axis_label": "Feature",
        "x_axis_units": "",
        "x_axis_range_min": "",
        "x_axis_range_max": "",
        "x_axis_scale_type": "categorical",
        "x_axis_tick_count": "4",
        "x_axis_tick_labels_json": json.dumps(["ADC", "D", "f", "D*"]),
        "y_axis_label": "Value",
        "y_axis_units": "mm²/s",
        "y_axis_range_min": "0",
        "y_axis_range_max": "0.003",
        "y_axis_scale_type": "linear",
        "y_axis_tick_count": "",
        "y_axis_tick_labels_json": "[]",
        "color_axis_label": "",
        "color_axis_units": "",
        "color_axis_range_min": "",
        "color_axis_range_max": "",
        "color_axis_scale_type": "",
        "color_axis_tick_count": "",
        "color_axis_tick_labels_json": "[]",
        "num_trends": "1",
        "trends_json": json.dumps([
            {"series": "ADC", "direction": "decreasing", "description": "ADC drops after denoising",
             "magnitude": None, "statistical_significance": None,
             "confidence_band": None, "start_value": None, "end_value": None}
        ]),
        "num_inflection_points": "0",
        "inflection_points_json": "[]",
        "num_statistical_tests": "0",
        "statistical_tests_json": "[]",
        "num_outliers": "1",
        "outliers_json": json.dumps([
            {"approximate_x": None, "approximate_y": 0.0028, "series": "ADC", "description": "extreme high ADC value"}
        ]),
        "num_reference_lines": "0",
        "reference_lines_json": "[]",
        "num_issues": "0",
        "issues_json": "[]",
        "summary": "Box plot showing p = 0.12 for ADC comparison.",
        "sample_size": "n=42",
        "data_series_count": "2",
        "error_bars": "IQR",
        "num_annotations": "0",
        "annotations_json": "[]",
        "clinical_relevance": "",
        "data_density": "moderate",
        "spatial_pattern": "",
        "num_legend_items": "2",
        "legend_items_json": json.dumps(["LF", "LC"]),
        "subpanel_count": "",
        "comparison_type": "unpaired",
        "figure_quality": "high",
    },
    # Row 3: IVIMnet longitudinal line plot — a different graph type ("line")
    # with two trends (LF decreasing, LC flat), one inflection point, and both
    # a significant p-value (0.02) and an r-squared value (0.45).  Exercises
    # multi-trend parsing, inflection point extraction, and longitudinal graph
    # grouping separately from the box plot rows above.
    {
        "file_path": "saved_files_20260301_120000/IVIMnet/Longitudinal_Mean_Metrics_IVIMnet.png",
        "graph_title": "Longitudinal Mean Metrics",
        "graph_type": "line",
        "x_axis_label": "Days",
        "x_axis_units": "days",
        "x_axis_range_min": "0",
        "x_axis_range_max": "180",
        "x_axis_scale_type": "linear",
        "x_axis_tick_count": "7",
        "x_axis_tick_labels_json": "[]",
        "y_axis_label": "ADC",
        "y_axis_units": "mm²/s",
        "y_axis_range_min": "0.0005",
        "y_axis_range_max": "0.002",
        "y_axis_scale_type": "linear",
        "y_axis_tick_count": "",
        "y_axis_tick_labels_json": "[]",
        "color_axis_label": "",
        "color_axis_units": "",
        "color_axis_range_min": "",
        "color_axis_range_max": "",
        "color_axis_scale_type": "",
        "color_axis_tick_count": "",
        "color_axis_tick_labels_json": "[]",
        "num_trends": "2",
        "trends_json": json.dumps([
            {"series": "LF", "direction": "decreasing", "description": "LF group shows decline",
             "magnitude": "~30% decline over 180 days", "statistical_significance": "p=0.02",
             "confidence_band": "95% CI shaded", "start_value": 0.0018, "end_value": 0.0012},
            {"series": "LC", "direction": "flat", "description": "LC group stays stable",
             "magnitude": None, "statistical_significance": None,
             "confidence_band": "95% CI shaded", "start_value": 0.0015, "end_value": 0.0014},
        ]),
        "num_inflection_points": "1",
        "inflection_points_json": json.dumps([
            {"approximate_x": 60.0, "approximate_y": 0.0012, "description": "LF diverges from LC at day 60",
             "magnitude": 0.0003}
        ]),
        "num_statistical_tests": "1",
        "statistical_tests_json": json.dumps([
            {"test_name": "log-rank", "statistic_value": 5.41, "p_value": 0.02,
             "comparison_groups": "LF vs LC"}
        ]),
        "num_outliers": "0",
        "outliers_json": "[]",
        "num_reference_lines": "1",
        "reference_lines_json": json.dumps([
            {"orientation": "horizontal", "value": 0.001, "label": "ADC threshold", "style": "dashed"}
        ]),
        "num_issues": "2",
        "issues_json": json.dumps(["Overlapping legend text", "X-axis tick labels truncated"]),
        "summary": "Longitudinal plot. p-value = 0.02 for group difference. r² = 0.45 fit quality.",
        "sample_size": "n=38",
        "data_series_count": "2",
        "error_bars": "SEM",
        "num_annotations": "2",
        "annotations_json": json.dumps(["*", "p=0.02"]),
        "clinical_relevance": "Declining ADC in LF group suggests progressive tumour cellularity associated with poor RT response.",
        "data_density": "moderate",
        "spatial_pattern": "",
        "num_legend_items": "2",
        "legend_items_json": json.dumps(["LF", "LC"]),
        "subpanel_count": "",
        "comparison_type": "longitudinal",
        "figure_quality": "medium",
    },
]

# Column order for the graph CSV (matches batch_graph_analysis.CSV_COLUMNS)
GRAPH_CSV_COLUMNS = [
    "file_path", "graph_title", "graph_type",
    "x_axis_label", "x_axis_units", "x_axis_range_min", "x_axis_range_max",
    "x_axis_scale_type", "x_axis_tick_count", "x_axis_tick_labels_json",
    "y_axis_label", "y_axis_units", "y_axis_range_min", "y_axis_range_max",
    "y_axis_scale_type", "y_axis_tick_count", "y_axis_tick_labels_json",
    "color_axis_label", "color_axis_units", "color_axis_range_min", "color_axis_range_max",
    "color_axis_scale_type", "color_axis_tick_count", "color_axis_tick_labels_json",
    "num_trends", "trends_json", "num_inflection_points", "inflection_points_json",
    "num_statistical_tests", "statistical_tests_json",
    "num_outliers", "outliers_json",
    "num_reference_lines", "reference_lines_json",
    "num_issues", "issues_json",
    "summary", "sample_size", "data_series_count", "error_bars",
    "num_annotations", "annotations_json",
    "clinical_relevance", "data_density", "spatial_pattern",
    "num_legend_items", "legend_items_json",
    "subpanel_count", "comparison_type", "figure_quality",
]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def saved_files_dir(tmp_path: Path) -> Path:
    """Create a minimal saved_files_* directory tree with DWI subfolders.

    Returns the path to the created saved_files_YYYYMMDD_HHMMSS directory.
    """
    folder = tmp_path / "saved_files_20260301_120000"
    for dwi_type in ("Standard", "dnCNN", "IVIMnet"):
        (folder / dwi_type).mkdir(parents=True)
    return folder


@pytest.fixture
def saved_files_with_graph_csv(saved_files_dir: Path) -> Path:
    """saved_files dir populated with a graph_analysis_results.csv.

    Writes all three SAMPLE_GRAPH_CSV_ROWS (Standard box, dnCNN box,
    IVIMnet line) into the CSV file at the root of the saved_files directory,
    mirroring the layout produced by batch_graph_analysis.py after a real
    pipeline run.
    """
    csv_path = saved_files_dir / "graph_analysis_results.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=GRAPH_CSV_COLUMNS)
        writer.writeheader()
        for row in SAMPLE_GRAPH_CSV_ROWS:
            writer.writerow(row)
    return saved_files_dir


@pytest.fixture
def saved_files_with_logs(saved_files_dir: Path) -> Path:
    """saved_files dir populated with synthetic MATLAB log files.

    Creates four log files under the Standard/ subfolder, each containing
    representative text that exercises a different parse_log_metrics parser:

    - metrics_stats_comparisons: GLME interaction p-values, per-metric
      details with adjusted alpha, FDR timepoint counts, and competing-risk
      patient exclusion.
    - metrics_stats_predictive: Elastic net feature selections at two
      timepoints, Firth refit confirmation, and two ROC analysis blocks
      with AUC/sensitivity/specificity.
    - metrics_survival: Two hazard ratio rows, a global LRT line, and
      IPCW weight range.
    - metrics_baseline: Per-metric outlier flags with group breakdowns,
      total outlier summary, and baseline exclusion with LF rate comparison.
    """
    # -- Standard logs --
    std = saved_files_dir / "Standard"

    # Stats comparisons log: contains GLME interaction p-values (decimal and
    # scientific notation), per-metric Wilcoxon results with FDR-adjusted
    # alpha, timepoint-level FDR summaries, and competing-risk exclusion count.
    (std / "metrics_stats_comparisons_output_Standard.txt").write_text(
        "Interaction P-Value for ADC: 0.023\n"
        "Interaction P-Value for D_star: 1.5e-04\n"
        "mean_adc: p=0.015, adj_alpha=0.025\n"
        "mean_d*: p=0.042, adj_alpha=0.0125\n"
        "Timepoint: BL — 3 significant\n"
        "Timepoint: W2 — 1 significant\n"
        "Excluded 5/42 (11.9%) competing-risk patients\n",
        encoding="utf-8",
    )

    # Predictive log: elastic net feature selections at BL and W2 with
    # different lambda values and feature counts, plus two PRIMARY ROC
    # ANALYSIS blocks with AUC, Youden cutoff, sensitivity, and specificity.
    (std / "metrics_stats_predictive_output_Standard.txt").write_text(
        "Elastic Net Selected Features for BL (Opt Lambda=0.0532): mean_adc, mean_d, f_ratio\n"
        "Elastic Net Selected Features for W2 (Opt Lambda=0.1200): delta_adc, mean_f\n"
        "Firth refit successful for BL (3 features)\n"
        "PRIMARY ROC ANALYSIS results for BL\n"
        "  AUC = 0.781\n"
        "  Youden Optimal Score Cutoff = 0.450\n"
        "  Sensitivity = 82.5% | Specificity = 68.2%\n"
        "PRIMARY ROC ANALYSIS results for W2\n"
        "  AUC = 0.692\n"
        "  Youden Optimal Score Cutoff = 0.520\n"
        "  Sensitivity = 75.0% | Specificity = 60.0%\n",
        encoding="utf-8",
    )

    # Survival log: bracket-format Cox PH table rows (matching MATLAB output),
    # Schoenfeld residuals PH test, TD Panel summary, a global likelihood
    # ratio test line, and an IPCW weight range annotation.
    (std / "metrics_survival_output_Standard.txt").write_text(
        "  [TD Panel] 35 patients \u2192 210 intervals (10 events of interest, 0 competing)\n"
        "\n"
        "  Feature       Coeff       HR       95% CI  p-value\n"
        "  ----------------------------------------------------\n"
        "  mean_adc      0.223    1.250 [ 0.98- 1.59]   0.0680\n"
        "  delta_d      -0.288    0.750 [ 0.55- 1.02]   0.0340\n"
        "Global LRT: chi2(2) = 7.82, p = 0.0200\n"
        "IPCW weights applied ... range [0.85, 1.42]\n"
        "\n"
        "  --- Schoenfeld Residuals: PH Assumption Test ---\n"
        "  Covariate        rho     chi2  p-value  PH_violated\n"
        "  --------------------------------------------------------\n"
        "  mean_adc      0.0909    0.0826    0.7737             \n"
        "  delta_d       0.7091    5.0281    0.0249          ***\n",
        encoding="utf-8",
    )

    # Baseline log: per-metric outlier flags with LF/LC/CR group breakdowns,
    # a total outlier summary line, and a baseline exclusion block including
    # LF rate comparison between included and excluded patients.
    (std / "metrics_baseline_output_Standard.txt").write_text(
        "Outlier flag (mean_adc): 3 flagged (LF=2, LC=1, CR=0)\n"
        "Outlier flag (mean_d): 1 flagged (LF=0, LC=0, CR=1)\n"
        "Total outliers removed: 4 / 42 (9.52%)\n"
        "Excluded 6/48 patients due to missing baseline\n"
        "LF rate: included=35.7%, excluded=50.0%\n",
        encoding="utf-8",
    )

    return saved_files_dir


@pytest.fixture
def saved_files_with_csvs(saved_files_dir: Path) -> Path:
    """saved_files dir populated with Significant_LF_Metrics.csv per DWI type.

    Creates CSVs with deliberately different metric sets across DWI types
    to exercise cross-DWI consistency analysis:

    - Standard: 3 rows (mean_adc@BL, mean_d@BL, mean_adc@W2)
    - dnCNN:    1 row  (mean_adc@BL only)
    - IVIMnet:  no CSV (tests missing-file graceful handling)

    This means mean_adc@BL appears in Standard and dnCNN but not IVIMnet
    (inconsistent), while mean_d@BL and mean_adc@W2 appear only in Standard
    (also inconsistent).  No metric is consistent across all three types.
    """
    for dwi_type in ("Standard", "dnCNN"):
        csv_path = saved_files_dir / dwi_type / "Significant_LF_Metrics.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["Metric", "Timepoint", "p_value"])
            writer.writeheader()
            if dwi_type == "Standard":
                writer.writerow({"Metric": "mean_adc", "Timepoint": "BL", "p_value": "0.003"})
                writer.writerow({"Metric": "mean_d", "Timepoint": "BL", "p_value": "0.01"})
                writer.writerow({"Metric": "mean_adc", "Timepoint": "W2", "p_value": "0.04"})
            else:
                # dnCNN: mean_adc@BL is significant, mean_d@BL is NOT (inconsistency)
                writer.writerow({"Metric": "mean_adc", "Timepoint": "BL", "p_value": "0.005"})

    return saved_files_dir


# ---------------------------------------------------------------------------
# Image helpers
# ---------------------------------------------------------------------------

def make_tiny_png() -> bytes:
    """Create a minimal valid 1x1 pixel PNG image for testing."""
    raw_data = b'\x00\x00\x00\x00'
    compressed = zlib.compress(raw_data)
    def chunk(ctype: bytes, data: bytes) -> bytes:
        c = ctype + data
        return struct.pack('>I', len(data)) + c + struct.pack('>I', zlib.crc32(c) & 0xffffffff)
    return (b'\x89PNG\r\n\x1a\n' +
            chunk(b'IHDR', struct.pack('>IIBBBBB', 1, 1, 8, 2, 0, 0, 0)) +
            chunk(b'IDAT', compressed) +
            chunk(b'IEND', b''))
