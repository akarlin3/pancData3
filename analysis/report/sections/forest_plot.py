"""Forest plot generation for hazard ratios across timepoints and DWI types.

Generates a proper matplotlib forest plot figure showing hazard ratios with
95% CIs, grouped by DWI type with visual separators.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from shared import DWI_TYPES  # type: ignore


def extract_hr_data(log_data: dict | None, dwi_types_present: list[str]) -> list[dict]:
    """Extract hazard ratio data from parsed log metrics.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.

    Returns
    -------
    list[dict]
        List of dicts with keys: dwi_type, covariate, hr, ci_lo, ci_hi, p
    """
    hr_entries: list[dict] = []
    if not log_data:
        return hr_entries

    for dwi_type in dwi_types_present:
        if dwi_type not in log_data:
            continue
        sv = log_data[dwi_type].get("survival", {})
        hrs = sv.get("hazard_ratios", [])
        for hr in hrs:
            entry = {
                "dwi_type": dwi_type,
                "covariate": hr.get("covariate", "Unknown"),
                "hr": hr.get("hr", 1.0),
                "ci_lo": hr.get("ci_lo", hr.get("hr", 1.0)),
                "ci_hi": hr.get("ci_hi", hr.get("hr", 1.0)),
                "p": hr.get("p", 1.0),
            }
            hr_entries.append(entry)

    return hr_entries


def generate_forest_plot(
    hr_data: list[dict],
    output_path: Path,
    title: str = "Forest Plot: Hazard Ratios",
) -> bool:
    """Generate a forest plot figure from hazard ratio data.

    Parameters
    ----------
    hr_data : list[dict]
        Hazard ratio entries from extract_hr_data.
    output_path : Path
        Path to save the PNG figure.
    title : str
        Plot title.

    Returns
    -------
    bool
        True if figure was saved successfully.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        print("  matplotlib not available; skipping forest plot.")
        return False

    if not hr_data:
        return False

    # Sort by DWI type, then by covariate
    type_order = {dt: i for i, dt in enumerate(DWI_TYPES)}
    hr_data_sorted = sorted(
        hr_data,
        key=lambda x: (type_order.get(x["dwi_type"], 99), x["covariate"]),
    )

    n_rows = len(hr_data_sorted)
    fig_height = max(4, 0.5 * n_rows + 2)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    # Colors by DWI type
    colors = {"Standard": "#2196F3", "dnCNN": "#FF9800", "IVIMnet": "#4CAF50"}
    y_positions = list(range(n_rows - 1, -1, -1))

    # Track DWI type transitions for separator lines
    prev_type = None
    separator_ys: list[float] = []

    labels: list[str] = []
    for i, entry in enumerate(hr_data_sorted):
        if prev_type is not None and entry["dwi_type"] != prev_type:
            separator_ys.append(y_positions[i] + 0.5)
        prev_type = entry["dwi_type"]

        hr_val = entry["hr"]
        ci_lo = entry["ci_lo"]
        ci_hi = entry["ci_hi"]
        p_val = entry["p"]
        color = colors.get(entry["dwi_type"], "#757575")

        # Log scale: plot on log10 axis
        if hr_val > 0 and ci_lo > 0 and ci_hi > 0:
            log_hr = math.log10(hr_val)
            log_ci_lo = math.log10(ci_lo)
            log_ci_hi = math.log10(ci_hi)
        else:
            continue

        # CI line
        ax.plot(
            [log_ci_lo, log_ci_hi],
            [y_positions[i], y_positions[i]],
            color=color,
            linewidth=1.5,
            solid_capstyle="round",
        )

        # Point estimate
        marker = "o" if p_val < 0.05 else "o"
        facecolor = color if p_val < 0.05 else "white"
        ax.plot(
            log_hr,
            y_positions[i],
            marker=marker,
            markersize=8,
            color=color,
            markerfacecolor=facecolor,
            markeredgecolor=color,
            markeredgewidth=1.5,
            zorder=5,
        )

        # Label
        label = f"{entry['dwi_type']}: {entry['covariate']}"
        labels.append(label)

    # Reference line at HR = 1 (log10(1) = 0)
    ax.axvline(x=0, color="black", linestyle="--", linewidth=0.8, alpha=0.7)

    # Separator lines between DWI types
    for sep_y in separator_ys:
        ax.axhline(y=sep_y, color="gray", linestyle=":", linewidth=0.5, alpha=0.5)

    # Axis formatting
    ax.set_yticks(y_positions[:len(labels)])
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Hazard Ratio (log₁₀ scale)", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold")

    # Custom x-axis ticks showing actual HR values
    import numpy as _np
    hr_ticks = [0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 10.0]
    log_ticks = [math.log10(t) for t in hr_ticks]
    ax.set_xticks(log_ticks)
    ax.set_xticklabels([str(t) for t in hr_ticks])

    # Legend
    legend_handles = []
    for dt in DWI_TYPES:
        if any(e["dwi_type"] == dt for e in hr_data_sorted):
            c = colors.get(dt, "#757575")
            legend_handles.append(mpatches.Patch(color=c, label=dt))

    filled = plt.Line2D([], [], marker="o", color="gray", markerfacecolor="gray",
                        linestyle="None", markersize=8, label="p < 0.05")
    hollow = plt.Line2D([], [], marker="o", color="gray", markerfacecolor="white",
                        markeredgecolor="gray", linestyle="None", markersize=8, label="p ≥ 0.05")
    legend_handles.extend([filled, hollow])

    ax.legend(handles=legend_handles, loc="lower right", fontsize=8, framealpha=0.9)
    ax.grid(axis="x", alpha=0.3)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=200, bbox_inches="tight")
    plt.close(fig)
    return True


def _section_forest_plot_figure(log_data, dwi_types_present, output_folder=None) -> list[str]:
    """Build an HTML section with the forest plot figure embedded.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types present.
    output_folder : Path or None
        If provided, generate and save the figure.

    Returns
    -------
    list[str]
        HTML chunks for the forest plot section.
    """
    from report.report_formatters import _h2, _figure_caption  # type: ignore

    h: list[str] = []

    hr_data = extract_hr_data(log_data, dwi_types_present)
    if not hr_data:
        return h

    if output_folder is not None:
        fig_path = Path(output_folder) / "forest_plot_hazard_ratios.png"
        success = generate_forest_plot(hr_data, fig_path)
        if success:
            h.append('<div class="figure-container">')
            h.append(f'<img src="forest_plot_hazard_ratios.png" '
                     f'alt="Forest Plot: Hazard Ratios" '
                     f'style="max-width: 100%;">')
            h.append(_figure_caption(
                "Forest plot of hazard ratios across DWI types",
                "Filled circles indicate p < 0.05; hollow circles indicate p ≥ 0.05. "
                "Vertical dashed line at HR = 1 (no effect). "
                "Error bars show 95% Wald confidence intervals.",
            ))
            h.append("</div>")

    return h
