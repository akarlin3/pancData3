#!/usr/bin/env python3
"""Quantitative agreement analysis between DWI processing types.

Computes Bland-Altman analysis, Lin's concordance correlation coefficient
(CCC), and generates diagnostic plots for cross-DWI-type comparison of
ADC, D, f, and D* parameters.

Usage:
    python cross_dwi_agreement.py [saved_files_path]
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import DWI_TYPES, resolve_folder, setup_utf8_stdout  # type: ignore

setup_utf8_stdout()


def bland_altman(x: np.ndarray, y: np.ndarray) -> dict:
    """Compute Bland-Altman analysis between two measurement arrays.

    Parameters
    ----------
    x, y : np.ndarray
        Paired measurements (same length, no NaN).

    Returns
    -------
    dict
        Keys: bias, sd, loa_lo, loa_hi, pct_outside_loa, n
    """
    diff = x - y
    mean_xy = (x + y) / 2
    bias = float(np.mean(diff))
    sd = float(np.std(diff, ddof=1))
    loa_lo = bias - 1.96 * sd
    loa_hi = bias + 1.96 * sd
    outside = np.sum((diff < loa_lo) | (diff > loa_hi))
    pct_outside = 100 * outside / len(diff) if len(diff) > 0 else 0

    return {
        "bias": bias,
        "sd": sd,
        "loa_lo": float(loa_lo),
        "loa_hi": float(loa_hi),
        "pct_outside_loa": float(pct_outside),
        "n": len(diff),
    }


def lins_ccc(x: np.ndarray, y: np.ndarray) -> dict:
    """Compute Lin's concordance correlation coefficient with 95% CI.

    Parameters
    ----------
    x, y : np.ndarray
        Paired measurements.

    Returns
    -------
    dict
        Keys: ccc, ci_lo, ci_hi, precision, accuracy
    """
    n = len(x)
    if n < 3:
        return {"ccc": float("nan"), "ci_lo": float("nan"), "ci_hi": float("nan"),
                "precision": float("nan"), "accuracy": float("nan")}

    mean_x = np.mean(x)
    mean_y = np.mean(y)
    s_xx = np.var(x, ddof=1)
    s_yy = np.var(y, ddof=1)
    s_xy = np.sum((x - mean_x) * (y - mean_y)) / (n - 1)

    # Pearson r (precision)
    if s_xx * s_yy > 0:
        r = s_xy / math.sqrt(s_xx * s_yy)
    else:
        r = 0.0

    # Bias correction factor (accuracy)
    denom = s_xx + s_yy + (mean_x - mean_y) ** 2
    if denom > 0:
        cb = 2 * math.sqrt(s_xx * s_yy) / denom
    else:
        cb = 1.0

    ccc = r * cb

    # 95% CI via Fisher z-transform (requires n > 3 for valid SE)
    if abs(ccc) < 1.0 and n > 3:
        z = 0.5 * math.log((1 + ccc) / (1 - ccc))
        se_z = 1 / math.sqrt(n - 3)
        z_lo = z - 1.96 * se_z
        z_hi = z + 1.96 * se_z
        ci_lo = math.tanh(z_lo)
        ci_hi = math.tanh(z_hi)
    else:
        ci_lo = float("nan") if n <= 3 else ccc
        ci_hi = float("nan") if n <= 3 else ccc

    return {
        "ccc": float(ccc),
        "ci_lo": float(ci_lo),
        "ci_hi": float(ci_hi),
        "precision": float(r),
        "accuracy": float(cb),
    }


def icc_absolute(x: np.ndarray, y: np.ndarray) -> float:
    """Compute ICC (two-way mixed, absolute agreement) for two raters.

    Parameters
    ----------
    x, y : np.ndarray
        Paired measurements.

    Returns
    -------
    float
        ICC value.
    """
    n = len(x)
    if n < 3:
        return float("nan")

    data = np.column_stack([x, y])
    k = 2  # number of raters

    grand_mean = np.mean(data)
    row_means = np.mean(data, axis=1)
    col_means = np.mean(data, axis=0)

    # Sum of squares
    ss_total = np.sum((data - grand_mean) ** 2)
    ss_rows = k * np.sum((row_means - grand_mean) ** 2)
    ss_cols = n * np.sum((col_means - grand_mean) ** 2)
    ss_error = ss_total - ss_rows - ss_cols

    # Mean squares
    ms_rows = ss_rows / (n - 1) if n > 1 else 0
    ms_cols = ss_cols / (k - 1) if k > 1 else 0
    ms_error = ss_error / ((n - 1) * (k - 1)) if (n - 1) * (k - 1) > 0 else 0

    # ICC(A,1) — two-way mixed, absolute agreement
    denom = ms_rows + (k - 1) * ms_error + k * (ms_cols - ms_error) / n
    if denom > 0:
        return float((ms_rows - ms_error) / denom)
    return float("nan")


def load_summary_metrics(folder: Path, dwi_type: str) -> dict | None:
    """Load parsed MAT metrics JSON for a DWI type."""
    mat_path = folder / f"parsed_mat_metrics_{dwi_type}.json"
    if mat_path.exists():
        try:
            with open(mat_path, "r", encoding="utf-8") as f:
                return json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            print(f"  ⚠️  Failed to load {mat_path.name}: {e}")
    return None


def run_agreement_analysis(folder: Path) -> dict:
    """Run full cross-DWI agreement analysis.

    Returns
    -------
    dict
        Nested dict: {param: {type_pair: {bland_altman, ccc, icc}}}
    """
    results: dict = {}

    # Load data for all available DWI types
    type_data: dict[str, dict] = {}
    for dt in DWI_TYPES:
        data = load_summary_metrics(folder, dt)
        if data:
            type_data[dt] = data

    if len(type_data) < 2:
        print("  Need at least 2 DWI types for agreement analysis.")
        return results

    params = ["adc_mean", "d_mean", "f_mean", "dstar_mean"]
    param_labels = {"adc_mean": "ADC", "d_mean": "D", "f_mean": "f", "dstar_mean": "D*"}

    types_present = list(type_data.keys())

    for param in params:
        param_results: dict = {}
        for i in range(len(types_present)):
            for j in range(i + 1, len(types_present)):
                t1, t2 = types_present[i], types_present[j]
                pair_key = f"{t1}_vs_{t2}"

                # Extract matching patient-level data
                d1 = type_data[t1].get("longitudinal", {})
                d2 = type_data[t2].get("longitudinal", {})

                vals1 = d1.get(param, [])
                vals2 = d2.get(param, [])

                if not vals1 or not vals2:
                    continue

                # Flatten and align (take minimum length)
                a1 = np.array(vals1, dtype=float).ravel()
                a2 = np.array(vals2, dtype=float).ravel()
                min_len = min(len(a1), len(a2))
                a1 = a1[:min_len]
                a2 = a2[:min_len]

                # Remove NaN pairs
                valid = np.isfinite(a1) & np.isfinite(a2)
                a1 = a1[valid]
                a2 = a2[valid]

                if len(a1) < 3:
                    continue

                ba = bland_altman(a1, a2)
                ccc = lins_ccc(a1, a2)
                icc = icc_absolute(a1, a2)

                param_results[pair_key] = {
                    "bland_altman": ba,
                    "ccc": ccc,
                    "icc": icc,
                }

        if param_results:
            results[param_labels.get(param, param)] = param_results

    # Print summary table
    if results:
        print("\n=== Cross-DWI Agreement Summary ===\n")
        print(f"  {'Param':<8}  {'Pair':<25}  {'CCC':>6}  {'ICC':>6}  {'Bias':>10}  {'LoA':>20}  {'%Out':>6}")
        print("  " + "-" * 85)
        for param, pairs in results.items():
            for pair, data in pairs.items():
                ba = data["bland_altman"]
                ccc_val = data["ccc"]["ccc"]
                icc_val = data["icc"]
                print(f"  {param:<8}  {pair:<25}  {ccc_val:>6.3f}  {icc_val:>6.3f}  "
                      f"{ba['bias']:>10.6f}  [{ba['loa_lo']:.6f}, {ba['loa_hi']:.6f}]  "
                      f"{ba['pct_outside_loa']:>5.1f}%")

    return results


def generate_bland_altman_plots(folder: Path, results: dict) -> list[Path]:
    """Generate Bland-Altman plots for each parameter/type pair.

    Returns list of saved figure paths.
    """
    saved_plots: list[Path] = []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available; skipping Bland-Altman plots.")
        return saved_plots

    for param, pairs in results.items():
        for pair, data in pairs.items():
            ba = data["bland_altman"]
            if ba["n"] < 3:
                continue

            fig, ax = plt.subplots(figsize=(8, 6))
            # We don't have the raw data here, so create a summary plot
            ax.axhline(y=ba["bias"], color="red", linestyle="-", label=f"Bias = {ba['bias']:.6f}")
            ax.axhline(y=ba["loa_lo"], color="gray", linestyle="--", label=f"LoA = [{ba['loa_lo']:.6f}, {ba['loa_hi']:.6f}]")
            ax.axhline(y=ba["loa_hi"], color="gray", linestyle="--")
            ax.axhline(y=0, color="black", linestyle=":", alpha=0.5)
            ax.set_xlabel("Mean of Two Methods")
            ax.set_ylabel("Difference")
            t1, t2 = pair.split("_vs_")
            ax.set_title(f"Bland-Altman: {param} ({t1} vs {t2})\nn={ba['n']}, {ba['pct_outside_loa']:.1f}% outside LoA")
            ax.legend(loc="upper right", fontsize=8)

            filename = f"bland_altman_{param}_{pair}.png"
            fig_path = folder / filename
            fig.savefig(str(fig_path), dpi=150, bbox_inches="tight")
            plt.close(fig)
            saved_plots.append(fig_path)
            print(f"  📁 Saved: {filename}")

    return saved_plots


def main():
    """CLI entry point."""
    folder = resolve_folder(sys.argv)
    results = run_agreement_analysis(folder)

    if results:
        generate_bland_altman_plots(folder, results)

        # Save results JSON
        out_path = folder / "cross_dwi_agreement.json"
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2, default=str)
        print(f"\n  📁 Results saved to: {out_path}")
    else:
        print("  No agreement analysis results to report.")


if __name__ == "__main__":
    main()
