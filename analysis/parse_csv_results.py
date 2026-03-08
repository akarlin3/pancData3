#!/usr/bin/env python3
"""Parse pipeline CSV exports (Significant_LF_Metrics.csv, FDR_Sig_Global.csv).

Reads CSV files directly from each DWI type subfolder and produces structured
cross-DWI comparisons without requiring the Claude Vision API.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

from shared import DWI_TYPES, resolve_folder, setup_utf8_stdout

setup_utf8_stdout()


def _read_csv(path: Path) -> list[dict]:
    """Read a CSV file, return empty list if missing."""
    if not path.exists():
        return []
    with open(path, encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f))


def parse_significant_metrics(folder: Path) -> dict[str, list[dict]]:
    """Parse Significant_LF_Metrics.csv from each DWI type subfolder."""
    results: dict[str, list[dict]] = {}
    for dwi_type in DWI_TYPES:
        csv_path = folder / dwi_type / "Significant_LF_Metrics.csv"
        rows = _read_csv(csv_path)
        if rows:
            results[dwi_type] = rows
    return results


def parse_fdr_global(folder: Path) -> dict[str, list[dict]]:
    """Parse FDR_Sig_Global.csv from each DWI type subfolder."""
    results: dict[str, list[dict]] = {}
    for dwi_type in DWI_TYPES:
        csv_path = folder / dwi_type / "FDR_Sig_Global.csv"
        rows = _read_csv(csv_path)
        if rows:
            results[dwi_type] = rows
    return results


def cross_reference_significance(sig_by_dwi: dict[str, list[dict]]) -> list[dict]:
    """Find metrics significant in one DWI type but not others.

    Returns list of dicts with metric, timepoint, and which DWI types
    show significance.
    """
    # Build a set of (metric_key) → set of DWI types where significant
    metric_dwi_map: dict[str, set[str]] = {}

    for dwi_type, rows in sig_by_dwi.items():
        for row in rows:
            # Try common column names for the metric identifier
            metric_key = row.get("Metric", row.get("metric", ""))
            tp = row.get("Timepoint", row.get("timepoint", ""))
            key = f"{metric_key}@{tp}" if tp else metric_key
            if key not in metric_dwi_map:
                metric_dwi_map[key] = set()
            metric_dwi_map[key].add(dwi_type)

    cross_ref = []
    for key, dwi_set in sorted(metric_dwi_map.items()):
        parts = key.split("@", 1)
        metric = parts[0]
        tp = parts[1] if len(parts) > 1 else ""
        missing = [d for d in DWI_TYPES if d not in dwi_set]
        cross_ref.append({
            "metric": metric,
            "timepoint": tp,
            "significant_in": sorted(dwi_set),
            "not_significant_in": missing,
            "consistent": len(missing) == 0 or len(dwi_set) == 0,
        })

    return cross_ref


def parse_all_csvs(folder: Path) -> dict:
    """Parse all pipeline CSV exports.

    Returns a dict with keys: significant_metrics, fdr_global, cross_reference.
    """
    sig = parse_significant_metrics(folder)
    fdr = parse_fdr_global(folder)
    cross_ref = cross_reference_significance(sig)

    return {
        "significant_metrics": sig,
        "fdr_global": fdr,
        "cross_reference": cross_ref,
    }


def main():
    folder = resolve_folder(sys.argv)
    results = parse_all_csvs(folder)

    sep = "=" * 80

    # ── Significant metrics per DWI type ──
    print(sep)
    print("  SIGNIFICANT LF METRICS (from pipeline CSV exports)")
    print(sep)

    sig = results["significant_metrics"]
    if not sig:
        print("\n  No Significant_LF_Metrics.csv files found.")
    else:
        for dwi_type in DWI_TYPES:
            if dwi_type not in sig:
                continue
            rows = sig[dwi_type]
            print(f"\n  [{dwi_type}] — {len(rows)} significant metric(s)")
            for r in rows:
                cols = list(r.values())
                # Print first few meaningful columns
                summary = " | ".join(f"{k}={v}" for k, v in list(r.items())[:5] if v)
                print(f"    {summary}")

    # ── Cross-DWI consistency ──
    print(f"\n{sep}")
    print("  CROSS-DWI SIGNIFICANCE CONSISTENCY")
    print(sep)

    cross_ref = results["cross_reference"]
    if not cross_ref:
        print("\n  No cross-DWI comparison possible (insufficient data).")
    else:
        inconsistent = [c for c in cross_ref if not c["consistent"]]
        consistent = [c for c in cross_ref if c["consistent"]]

        if inconsistent:
            print(f"\n  Inconsistent across DWI types ({len(inconsistent)}):")
            for c in inconsistent:
                sig_in = ", ".join(c["significant_in"])
                not_in = ", ".join(c["not_significant_in"]) or "none"
                tp = f" @ {c['timepoint']}" if c["timepoint"] else ""
                print(f"    {c['metric']}{tp}")
                print(f"      Significant in: {sig_in}")
                print(f"      NOT significant in: {not_in}")

        if consistent:
            print(f"\n  Consistent across all DWI types ({len(consistent)}):")
            for c in consistent:
                tp = f" @ {c['timepoint']}" if c["timepoint"] else ""
                print(f"    {c['metric']}{tp}: {', '.join(c['significant_in'])}")

    # ── FDR summary ──
    fdr = results["fdr_global"]
    if fdr:
        print(f"\n{sep}")
        print("  FDR GLOBAL CORRECTION SUMMARY")
        print(sep)
        for dwi_type in DWI_TYPES:
            if dwi_type not in fdr:
                continue
            rows = fdr[dwi_type]
            n_sig = sum(1 for r in rows if any(
                v and v not in ("", "NaN", "nan") and "sig" in k.lower()
                for k, v in r.items()
            ))
            print(f"  [{dwi_type}] {len(rows)} tests, {n_sig} FDR-significant")

    print()


if __name__ == "__main__":
    main()
