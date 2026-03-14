#!/usr/bin/env python3
"""Parse pipeline CSV exports (Significant_LF_Metrics.csv, FDR_Sig_Global.csv).

The MATLAB pipeline exports per-DWI-type CSV files containing metrics that
reached statistical significance in Wilcoxon rank-sum tests (before and
after BH-FDR correction).  This script reads those CSVs directly and
produces a structured cross-DWI comparison to identify:

- Which metrics are significant in **all** DWI types (consistent).
- Which metrics are significant in **some but not all** types (inconsistent).

The results are also consumed by :func:`generate_report.generate_report` to
populate the "Cross-DWI Significance Inconsistencies" section of the HTML
report.

Usage:
    python parse_csv_results.py [saved_files_path]
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Optional

from tqdm import tqdm  # type: ignore

# Ensure analysis/ root is on sys.path so 'shared' is importable when run as subprocess.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import DWI_TYPES, resolve_folder, setup_utf8_stdout  # type: ignore

setup_utf8_stdout()


def _read_csv(path: Path) -> list[dict]:
    """Read a CSV file into a list of row dicts.

    Parameters
    ----------
    path : Path
        Path to the CSV file.

    Returns
    -------
    list[dict]
        Rows as dicts keyed by column header.  Returns an empty list if
        the file does not exist (graceful degradation).
    """
    if not path.exists():
        return []
    with open(path, encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f))


def _find_column(row: dict, *candidates: str) -> Optional[str]:
    """Return the first column name from *candidates* found in *row* (case-insensitive).

    Parameters
    ----------
    row : dict
        A CSV row dict whose keys are the column headers.
    *candidates : str
        Column name candidates to search for, in priority order.

    Returns
    -------
    str or None
        The actual column key as it appears in the row, or ``None`` if no
        candidate matches any column.
    """
    lower_map = {k.lower(): k for k in row}
    for c in candidates:
        actual = lower_map.get(c.lower())
        if actual is not None:
            return actual
    return None


def _normalize_timepoint(tp: str) -> str:
    """Normalise a timepoint string for consistent cross-DWI key comparison.

    Strips whitespace and lowercases so that ``"Fx5"``, ``"fx5"``, and
    ``"FX5"`` all map to the same key.

    Parameters
    ----------
    tp : str
        Raw timepoint string from a CSV cell.

    Returns
    -------
    str
        Lowercased, whitespace-stripped timepoint string.
    """
    return tp.strip().lower()


def _extract_record(row: dict) -> dict:
    """Extract a structured record from a CSV row with case-insensitive columns.

    Reads metric, timepoint, p-value, and effect size from a row using
    case-insensitive column matching.  Missing columns gracefully yield
    ``None`` rather than raising.

    Parameters
    ----------
    row : dict
        A single row dict from a Significant_LF_Metrics.csv file.

    Returns
    -------
    dict
        Keys: ``metric``, ``timepoint``, ``p_value`` (float or None),
        ``effect_size`` (float or None), ``significant`` (bool).
    """
    metric_col = _find_column(row, "Metric", "metric", "METRIC")
    tp_col = _find_column(row, "Timepoint", "timepoint", "TIMEPOINT")
    pval_col = _find_column(row, "p_value", "p value", "p-value", "pval", "P", "p")
    es_col = _find_column(row, "effect_size", "cohens_d", "auc", "d", "AUC",
                          "Effect_Size", "EffectSize")

    metric = row.get(metric_col, "") if metric_col else ""
    tp_raw = row.get(tp_col, "") if tp_col else ""
    tp = _normalize_timepoint(tp_raw) if tp_raw else tp_raw

    # Parse p-value as float; None if column absent or value non-numeric.
    p_value: Optional[float] = None
    if pval_col:
        try:
            p_value = float(row[pval_col])
        except (ValueError, TypeError):
            pass

    # Parse effect size as float; None if column absent or value non-numeric.
    effect_size: Optional[float] = None
    if es_col:
        try:
            effect_size = float(row[es_col])
        except (ValueError, TypeError):
            pass

    significant = p_value is not None and p_value < 0.05

    return {
        "metric": metric,
        "timepoint": tp,
        "p_value": p_value,
        "effect_size": effect_size,
        "significant": significant,
    }


def _temporal_pattern(records: list[dict]) -> dict[str, list[str]]:
    """Build a mapping of metric → timepoints where it is significant.

    Used by report sections to show temporal stability: a metric that is
    significant at every timepoint is more clinically robust than one that
    is only significant at a single timepoint.

    Parameters
    ----------
    records : list[dict]
        Output records from :func:`_extract_record`.

    Returns
    -------
    dict[str, list[str]]
        Mapping of ``metric_name`` → sorted list of timepoints at which that
        metric is significant (``p_value < 0.05`` or ``significant == True``).
    """
    pattern: dict[str, list[str]] = {}
    for rec in records:
        metric = rec.get("metric", "")
        tp = rec.get("timepoint", "")
        sig = rec.get("significant", False) or (
            rec.get("p_value") is not None and rec["p_value"] < 0.05
        )
        if sig and metric:
            pattern.setdefault(metric, [])
            if tp and tp not in pattern[metric]:
                pattern[metric].append(tp)
    # Sort timepoints within each metric for deterministic output.
    for metric in pattern:
        pattern[metric] = sorted(pattern[metric])
    return pattern


def parse_significant_metrics(folder: Path) -> dict[str, list[dict]]:
    """Parse ``Significant_LF_Metrics.csv`` from each DWI type subfolder.

    This CSV is exported by ``metrics_stats_comparisons.m`` and lists
    metrics with Wilcoxon rank-sum p < 0.05 for local failure vs control.
    Columns are matched case-insensitively; p-value and effect-size columns
    are extracted as floats when present.

    Parameters
    ----------
    folder : Path
        Path to the ``saved_files_*`` output folder.

    Returns
    -------
    dict[str, list[dict]]
        Mapping of DWI type name to its list of structured metric records.
        Each record has keys: ``metric``, ``timepoint``, ``p_value``,
        ``effect_size``, ``significant``.
    """
    results: dict[str, list[dict]] = {}
    for dwi_type in DWI_TYPES:
        csv_path = folder / dwi_type / "Significant_LF_Metrics.csv"
        rows = _read_csv(csv_path)
        if rows:
            results[dwi_type] = [_extract_record(r) for r in rows]
    return results


def parse_fdr_global(folder: Path) -> dict[str, list[dict]]:
    """Parse ``FDR_Sig_Global.csv`` from each DWI type subfolder.

    This CSV contains metrics surviving Benjamini-Hochberg FDR correction
    across the full feature set.

    Parameters
    ----------
    folder : Path
        Path to the ``saved_files_*`` output folder.

    Returns
    -------
    dict[str, list[dict]]
        Mapping of DWI type name to its list of FDR-significant rows.
    """
    results: dict[str, list[dict]] = {}
    for dwi_type in DWI_TYPES:
        csv_path = folder / dwi_type / "FDR_Sig_Global.csv"
        rows = _read_csv(csv_path)
        if rows:
            results[dwi_type] = rows
    return results


def cross_reference_significance(sig_by_dwi: dict[str, list[dict]]) -> list[dict]:
    """Cross-reference significant metrics across DWI types.

    For each unique ``(metric, timepoint)`` pair found in any DWI type's
    significance CSV, records which DWI types show significance and which
    do not.  A metric is "consistent" if it is significant in all three
    DWI types (or in none -- though the latter would not appear in the CSV).

    Timepoint values are normalised (lowercased, stripped) so that ``"Fx5"``
    and ``"fx5"`` are treated as the same timepoint.

    Parameters
    ----------
    sig_by_dwi : dict[str, list[dict]]
        Output of :func:`parse_significant_metrics`.

    Returns
    -------
    list[dict]
        Each dict has keys: ``metric``, ``timepoint``, ``significant_in``,
        ``not_significant_in``, ``consistent``.
    """
    # Build mapping: "metric@timepoint" -> set of DWI types.
    metric_dwi_map: dict[str, set[str]] = {}

    for dwi_type, records in sig_by_dwi.items():
        for rec in records:
            # Support both structured records (with "metric"/"timepoint" keys)
            # and raw CSV row dicts (with "Metric"/"Timepoint" keys).
            if "metric" in rec:
                metric_key = str(rec.get("metric") or "")
                tp_raw = str(rec.get("timepoint") or "")
            else:
                metric_key = str(rec.get("Metric") or rec.get("metric") or "")
                tp_raw = str(rec.get("Timepoint") or rec.get("timepoint") or "")

            # Normalise timepoint for consistent cross-DWI keying.
            tp: str = _normalize_timepoint(tp_raw) if tp_raw else ""
            key: str = f"{metric_key}@{tp}" if tp else metric_key
            if key not in metric_dwi_map:  # type: ignore
                metric_dwi_map[key] = set()  # type: ignore
            metric_dwi_map[key].add(dwi_type)  # type: ignore

    cross_ref = []
    for key, dwi_set in sorted(metric_dwi_map.items()):
        parts = key.split("@", 1)
        metric = parts[0]
        tp = parts[1] if len(parts) > 1 else ""
        # Determine which DWI types did NOT reach significance.
        missing = [d for d in DWI_TYPES if d not in dwi_set]
        cross_ref.append({
            "metric": metric,
            "timepoint": tp,
            "significant_in": sorted(dwi_set),
            "not_significant_in": missing,
            # Consistent if ALL known DWI types show significance (no missing).
            "consistent": len(missing) == 0,
        })

    return cross_ref


def parse_all_csvs(folder: Path) -> dict:
    """Parse all pipeline CSV exports and cross-reference them.

    Parameters
    ----------
    folder : Path
        Path to the ``saved_files_*`` output folder.

    Returns
    -------
    dict
        Keys: ``significant_metrics``, ``fdr_global``, ``cross_reference``,
        ``temporal_patterns``.
    """
    steps = [
        ("Significant metrics", lambda: parse_significant_metrics(folder)),
        ("FDR global", lambda: parse_fdr_global(folder)),
    ]
    pbar = tqdm(
        steps,
        desc="Parsing CSVs",
        unit="file",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}] {postfix}",
    )
    results_list = []
    for name, fn in pbar:
        pbar.set_postfix_str(name, refresh=True)
        results_list.append(fn())

    sig, fdr = results_list
    cross_ref = cross_reference_significance(sig)

    # Build temporal patterns per DWI type.
    temporal_patterns: dict[str, dict] = {}
    for dwi_type, records in sig.items():
        temporal_patterns[dwi_type] = _temporal_pattern(records)

    return {
        "significant_metrics": sig,
        "fdr_global": fdr,
        "cross_reference": cross_ref,
        "temporal_patterns": temporal_patterns,
    }


def main():
    """CLI entry point: parse CSVs and print a human-readable summary."""
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
            records = sig[dwi_type]
            print(f"\n  [{dwi_type}] \u2014 {len(records)} significant metric(s)")
            for rec in records:
                pval_str = f"p={rec['p_value']:.4f}" if rec.get("p_value") is not None else "p=N/A"
                es_str = f", effect={rec['effect_size']:.3f}" if rec.get("effect_size") is not None else ""
                print(f"    {rec['metric']} @ {rec['timepoint']}: {pval_str}{es_str}")

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

    # ── Temporal patterns ──
    tp_patterns = results.get("temporal_patterns", {})
    if tp_patterns:
        print(f"\n{sep}")
        print("  TEMPORAL PATTERNS (metric → timepoints where significant)")
        print(sep)
        for dwi_type in DWI_TYPES:
            if dwi_type not in tp_patterns:
                continue
            pat = tp_patterns[dwi_type]
            if pat:
                print(f"\n  [{dwi_type}]")
                for metric, tps in sorted(pat.items()):
                    print(f"    {metric}: {', '.join(tps)}")

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
            # Count rows where any column containing "sig" has a non-empty value.
            n_sig = sum(1 for r in rows if any(
                v and v not in ("", "NaN", "nan") and "sig" in k.lower()
                for k, v in r.items()
            ))
            print(f"  [{dwi_type}] {len(rows)} tests, {n_sig} FDR-significant")

    print()


if __name__ == "__main__":
    main()
