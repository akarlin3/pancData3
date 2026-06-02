#!/usr/bin/env python3
"""Dosimetry and longitudinal MAT-file parsers for the analysis suite.

These per-block parsers were extracted from
:mod:`parsers.parse_mat_metrics` to keep that module under the project
file-length limit.  Each function reads one ``.mat`` file for a single DWI
type and mutates the shared ``out_data`` dict in place.  They are
re-imported back into ``parse_mat_metrics`` so its public API is unchanged.

Covers: sub-volume dosimetry, whole-GTV dosimetry from baseline results,
and the summary-metrics longitudinal scope (cohort dimensions, Fx1 repeat
Dice, whole-GTV mean dose, and sub-volume sizes).
"""

import sys
import typing
from pathlib import Path

# Ensure analysis/ root is on sys.path so 'parsers' is importable when this
# module is imported from a subprocess-executed sibling.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from parsers.parse_mat_helpers import (  # type: ignore
    METHOD_DESC,
    _array_to_list,
    _nested_safe,
    _safe_float,
    numpy_np,
    scipy_io,
)


def _parse_dosimetry(folder, dwi, out_data):
    # ── Dosimetry ──
    # The dosimetry MAT file contains 1-D arrays of per-patient dose
    # coverage metrics.  We compute the cohort mean and std for each.
    dosimetry_mat = folder / f"{dwi}" / f"metrics_dosimetry_results_{dwi}.mat"
    if dosimetry_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(dosimetry_mat), squeeze_me=True)  # type: ignore

            def _safe_nanmean_std(arr: typing.Any) -> dict:
                """Compute nanmean and nanstd; return dict with mean and std.

                NaN results are returned as None so JSON serialisation does
                not produce invalid output.
                """
                a = numpy_np.asarray(arr, dtype=float).ravel()  # type: ignore
                if a.size == 0 or numpy_np.all(numpy_np.isnan(a)):  # type: ignore
                    return {"mean": None, "std": None}
                mean_val = _safe_float(numpy_np.nanmean(a))  # type: ignore
                std_val = _safe_float(numpy_np.nanstd(a))  # type: ignore
                return {"mean": mean_val, "std": std_val}

            nan_scalar = numpy_np.array([float("nan")])  # type: ignore
            out_data["dosimetry"] = {  # type: ignore
                "d95_adc_mean": _safe_nanmean_std(mat.get("d95_adc_sub", nan_scalar)),  # type: ignore
                "v50_adc_mean": _safe_nanmean_std(mat.get("v50_adc_sub", nan_scalar)),  # type: ignore
                "dmean_adc_mean": _safe_nanmean_std(mat.get("dmean_adc_sub", nan_scalar)),  # type: ignore
                "d95_d_mean": _safe_nanmean_std(mat.get("d95_d_sub", nan_scalar)),  # type: ignore
                "v50_d_mean": _safe_nanmean_std(mat.get("v50_d_sub", nan_scalar)),  # type: ignore
                "dmean_d_mean": _safe_nanmean_std(mat.get("dmean_d_sub", nan_scalar)),  # type: ignore
            }
        except Exception as e:
            print(f"Error parsing {dosimetry_mat.name}: {e}")



def _parse_baseline_dosimetry(folder, dwi, out_data):
    # ── Whole-GTV dosimetry (from baseline_results) ──
    # The baseline results MAT file stores m_d95_gtvp and m_v50gy_gtvp
    # (per-patient whole-GTV D95 / V50 arrays).  We compute cohort
    # nanmean / nanstd and add them to out_data["dosimetry"] alongside
    # the sub-volume metrics.  dmean_gtvp lives on summary_metrics.
    baseline_mat = folder / f"{dwi}" / f"metrics_baseline_results_{dwi}.mat"
    if baseline_mat.exists() and scipy_io and numpy_np:
        try:
            mat: typing.Any = scipy_io.loadmat(str(baseline_mat), squeeze_me=True, struct_as_record=False)  # type: ignore

            def _whole_mean_std(arr_in: typing.Any) -> typing.Tuple[typing.Optional[float], typing.Optional[float]]:
                if arr_in is None:
                    return None, None
                try:
                    a = numpy_np.asarray(arr_in, dtype=float)  # type: ignore
                    # Use baseline (column 0) if 2-D; else flatten.
                    if a.ndim >= 2:
                        col = a[:, 0]
                    else:
                        col = a
                    return (
                        _safe_float(numpy_np.nanmean(col)),  # type: ignore
                        _safe_float(numpy_np.nanstd(col)),   # type: ignore
                    )
                except Exception:
                    return None, None

            dosi_entry = out_data.get("dosimetry") or {}
            d95_m, d95_s = _whole_mean_std(mat.get("m_d95_gtvp"))
            v50_m, v50_s = _whole_mean_std(mat.get("m_v50gy_gtvp"))
            if d95_m is not None:
                dosi_entry["whole_gtv_d95_mean"] = d95_m
                dosi_entry["whole_gtv_d95_std"] = d95_s
            if v50_m is not None:
                dosi_entry["whole_gtv_v50_mean"] = v50_m
                dosi_entry["whole_gtv_v50_std"] = v50_s

            # m_lf may live in baseline file (used by subvolume_sizes).
            m_lf_arr = mat.get("m_lf")
            if m_lf_arr is not None:
                try:
                    out_data["_baseline_m_lf"] = _array_to_list(numpy_np.asarray(m_lf_arr).ravel())  # type: ignore
                except Exception:
                    pass

            out_data["dosimetry"] = dosi_entry
        except Exception as e:
            print(f"Error parsing {baseline_mat.name}: {e}")



def _parse_summary(folder, dwi, out_data):
    # ── Longitudinal scope (from summary_metrics) ──
    # We only extract the array shape to determine how many patients and
    # timepoints were processed.
    # Expected shape: (N_patients, N_params, N_timepoints) or (N_patients, N_timepoints).
    # num_patients = shape[0]; num_timepoints = shape[-1] (last axis).
    summary_mat = folder / f"{dwi}" / f"summary_metrics_{dwi}.mat"
    if summary_mat.exists():
        try:
            mat: typing.Any = scipy_io.loadmat(str(summary_mat), squeeze_me=True, struct_as_record=False)  # type: ignore
            summary = mat.get("summary_metrics")
            if summary and hasattr(summary, "ADC_abs"):
                adc = summary.ADC_abs
                if hasattr(adc, "shape") and len(adc.shape) >= 2:
                    out_data["longitudinal"]["num_patients"] = int(adc.shape[0])
                    out_data["longitudinal"]["num_timepoints"] = int(adc.shape[-1])

            # ── Fx1 Repeat Dice (spatial repeatability) ──
            # summary_metrics stores dice_rpt_{adc,d,f,dstar} as
            # [nPatients x 3] arrays (3 DWI types: Standard=0, DnCNN=1, IVIMnet=2).
            # Extract the column for the current DWI type and compute
            # mean/std across patients with non-NaN entries.
            if summary:
                dwi_col_map = {"Standard": 0, "dnCNN": 1, "DnCNN": 1, "IVIMnet": 2}
                dwi_col = dwi_col_map.get(dwi)
                if dwi_col is not None:
                    rpt_fields = {
                        "adc": "dice_rpt_adc",
                        "d": "dice_rpt_d",
                        "f": "dice_rpt_f",
                        "dstar": "dice_rpt_dstar",
                    }
                    rpt_entry: dict = {}
                    for param, field_name in rpt_fields.items():
                        if not hasattr(summary, field_name):
                            continue
                        arr = getattr(summary, field_name)
                        if not hasattr(arr, "shape"):
                            continue
                        # Select column for this DWI type (may be 1-D if single DWI type).
                        try:
                            if arr.ndim >= 2 and arr.shape[-1] > dwi_col:
                                col = arr[:, dwi_col]
                            elif arr.ndim == 1:
                                col = arr
                            else:
                                continue
                        except Exception:
                            continue
                        col = numpy_np.asarray(col, dtype=float)  # type: ignore
                        valid = col[~numpy_np.isnan(col)]  # type: ignore
                        n_valid = int(valid.size)
                        if n_valid == 0:
                            rpt_entry[param] = {"mean": None, "std": None, "n": 0}
                        else:
                            rpt_entry[param] = {
                                "mean": _safe_float(numpy_np.mean(valid)),  # type: ignore
                                "std": _safe_float(numpy_np.std(valid)),  # type: ignore
                                "n": n_valid,
                            }
                    if rpt_entry:
                        out_data["repeatability_dice"] = rpt_entry

            # ── Whole-GTV mean dose (dmean_gtvp on summary_metrics) ──
            if summary and hasattr(summary, "dmean_gtvp"):
                try:
                    dm = numpy_np.asarray(summary.dmean_gtvp, dtype=float)  # type: ignore
                    if dm.ndim >= 2:
                        dm_col = dm[:, 0]
                    else:
                        dm_col = dm
                    dosi_entry = out_data.get("dosimetry") or {}
                    dosi_entry["whole_gtv_dmean_mean"] = _safe_float(numpy_np.nanmean(dm_col))  # type: ignore
                    dosi_entry["whole_gtv_dmean_std"] = _safe_float(numpy_np.nanstd(dm_col))  # type: ignore
                    out_data["dosimetry"] = dosi_entry
                except Exception:
                    pass

            # ── Sub-Volume Sizes (adc_sub_vol / adc_sub_vol_pc) ──
            # Shapes: [nPatients x nTimepoints x 3].  Extract the column for
            # this DWI type and stratify by outcome (m_lf from baseline).
            if summary:
                dwi_col_map = {"Standard": 0, "dnCNN": 1, "DnCNN": 1, "IVIMnet": 2}
                dwi_col = dwi_col_map.get(dwi)
                if dwi_col is not None and hasattr(summary, "adc_sub_vol"):
                    try:
                        vol = numpy_np.asarray(summary.adc_sub_vol, dtype=float)  # type: ignore
                        if vol.ndim == 3 and vol.shape[-1] > dwi_col:
                            vol2d = vol[:, :, dwi_col]  # [nPat x nTp]
                        elif vol.ndim == 2:
                            vol2d = vol
                        else:
                            vol2d = None

                        pc2d = None
                        if hasattr(summary, "adc_sub_vol_pc"):
                            pc = numpy_np.asarray(summary.adc_sub_vol_pc, dtype=float)  # type: ignore
                            if pc.ndim == 3 and pc.shape[-1] > dwi_col:
                                pc2d = pc[:, :, dwi_col]
                            elif pc.ndim == 2:
                                pc2d = pc

                        if vol2d is not None and vol2d.ndim == 2 and vol2d.size > 0:
                            n_tp = int(vol2d.shape[1])
                            # Outcome vector — prefer the raw clinical
                            # encoding from summary_metrics.lf (binary
                            # 0 = LC / 1 = LF, read straight from the
                            # clinical spreadsheet) over the baseline
                            # m_lf variant, which may re-encode competing-
                            # risk deaths as 2 and therefore empty the LF
                            # group for cohorts with few pure local-failure
                            # events.  Fall back to baseline m_lf only when
                            # summary.lf is unavailable or sized wrong.
                            n_pat_summary = int(vol2d.shape[0])
                            m_lf = None
                            if hasattr(summary, "lf"):
                                try:
                                    cand = numpy_np.asarray(summary.lf, dtype=float).ravel()  # type: ignore
                                    if cand.size == n_pat_summary:
                                        m_lf = cand
                                except Exception:
                                    m_lf = None
                            if m_lf is None:
                                m_lf_list = out_data.get("_baseline_m_lf") or []
                                if m_lf_list:
                                    cand = numpy_np.asarray(m_lf_list, dtype=float)  # type: ignore
                                    if cand.size == n_pat_summary:
                                        m_lf = cand

                            # Attempt scipy Mann-Whitney (optional dependency).
                            try:
                                from scipy.stats import mannwhitneyu  # type: ignore
                                have_mwu = True
                            except Exception:
                                mannwhitneyu = None  # type: ignore
                                have_mwu = False

                            timepoints: list[int] = []
                            lc_mean_vol: list[typing.Optional[float]] = []
                            lc_std_vol: list[typing.Optional[float]] = []
                            lf_mean_vol: list[typing.Optional[float]] = []
                            lf_std_vol: list[typing.Optional[float]] = []
                            lc_mean_pct: list[typing.Optional[float]] = []
                            lf_mean_pct: list[typing.Optional[float]] = []
                            wilcoxon_p: list[typing.Optional[float]] = []
                            median_vox: list[typing.Optional[float]] = []
                            overall_mean_vol: list[typing.Optional[float]] = []
                            overall_std_vol: list[typing.Optional[float]] = []
                            overall_mean_pct: list[typing.Optional[float]] = []

                            for k in range(n_tp):
                                timepoints.append(k + 1)
                                col_vol = vol2d[:, k]
                                col_pc = pc2d[:, k] if pc2d is not None else None

                                # Overall median voxel count (vol is cm^3; use as-is).
                                try:
                                    median_vox.append(_safe_float(numpy_np.nanmedian(col_vol)))  # type: ignore
                                except Exception:
                                    median_vox.append(None)

                                # Overall mean/std across all patients
                                # (not stratified by outcome).
                                try:
                                    finite_vol = col_vol[~numpy_np.isnan(col_vol)]  # type: ignore
                                    overall_mean_vol.append(
                                        _safe_float(numpy_np.nanmean(finite_vol)) if finite_vol.size else None)  # type: ignore
                                    overall_std_vol.append(
                                        _safe_float(numpy_np.nanstd(finite_vol)) if finite_vol.size else None)  # type: ignore
                                except Exception:
                                    overall_mean_vol.append(None)
                                    overall_std_vol.append(None)
                                try:
                                    if col_pc is not None:
                                        finite_pc = col_pc[~numpy_np.isnan(col_pc)]  # type: ignore
                                        overall_mean_pct.append(
                                            _safe_float(numpy_np.nanmean(finite_pc)) if finite_pc.size else None)  # type: ignore
                                    else:
                                        overall_mean_pct.append(None)
                                except Exception:
                                    overall_mean_pct.append(None)

                                if m_lf is not None and m_lf.size == col_vol.size:
                                    lf_mask = (m_lf == 1) & ~numpy_np.isnan(col_vol)  # type: ignore
                                    lc_mask = (m_lf == 0) & ~numpy_np.isnan(col_vol)  # type: ignore
                                    lf_vals = col_vol[lf_mask]
                                    lc_vals = col_vol[lc_mask]
                                else:
                                    lf_vals = numpy_np.array([])  # type: ignore
                                    lc_vals = numpy_np.array([])  # type: ignore

                                lc_mean_vol.append(_safe_float(numpy_np.nanmean(lc_vals)) if lc_vals.size else None)  # type: ignore
                                lc_std_vol.append(_safe_float(numpy_np.nanstd(lc_vals)) if lc_vals.size else None)  # type: ignore
                                lf_mean_vol.append(_safe_float(numpy_np.nanmean(lf_vals)) if lf_vals.size else None)  # type: ignore
                                lf_std_vol.append(_safe_float(numpy_np.nanstd(lf_vals)) if lf_vals.size else None)  # type: ignore

                                if col_pc is not None and m_lf is not None and m_lf.size == col_pc.size:
                                    lf_pc_mask = (m_lf == 1) & ~numpy_np.isnan(col_pc)  # type: ignore
                                    lc_pc_mask = (m_lf == 0) & ~numpy_np.isnan(col_pc)  # type: ignore
                                    lf_pc_vals = col_pc[lf_pc_mask]
                                    lc_pc_vals = col_pc[lc_pc_mask]
                                    lc_mean_pct.append(_safe_float(numpy_np.nanmean(lc_pc_vals)) if lc_pc_vals.size else None)  # type: ignore
                                    lf_mean_pct.append(_safe_float(numpy_np.nanmean(lf_pc_vals)) if lf_pc_vals.size else None)  # type: ignore
                                else:
                                    lc_mean_pct.append(None)
                                    lf_mean_pct.append(None)

                                p_val = None
                                if have_mwu and lf_vals.size > 0 and lc_vals.size > 0:
                                    try:
                                        _, p_val = mannwhitneyu(lf_vals, lc_vals, alternative="two-sided")  # type: ignore
                                        p_val = _safe_float(p_val)
                                    except Exception:
                                        p_val = None
                                wilcoxon_p.append(p_val)

                            out_data["subvolume_sizes"] = {
                                "timepoints": timepoints,
                                "overall": {
                                    "mean_vol_cm3": overall_mean_vol,
                                    "std_vol_cm3": overall_std_vol,
                                    "mean_frac_pct": overall_mean_pct,
                                },
                                "lc": {
                                    "mean_vol_cm3": lc_mean_vol,
                                    "std_vol_cm3": lc_std_vol,
                                    "mean_frac_pct": lc_mean_pct,
                                },
                                "lf": {
                                    "mean_vol_cm3": lf_mean_vol,
                                    "std_vol_cm3": lf_std_vol,
                                    "mean_frac_pct": lf_mean_pct,
                                },
                                "wilcoxon_p": wilcoxon_p,
                                "median_voxel_count": median_vox,
                            }
                    except Exception:
                        pass
        except Exception as e:
            print(f"Error parsing {summary_mat.name}: {e}")

