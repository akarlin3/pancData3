"""Report sections: statistical power commentary."""

from __future__ import annotations

import math

from shared import (  # type: ignore
    DWI_TYPES,
)
from report.report_formatters import (  # type: ignore
    _dwi_badge,
    _esc,
    _h2,
)


def _section_power_analysis(log_data, dwi_types_present, mat_data) -> list[str]:
    """Build the Statistical Power Commentary section.

    Provides context on statistical power based on the observed cohort
    size and effect sizes, helping readers interpret non-significant
    results appropriately.

    Parameters
    ----------
    log_data : dict or None
        Parsed log metrics.
    dwi_types_present : list[str]
        DWI types found in this pipeline run.
    mat_data : dict
        Parsed MAT file metrics.

    Returns
    -------
    list[str]
        HTML chunks for the power analysis section.
    """
    h: list[str] = []

    # Determine cohort size
    n_patients = 0
    if mat_data:
        for dt in DWI_TYPES:
            if dt in mat_data and "longitudinal" in mat_data[dt]:  # type: ignore
                n = mat_data[dt]["longitudinal"].get("num_patients", 0)  # type: ignore
                if n > n_patients:
                    n_patients = n

    if n_patients == 0:
        return h

    h.append(_h2("Statistical Power Context", "power"))
    h.append(
        '<p class="meta">Post-hoc power commentary based on observed cohort '
        'size and effect sizes. These are approximate guidelines, not formal '
        'power calculations (which require pre-specified effect sizes and '
        'significance levels).</p>'
    )

    # Compute minimum detectable effect sizes for common tests
    # Using rule-of-thumb: for two-sided test at alpha=0.05, 80% power,
    # Wilcoxon requires ~n per group; Cohen's d detectable ≈ 2.8/sqrt(n/2)
    n_per_group = n_patients // 2  # approximate split
    if n_per_group > 0:
        min_cohen_d = 2.8 / math.sqrt(n_per_group)
        # Minimum detectable HR (approximate): HR = exp(d * pi/sqrt(3))
        # where d is the standardised effect size on the log-hazard scale
        min_hr = math.exp(min_cohen_d * 0.5)

        h.append('<div class="methods-box">')
        h.append(f"<h3>Approximate Detectable Effect Sizes (n = {n_patients})</h3>")
        h.append(
            f"<p>With approximately {n_per_group} patients per outcome group "
            f"and conventional settings (\u03b1 = 0.05, 80% power):</p>"
        )
        h.append("<ul>")
        h.append(
            f"<li><strong>Group comparison (Wilcoxon):</strong> Minimum "
            f"detectable standardised effect size d \u2248 {min_cohen_d:.2f} "
            f"({'large' if min_cohen_d >= 0.8 else 'medium' if min_cohen_d >= 0.5 else 'small'} "
            f"effect). "
            f"{'The study is underpowered for detecting small to medium effects.' if min_cohen_d > 0.5 else 'The study has adequate power for medium effects.'}"
            f"</li>"
        )
        h.append(
            f"<li><strong>Survival analysis (Cox PH):</strong> Minimum "
            f"detectable hazard ratio HR \u2248 {min_hr:.2f} (or 1/{min_hr:.2f} = "
            f"{1/min_hr:.2f} for protective effects). Hazard ratios closer to 1.0 "
            f"than this threshold cannot be reliably detected.</li>"
        )

        n_tests_approx: int = 0
        if log_data:
            for dt in dwi_types_present:
                if dt in log_data:
                    sc = log_data[dt].get("stats_comparisons", {})
                    n_tests_approx = int(n_tests_approx + len(sc.get("glme_details", [])))  # type: ignore
        if n_tests_approx > 0:
            # After BH-FDR, effective alpha is approximately alpha * k / m
            # where k is rank; for the median test, effective alpha ≈ 0.025
            adj_min_d = 2.8 / math.sqrt(n_per_group) * 1.15  # ~15% penalty
            h.append(
                f"<li><strong>After FDR correction ({n_tests_approx} tests):</strong> "
                f"The effective significance threshold is reduced, requiring "
                f"slightly larger effects (d \u2248 {adj_min_d:.2f}) for detection. "
                f"Non-significant results after FDR correction should not be "
                f"interpreted as evidence of no effect.</li>"
            )
        h.append("</ul>")

        # Classify observed effect sizes
        if log_data:
            detectable_hrs = []
            undetectable_hrs = []
            for dt in dwi_types_present:
                if dt not in log_data:
                    continue
                hrs = log_data[dt].get("survival", {}).get("hazard_ratios", [])
                for hr_item in hrs:
                    hr_val = hr_item.get("hr", 1.0)
                    log_hr = abs(math.log(hr_val)) if hr_val > 0 else 0
                    if log_hr >= math.log(min_hr):
                        detectable_hrs.append(hr_item)
                    else:
                        undetectable_hrs.append(hr_item)

            if detectable_hrs or undetectable_hrs:
                h.append("<h4>Observed vs Detectable Effect Sizes</h4>")
                if detectable_hrs:
                    h.append(
                        f'<div class="info-box">'
                        f'{len(detectable_hrs)} covariate(s) have effect sizes '
                        f'within the study\'s detection range (|log HR| \u2265 '
                        f'{math.log(min_hr):.2f}).</div>'
                    )
                if undetectable_hrs:
                    h.append(
                        f'<div class="warn-box">'
                        f'{len(undetectable_hrs)} covariate(s) have effect sizes '
                        f'below the study\'s detection threshold. Non-significant '
                        f'p-values for these covariates may reflect insufficient '
                        f'power rather than absence of effect.</div>'
                    )

        h.append("</div>")

    return h
