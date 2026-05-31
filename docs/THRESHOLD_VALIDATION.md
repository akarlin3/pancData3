# ADC Threshold Selection — Validation & Selection-Bias Correction

> **Status:** Bias measured (permutation test), corrected (nested CV: LOO +
> repeated K-fold), and surfaced in the report; a config selector lets the
> sub-volume cut be pre-specified or tactic-derived (bias-corrected). The
> headline real-cohort numbers must be produced by running the instruments on
> the secure PHI-holding machine — see §4.

This document explains why **Tactic 3** of the ADC-threshold optimizer
(`pipeline/utils/optimize_adc_threshold.m`) reports an optimistically biased
p-value, how large that bias is, and the instruments built to measure and
correct it. It is written to seed a methods paragraph for the manuscript.

---

## 1. The selection-on-outcome problem (in plain terms)

`optimize_adc_threshold` sweeps **13 candidate ADC thresholds**
(0.8×10⁻³ → 2.0×10⁻³ mm²/s, step 0.1×10⁻³). For each candidate it computes,
per patient, the **sub-volume fraction** (fraction of GTV voxels with
ADC < threshold) and runs a **Wilcoxon rank-sum test** of that fraction between
local-control (LC, `LF==0`) and local-failure (LF, `LF==1`) patients,
requiring ≥3 patients per group.

**Tactic 3 then reports the threshold with the _smallest_ p-value** and that
p-value (`significance_pvalue`).

This is **selection-on-outcome**: the *same* outcome data both **chooses** the
threshold (argmin over 13 candidates) and is then **tested** at it. The
reported p-value is therefore the **minimum of up to 13 correlated tests,
reported as if it were a single pre-specified test**. Under the null hypothesis
(no true LC/LF difference at any threshold), the minimum of several p-values is
systematically smaller than any single p-value would be — so
`significance_pvalue` is **optimistically biased**.

### Ground-truth code references
- Sweep definition: `optimize_adc_threshold.m:81` (`thresholds = 0.8e-3 : 0.1e-3 : 2.0e-3`).
- Per-threshold rank-sum + **argmin selection**: now factored into
  `pipeline/utils/select_significance_threshold.m`
  (`[best_p, best_idx] = min(pvalues)`), called from `optimize_adc_threshold.m`.
- The biased numbers are `significance_pvalue`, `significance_thresh`, and the
  per-threshold `significance_pvalues`.

### What inherits the bias (scope of contamination)
Traced end-to-end (`dispatch_pipeline_steps.m:254`, the only caller):

- The optimizer's output `opt_results` is **saved to disk and read only by the
  Python report** (`analysis/parsers/parse_mat_metrics.py:800–844`,
  `analysis/report/sections/threshold_optimization.py`).
- **The three tactic thresholds are _never_ fed back into sub-volume
  creation.** Production sub-volumes (and every downstream dose metric —
  `d95_adc_sub`, `v50_adc_sub`, and hence `compute_dose_response_roc.m`) are
  built from the **fixed, pre-specified** `config.adc_thresh` (default
  `0.001`), set once at pipeline start
  (`extract_tumor_core.m:148,256`; `compute_adc_metrics`).
- The `optimal_threshold` at `compute_dose_response_roc.m:127` is a **dose**
  cutoff from `perfcurve` (Youden index on D95/V50) — **unrelated** to the
  Tactic-3 ADC threshold.

**Conclusion:** the bias is confined to the *reported* Tactic-3 numbers in the
HTML report. No sub-volume mask, dose metric, or ROC result inherits it. The
production pipeline already uses a pre-specified threshold — itself a valid
escape from the selection bias (see §5).

---

## 2. Why not a Bonferroni correction?

The 13 thresholds are **nested / monotone**: the sub-volume at a higher
threshold is roughly a superset of the one below it, so the 13 sub-volume
fractions — and the 13 p-values — are **highly correlated**. In simulation
(below) the mean cross-threshold correlation of the per-patient sub-volume
fraction is ≈ **0.64**, so the **effective number of independent tests is
between 1 and 13**. A naive Bonferroni (×13) assumes 13 *independent* tests and
would therefore be **too conservative**.

This is precisely why a **resampling (permutation) estimate** is the right tool:
it captures the cohort's *actual* correlation structure and neither over- nor
under-corrects.

---

## 3. The measurement instrument

`pipeline/utils/optimize_adc_threshold_permutation_test.m` wraps the **entire**
Tactic-3 selection in a **label-permutation test** (Westfall–Young single-step
minP):

1. Compute the observed naive min-p on the real labels (via the identical
   production selection, `select_significance_threshold`).
2. For each of `n_perm` random permutations of the LF labels, **re-run the full
   "sweep 13 thresholds, take the minimum p" procedure** and record the
   permuted minimum p. This builds the **null distribution of the _minimum_
   statistic**.
3. **Selection-adjusted p** = fraction of permutations whose minimum p ≤ the
   observed minimum p (with the +1 plug-in estimator). The **gap between the
   naive min-p and this adjusted p is the bias, quantified.**
4. A per-threshold Westfall–Young minP-adjusted p is also returned, giving the
   count of thresholds that remain notable *after* adjusting for the 13-way
   selection.

Folds/thresholds that are inestimable (a group falls below the ≥3 floor) are
**NaN-propagated and reported**, never silently dropped. The instrument
returns NaN honestly when nothing is estimable.

---

## 4. Measured bias (simulation)

> ⚠️ **The real-cohort headline numbers must be produced by running the
> instrument on the actual cohort**, which lives only on the secure
> (PHI-holding) machine — patient data never enters the cloud/CI container per
> the project safety rules. The numbers below are from a **faithful Monte-Carlo
> model** (N=42, 14 LF / 28 LC, per-patient 2-component ADC voxel histograms so
> the per-threshold tests genuinely decorrelate) that quantifies the bias
> *mechanism* and validates the instrument. The MATLAB test
> `test_optimize_adc_threshold_permutation_test.m` reproduces the same finding.

Under **pure null** (sub-volume fraction independent of LC/LF label):

| Quantity | Naive Tactic-3 (min over 13) | After permutation adjustment |
|---|---|---|
| P(reported p < 0.05) | **0.181** | **0.036** |
| P(reported p < 0.01) | 0.041 | — |
| median reported p | **0.182** | **0.497** |

- The naive procedure calls a *non-existent* effect "significant" **~18% of the
  time** — a **≈3.6× inflation** of the nominal 5% false-positive rate. Its
  median p (0.18) sits far below the 0.5 a calibrated test averages to.
- The permutation adjustment **restores calibration**: false-positive rate back
  to ≈0.05 and median adjusted p ≈ 0.50. **The instrument detects and removes
  the bias.**
- The exact magnitude on the real cohort depends on that cohort's histogram
  heterogeneity (how much the patient ranking reorders across thresholds) — an
  empirical property the permutation test reads directly from the data.

---

## 5. What this means for the reported claim (decision pending review)

The size of the measured bias bears on whether Tactic 3 survives as a headline:

- If, on the real cohort, the **permutation-adjusted p remains small**, Tactic 3
  has genuine outcome-discriminative value beyond the selection artifact.
- If the **adjusted p washes out**, the naive `significance_pvalue` was an
  artifact of selecting over 13 thresholds, and the honest report is a
  cautionary/negative result.

Because production sub-volumes already use the **pre-specified**
`config.adc_thresh`, a clean manuscript story is available either way: the
operating threshold is pre-specified (0.001 default / 0.0016 proposed), and
Tactic 3 is reported as an *exploratory* selection with its bias measured and,
in CP1, corrected by nested LOOCV.

---

## 6. The cross-validated correction (implemented)

`pipeline/utils/optimize_adc_threshold_nested_cv.m` selects the threshold
**inside each training split** (reusing `select_significance_threshold`, so the
inner selection is identical to production) and evaluates the LC-vs-LF
association on **out-of-fold** patients whose threshold was never chosen using
their own data. One entry point, two schemes:

- **Leave-one-patient-out (LOO).** Deterministic, but a **weak corrector at
  N≈42**: removing one patient rarely changes the argmin threshold, so the
  out-of-fold p only *partially* relaxes toward the in-sample min-p. In null
  simulation it moves the mean p from ≈0.28 (in-sample) to ≈0.40 and roughly
  halves the false-positive rate — better, but not calibrated. **This is itself
  a finding:** the cut is not robustly identifiable by leave-one-out at this N.
- **Repeated stratified grouped K-fold** (`make_grouped_folds`, e.g. 5-fold ×20).
  Leaving out ~N/5 patients perturbs the selection enough to be genuinely
  out-of-fold. In simulation it is well-calibrated on null data
  (false-positive rate ≈ nominal, mean p ≈ 0.5–0.6) and recovers a true effect
  (out-of-fold significant; recommended cut within ±0.0004 of truth in 12/12
  trials). This is the **stronger** corrector.

**Hierarchy of honesty for the reported claim:** the **permutation
selection-adjusted p (§3–4) is the primary bias-corrected number**; the
repeated-K-fold out-of-fold association is the cross-validated confirmation; the
LOO estimate is reported with its weakness stated; the naive in-sample min-p is
shown only labeled "optimistic / not for headline." All four, plus the
fold-to-fold **threshold stability** (modal cut, fraction of folds choosing it,
number of distinct cuts), render in
`analysis/report/sections/threshold_optimization.py`.

## 7. Sub-volume threshold selector (config)

`config.adc_thresh_source` chooses how the sub-volume's ADC cut is set:

| value | meaning |
|---|---|
| `preset` (default) | the pre-specified `config.adc_thresh` — current behaviour, unchanged |
| `tactic1` | reproducibility optimum (max Fx1-repeat Dice) |
| `tactic2` | volume-inflection knee |
| `tactic3` | outcome-significance cut, **bias-corrected** |

`resolve_adc_thresh.m` enforces the safety rule that **Tactic 3 uses only the
nested-CV out-of-fold recommendation (repeated K-fold, else LOO), never the
biased in-sample `significance_thresh`**; any unavailable tactic falls back to
the preset. When a non-`preset` source is set, `dispatch_pipeline_steps.m`
resolves the threshold *before* sub-volumes are built. Default `preset` is a
complete no-op, so existing configs and every downstream dose/ROC metric are
unaffected unless a tactic is explicitly opted into.

## 8. What N≈42 can and cannot support

- **Can:** detect whether a *pre-specified* ADC cut separates LC/LF, and
  quantify/correct the inflation from selecting among candidate cuts.
- **Cannot:** robustly *identify* the single best cut by in-sample selection —
  the threshold is unstable fold-to-fold (LOO barely moves it; repeated K-fold
  shows multiple distinct cuts). The defensible operating threshold is therefore
  **pre-specified** (0.001 default / 0.0016 proposed), with Tactic 3 reported as
  an *exploratory, bias-corrected* association rather than a fitted cutoff.

This is publishable either way: a validated bias-corrected association if it
survives, or a rigorous cautionary result ("in-sample threshold selection on
small radiomic cohorts inflates significance; here is the corrected estimate")
if it does not.
