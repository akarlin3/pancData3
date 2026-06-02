# Synthetic IVIM demo (PHI-free)

A self-contained, **zero-PHI** demo of the pancData3 diffusion pipeline. Every
"patient" is simulated from a known IVIM ground truth, so a reviewer can clone
and run the project end-to-end on data they can trust is fake.

> **All artifacts here are synthetic phantom data, not clinical results.**

## Run it

```bash
# MATLAB:        run('demo/run_demo.m')
# Octave (CLI):  octave demo/run_demo.m
```

Outputs land in `demo/output/` (git-ignored, regenerable):

| File | What |
|---|---|
| `synthetic_longitudinal_D.png` | Recovered D vs fraction, local control vs local failure |
| `synthetic_longitudinal_f.png` | Recovered f vs fraction, by outcome |
| `synthetic_recovery_vs_snr.png` | Recovery error vs SNR for D, f, D\* |
| `synthetic_recovery_scatter.png` | Truth-vs-recovered D at best/worst SNR |
| `synthetic_ground_truth.csv` | The answer key: per-region true (D, f, D\*, ADC) |
| `synthetic_recovered_longitudinal.csv` | Per-patient/fraction recovered means |

## Files

**Physics primitives** (standalone, reusable, unit-tested):

| File | Role |
|---|---|
| [`ivim_signal.m`](ivim_signal.m) | IVIM bi-exponential forward model `S(b)=S0·(f·e^{-bD*}+(1-f)·e^{-bD})`. |
| [`add_rician_noise.m`](add_rician_noise.m) | Magnitude-MR (Rician) noise with the noise floor that destabilises D*/f. |
| [`adc_from_signal.m`](adc_from_signal.m) | Weighted log-linear ADC estimator (matches the pipeline). |

**Cohort / driver:**

| File | Role |
|---|---|
| [`synthetic_ivim.m`](synthetic_ivim.m) | Builds the phantom cohort: per-voxel ground truth, resistant low-ADC core, longitudinal D/f/D* response model. Heavily documented. |
| [`validate_fitting.m`](validate_fitting.m) | Recover-known-truth check: runs the real fitter across an SNR sweep and scores recovery per parameter. |
| [`run_demo.m`](run_demo.m) | End-to-end entry: generate → real fit → validate → headline figures. |
| [`run_synthetic_demo.m`](run_synthetic_demo.m) | The function that `run_demo.m` calls (does the actual work). |
| [`fit_synthetic_scan.m`](fit_synthetic_scan.m) | Thin wrapper that calls the **real** `pipeline/core/fit_models.m` on one phantom scan. |
| [`demo_plot.m`](demo_plot.m) | Figure renderer (native MATLAB/Octave; gnuplot fallback for headless/CI). |
| [`add_demo_paths.m`](add_demo_paths.m) | Puts the pipeline modules (+ Octave shims) on the path. |
| [`tests/`](tests) | Unit tests for the physics primitives and the cohort generator. |

## What's real vs. synthetic

- **Synthetic:** the input DWI volumes and the (D, f, D\*) ground truth — generated here from documented distributions with a fixed seed.
- **Real (unmodified):** the fitting — `pipeline/core/fit_models.m` → `IVIMmodelfit` → `IVIM_seg`. The demo exercises the actual pipeline physics code; it does not re-implement it.

## The physics

The IVIM forward model and the Rician noise model are implemented as small,
standalone, documented functions (`ivim_signal.m`, `add_rician_noise.m`,
`adc_from_signal.m`) with unit tests in `tests/`. They are the load-bearing
science of the demo and can be reused or extended directly. The comments
explain every equation and parameter range so the model is fully transparent,
not a black box.

## Tests

```bash
octave demo/tests/run_demo_tests.m
```

Covers: forward-model correctness (analytic checks at b=0 and limits),
Rician-noise statistics (non-negativity, the Rayleigh noise floor),
ADC round-trip on noiseless signal, cohort reproducibility under a fixed seed,
and that the cohort feeds the real fitter with recovery within tolerance.
