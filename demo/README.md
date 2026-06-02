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

| File | Role |
|---|---|
| [`synthetic_ivim.m`](synthetic_ivim.m) | **Teaching reference.** IVIM forward model `S(b)=S0·(f·e^{-bD*}+(1-f)·e^{-bD})`, Rician noise, ground-truth cohort + resistant low-ADC sub-volume. Over-commented for study/rewrite. |
| [`validate_fitting.m`](validate_fitting.m) | Recover-known-truth check: runs the real fitter across an SNR sweep and scores recovery per parameter. |
| [`run_demo.m`](run_demo.m) | End-to-end entry: generate → real fit → validate → headline figures. |
| [`fit_synthetic_scan.m`](fit_synthetic_scan.m) | Thin wrapper that calls the **real** `pipeline/core/fit_models.m` on one phantom scan. |
| [`demo_plot.m`](demo_plot.m) | Figure renderer (native MATLAB/Octave; gnuplot fallback for headless/CI). |
| [`add_demo_paths.m`](add_demo_paths.m) | Puts the pipeline modules (+ Octave shims) on the path. |

## What's real vs. synthetic

- **Synthetic:** the input DWI volumes and the (D, f, D\*) ground truth — generated here from documented distributions with a fixed seed.
- **Real (unmodified):** the fitting — `pipeline/core/fit_models.m` → `IVIMmodelfit` → `IVIM_seg`. The demo exercises the actual pipeline physics code; it does not re-implement it.

## Note on the forward model

`synthetic_ivim.m` is intentionally scaffolding: the comments are a tutorial
for re-deriving the IVIM forward model and the Rician noise floor from scratch.
It is meant to be **rewritten by hand** by the project owner, not treated as a
finished module.
