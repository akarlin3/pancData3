# pancData3 — For Reviewers

A one-page tour for PIs, committee reviewers, and successor researchers.
If you have 5 minutes, read this. If you have 30 minutes, read this then ONBOARDING.md.

## What this pipeline does (3 sentences)

pancData3 processes diffusion-weighted MRI from pancreatic cancer patients treated on the Elekta Unity 1.5T MR-Linac, fitting IVIM and ADC biophysical models (with optional DnCNN / IVIMnet deep-learning denoising) to identify treatment-resistant tumor sub-volumes. It then overlays radiotherapy dose maps onto those sub-volumes and runs competing-risks survival analysis correlating dose coverage with local failure. The clinical question driving the project: do ADC-defined radioresistant sub-volumes receive adequate dose under the standard 33/50 Gy SIB schedule?

## Headline result (v2.4, N=42)

| Metric | Value | Interpretation |
|---|---|---|
| D95 (ADC sub-volume) | ≈ 36.9 Gy | 26% below the 50 Gy SIB prescription |
| V50 (ADC sub-volume) | ≈ 69% | Below the 90% radical-intent target |
| Spearman ρ (D95 vs LF proxy) | −0.37 (p ≈ 0.016) | Lower dose coverage → higher local-failure signal |
| GTV-volume confounding | <5% HR change | Signal is not driven by tumor size |
| Cohort | 42 patients (12 LF, 30 LC) | Underpowered; nominally significant |

Clinical interpretation: ADC-defined resistant regions receive roughly 26%
less dose than the radical-intent threshold. Motivates daily adaptive
replanning targeted at the diffusion-defined sub-volume.

## Reproducing it

1. Clone `main`, follow ONBOARDING.md §1.2 to configure.
2. `cd pipeline; execute_all_workflows` — runs Standard, dnCNN, IVIMnet sequentially (~hours).
3. `py -3.12 analysis/run_analysis.py` — produces the PDF report.
4. The numbers above appear in the dose_context and data_supplemental sections of the report.

## What to be skeptical of

- N=42 is underpowered. The dose-coverage finding is the most defensible
  result; IVIM biomarker survival results are directionally consistent but
  nominally significant at best.
- The earlier Cox PH "D under DnCNN HR=0.14, p=0.036" headline from v2.3
  did not replicate after the cohort grew from 40 → 42 patients. v2.4 leads
  with dosimetry for this reason.
- The pruning cutoff `max_core_failure_rate = 0.25` keeps 6 of 11 tumor-core
  methods. `adc_threshold` is retained regardless of cutoff (hardcoded
  safety guard in `pipeline/utils/filter_core_methods.m`).
- One of the 11 methods (`fdm`, functional diffusion map) requires
  longitudinal data and is excluded from Fx1-only analyses.

## Where the science lives in the code

| Question | File |
|---|---|
| How is the resistant sub-volume defined? | `pipeline/utils/extract_tumor_core.m` (11 methods) |
| How is D95 computed? | `pipeline/core/metrics_dosimetry.m` |
| How is the Cox PH model fit? | `pipeline/core/metrics_survival.m` |
| How is cross-pipeline agreement computed? | `pipeline/utils/compute_cross_pipeline_dice.m` |
| Where do the report numbers come from? | `analysis/report/sections/*.py` |

## Document map

- **You are here:** FOR_REVIEWERS.md — 5-minute tour
- ONBOARDING.md — 30-minute hands-on guide for someone who will rerun it
- README.md — full feature list, install, Docker, CI
- GRAPH_FAQ.md — every figure in the report, with caveats
- CHANGELOG.md — what changed and when
- CLAUDE.md / CLAUDE_REFERENCE.md / CLAUDE_WORKFLOWS.md — *internal AI-assistant scaffolding, not user docs*
