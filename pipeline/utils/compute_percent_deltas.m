function [ADC_pct, D_pct, f_delta, Dstar_pct] = compute_percent_deltas(ADC_abs, D_abs, f_abs, Dstar_abs)
% COMPUTE_PERCENT_DELTAS — Treatment-induced changes from baseline.
%
% Computes percent change = (value_at_timepoint - baseline) / baseline * 100
% for ADC, D, and D*. Uses absolute delta for f because f values near zero
% (typical range 0.05-0.15) make percent change numerically unstable.
%
% Fixed, physiologically motivated epsilon values prevent inflated percent
% changes when baseline values are near zero. Percent changes are winsorized
% at +/-500% to cap artifacts from near-zero baselines.
%
% Inputs:
%   ADC_abs   - [nPatients x nTimepoints] ADC absolute values
%   D_abs     - [nPatients x nTimepoints] D absolute values
%   f_abs     - [nPatients x nTimepoints] f absolute values
%   Dstar_abs - [nPatients x nTimepoints] D* absolute values
%
% Outputs:
%   ADC_pct   - [nPatients x nTimepoints] ADC percent change from baseline
%   D_pct     - [nPatients x nTimepoints] D percent change from baseline
%   f_delta   - [nPatients x nTimepoints] f absolute change from baseline
%   Dstar_pct - [nPatients x nTimepoints] D* percent change from baseline

% Fixed epsilon values (~1% of typical physiological range):
%   ADC: 0.001-0.003 mm^2/s -> eps = 1e-5
%   D:   0.001-0.003 mm^2/s -> eps = 1e-5
%   D*:  0.005-0.050 mm^2/s -> eps = 5e-5
adc_eps  = 1e-5;
d_eps    = 1e-5;
dstar_eps = 5e-5;

% Exclude patients with near-zero or negative baselines from percent change
% computation to avoid sign-flipped or inflated ratios.
adc_bl = ADC_abs(:,1);  adc_bl(adc_bl < adc_eps) = NaN;
d_bl   = D_abs(:,1);    d_bl(d_bl < d_eps) = NaN;
dstar_bl = Dstar_abs(:,1); dstar_bl(dstar_bl < dstar_eps) = NaN;
ADC_pct = ((ADC_abs - ADC_abs(:,1)) ./ adc_bl) * 100;
D_pct   = ((D_abs - D_abs(:,1)) ./ d_bl) * 100;
% f uses ABSOLUTE delta instead of percent change because:
%   1. Baseline f values are often near zero (0.05-0.15)
%   2. Percent change from 0.05 can be misleadingly large
%   3. The physiological range of f is bounded [0,1]
f_delta = (f_abs - f_abs(:,1));
Dstar_pct = ((Dstar_abs - Dstar_abs(:,1)) ./ dstar_bl) * 100;

% Warn if a significant fraction of the cohort lost percent-change data
n_total = size(ADC_abs, 1);
n_nan_adc = sum(isnan(adc_bl));
n_nan_d = sum(isnan(d_bl));
n_nan_dstar = sum(isnan(dstar_bl));
if n_nan_adc > 0.2 * n_total || n_nan_d > 0.2 * n_total || n_nan_dstar > 0.2 * n_total
    fprintf('  ⚠️  Near-zero baseline exclusions: ADC=%d/%d, D=%d/%d, D*=%d/%d patients — percent changes set to NaN.\n', ...
        n_nan_adc, n_total, n_nan_d, n_total, n_nan_dstar, n_total);
end

% Winsorize percent changes at +/-500% to limit influence of near-zero baselines.
pct_clip = 500;
ADC_pct(ADC_pct < -pct_clip) = -pct_clip;  ADC_pct(ADC_pct > pct_clip) = pct_clip;
D_pct(D_pct < -pct_clip) = -pct_clip;      D_pct(D_pct > pct_clip) = pct_clip;
Dstar_pct(Dstar_pct < -pct_clip) = -pct_clip;  Dstar_pct(Dstar_pct > pct_clip) = pct_clip;
end
