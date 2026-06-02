function adc = adc_from_signal(b, S)
% ADC_FROM_SIGNAL  Mono-exponential ADC by weighted log-linear least squares.
%
%   adc = adc_from_signal(b, S)
%
% Estimates the apparent diffusion coefficient from a single voxel's multi-b
% signal, using the SAME weighted estimator the pipeline uses
% (pipeline/core/fit_models.m). Computing the phantom's "true" ADC with this
% estimator keeps ground truth defined consistently with what the pipeline
% recovers, so ADC recovery can be scored fairly.
%
% MODEL & DERIVATION
%   Mono-exponential:        S(b) = S0 * exp(-b * ADC)
%   Linearised:              ln(S/S0) = -b * ADC
%   The log transform turns (approximately Gaussian) signal noise into
%   heteroscedastic noise with variance ~ 1/S^2 (delta method). Weighted least
%   squares with weights w = S^2 corrects for this, up-weighting high-SNR
%   points. Closed form for the single predictor A = -b:
%       ADC = sum(w .* A .* y) / sum(w .* A.^2),   y = ln(S/S0)
%   b=0 is the S0 reference and is excluded from the regression (ln(S0/S0)=0
%   carries no information about ADC).
%
% INPUTS
%   b   vector of b-values (s/mm^2); b(1) must be 0 (the S0 reference).
%   S   matching vector of signal magnitudes.
%
% OUTPUT
%   adc apparent diffusion coefficient (mm^2/s). Returns NaN when the signal is
%       non-positive (log undefined) or the fit is non-physical (adc < 0).

    if nargin ~= 2
        error('adc_from_signal:nargin', 'Usage: adc = adc_from_signal(b, S)');
    end
    b = b(:); S = S(:);
    if b(1) ~= 0
        error('adc_from_signal:noB0', 'b(1) must be 0 (the S0 reference).');
    end
    S0 = S(1);
    if S0 <= 0 || any(S(2:end) <= 0)
        adc = NaN; return;
    end
    A = -b(2:end);
    y = log(S(2:end) ./ S0);
    w = S(2:end).^2;
    adc = sum(w .* y .* A) / sum(w .* (A.^2));
    if ~isfinite(adc) || adc < 0
        adc = NaN;
    end
end
