function out = add_rician_noise(S, sigma)
% ADD_RICIAN_NOISE  Corrupt a clean magnitude-MR signal with Rician noise.
%
%   out = add_rician_noise(S, sigma)
%
% MRI reconstructs MAGNITUDE images from a complex signal whose real and
% imaginary channels each carry independent zero-mean Gaussian noise with the
% same standard deviation sigma. Taking the magnitude makes the measured value
% Rician-distributed:
%
%       S_meas = sqrt( (S + eps_r)^2 + eps_i^2 ),   eps_r, eps_i ~ N(0, sigma)
%
% WHY RICIAN, NOT GAUSSIAN (this is the physically important part)
%   At high SNR the magnitude is approximately Gaussian about the true signal
%   S. But as S -> 0 (high b-values, and the dense resistant tumour core) the
%   distribution collapses to a Rayleigh law with a strictly POSITIVE mean,
%       E[S_meas | S=0] = sigma * sqrt(pi/2),
%   the "noise floor". Plain additive Gaussian noise would (wrongly) let the
%   high-b signal go negative and would HIDE this floor. The floor is precisely
%   what destabilises IVIM fitting: the perfusion signal that should decay
%   toward zero instead rests on a noise pedestal, and the fitter mistakes that
%   pedestal for perfusion — inflating/whitening the f and D* estimates. Using
%   the correct Rician model reproduces that instability faithfully.
%
% INPUTS
%   S      array of clean (noise-free) magnitude signal values.
%   sigma  noise standard deviation per quadrature channel (same units as S).
%          A common convention is sigma = S0 / SNR, i.e. SNR referenced to the
%          b=0 signal.
%
% OUTPUT
%   out    Rician-corrupted signal, same size as S, strictly non-negative.

    if nargin ~= 2
        error('add_rician_noise:nargin', 'Usage: out = add_rician_noise(S, sigma)');
    end
    if sigma < 0
        error('add_rician_noise:sigma', 'sigma must be non-negative (got %g).', sigma);
    end
    if sigma == 0
        out = abs(S);   % degenerate noiseless case
        return;
    end
    eps_r = sigma .* randn(size(S));
    eps_i = sigma .* randn(size(S));
    out = sqrt((S + eps_r).^2 + eps_i.^2);
end
