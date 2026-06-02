function S = ivim_signal(b, S0, D, f, Dstar)
% IVIM_SIGNAL  Bi-exponential intravoxel-incoherent-motion (IVIM) forward model.
%
%   S = ivim_signal(b, S0, D, f, Dstar)
%
% Computes the diffusion-weighted signal predicted by the IVIM model for one
% or more b-values:
%
%       S(b) = S0 * [ f * exp(-b * D*) + (1 - f) * exp(-b * D) ]
%
% This is the canonical two-compartment model the pancData3 pipeline fits
% (see pipeline/core/fit_models.m and pipeline/dependencies/IVIM_seg.m). It is
% a complete, self-contained implementation used to synthesise phantom data
% and as the reference signal in the demo's validation.
%
% TERM-BY-TERM PHYSICS
%   S0                 Signal with no diffusion weighting (b=0). Physically the
%                      T2-weighted magnitude; sets overall voxel brightness.
%   f * exp(-b*D*)     PERFUSION compartment. f is the perfusion fraction
%                      (0..1): the share of the b=0 signal coming from blood in
%                      the capillary bed. D* is the pseudo-diffusion
%                      coefficient — the apparent "diffusion" of blood water as
%                      it follows the randomly oriented capillary network.
%                      Because D* >> D, this term decays away by b ~ 150 s/mm^2,
%                      so only LOW b-values constrain (f, D*).
%   (1-f)*exp(-b*D)    TISSUE-DIFFUSION compartment. D is the true diffusion
%                      coefficient — Brownian motion of water hindered by cell
%                      membranes. It survives to high b and reports on
%                      cellularity (low D = densely packed cells).
%
% WHY THE MODEL IS SEPARABLE
%   D* is an order of magnitude larger than D, so the two exponentials live on
%   well-separated b-value scales. That is exactly what makes the segmented
%   fit possible: estimate D from the high-b tail (perfusion already gone),
%   then peel it off to recover f and D*.
%
% INPUTS
%   b      vector of b-values (s/mm^2). Any shape; output matches its shape.
%   S0     scalar signal at b=0 (arbitrary units).
%   D      true diffusion coefficient (mm^2/s),  e.g. ~1.1e-3 in pancreatic tumour.
%   f      perfusion fraction (dimensionless 0..1), e.g. ~0.12.
%   Dstar  pseudo-diffusion coefficient (mm^2/s), e.g. ~20e-3.
%
% OUTPUT
%   S      signal at each b, same shape as b.
%
% Vectorised over b; (S0, D, f, Dstar) are scalars for a single voxel.

    if nargin ~= 5
        error('ivim_signal:nargin', 'Usage: S = ivim_signal(b, S0, D, f, Dstar)');
    end
    if f < 0 || f > 1
        error('ivim_signal:fRange', 'Perfusion fraction f must lie in [0, 1] (got %g).', f);
    end
    S = S0 .* ( f .* exp(-b .* Dstar) + (1 - f) .* exp(-b .* D) );
end
