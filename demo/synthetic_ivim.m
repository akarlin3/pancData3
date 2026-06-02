function cohort = synthetic_ivim(cfg)
% TEACHING REFERENCE — synthetic phantom data, not clinical.
% =========================================================================
% synthetic_ivim.m  —  IVIM forward model + phantom cohort generator
% =========================================================================
%
% ##  READ THIS FIRST  ######################################################
% #                                                                         #
% #  This module is SCAFFOLDING the user intends to rewrite by hand to own  #
% #  the physics. It is deliberately *over-commented*: every equation,      #
% #  parameter range, and noise choice is annotated with the physical and   #
% #  clinical reasoning so the whole module can be re-derived from the      #
% #  comments alone. Do not treat it as a finished/owned implementation.    #
% #                                                                         #
% #  It generates PURELY SYNTHETIC pancreatic-DWI phantoms from a KNOWN     #
% #  IVIM ground truth. There is NO patient data here — every "patient" is  #
% #  simulated, every voxel's true (D, f, D*) is something we chose. That   #
% #  is the point: the downstream demo feeds these phantoms through the     #
% #  *real* pipeline fitter (fit_models.m) and checks whether it recovers   #
% #  the numbers we put in. A reviewer can run the whole thing with zero    #
% #  PHI risk.                                                              #
% #                                                                         #
% ###########################################################################
%
% WHAT THIS FILE PRODUCES
%   A `cohort` struct holding, for every phantom patient x treatment
%   fraction:
%     - a 4-D DWI volume  dwi(x, y, z, b)  in the *exact* shape the real
%       pipeline's fit_models() consumes (see CHECKPOINT 0 enumeration),
%     - the GTV mask that selects the tumour voxels,
%     - the per-voxel GROUND-TRUTH (D, f, D*, ADC) we used to synthesise it,
%     - the region label (1 = ordinary tumour, 2 = resistant sub-volume),
%     - the patient-level outcome label (local control vs local failure).
%
% USAGE
%   cohort = synthetic_ivim();                 % all defaults
%   cohort = synthetic_ivim(struct('snr',20)); % override any field below
%
% CONFIG FIELDS (all optional; sensible defaults provided)
%   .seed              RNG seed for reproducibility            (default 20260602)
%   .n_patients        number of phantom patients              (default 24)
%   .n_fractions       DWI timepoints (Fx1..FxN)               (default 5)
%   .snr               signal-to-noise ratio at b=0            (default 25)
%   .grid              [nx ny] phantom slice size              (default [24 24])
%   .gtv_radius        GTV disk radius in voxels               (default 9)
%   .resistant_frac    fraction of GTV that is resistant core  (default 0.30)
%   .bvalues           column vector of b-values (s/mm^2)      (default [0;30;150;550])
%
% =========================================================================

% -------------------------------------------------------------------------
% 0. Configuration & reproducibility
% -------------------------------------------------------------------------
% A FIXED SEED is non-negotiable for a teaching/demo asset: the published
% ground-truth parameters and every figure must be byte-for-byte
% reproducible so a reviewer can confirm the phantom is exactly that — a
% phantom — and not a back-door for real data.
if nargin < 1 || isempty(cfg); cfg = struct(); end
cfg = set_default(cfg, 'seed',           20260602);
cfg = set_default(cfg, 'n_patients',     24);
cfg = set_default(cfg, 'n_fractions',    5);
cfg = set_default(cfg, 'snr',            25);
cfg = set_default(cfg, 'grid',           [24 24]);
cfg = set_default(cfg, 'gtv_radius',     9);
cfg = set_default(cfg, 'resistant_frac', 0.30);
cfg = set_default(cfg, 'bvalues',        [0; 30; 150; 550]);

% rng() with a named generator gives identical draws across MATLAB and
% Octave. We draw the GROUND TRUTH before any noise, so sweeping SNR with a
% fixed seed leaves the truth unchanged and isolates the noise effect.
rng(cfg.seed, 'twister');

b = cfg.bvalues(:);          % force column — fit_models divides by S(:,1)=S(b=0)
if b(1) ~= 0
    error('synthetic_ivim:noB0', 'b-values must start at b=0 (the S0 reference).');
end

% -------------------------------------------------------------------------
% 1. GROUND-TRUTH PARAMETER RANGES  (pancreatic tumour IVIM, units mm^2/s)
% -------------------------------------------------------------------------
% These ranges are grounded in this repo's own conventions
% (config.example.json thresholds and the IVIMmodelfit physiological bounds
% lim = [0 0 0 0; 3e-3 2*max(S) 0.4 0.1]):
%
%   D   (true tissue diffusion)  ~ 1.1e-3   bound [0, 3e-3]
%       Brownian motion of water hindered by cell membranes. LOW D = densely
%       packed viable tumour cells (restricted diffusion). This is the
%       cellularity surrogate and the most clinically trusted IVIM parameter.
%
%   f   (perfusion fraction)     ~ 0.12     bound [0, 0.4]   (dimensionless)
%       Fraction of the b=0 signal coming from blood in the capillary bed
%       rather than tissue water. Pancreatic adenocarcinoma is famously
%       hypovascular (dense desmoplastic stroma), hence a LOW f. Above ~0.4
%       a voxel would be modelled as mostly blood, which is unphysical here.
%
%   D*  (pseudo-diffusion)       ~ 20e-3    bound [0, 0.1]
%       "Diffusion" of blood water as it perfuses through the randomly
%       oriented capillary network — an order of magnitude faster than true
%       diffusion (D* >> D), which is exactly WHY the bi-exponential can be
%       separated at all: the D* component decays away by b~150 s/mm^2,
%       leaving only the D component at high b. D* is the LEAST stable
%       parameter (tiny, fast-decaying signal, swamped by noise) — the
%       phantom is built to expose that instability on purpose.
%
%   S0  (b=0 signal)             = 1000 a.u.
%       The T2-weighted signal with NO diffusion weighting. Sets the overall
%       brightness; arbitrary units. SNR below is referenced to this.
%
% Per-PATIENT baseline means are drawn from these distributions; per-VOXEL
% values then scatter around the patient mean to mimic intratumoural
% heterogeneity. Clamps keep every draw inside the physiological bounds the
% real fitter enforces.
gt.D_mean      = 1.10e-3;  gt.D_sd_pat   = 0.10e-3;  gt.D_clamp   = [0.70e-3 1.60e-3];
gt.f_mean      = 0.12;     gt.f_sd_pat   = 0.030;    gt.f_clamp   = [0.04 0.25];
gt.Dstar_mean  = 20e-3;    gt.Dstar_sd_pat = 5e-3;   gt.Dstar_clamp = [8e-3 40e-3];
gt.S0          = 1000;

% Per-voxel (within-tumour) spread around the patient/region mean.
gt.D_sd_vox     = 0.12e-3;
gt.f_sd_vox     = 0.025;
gt.Dstar_sd_vox = 6e-3;     % D* scatters most — deliberately wide

% RESISTANT SUB-VOLUME multipliers.
% Clinical rationale: the treatment-resistant core of a pancreatic tumour is
% the densest, most hypoxic tissue. Dense cell packing => more restricted
% diffusion => LOWER D (and lower ADC). Hypovascular => LOWER f. The
% downstream pipeline's sub-volume/thresholding logic exists precisely to
% find this low-ADC core, so the phantom must contain one for that logic to
% have something to detect.
gt.D_resist_mult = 0.70;    % ~0.8e-3 mm^2/s core diffusion
gt.f_resist_mult = 0.60;    % lower perfusion in the dense core

% LONGITUDINAL response model.
% Over a course of radiotherapy, RESPONDING tumour loses cellularity (cells
% die, packing loosens) => water diffuses more freely => D RISES fraction to
% fraction. NON-RESPONDERS stay densely packed => D ~ flat. We encode the
% patient outcome (local control vs local failure) as this divergence, so
% the headline longitudinal figure shows the two groups separating — on data
% where we *know* the answer.
gt.D_response_per_fx = 0.060e-3;   % responders: +0.06e-3 mm^2/s per fraction
gt.lf_prevalence     = 0.40;       % ~40% local-failure (non-responder) cohort

% -------------------------------------------------------------------------
% 2. Assemble the cohort
% -------------------------------------------------------------------------
nx = cfg.grid(1); ny = cfg.grid(2); nb = numel(b);

% Build the GTV mask once: a centred disk. A disk (not the whole frame)
% means fit_models actually exercises its masking path, just like a real
% physician-drawn contour selecting a few hundred tumour voxels out of the
% full image.
[XX, YY] = meshgrid(1:ny, 1:nx);
cx = (nx + 1) / 2; cy = (ny + 1) / 2;
gtv_mask = sqrt((XX - cy).^2 + (YY - cx).^2) <= cfg.gtv_radius;   % logical nx-by-ny
n_gtv = sum(gtv_mask(:));

% Sanity: IVIM_seg has an orientation ambiguity when #voxels == #b-values,
% so the GTV must hold many more voxels than b-values (it does: ~250 vs 4).
if n_gtv <= nb + 4
    error('synthetic_ivim:gtvTooSmall', ...
        'GTV has only %d voxels; need many more than %d b-values for a stable fit. Increase gtv_radius.', n_gtv, nb);
end

% Designate the resistant sub-volume: the innermost voxels (smallest radius)
% of the GTV, i.e. the tumour core. Picking by radius (not at random) makes
% the resistant region spatially contiguous, like a real dense core.
gtv_lin = find(gtv_mask);
rad = sqrt((XX(gtv_lin) - cy).^2 + (YY(gtv_lin) - cx).^2);
[~, order] = sort(rad, 'ascend');
n_resist = round(cfg.resistant_frac * n_gtv);
resist_lin = gtv_lin(order(1:n_resist));         % linear idx of core voxels
region = ones(nx, ny);                            % 1 = ordinary tumour
region(resist_lin) = 2;                           % 2 = resistant core

% Assign patient outcomes deterministically given the seed.
is_lf = rand(cfg.n_patients, 1) < gt.lf_prevalence;   % true => local failure (non-responder)

cohort = struct();
cohort.meta = struct('seed', cfg.seed, 'snr', cfg.snr, 'bvalues', b, ...
    'n_patients', cfg.n_patients, 'n_fractions', cfg.n_fractions, ...
    'grid', cfg.grid, 'gtv_radius', cfg.gtv_radius, 'n_gtv_voxels', n_gtv, ...
    'n_resistant_voxels', n_resist, 'ground_truth', gt, ...
    'synthetic', true, 'note', 'SYNTHETIC PHANTOM — no patient data');
cohort.gtv_mask = gtv_mask;
cohort.region   = region;
cohort.is_lf    = is_lf;

% Pre-allocate a flat ground-truth table (one row per patient x fraction x
% region) for saving alongside the figures, so the "answer key" is published.
gt_rows = {};

for p = 1:cfg.n_patients
    % Draw this patient's baseline (Fx1) tumour means, clamped to physiology.
    Dp     = clamp(gt.D_mean     + gt.D_sd_pat     * randn(), gt.D_clamp);
    fp     = clamp(gt.f_mean     + gt.f_sd_pat     * randn(), gt.f_clamp);
    Dstarp = clamp(gt.Dstar_mean + gt.Dstar_sd_pat * randn(), gt.Dstar_clamp);
    responder = ~is_lf(p);

    for k = 1:cfg.n_fractions
        % Longitudinal drift of the diffusion mean for this fraction.
        % (k-1) so Fx1 is the untreated baseline.
        D_drift = gt.D_response_per_fx * (k - 1) * double(responder);
        D_fx_mean = clamp(Dp + D_drift, gt.D_clamp);

        % --- Build per-voxel ground-truth maps over the GTV ---------------
        Dtrue     = nan(nx, ny);
        ftrue     = nan(nx, ny);
        Dstartrue = nan(nx, ny);
        for idx = gtv_lin'
            if region(idx) == 2     % resistant core: lower D and f
                mu_D = D_fx_mean * gt.D_resist_mult;
                mu_f = fp        * gt.f_resist_mult;
            else                    % ordinary tumour
                mu_D = D_fx_mean;
                mu_f = fp;
            end
            Dtrue(idx)     = clamp(mu_D        + gt.D_sd_vox     * randn(), gt.D_clamp);
            ftrue(idx)     = clamp(mu_f        + gt.f_sd_vox     * randn(), gt.f_clamp);
            Dstartrue(idx) = clamp(Dstarp      + gt.Dstar_sd_vox * randn(), gt.Dstar_clamp);
        end

        % --- Forward model: synthesise the noiseless multi-b signal -------
        % S(b) = S0 * [ f*exp(-b*D*) + (1-f)*exp(-b*D) ]  for every GTV voxel.
        dwi_clean = zeros(nx, ny, 1, nb);
        for vi = gtv_lin'
            [ix, iy] = ind2sub([nx ny], vi);
            S = ivim_forward(b, gt.S0, Dtrue(vi), ftrue(vi), Dstartrue(vi));
            dwi_clean(ix, iy, 1, :) = reshape(S, [1 1 1 nb]);
        end
        % Background (outside GTV): a low DC level so the volume looks like a
        % real image. It is never fitted (the mask excludes it) but keeps the
        % array shape honest.
        bg = 0.02 * gt.S0;
        for kk = 1:nb
            slice = dwi_clean(:, :, 1, kk);
            slice(~gtv_mask) = bg;
            dwi_clean(:, :, 1, kk) = slice;
        end

        % --- Rician noise -------------------------------------------------
        % MR magnitude images are |complex Gaussian|. Real and imaginary
        % channels each carry independent zero-mean Gaussian noise with the
        % same sigma; the scanner reports the MAGNITUDE. So the *correct*
        % noise model is Rician, NOT additive Gaussian. The distinction
        % matters most at LOW signal (high b-values, and the resistant core):
        % there the magnitude operation rectifies noise to a positive bias
        % (the "noise floor"), which is precisely what destabilises the D*/f
        % estimate — the signal that should decay toward zero instead piles
        % up on a noise pedestal, and the fitter mistakes that pedestal for
        % perfusion. This is the heart of IVIM ill-conditioning, and the
        % phantom reproduces it faithfully.
        sigma = gt.S0 / cfg.snr;            % SNR referenced to the b=0 signal
        dwi_noisy = add_rician_noise(dwi_clean, sigma);

        % --- Ground-truth ADC --------------------------------------------
        % ADC is not an independent input — it is an EMERGENT property of the
        % (D, f, D*) signal. The "true" ADC of a voxel is what a clean
        % mono-exponential fit of its NOISELESS signal would return, using
        % the same weighted-log-linear estimator the pipeline uses. We
        % compute it here so the validation can score ADC recovery too.
        ADCtrue = nan(nx, ny);
        for vi = gtv_lin'
            [ix, iy] = ind2sub([nx ny], vi);
            S = squeeze(dwi_clean(ix, iy, 1, :));
            ADCtrue(vi) = adc_truth_from_signal(b, S);
        end

        % --- Stash this scan ---------------------------------------------
        cohort.patients(p).responder           = responder;
        cohort.patients(p).is_lf               = is_lf(p);
        cohort.patients(p).fraction(k).dwi     = dwi_noisy;     % feeds fit_models
        cohort.patients(p).fraction(k).mask    = logical(gtv_mask);
        cohort.patients(p).fraction(k).D_true  = Dtrue;
        cohort.patients(p).fraction(k).f_true  = ftrue;
        cohort.patients(p).fraction(k).Dstar_true = Dstartrue;
        cohort.patients(p).fraction(k).ADC_true   = ADCtrue;
        cohort.patients(p).fraction(k).region     = region;

        % Append region-mean rows to the flat answer key.
        for rlabel = 1:2
            rmask = (region == rlabel) & gtv_mask;
            if ~any(rmask(:)); continue; end
            gt_rows(end+1, :) = { p, k, responder, rlabel, ...
                mean(Dtrue(rmask)), mean(ftrue(rmask)), ...
                mean(Dstartrue(rmask)), mean(ADCtrue(rmask)) }; %#ok<AGROW>
        end
    end
end

cohort.ground_truth_table = gt_rows;   % cols: pat fx responder region D f Dstar ADC
cohort.ground_truth_columns = {'patient','fraction','responder','region', ...
    'D_true','f_true','Dstar_true','ADC_true'};
end

% =========================================================================
% LOCAL FUNCTIONS  (the physics primitives — study these)
% =========================================================================

function S = ivim_forward(b, S0, D, f, Dstar)
% IVIM_FORWARD  Bi-exponential intravoxel-incoherent-motion signal model.
%
%   S(b) = S0 * [ f * exp(-b * D*) + (1 - f) * exp(-b * D) ]
%
% Term by term:
%   S0                 the undiffused (b=0) signal — overall voxel brightness.
%   f * exp(-b*D*)     the PERFUSION compartment. Weight f (perfusion
%                      fraction); decays at the fast pseudo-diffusion rate D*.
%                      Because D* is large, this term is essentially gone by
%                      b ~ 150 s/mm^2 — only the low-b points constrain it.
%   (1-f)*exp(-b*D)    the TISSUE-DIFFUSION compartment. Weight (1-f); decays
%                      at the slow true-diffusion rate D. This is what
%                      survives to high b and reports on cellularity.
%
% The whole reason the segmented fit works: D* >> D means the two
% exponentials live on separated b-value scales, so you can estimate D from
% the high-b tail first, then peel it off to get f and D*.
S = S0 .* ( f .* exp(-b .* Dstar) + (1 - f) .* exp(-b .* D) );
end

function out = add_rician_noise(S, sigma)
% ADD_RICIAN_NOISE  Corrupt a clean magnitude signal with Rician noise.
%
% Physics: the scanner measures a complex signal S + (eps_r + i*eps_i), with
% eps_r, eps_i ~ N(0, sigma) independent, then takes the MAGNITUDE:
%
%   S_meas = sqrt( (S + eps_r)^2 + eps_i^2 )
%
% That magnitude is Rician-distributed. At high SNR it looks Gaussian about
% S; at low SNR (S -> 0) it collapses to a Rayleigh distribution with a
% strictly POSITIVE mean ~ sigma*sqrt(pi/2) — the noise floor. Modelling this
% correctly is essential: using plain additive Gaussian noise would let the
% high-b signal go negative and would HIDE the floor-induced bias that makes
% real IVIM D*/f estimation unstable.
eps_r = sigma .* randn(size(S));
eps_i = sigma .* randn(size(S));
out = sqrt((S + eps_r).^2 + eps_i.^2);
end

function adc = adc_truth_from_signal(b, S)
% ADC_TRUTH_FROM_SIGNAL  Weighted log-linear mono-exponential ADC.
%
% Mirrors the pipeline's ADC estimator (fit_models.m) so the "true" ADC is
% defined consistently with what the pipeline recovers. Model S = S0*exp(-b*ADC);
% linearise to ln(S/S0) = -b*ADC and solve by weighted least squares with
% weights w = S^2 (the delta-method correction for log-transformed Gaussian
% noise, which is heteroscedastic: var(ln S) ~ 1/S^2). b=0 is the reference
% only and is excluded from the regression.
S  = S(:);
b  = b(:);
S0 = S(1);
if S0 <= 0 || any(S(2:end) <= 0)
    adc = NaN; return;
end
A = -b(2:end);
y = log(S(2:end) ./ S0);
w = S(2:end).^2;
adc = sum(w .* y .* A) / sum(w .* (A.^2));
if ~isfinite(adc) || adc < 0; adc = NaN; end
end

function v = clamp(v, lim)
% CLAMP  Keep a draw inside [lo hi] so synthetic truth never violates the
% physiological bounds the real fitter enforces.
v = min(max(v, lim(1)), lim(2));
end

function s = set_default(s, field, val)
if ~isfield(s, field) || isempty(s.(field)); s.(field) = val; end
end
