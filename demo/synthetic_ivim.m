function cohort = synthetic_ivim(cfg)
% SYNTHETIC_IVIM  Generate a synthetic pancreatic-DWI cohort from known IVIM truth.
% SYNTHETIC PHANTOM DATA — not clinical. Contains NO patient data.
% =========================================================================
%
% Produces a fully-specified, reproducible phantom cohort whose per-voxel
% IVIM ground truth (D, f, D*, ADC) is known exactly, so the downstream demo
% can feed the phantoms through the real pipeline fitter (fit_models.m) and
% score how well it recovers the numbers we put in. Every "patient" is
% simulated — there is no PHI anywhere in this module.
%
% The physics is implemented in three reusable, separately-tested functions:
%     ivim_signal.m       bi-exponential forward model S(b)
%     add_rician_noise.m  magnitude-MR (Rician) noise with the noise floor
%     adc_from_signal.m   weighted log-linear ADC (matches the pipeline)
% This file composes them into a cohort with realistic spatial structure
% (a GTV with a resistant low-ADC core) and a longitudinal treatment-response
% model.
%
% OUTPUT  `cohort` struct with, for every phantom patient x fraction:
%     .patients(p).fraction(k).dwi          4-D volume (x,y,z,b) for fit_models
%     .patients(p).fraction(k).mask         GTV logical mask
%     .patients(p).fraction(k).D_true/f_true/Dstar_true/ADC_true   truth maps
%     .patients(p).fraction(k).region       1 = ordinary tumour, 2 = resistant core
%     .patients(p).responder / .is_lf       simulated outcome
%     .gtv_mask, .region, .is_lf, .meta
%     .ground_truth_table / .ground_truth_columns   flat answer key
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
% A FIXED SEED is required for a reproducible phantom: the published ground
% truth and every figure must be byte-for-byte reproducible so the data is
% provably synthetic. rng() with a named generator gives identical draws on
% MATLAB and Octave.
if nargin < 1 || isempty(cfg); cfg = struct(); end
cfg = set_default(cfg, 'seed',           20260602);
cfg = set_default(cfg, 'n_patients',     24);
cfg = set_default(cfg, 'n_fractions',    5);
cfg = set_default(cfg, 'snr',            25);
cfg = set_default(cfg, 'grid',           [24 24]);
cfg = set_default(cfg, 'gtv_radius',     9);
cfg = set_default(cfg, 'resistant_frac', 0.30);
cfg = set_default(cfg, 'bvalues',        [0; 30; 150; 550]);

rng(cfg.seed, 'twister');

b = cfg.bvalues(:);          % force column — fit_models divides by S(:,1)=S(b=0)
if b(1) ~= 0
    error('synthetic_ivim:noB0', 'b-values must start at b=0 (the S0 reference).');
end

% -------------------------------------------------------------------------
% 1. GROUND-TRUTH PARAMETER RANGES  (pancreatic tumour IVIM, units mm^2/s)
% -------------------------------------------------------------------------
% Grounded in this repo's own conventions (config.example.json thresholds and
% the IVIMmodelfit physiological bounds lim = [0 0 0 0; 3e-3 2*max(S) 0.4 0.1]):
%
%   D   (true tissue diffusion)  ~ 1.1e-3   bound [0, 3e-3]
%       Cellularity surrogate; LOW D = densely packed viable tumour.
%   f   (perfusion fraction)     ~ 0.12     bound [0, 0.4]   (dimensionless)
%       Pancreatic adenocarcinoma is hypovascular (dense desmoplastic stroma)
%       => LOW f.
%   D*  (pseudo-diffusion)       ~ 20e-3    bound [0, 0.1]
%       D* >> D (separable model); the LEAST stable parameter (tiny,
%       fast-decaying signal swamped by noise).
%   S0  (b=0 signal)             = 1000 a.u.   SNR below is referenced to this.
%
% Per-PATIENT baseline means are drawn from these distributions; per-VOXEL
% values then scatter around the patient mean (intratumoural heterogeneity).
% Clamps keep every draw inside the physiological bounds the fitter enforces.
gt.D_mean      = 1.10e-3;  gt.D_sd_pat   = 0.10e-3;  gt.D_clamp   = [0.70e-3 1.60e-3];
gt.f_mean      = 0.12;     gt.f_sd_pat   = 0.030;    gt.f_clamp   = [0.04 0.30];
gt.Dstar_mean  = 20e-3;    gt.Dstar_sd_pat = 5e-3;   gt.Dstar_clamp = [8e-3 45e-3];
gt.S0          = 1000;

% Per-voxel (within-tumour) spread around the patient/region mean.
gt.D_sd_vox     = 0.12e-3;
gt.f_sd_vox     = 0.025;
gt.Dstar_sd_vox = 6e-3;     % D* scatters most

% RESISTANT SUB-VOLUME multipliers. Clinical rationale: the treatment-resistant
% core is the densest, most hypoxic tissue — dense packing => more restricted
% diffusion => LOWER D and ADC; hypovascular => LOWER f. The pipeline's
% sub-volume/thresholding logic exists to find this low-ADC core, so the
% phantom must contain one.
gt.D_resist_mult = 0.70;    % ~0.8e-3 mm^2/s core diffusion
gt.f_resist_mult = 0.60;    % lower perfusion in the dense core

% LONGITUDINAL RESPONSE MODEL (per-fraction drift of the tumour means).
% Over a course of radiotherapy a RESPONDING tumour changes; a NON-RESPONDER
% stays roughly fixed. We encode the simulated outcome (local control vs local
% failure) as this divergence so the headline trajectories separate the groups
% on data where we know the answer.
%   D  : responders' diffusion RISES — cells die, packing loosens, water moves
%        more freely. The strongest, most clinically trusted response signal.
%   f  : responders' perfusion fraction RISES modestly — as dense desmoplastic
%        stroma loosens and microvasculature normalises, the perfused fraction
%        grows. (Direction is a modelling choice; the response literature is
%        mixed, hence a small effect.)
%   D* : a small upward drift is encoded for completeness, but note D* is
%        nearly unrecoverable at clinical SNR — the demo's validation shows the
%        fitter will NOT reliably track this trajectory, which is itself the
%        lesson about D* instability.
gt.D_response_per_fx     = 0.060e-3;   % responders: +0.06e-3 mm^2/s per fraction
gt.f_response_per_fx     = 0.010;      % responders: +0.010 per fraction
gt.Dstar_response_per_fx = 0.5e-3;     % responders: +0.5e-3 mm^2/s per fraction
gt.lf_prevalence         = 0.40;       % ~40% local-failure (non-responder) cohort

% -------------------------------------------------------------------------
% 2. Assemble the cohort
% -------------------------------------------------------------------------
nx = cfg.grid(1); ny = cfg.grid(2); nb = numel(b);

% Build the GTV mask once: a centred disk, so fit_models exercises its masking
% path just like a physician-drawn contour selecting a few hundred tumour
% voxels out of the full image.
[XX, YY] = meshgrid(1:ny, 1:nx);
cx = (nx + 1) / 2; cy = (ny + 1) / 2;
gtv_mask = sqrt((XX - cy).^2 + (YY - cx).^2) <= cfg.gtv_radius;   % logical nx-by-ny
n_gtv = sum(gtv_mask(:));

% IVIM_seg has an orientation ambiguity when #voxels == #b-values, so the GTV
% must hold many more voxels than b-values (it does: ~250 vs 4).
if n_gtv <= nb + 4
    error('synthetic_ivim:gtvTooSmall', ...
        'GTV has only %d voxels; need many more than %d b-values for a stable fit. Increase gtv_radius.', n_gtv, nb);
end

% Designate the resistant sub-volume: the innermost (smallest-radius) GTV
% voxels — the tumour core. Picking by radius makes it spatially contiguous,
% like a real dense core.
gtv_lin = find(gtv_mask);
rad = sqrt((XX(gtv_lin) - cy).^2 + (YY(gtv_lin) - cx).^2);
[~, order] = sort(rad, 'ascend');
n_resist = round(cfg.resistant_frac * n_gtv);
resist_lin = gtv_lin(order(1:n_resist));
region = ones(nx, ny);                            % 1 = ordinary tumour
region(resist_lin) = 2;                           % 2 = resistant core

% Assign patient outcomes deterministically given the seed.
is_lf = rand(cfg.n_patients, 1) < gt.lf_prevalence;   % true => local failure

cohort = struct();
cohort.meta = struct('seed', cfg.seed, 'snr', cfg.snr, 'bvalues', b, ...
    'n_patients', cfg.n_patients, 'n_fractions', cfg.n_fractions, ...
    'grid', cfg.grid, 'gtv_radius', cfg.gtv_radius, 'n_gtv_voxels', n_gtv, ...
    'n_resistant_voxels', n_resist, 'ground_truth', gt, ...
    'synthetic', true, 'note', 'SYNTHETIC PHANTOM — no patient data');
cohort.gtv_mask = gtv_mask;
cohort.region   = region;
cohort.is_lf    = is_lf;

% Flat ground-truth table (one row per patient x fraction x region) for the
% published answer key.
gt_rows = {};

for p = 1:cfg.n_patients
    % This patient's baseline (Fx1) tumour means, clamped to physiology.
    Dp     = clamp(gt.D_mean     + gt.D_sd_pat     * randn(), gt.D_clamp);
    fp     = clamp(gt.f_mean     + gt.f_sd_pat     * randn(), gt.f_clamp);
    Dstarp = clamp(gt.Dstar_mean + gt.Dstar_sd_pat * randn(), gt.Dstar_clamp);
    responder = ~is_lf(p);

    for k = 1:cfg.n_fractions
        % Longitudinal drift of each tumour mean for this fraction. (k-1) so
        % Fx1 is the untreated baseline; drift applies only to responders.
        rr = double(responder);
        D_fx_mean     = clamp(Dp     + gt.D_response_per_fx     * (k-1) * rr, gt.D_clamp);
        f_fx_mean     = clamp(fp     + gt.f_response_per_fx     * (k-1) * rr, gt.f_clamp);
        Dstar_fx_mean = clamp(Dstarp + gt.Dstar_response_per_fx * (k-1) * rr, gt.Dstar_clamp);

        % --- Per-voxel ground-truth maps over the GTV ---------------------
        Dtrue     = nan(nx, ny);
        ftrue     = nan(nx, ny);
        Dstartrue = nan(nx, ny);
        for idx = gtv_lin'
            if region(idx) == 2     % resistant core: lower D and f
                mu_D = D_fx_mean * gt.D_resist_mult;
                mu_f = f_fx_mean * gt.f_resist_mult;
            else                    % ordinary tumour
                mu_D = D_fx_mean;
                mu_f = f_fx_mean;
            end
            Dtrue(idx)     = clamp(mu_D          + gt.D_sd_vox     * randn(), gt.D_clamp);
            ftrue(idx)     = clamp(mu_f          + gt.f_sd_vox     * randn(), gt.f_clamp);
            Dstartrue(idx) = clamp(Dstar_fx_mean + gt.Dstar_sd_vox * randn(), gt.Dstar_clamp);
        end

        % --- Forward model: synthesise the noiseless multi-b signal -------
        dwi_clean = zeros(nx, ny, 1, nb);
        for vi = gtv_lin'
            [ix, iy] = ind2sub([nx ny], vi);
            S = ivim_signal(b, gt.S0, Dtrue(vi), ftrue(vi), Dstartrue(vi));
            dwi_clean(ix, iy, 1, :) = reshape(S, [1 1 1 nb]);
        end
        % Background (outside GTV): a low DC level so the volume looks like a
        % real image. Never fitted (the mask excludes it).
        bg = 0.02 * gt.S0;
        for kk = 1:nb
            slice = dwi_clean(:, :, 1, kk);
            slice(~gtv_mask) = bg;
            dwi_clean(:, :, 1, kk) = slice;
        end

        % --- Rician noise -------------------------------------------------
        % SNR referenced to the b=0 signal (sigma = S0/SNR). add_rician_noise
        % documents why magnitude-MR noise is Rician and why the noise floor
        % is what destabilises the D*/f fit at high b.
        sigma = gt.S0 / cfg.snr;
        dwi_noisy = add_rician_noise(dwi_clean, sigma);

        % --- Ground-truth ADC --------------------------------------------
        % ADC is emergent, not an input: the "true" ADC is what a clean
        % mono-exponential fit of the NOISELESS signal returns, using the same
        % estimator as the pipeline (adc_from_signal).
        ADCtrue = nan(nx, ny);
        for vi = gtv_lin'
            [ix, iy] = ind2sub([nx ny], vi);
            S = squeeze(dwi_clean(ix, iy, 1, :));
            ADCtrue(vi) = adc_from_signal(b, S);
        end

        % --- Stash this scan ---------------------------------------------
        cohort.patients(p).responder           = responder;
        cohort.patients(p).is_lf               = is_lf(p);
        cohort.patients(p).fraction(k).dwi     = dwi_noisy;
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

cohort.ground_truth_table = gt_rows;
cohort.ground_truth_columns = {'patient','fraction','responder','region', ...
    'D_true','f_true','Dstar_true','ADC_true'};
end

% =========================================================================
% LOCAL HELPERS
% =========================================================================
function v = clamp(v, lim)
% CLAMP  Keep a draw inside [lo hi] so synthetic truth never violates the
% physiological bounds the fitter enforces.
v = min(max(v, lim(1)), lim(2));
end

function s = set_default(s, field, val)
if ~isfield(s, field) || isempty(s.(field)); s.(field) = val; end
end
