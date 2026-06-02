function run_synthetic_demo(varargin)
% RUN_SYNTHETIC_DEMO  End-to-end synthetic IVIM demo for pancData3 (PHI-free).
% SYNTHETIC PHANTOM DATA — not clinical.
% =========================================================================
%
% One command, no patient data (via the run_demo.m script entry point):
%     >> run('demo/run_demo.m')        % or:  octave demo/run_demo.m
%
% Pipeline of the demo itself:
%   1. synthetic_ivim.m   builds a phantom cohort from a KNOWN IVIM ground
%                         truth (heavily-commented forward model + Rician
%                         noise + a resistant low-ADC sub-volume).
%   2. fit_models.m       the REAL pipeline fitter runs on the phantom scans,
%                         unmodified — segmented IVIM + weighted ADC.
%   3. validate_fitting.m scores recovered (D, f, D*) against the truth across
%                         SNR (the recover-known-truth correctness check).
%   4. Headline figures   longitudinal D and f trajectories, recovered by the
%                         real fitter, split by simulated outcome (local
%                         control vs local failure) — every artefact labelled
%                         SYNTHETIC.
%
% Optional name/value args (passed through to synthetic_ivim):
%   'snr' (default 25), 'n_patients' (24), 'n_fractions' (5), 'seed'.
%
% Outputs land in demo/output/ (figures + CSV ground-truth/recovery tables).
% =========================================================================

    add_demo_paths();
    here = fileparts(mfilename('fullpath'));
    outdir = fullfile(here, 'output');
    if ~exist(outdir, 'dir'); mkdir(outdir); end

    % ---- config ----
    cfg = struct('snr', 25, 'n_patients', 24, 'n_fractions', 5, 'seed', 20260602);
    for i = 1:2:numel(varargin); cfg.(varargin{i}) = varargin{i+1}; end

    fprintf('\n');
    fprintf('############################################################\n');
    fprintf('#  pancData3 SYNTHETIC DEMO  —  NO PHI, every patient fake  #\n');
    fprintf('############################################################\n');
    fprintf('seed=%d  patients=%d  fractions=%d  SNR=%g  b=[%s]\n', ...
        cfg.seed, cfg.n_patients, cfg.n_fractions, cfg.snr, ...
        strtrim(sprintf('%g ', [0;30;150;550])));

    % ---- 1. generate phantom cohort ----
    fprintf('\n[1/4] Generating synthetic phantom cohort...\n');
    cohort = synthetic_ivim(cfg);
    fprintf('      %d patients x %d fractions, GTV=%d voxels (%d resistant core)\n', ...
        cfg.n_patients, cfg.n_fractions, cohort.meta.n_gtv_voxels, cohort.meta.n_resistant_voxels);
    write_ground_truth_csv(cohort, fullfile(outdir, 'synthetic_ground_truth.csv'));

    % ---- 2+3. run REAL fitter over every scan, aggregate longitudinally ----
    fprintf('\n[2/4] Running REAL pipeline fitter (fit_models.m) on every scan...\n');
    b = cohort.meta.bvalues; bthr = 100;
    np = cfg.n_patients; nk = cfg.n_fractions;
    D_rec = nan(np, nk); f_rec = nan(np, nk);   % per-patient/fraction mean recovered
    for p = 1:np
        for k = 1:nk
            sc = cohort.patients(p).fraction(k);
            rec = fit_synthetic_scan(sc.dwi, b, sc.mask, bthr);
            D_rec(p, k) = nanmean(rec.D);
            f_rec(p, k) = nanmean(rec.f);
        end
        fprintf('      patient %2d/%d fitted\n', p, np);
    end
    is_lf = logical(cohort.is_lf);
    write_recovery_csv(D_rec, f_rec, is_lf, fullfile(outdir, 'synthetic_recovered_longitudinal.csv'));

    % ---- 4. headline longitudinal figures (D and f), by outcome ----
    fprintf('\n[3/4] Building headline longitudinal figures...\n');
    fx = 1:nk;
    make_longitudinal_fig(fx, D_rec, is_lf, 'D (true diffusion, mm^2/s)', ...
        'Synthetic longitudinal D by outcome (phantom)', ...
        fullfile(outdir, 'synthetic_longitudinal_D.png'));
    make_longitudinal_fig(fx, f_rec, is_lf, 'f (perfusion fraction)', ...
        'Synthetic longitudinal f by outcome (phantom)', ...
        fullfile(outdir, 'synthetic_longitudinal_f.png'));

    % ---- recover-known-truth validation + its figures ----
    fprintf('\n[4/4] Running fitting validation across SNR...\n');
    validate_fitting([10 20 40 80], outdir);

    % ---- summary ----
    [mD_lc, mD_lf] = group_means(D_rec, is_lf);
    fprintf('\n=== DEMO SUMMARY (synthetic) ===\n');
    fprintf('Local-control (responder) mean recovered D: Fx1=%.3e -> Fx%d=%.3e\n', ...
        mD_lc(1), nk, mD_lc(end));
    fprintf('Local-failure (non-responder) mean recovered D: Fx1=%.3e -> Fx%d=%.3e\n', ...
        mD_lf(1), nk, mD_lf(end));
    fprintf('Responders show rising D (cell kill); non-responders ~flat — as designed.\n');
    fprintf('Artefacts (figures + CSV) in: %s\n', outdir);
    fprintf('Ground-truth answer key: synthetic_ground_truth.csv\n');
    fprintf('================================\n\n');
end

% -------------------------------------------------------------------------
function make_longitudinal_fig(fx, M, is_lf, ylab, ttl, outpath)
% Group-mean +/- SEM trajectory, local control vs local failure.
    [m_lc, ~, se_lc] = stats_group(M(~is_lf, :));
    [m_lf, ~, se_lf] = stats_group(M(is_lf, :));
    spec = struct();
    spec.title  = ttl;
    spec.xlabel = 'treatment fraction (Fx)';
    spec.ylabel = ylab;
    spec.series = struct( ...
        'x', {fx, fx}, ...
        'y', {m_lc, m_lf}, ...
        'yerr', {se_lc, se_lf}, ...
        'label', {'local control (responder)', 'local failure (non-responder)'}, ...
        'kind', {'linepoints', 'linepoints'});
    demo_plot(spec, outpath);
end

function [mu, sd, se] = stats_group(X)
    mu = nanmean(X, 1);
    sd = nanstd(X, 0, 1);
    n  = sum(~isnan(X), 1);
    se = sd ./ sqrt(max(n, 1));
end

function [m_lc, m_lf] = group_means(M, is_lf)
    m_lc = nanmean(M(~is_lf, :), 1);
    m_lf = nanmean(M(is_lf, :), 1);
end

% -------------------------------------------------------------------------
function write_ground_truth_csv(cohort, path)
    fid = fopen(path, 'w');
    fprintf(fid, '# SYNTHETIC PHANTOM GROUND TRUTH — no patient data. seed=%d snr=%g\n', ...
        cohort.meta.seed, cohort.meta.snr);
    fprintf(fid, '%s\n', strjoin(cohort.ground_truth_columns, ','));
    T = cohort.ground_truth_table;
    for r = 1:size(T, 1)
        fprintf(fid, '%d,%d,%d,%d,%.6e,%.6f,%.6e,%.6e\n', ...
            T{r,1}, T{r,2}, T{r,3}, T{r,4}, T{r,5}, T{r,6}, T{r,7}, T{r,8});
    end
    fclose(fid);
end

function write_recovery_csv(D_rec, f_rec, is_lf, path)
    fid = fopen(path, 'w');
    fprintf(fid, '# SYNTHETIC recovered longitudinal means (real fitter). is_lf: 1=local failure\n');
    fprintf(fid, 'patient,is_lf,fraction,D_recovered,f_recovered\n');
    [np, nk] = size(D_rec);
    for p = 1:np
        for k = 1:nk
            fprintf(fid, '%d,%d,%d,%.6e,%.6f\n', p, is_lf(p), k, D_rec(p,k), f_rec(p,k));
        end
    end
    fclose(fid);
end
