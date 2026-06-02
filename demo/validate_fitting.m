function results = validate_fitting(snr_list, outdir)
% VALIDATE_FITTING  Recover-known-truth check for the IVIM pipeline fitter.
% TEACHING REFERENCE — synthetic phantom data, not clinical.
%
% Generates phantom scans with a KNOWN (D, f, D*) ground truth across a range
% of SNRs, runs the REAL pipeline fitter (fit_models.m) on them, and scores
% how well it recovers each parameter. This is simultaneously:
%   (a) a demo asset — it visualises where IVIM fitting is trustworthy, and
%   (b) a correctness check on the fitting code — on data where we know the
%       answer, D should recover tightly while f and especially D* degrade as
%       SNR falls. That degradation is the documented IVIM ill-conditioning,
%       now demonstrated quantitatively rather than asserted.
%
% Usage:
%   results = validate_fitting();                 % default SNR sweep
%   results = validate_fitting([10 20 40 80], 'demo/output');
%
% Returns a struct array `results`, one entry per SNR, with per-parameter
% recovery metrics. Also writes recovery figures to OUTDIR.

    if nargin < 1 || isempty(snr_list); snr_list = [10 20 40 80]; end
    if nargin < 2 || isempty(outdir);   outdir = fullfile('demo', 'output'); end
    if ~exist(outdir, 'dir'); mkdir(outdir); end

    params = {'D', 'f', 'Dstar', 'ADC'};
    fprintf('\n=== CHECKPOINT 2: FITTING VALIDATION (recover-known-truth) ===\n');
    fprintf('Real fitter: pipeline/core/fit_models.m (unmodified)\n');
    fprintf('Phantom: fixed seed, baseline fraction only, %d SNR levels\n\n', numel(snr_list));

    % Header for the recovery table.
    fprintf('%-6s %-7s %10s %12s %10s %10s\n', ...
        'SNR', 'param', 'median%err', 'bias', 'corr', 'fit%ok');
    fprintf('%s\n', repmat('-', 1, 62));

    results = struct('snr', {}, 'metrics', {});
    % Keep per-voxel D recovery at the extreme SNRs for the scatter figure.
    scatter_cache = struct();

    for si = 1:numel(snr_list)
        snr = snr_list(si);
        % Fixed seed => identical ground truth at every SNR, so differences
        % are purely the noise effect. One fraction keeps it a clean
        % noise-vs-recovery study (no longitudinal drift).
        cohort = synthetic_ivim(struct('snr', snr, 'n_fractions', 1, ...
            'n_patients', 8, 'seed', 20260602));
        b = cohort.meta.bvalues;
        bthr = 100;   % matches config.example.json ivim_bthr

        truth = struct('D', [], 'f', [], 'Dstar', [], 'ADC', []);
        recov = struct('D', [], 'f', [], 'Dstar', [], 'ADC', []);
        for p = 1:numel(cohort.patients)
            sc = cohort.patients(p).fraction(1);
            m = sc.mask;
            rec = fit_synthetic_scan(sc.dwi, b, m, bthr);
            truth.D     = [truth.D;     sc.D_true(m)];
            truth.f     = [truth.f;     sc.f_true(m)];
            truth.Dstar = [truth.Dstar; sc.Dstar_true(m)];
            truth.ADC   = [truth.ADC;   sc.ADC_true(m)];
            recov.D     = [recov.D;     rec.D];
            recov.f     = [recov.f;     rec.f];
            recov.Dstar = [recov.Dstar; rec.Dstar];
            recov.ADC   = [recov.ADC;   rec.ADC];
        end

        m_struct = struct();
        for pj = 1:numel(params)
            pn = params{pj};
            t = truth.(pn); r = recov.(pn);
            mm = compute_recovery_metrics(t, r);
            m_struct.(pn) = mm;
            fprintf('%-6g %-7s %9.1f%% %12.3g %10.3f %9.0f%%\n', ...
                snr, pn, 100*mm.median_rel_err, mm.bias, mm.corr, 100*mm.frac_ok);
        end
        fprintf('%s\n', repmat('-', 1, 62));
        results(si).snr = snr; %#ok<AGROW>
        results(si).metrics = m_struct; %#ok<AGROW>

        if snr == min(snr_list); scatter_cache.low  = struct('snr', snr, 't', truth.D, 'r', recov.D); end
        if snr == max(snr_list); scatter_cache.high = struct('snr', snr, 't', truth.D, 'r', recov.D); end
    end

    % --- Figure 1: median relative error vs SNR (the ill-conditioning curve)
    relerr = @(field) arrayfun(@(rr) 100*rr.metrics.(field).median_rel_err, results);
    spec1 = struct();
    spec1.title  = 'Synthetic IVIM recovery error vs SNR (phantom)';
    spec1.xlabel = 'SNR at b=0';
    spec1.ylabel = 'median |relative error| (%)';
    spec1.series = struct( ...
        'x', {[results.snr], [results.snr], [results.snr]}, ...
        'y', {relerr('D'), relerr('f'), relerr('Dstar')}, ...
        'label', {'D (true diffusion)', 'f (perfusion fraction)', 'D* (pseudo-diffusion)'}, ...
        'kind', {'linepoints', 'linepoints', 'linepoints'});
    demo_plot(spec1, fullfile(outdir, 'synthetic_recovery_vs_snr.png'));

    % --- Figure 2: D recovery scatter at the best vs worst SNR
    if isfield(scatter_cache, 'low') && isfield(scatter_cache, 'high')
        spec2 = struct();
        spec2.title  = 'Synthetic D recovery: truth vs fitted (phantom)';
        spec2.xlabel = 'ground-truth D (mm^2/s)';
        spec2.ylabel = 'pipeline-recovered D (mm^2/s)';
        spec2.identity = true;
        lo = scatter_cache.low; hi = scatter_cache.high;
        spec2.series = struct( ...
            'x', {hi.t(:)', lo.t(:)'}, ...
            'y', {hi.r(:)', lo.r(:)'}, ...
            'label', {sprintf('SNR %g', hi.snr), sprintf('SNR %g', lo.snr)}, ...
            'kind', {'points', 'points'});
        demo_plot(spec2, fullfile(outdir, 'synthetic_recovery_scatter.png'));
    end

    fprintf('Recovery figures written to %s\n', outdir);
end

% -------------------------------------------------------------------------
function mm = compute_recovery_metrics(truth, recov)
% Per-parameter recovery summary over voxels where the fit returned a finite
% value. frac_ok exposes how often the fit simply failed (NaN) — which is
% itself a key result: D* fits fail far more often than D at low SNR.
    truth = truth(:); recov = recov(:);
    ok = isfinite(recov) & isfinite(truth) & (truth ~= 0);
    mm = struct('median_rel_err', NaN, 'bias', NaN, 'corr', NaN, 'frac_ok', 0);
    mm.frac_ok = mean(isfinite(recov));
    if sum(ok) < 3; return; end
    t = truth(ok); r = recov(ok);
    mm.median_rel_err = median(abs(r - t) ./ abs(t));
    mm.bias = median(r - t);
    cc = corrcoef(t, r);
    mm.corr = cc(1, 2);
end
