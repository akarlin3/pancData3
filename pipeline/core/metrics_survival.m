function metrics_survival(valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, m_lf, m_total_time, m_total_follow_up_time, nTp, fx_label, dtype_label, m_gtv_vol, output_folder, actual_scan_days, config_struct_in)
% METRICS_SURVIVAL — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 5/5 of the metrics step. Fits a Time-Dependent Cox Proportional Hazards
% model with dynamic covariate updating.
%
% ANALYTICAL OVERVIEW:
%   This module fits a time-dependent Cox Proportional Hazards (Cox PH) model
%   using the Andersen-Gill counting process formulation.  This is the most
%   sophisticated analysis in the pipeline because it addresses several
%   challenges unique to longitudinal imaging biomarker studies:
%
%   1. TIME-VARYING COVARIATES — Standard Cox regression assumes covariates
%      are fixed at baseline.  In this study, diffusion parameters are
%      measured at each treatment fraction and evolve over time.  The
%      counting-process formulation splits each patient's follow-up into
%      intervals [t_start, t_stop], with covariates updated at each scan.
%      This allows the model to use the most recent biomarker values when
%      estimating instantaneous hazard of local failure.
%
%   2. COMPETING RISKS (CAUSE-SPECIFIC HAZARDS) — Pancreatic cancer patients
%      may die from systemic disease, treatment toxicity, or unrelated causes
%      before local failure can be observed.  The Cause-Specific Hazard (CSH)
%      approach treats these competing events as censored observations when
%      modelling the hazard of local failure.  This is preferred over the
%      subdistribution hazard (Fine-Gray) for aetiological questions because
%      CSH estimates the actual biological hazard rate, not the cumulative
%      incidence accounting for competing events.
%
%   3. INVERSE PROBABILITY OF CENSORING WEIGHTING (IPCW) — Administrative
%      censoring (patients lost to follow-up) may be informative if sicker
%      patients are more likely to drop out.  IPCW models the censoring
%      probability as a function of covariates and upweights observations
%      from patients similar to those who were censored, reducing bias.
%      Competing events are NOT included in the IPCW model because they
%      are not uninformative censoring — they are a distinct endpoint.
%
%   4. LANDMARK ANALYSIS — Only intervals AFTER the end of radiotherapy
%      are included in the primary analysis.  Using covariates from scans
%      that occurred after the event they predict would violate the temporal
%      ordering assumption.  The landmark ensures all patients in the
%      analysis were still at risk at the time of the last treatment scan.
%
%   5. IMPUTATION HALF-LIFE SENSITIVITY — The time-dependent panel uses
%      exponential decay imputation for missing inter-scan values (via
%      build_td_panel).  The sensitivity analysis varies the decay half-life
%      to assess robustness of hazard ratio estimates to imputation assumptions.
%
% Inputs:
%   valid_pts         - Logical mask of patients mapped to LF/LC groups
%   ADC_abs, D_abs, f_abs, Dstar_abs - Baseline covariate values matrices
%   m_lf              - Local failure clinical grouping statuses
%   m_total_time      - Time-to-local-failure (events)
%   m_total_follow_up_time - Time-to-censoring for event-free subjects
%   nTp               - Counter of number of timepoints max
%   fx_label          - Fraction labels used in logging
%   dtype_label       - DWI type name used in output
%   m_gtv_vol         - (Optional) GTV volume matrix (patients x fractions)
%   output_folder     - (Optional) Directory for diary output
%   actual_scan_days  - (Optional) Vector of actual scan days from DICOM
%                       headers or clinical records.  When provided,
%                       replaces the default [0,5,10,15,20,90] to avoid
%                       immortal time bias.
%
% Outputs:
%   None. Outputs printed to console (HR and p-value tables).
%

% Handle optional arguments for backward compatibility
if nargin < 12, m_gtv_vol = []; end
if nargin < 13, output_folder = ''; end
if nargin < 14, actual_scan_days = []; end
if nargin < 15 || isempty(config_struct_in)
    config_struct_internal = struct();
else
    config_struct_internal = config_struct_in;
end

fprintf('\n--- TIME-DEPENDENT COX PH MODEL (Counting Process) ---\n');

% Diary: capture console output to output_folder
if ~isempty(output_folder)
    diary_file = fullfile(output_folder, ['metrics_survival_output_' dtype_label '.txt']);
    if exist(diary_file, 'file'), delete(diary_file); end
    diary(diary_file);
end

% Scan day timing for the time-dependent Cox model.
if ~isempty(actual_scan_days)
    td_scan_days = actual_scan_days;
    fprintf('  Using provided scan days [%s].\n', num2str(td_scan_days));
else
    % Default scan days assume 5 on-treatment fractions + 1 post-treatment scan.
    td_scan_days = [0, 5, 10, 15, 20, 90];
    fprintf('  ⚠️  CAUTION: Using default scan days [%s].\n', num2str(td_scan_days));
    fprintf('      Pass actual DICOM-derived scan days via config.json td_scan_days field\n');
    fprintf('      or as the 14th argument to avoid immortal time bias.\n');
end

% Covariates: all four diffusion/IVIM parameters (absolute values across
% all treatment fractions).  These are the time-varying covariates whose
% trajectories are hypothesised to predict local failure hazard.
% Interpretation of hazard ratios:
%   HR(ADC) > 1 → higher ADC associated with INCREASED failure risk
%                 (counterintuitive: high ADC usually = good response)
%   HR(ADC) < 1 → higher ADC associated with DECREASED failure risk
%                 (expected: radiation-induced cell death → higher ADC)
%   HR(D)   — same interpretation as ADC (reflects cellularity)
%   HR(f)   — higher perfusion fraction may indicate better oxygenation
%             and thus better radiation response (HR < 1 expected)
%   HR(D*)  — physiologically noisy, interpretation less reliable
td_feat_arrays = { ADC_abs(valid_pts,:), D_abs(valid_pts,:), ...
                   f_abs(valid_pts,:),   Dstar_abs(valid_pts,:) };
td_feat_names  = {'ADC', 'D', 'f', 'D*'};

% Include baseline GTV volume as a time-constant confounder when available.
% Tumor volume is a known prognostic factor; omitting it risks confounding
% imaging biomarker effect estimates.
has_vol = ~isempty(m_gtv_vol) && any(isfinite(m_gtv_vol(valid_pts, 1)));
if has_vol
    vol_baseline = m_gtv_vol(valid_pts, 1);
    % Replicate baseline volume across timepoints (time-constant covariate)
    vol_rep = repmat(vol_baseline, 1, nTp);
    td_feat_arrays{end+1} = vol_rep;
    td_feat_names{end+1}  = 'GTVvol';
    fprintf('  Including baseline GTV volume as a covariate.\n');
end

td_n_feat      = numel(td_feat_arrays);

td_lf       = m_lf(valid_pts);
td_tot_time = m_total_time(valid_pts);

% Censored and competing-risk patients use follow-up time; events use
% time-to-event.  In a Cause-Specific Hazard (CSH) model, competing risk
% patients (lf==2) are censored at their last follow-up — they did not
% experience the event of interest (local failure).
follow_up_valid = m_total_follow_up_time(valid_pts);
cens_mask_td = (td_lf == 0 | td_lf == 2) & ~isnan(follow_up_valid);
td_tot_time(cens_mask_td) = follow_up_valid(cens_mask_td);

% Build the counting-process panel: each patient's follow-up is split into
% intervals [t_start, t_stop] aligned to scan times.  Covariates are
% constant within each interval and updated at each scan boundary.
% The default half-life of 18 months controls exponential decay imputation
% for covariate values between scans (e.g., if a scan is missing at Fx3,
% the Fx2 value decays towards the population mean with this half-life).
[X_td_def, t_start_td_def, t_stop_td_def, event_td_def, pat_id_td_def, frac_td_def] = ...
    build_td_panel(td_feat_arrays, td_feat_names, td_lf, td_tot_time, nTp, td_scan_days, 18);

% NOTE: Scaling is deferred until after landmark subsetting (below) to
% prevent pre-landmark covariate distributions from contaminating the
% statistics used to standardise the post-landmark analysis set.
td_ok = (sum(event_td_def == 1) >= 3) && (size(X_td_def, 1) > td_n_feat + 1);

% Sensitivity analysis: vary the imputation half-life to assess robustness.
% Short half-lives (3 months) assume biomarker values revert quickly to
% the population mean between scans (conservative — discounts old data).
% Long half-lives (24 months) assume biomarker values persist over time
% (aggressive — may carry forward outdated measurements).
% If hazard ratios are stable across the grid, the results are robust to
% imputation assumptions.  Large variation indicates sensitivity to
% missing data handling and warrants caution in interpretation.
half_life_grid = [3, 6, 12, 18, 24];
n_halflife = length(half_life_grid);
td_panels = cell(n_halflife, 1);
for hl_idx = 1:n_halflife
    text_progress_bar(hl_idx, n_halflife, 'Half-life sensitivity');
    [X_i, t_start_i, t_stop_i, ev_i, pid_i, frac_i] = ...
        build_td_panel(td_feat_arrays, td_feat_names, td_lf, td_tot_time, nTp, td_scan_days, half_life_grid(hl_idx));
    td_panels{hl_idx} = struct('X', X_i, 't_start', t_start_i, 't_stop', t_stop_i, 'event', ev_i, 'pat_id', pid_i, 'frac', frac_i);
end

if ~td_ok
    fprintf('  Insufficient events (%d) or intervals for time-dependent Cox model.\n', sum(event_td_def == 1));
    if ~isempty(output_folder), diary off; end
    return;
end

% Copy the default half-life (18-month) panel into working variables for
% landmark subsetting and Cox fitting below.
X_td = X_td_def; t_start_td = t_start_td_def; t_stop_td = t_stop_td_def; event_td = event_td_def; pat_id_td = pat_id_td_def; frac_td = frac_td_def;

% ---- Landmark analysis: discard intervals before end-of-RT -----------
% Using covariates from scans that occurred AFTER the event they predict
% violates the time-dependent Cox assumption.  Landmark subsetting keeps
% only patients still at risk after the last treatment fraction, using
% covariates measured at or before the landmark.
% Select the landmark as the last scan before the largest inter-scan gap,
% which marks the transition from on-treatment to post-treatment scans.
% This is more robust than assuming the last element is always a single
% post-treatment scan (e.g., handles protocols with multiple post-RT scans).
scan_gaps = diff(td_scan_days);
[~, gap_idx] = max(scan_gaps);
landmark_idx = gap_idx;  % last on-treatment scan is before the largest gap
landmark_day = td_scan_days(landmark_idx);
lm_keep = (t_start_td >= landmark_day);
if any(lm_keep)
    n_before = length(t_start_td);
    X_td        = X_td(lm_keep, :);
    t_start_td  = t_start_td(lm_keep);
    t_stop_td   = t_stop_td(lm_keep);
    event_td    = event_td(lm_keep);
    pat_id_td   = pat_id_td(lm_keep);
    frac_td     = frac_td(lm_keep);
    fprintf('  Landmark at day %d: %d → %d intervals (%d patients, %d events)\n', ...
        landmark_day, n_before, sum(lm_keep), numel(unique(pat_id_td)), sum(event_td == 1));

    % Re-validate event count after landmark subsetting — the pre-landmark
    % check (line 76) may have passed, but removing early intervals can
    % reduce event count below the minimum for reliable Cox estimation.
    n_events_post_lm = sum(event_td == 1);
    if n_events_post_lm < 3 || size(X_td, 1) <= td_n_feat + 1
        fprintf('  Insufficient events (%d) or intervals after landmark subsetting for Cox model.\n', n_events_post_lm);
        if ~isempty(output_folder), diary off; end
        return;
    end
end

% Z-SCORE STANDARDISATION (post-landmark only)
% Scale AFTER landmark subsetting so that standardisation statistics are
% computed exclusively on post-landmark intervals, preventing pre-landmark
% covariate distributions from leaking into the primary analysis.
% Using 'baseline' mode: z-scores are computed from each patient's first
% post-landmark observation, so that coefficients represent hazard per
% SD change from the patient's own post-treatment baseline.  This makes
% hazard ratios interpretable in clinical units: HR=2.0 means a 1-SD
% increase in the biomarker doubles the instantaneous failure hazard.
X_td_global = scale_td_panel(X_td, td_feat_names, pat_id_td, t_start_td, unique(pat_id_td), 'baseline');

% ---- Cause-Specific Hazard (CSH) Event Recoding ----------------------
% In the CSH framework for local failure analysis:
%   event=1 → local/regional failure (the event of interest)
%   event=2 → competing risk (non-cancer death without prior LF)
%   event=0 → administratively censored (still alive, no LF, lost to follow-up)
%
% For the CSH model, competing events are recoded as censored (event=0).
% This means the CSH estimates the hazard of local failure among patients
% who have NOT yet experienced any event — the "cause-specific" hazard
% rate.  This is the correct estimand when asking "does this biomarker
% predict local failure?" (aetiological question).
event_td_csh = event_td;
event_td_csh(event_td_csh == 2) = 0;  % CSH: competing risks → censored

% IPCW corrects for informative *administrative* censoring only.
% Competing events are excluded from the censoring model to preserve
% the IPCW independence assumption (Robins & Finkelstein 2000).
ipcw_weights = compute_ipcw_weights(event_td, t_start_td, t_stop_td, X_td_global, pat_id_td);

% ---- Detect tied event times for special handling --------------------
% Check for tied event times which can affect Schoenfeld residuals
% and require specialized Cox regression methods
event_times = t_stop_td(event_td_csh == 1);
[unique_times, ~, time_idx] = unique(event_times);
tied_times = false(length(unique_times), 1);
for i = 1:length(unique_times)
    tied_times(i) = sum(time_idx == i) > 1;
end
has_tied_times = any(tied_times);

if has_tied_times
    n_tied = sum(tied_times);
    fprintf('  ⚠️  Detected tied event times at %d time point(s). Using Efron approximation.\n', n_tied);
    ties_method = 'efron';
else
    ties_method = 'breslow';
end

% ---- Fit the time-dependent Cox PH model ------------------------------
% The Cox PH model estimates the hazard function:
%   h(t|X) = h0(t) * exp(beta * X(t))
% where h0(t) is the unspecified baseline hazard (semi-parametric) and
% beta * X(t) is the linear predictor from time-varying covariates.
%
% coxphfit accepts a two-column [t_start t_stop] matrix for counting-
% process (start-stop) data, with each row representing one interval
% for one patient.  The Efron method provides better handling of tied
% event times compared to Breslow when there are many ties.
warning('error', 'stats:coxphfit:FitWarning');
warning('error', 'stats:coxphfit:IterationLimit');
try
    T_td = [t_start_td, t_stop_td];
    is_censored = (event_td_csh == 0);

    % Remove constant columns to prevent DisallowedConstantTerm warning
    [X_td_clean, keep_main] = remove_constant_columns(X_td_global);
    if sum(~keep_main) > 0
        fprintf('  💡 Removed %d constant column(s) before Cox fit: %s\n', ...
            sum(~keep_main), strjoin(td_feat_names(~keep_main), ', '));
    end
    if size(X_td_clean, 2) == 0
        error('Cox:NoVariableColumns', 'All covariate columns are constant after scaling.');
    end

    % Suppress trivial warnings, but KEEP error triggers for convergence failures
    w_temp = warning('off', 'all');
    warning('error', 'stats:coxphfit:FitWarning');
    warning('error', 'stats:coxphfit:IterationLimit');
    % IPCW weighting via coxphfit's 'Frequency' parameter.  'Frequency'
    % expects integer counts, so we scale weights by 100 before rounding to
    % preserve relative differences.  A scale factor of 1 (i.e., simple
    % rounding of mean-stabilised ~1 weights) collapses most weights to 1,
    % losing the IPCW correction.  Scaling by 100 retains two decimal
    % places of precision.  The ideal solution would be a custom weighted
    % Cox partial likelihood, but MATLAB's coxphfit does not support
    % continuous probability weights natively.
    ipcw_scale = 100;
    ipcw_freq = max(1, round(ipcw_weights * ipcw_scale));
    fprintf('  IPCW→Frequency: scale=%d, effective N=%d (actual N=%d)\n', ...
        ipcw_scale, sum(ipcw_freq), length(ipcw_weights));
    [b_td_short, logl_td, ~, stats_td_short] = coxphfit(X_td_clean, T_td, ...
        'Censoring', is_censored, 'Ties', ties_method, ...
        'Frequency', ipcw_freq);
    warning(w_temp);

    % Correct variance inflation from the Frequency workaround.
    % coxphfit's Frequency parameter treats each "replicate" as an
    % independent observation, inflating effective N by ~ipcw_scale and
    % deflating SEs by ~sqrt(ipcw_scale).  Rescale SEs and recompute
    % p-values to recover correct inference.
    % Use actual mean frequency weight (accounts for rounding/truncation)
    % rather than the nominal ipcw_scale constant.
    eff_ipcw_scale = mean(ipcw_freq);
    stats_td_short.se = stats_td_short.se * sqrt(eff_ipcw_scale);
    stats_td_short.p  = 2 * (1 - normcdf(abs(b_td_short ./ stats_td_short.se)));

    % Map back to full feature space (removed columns get coef=0, SE/p=NaN)
    % Map coefficients from the reduced (non-constant) column space back to
    % the full feature space. Removed columns get beta=0 (no effect),
    % SE=NaN, p=NaN (indeterminate), so they print as NaN in the results
    % table and do not contribute to the hazard ratio.
    b_td = zeros(td_n_feat, 1);
    b_td(keep_main) = b_td_short;
    stats_td.se = nan(td_n_feat, 1);
    stats_td.p  = nan(td_n_feat, 1);
    stats_td.se(keep_main) = stats_td_short.se;
    stats_td.p(keep_main)  = stats_td_short.p;
catch ME_td
    if exist('w_temp', 'var'), warning(w_temp); end
    if ~isempty(strfind(ME_td.identifier, 'FitWarning')) || ...
       ~isempty(strfind(ME_td.identifier, 'IterationLimit')) || ...
       ~isempty(strfind(lower(ME_td.message), 'perfect'))
        % Cox model failed to converge (e.g., perfect separation).
        % Set outputs to NaN rather than falling back to a different model
        % to avoid mixing Hazard Ratios with Odds Ratios.
        fprintf('  ⚠️  Cox model did not converge: %s\n', ME_td.message);
        b_td        = nan(td_n_feat, 1);
        stats_td.se = nan(td_n_feat, 1);
        stats_td.p  = nan(td_n_feat, 1);
        logl_td     = NaN;
    else
        warning('on', 'stats:coxphfit:FitWarning');
        warning('on', 'stats:coxphfit:IterationLimit');
        rethrow(ME_td);
    end
end
warning('on', 'stats:coxphfit:FitWarning');
warning('on', 'stats:coxphfit:IterationLimit');

% ---- Likelihood-ratio test (LRT) vs. null model (no covariates) ------
% The LRT tests the global null hypothesis H0: all beta = 0 (no covariate
% has any association with failure hazard).  If the LRT p-value is
% significant, at least one diffusion biomarker is associated with local
% failure risk.  This is a more powerful omnibus test than examining
% individual covariate p-values, especially with correlated predictors.
%
% Compute both log-likelihoods manually using the original (unscaled)
% IPCW weights to avoid the inflation/deflation approximation inherent
% in dividing by mean(ipcw_freq).  The Frequency workaround inflates
% log-likelihoods non-uniformly across risk sets, so a scalar deflation
% factor can produce incorrect chi-squared statistics.
try
    is_censored_null = (event_td_csh == 0);
    w_temp_null = warning('off', 'all');
    cleanupObj = onCleanup(@() warning(w_temp_null));
    ev_lrt = ~is_censored_null;
    w_lrt  = ipcw_weights;  % original continuous weights, not inflated
    unique_fail = unique(t_stop_td(ev_lrt));

    % Null partial log-likelihood (no covariates):
    %   sum_events [ w_i * ( 0 - log(sum_{j in R_i} w_j) ) ]
    logl_null_lrt = 0;
    for uf_i = 1:length(unique_fail)
        tf = unique_fail(uf_i);
        ev_at_t = (t_stop_td == tf) & ev_lrt;
        at_risk = (t_start_td < tf) & (t_stop_td >= tf);
        R_t = sum(w_lrt(at_risk));
        if R_t > 0
            logl_null_lrt = logl_null_lrt - sum(w_lrt(ev_at_t)) * log(R_t);
        end
    end

    % Model partial log-likelihood using fitted coefficients b_td_short:
    %   sum_events [ w_i * ( X_i*b - log(sum_{j in R_i} w_j*exp(X_j*b)) ) ]
    eta = X_td_clean * b_td_short;  % linear predictor
    logl_model_lrt = 0;
    for uf_i = 1:length(unique_fail)
        tf = unique_fail(uf_i);
        ev_at_t = (t_stop_td == tf) & ev_lrt;
        at_risk = (t_start_td < tf) & (t_stop_td >= tf);
        R_t = sum(w_lrt(at_risk) .* exp(eta(at_risk)));
        if R_t > 0
            logl_model_lrt = logl_model_lrt + ...
                sum(w_lrt(ev_at_t) .* eta(ev_at_t)) - ...
                sum(w_lrt(ev_at_t)) * log(R_t);
        end
    end

    % LRT statistic follows chi-squared distribution under the null.
    % max(0, ...) guards against numerical rounding producing negative values.
    LRT_stat = 2 * (logl_model_lrt - logl_null_lrt);
    LRT_df   = sum(keep_main);  % degrees of freedom = number of non-constant features actually fit
    LRT_p    = 1 - chi2cdf(max(0, LRT_stat), LRT_df);
catch
    LRT_stat = NaN; LRT_p = NaN;
end

% ---- Print results table -------------------------------------------
% Hazard Ratio (HR) interpretation:
%   HR > 1: higher covariate value → increased failure hazard
%   HR < 1: higher covariate value → decreased failure hazard
%   HR = 1: no association
%   95% CI not crossing 1.0 → significant at alpha=0.05
% Clinical expectation: HR(ADC) < 1 and HR(D) < 1 because higher
% diffusion (less cellularity) should be protective against local failure.
fprintf('  %-10s  %6s  %6s  %6s  %6s\n', 'Covariate', 'HR ', 'CI_lo', 'CI_hi', 'p');
fprintf('  %s\n', repmat('-', 1, 52));
for fi = 1:td_n_feat
    % Convert log-hazard coefficients to hazard ratios via exponentiation.
    % 95% CI uses the Wald interval: exp(beta +/- 1.96 * SE(beta)).
    hr_i  = exp(b_td(fi));
    ci_lo = exp(b_td(fi) - 1.96*stats_td.se(fi));
    ci_hi = exp(b_td(fi) + 1.96*stats_td.se(fi));
    fprintf('  %-10s  %6.3f  %6.3f  %6.3f  %6.4f\n', ...
        td_feat_names{fi}, hr_i, ci_lo, ci_hi, stats_td.p(fi));
end

% --- Bootstrap CIs for hazard ratios ---
% BCa bootstrap provides non-parametric CI estimates that do not rely on
% the Wald normal approximation, which can be inaccurate for small samples.
try
    n_boot_hr = 1000;
    fprintf('\n  --- Bootstrap 95%% CIs for HR (B=%d) ---\n', n_boot_hr);
    fprintf('  %-10s  %6s  %6s  %6s\n', 'Covariate', 'HR', 'BCa_lo', 'BCa_hi');
    fprintf('  %s\n', repmat('-', 1, 44));
    boot_data_hr = [X_td_global, T_td(:), is_censored(:), ipcw_freq(:)];
    for fi = 1:td_n_feat
        hr_fn = @(d) exp(local_coxph_coef(d, fi, keep_main, td_n_feat));
        [bci_lo, bci_hi] = bootstrap_ci(boot_data_hr, hr_fn, n_boot_hr, 0.05);
        fprintf('  %-10s  %6.3f  %6.3f  %6.3f\n', ...
            td_feat_names{fi}, exp(b_td(fi)), bci_lo, bci_hi);
    end
catch boot_err
    fprintf('  ⚠️  Bootstrap CI skipped: %s\n', boot_err.message);
end
if ~isnan(LRT_p)
    fprintf('  Global LRT: chi2(%d) = %.2f, p = %.4f\n', LRT_df, LRT_stat, LRT_p);
end

% ---- Imputation half-life sensitivity analysis -----------------------
% Report how HRs change across the half-life grid so users can assess
% robustness of the exponential decay imputation assumption.
fprintf('\n  --- Imputation Sensitivity (half-life months → HR) ---\n');
fprintf('  %-6s', 'HL');
for fi = 1:td_n_feat
    fprintf('  %10s', td_feat_names{fi});
end
fprintf('\n');
for hl_idx = 1:length(half_life_grid)
    pnl = td_panels{hl_idx};
    try
        % For each half-life, repeat the full analysis chain: landmark
        % subset, z-score scaling, IPCW weighting, and Cox PH fitting.
        % This ensures the sensitivity analysis is self-consistent and
        % not contaminated by the default half-life's statistics.
        % Apply landmark subsetting to match the primary analysis
        lm_keep_hl = (pnl.t_start >= landmark_day);
        if ~any(lm_keep_hl), error('sensitivity:noData', 'No post-landmark intervals.'); end
        pnl_X = pnl.X(lm_keep_hl, :);
        pnl_tstart = pnl.t_start(lm_keep_hl);
        pnl_tstop  = pnl.t_stop(lm_keep_hl);
        pnl_event  = pnl.event(lm_keep_hl);
        pnl_pid    = pnl.pat_id(lm_keep_hl);
        ev_csh = p