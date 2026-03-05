function [m_lf, m_total_time, m_total_follow_up_time, m_gtv_vol, m_adc_mean, m_d_mean, m_f_mean, m_dstar_mean, m_id_list, m_mrn_list, m_d95_gtvp, m_v50gy_gtvp, m_data_vectors_gtvp, lf_group, valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, nTp, metric_sets, set_names, time_labels, dtype_label, dl_provenance] = metrics_baseline(data_vectors_gtvp, data_vectors_gtvn, summary_metrics, config_struct)
% METRICS_BASELINE — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 1/5 of the metrics step. Compiles baseline measures, cleans outliers,
% computes relative changes (percent delta), and groups metric sets for later steps.
%
% Inputs:
%   data_vectors_gtvp - Struct array containing primary GTV parameter maps (by iter)
%   data_vectors_gtvn - Struct array containing nodal GTV parameter maps
%   summary_metrics   - Pre-computed overall summary values across timepoints
%   config_struct     - Configuration struct
%
% Outputs:
%   [Multiple Arrays] - Includes logical masks for valid_pts, clean arrays for 
%                       ADC_abs, D_abs, f_abs, and their delta percent variations,
%                       as well as organized sets (metric_sets) for downstream analysis.
%

% Extract parameters
dataloc = config_struct.dataloc;
adc_mean = summary_metrics.adc_mean;
adc_sd = summary_metrics.adc_sd;
d_mean = summary_metrics.d_mean;
f_mean = summary_metrics.f_mean;
dstar_mean = summary_metrics.dstar_mean;
adc_mean_rpt = summary_metrics.adc_mean_rpt;
adc_sub_rpt = summary_metrics.adc_sub_rpt;
d_mean_rpt = summary_metrics.d_mean_rpt;
f_mean_rpt = summary_metrics.f_mean_rpt;
dstar_mean_rpt = summary_metrics.dstar_mean_rpt;
n_rpt = summary_metrics.n_rpt;
id_list = summary_metrics.id_list;
mrn_list = summary_metrics.mrn_list;
gtv_vol = summary_metrics.gtv_vol;
d95_gtvp = summary_metrics.d95_gtvp;
v50gy_gtvp = summary_metrics.v50gy_gtvp;
dmean_gtvp = summary_metrics.dmean_gtvp;
gtv_locations = summary_metrics.gtv_locations;
dwi_locations = summary_metrics.dwi_locations;

fprintf('  --- SECTION 1: Repeatability Analysis ---\n');
% Within-patient coefficient of variation (wCV = SD / mean).
% Guard against division by near-zero denominators (especially for f, D*).
wcv_denom_floor = 1e-10;
if exist('OCTAVE_VERSION', 'builtin')
    % Use reshape instead of squeeze to preserve [nPatients x nDwiTypes]
    % shape.  squeeze() collapses a [1 x 1 x 3] result to [3 x 1] when
    % nPatients==1, causing dimension mismatches with the [nPatients x 1]
    % masking vectors (n_rpt, denom floors).
    nPat_rpt_oct = size(adc_mean_rpt, 1);
    nDwi_rpt_oct = size(adc_mean_rpt, 3);
    adc_denom = reshape(nanmean(adc_mean_rpt,2), [nPat_rpt_oct, nDwi_rpt_oct]);
    adc_wCV = reshape(nanstd(adc_mean_rpt,0,2), [nPat_rpt_oct, nDwi_rpt_oct]) ./ adc_denom;
    adc_wCV(n_rpt<2 | abs(adc_denom) < wcv_denom_floor, :) = nan;
    adc_sub_denom = reshape(nanmean(adc_sub_rpt,2), [nPat_rpt_oct, nDwi_rpt_oct]);
    adc_wCV_sub = reshape(nanstd(adc_sub_rpt,0,2), [nPat_rpt_oct, nDwi_rpt_oct]) ./ adc_sub_denom;
    adc_wCV_sub(n_rpt<2 | abs(adc_sub_denom) < wcv_denom_floor, :) = nan;
    d_denom = reshape(nanmean(d_mean_rpt,2), [nPat_rpt_oct, nDwi_rpt_oct]);
    d_wCV = reshape(nanstd(d_mean_rpt,0,2), [nPat_rpt_oct, nDwi_rpt_oct]) ./ d_denom;
    d_wCV(n_rpt<2 | abs(d_denom) < wcv_denom_floor, :) = nan;
    f_denom = reshape(nanmean(f_mean_rpt,2), [nPat_rpt_oct, nDwi_rpt_oct]);
    f_wCV = reshape(nanstd(f_mean_rpt,0,2), [nPat_rpt_oct, nDwi_rpt_oct]) ./ f_denom;
    f_wCV(n_rpt<2 | abs(f_denom) < wcv_denom_floor, :) = nan;
    dstar_denom = reshape(nanmean(dstar_mean_rpt,2), [nPat_rpt_oct, nDwi_rpt_oct]);
    dstar_wCV = reshape(nanstd(dstar_mean_rpt,0,2), [nPat_rpt_oct, nDwi_rpt_oct]) ./ dstar_denom;
    dstar_wCV(n_rpt<2 | abs(dstar_denom) < wcv_denom_floor, :) = nan;
else
    % Use reshape instead of squeeze to preserve [nPatients x nDwiTypes]
    % shape.  squeeze() collapses a [1 x 1 x 3] result to [3 x 1] when
    % nPatients==1, causing dimension mismatches with the [nPatients x 1]
    % masking vectors (n_rpt, denom floors).
    nPat_rpt = size(adc_mean_rpt, 1);
    nDwi_rpt = size(adc_mean_rpt, 3);
    adc_denom = reshape(mean(adc_mean_rpt,2,'omitnan'), [nPat_rpt, nDwi_rpt]);
    adc_wCV = reshape(std(adc_mean_rpt,0,2,'omitnan'), [nPat_rpt, nDwi_rpt]) ./ adc_denom;
    adc_wCV(n_rpt<2 | abs(adc_denom) < wcv_denom_floor, :) = nan;
    adc_sub_denom = reshape(mean(adc_sub_rpt,2,'omitnan'), [nPat_rpt, nDwi_rpt]);
    adc_wCV_sub = reshape(std(adc_sub_rpt,0,2,'omitnan'), [nPat_rpt, nDwi_rpt]) ./ adc_sub_denom;
    adc_wCV_sub(n_rpt<2 | abs(adc_sub_denom) < wcv_denom_floor, :) = nan;
    d_denom = reshape(mean(d_mean_rpt,2,'omitnan'), [nPat_rpt, nDwi_rpt]);
    d_wCV = reshape(std(d_mean_rpt,0,2,'omitnan'), [nPat_rpt, nDwi_rpt]) ./ d_denom;
    d_wCV(n_rpt<2 | abs(d_denom) < wcv_denom_floor, :) = nan;
    f_denom = reshape(mean(f_mean_rpt,2,'omitnan'), [nPat_rpt, nDwi_rpt]);
    f_wCV = reshape(std(f_mean_rpt,0,2,'omitnan'), [nPat_rpt, nDwi_rpt]) ./ f_denom;
    f_wCV(n_rpt<2 | abs(f_denom) < wcv_denom_floor, :) = nan;
    dstar_denom = reshape(mean(dstar_mean_rpt,2,'omitnan'), [nPat_rpt, nDwi_rpt]);
    dstar_wCV = reshape(std(dstar_mean_rpt,0,2,'omitnan'), [nPat_rpt, nDwi_rpt]) ./ dstar_denom;
    dstar_wCV(n_rpt<2 | abs(dstar_denom) < wcv_denom_floor, :) = nan;
end

fprintf('  --- SECTION 2: Load Clinical Outcome Data ---\n');
clinical_data_sheet = fullfile(dataloc, config_struct.clinical_data_sheet);

if exist('OCTAVE_VERSION', 'builtin')
    warning('Skipping actual readtable for clinical data in Octave mock environment.');
    T = struct();
    T.Pat = id_list;
    T.MRN = mrn_list;
    T.LocalOrRegionalFailure = zeros(size(id_list));
    T.LocoregionalFailureDateOfLocalOrRegionalFailure = zeros(size(id_list));
    T.LocalFailureDateOfLocalFailureOrCensor = zeros(size(id_list));
    T.RegionalFailureDateOfRegionalFailureOrCensor = zeros(size(id_list));
    T.RTStartDate = zeros(size(id_list));
    T.RTStopDate = zeros(size(id_list));

    T_prop = struct();
    T_prop.VariableNames = {'Pat', 'MRN', 'LocalOrRegionalFailure', 'LocoregionalFailureDateOfLocalOrRegionalFailure', 'LocalFailureDateOfLocalFailureOrCensor', 'RegionalFailureDateOfRegionalFailureOrCensor', 'RTStartDate', 'RTStopDate'};
    T.Properties = T_prop;
else
    if isfield(config_struct, 'clinical_sheet_name')
        T = readtable(clinical_data_sheet,'Sheet', config_struct.clinical_sheet_name);
    else
        T = readtable(clinical_data_sheet);
    end
end

if exist('OCTAVE_VERSION', 'builtin')
    % Create mock outcome data
    lf = zeros(length(id_list), 1);
    m_total_time = zeros(length(id_list), 1);
    m_total_follow_up_time = zeros(length(id_list), 1);
    lf_date = zeros(length(id_list), 1);
    censor_date = zeros(length(id_list), 1);
    rtstartdate = zeros(length(id_list), 1);
    rtenddate = zeros(length(id_list), 1);
else
    lf = nan(length(id_list), 1);
    lf_date = repmat(datetime(NaT), length(id_list), 1);
    censor_date = repmat(datetime(NaT), length(id_list), 1);
    rtstartdate = repmat(datetime(NaT), length(id_list), 1);
    rtenddate = repmat(datetime(NaT), length(id_list), 1);

    pat_normalized = strrep(T.Pat,'_','-');
    id_list_normalized = strrep(id_list,'_','-');

    n_ids_clinical = length(id_list);
    for j = 1:n_ids_clinical
        text_progress_bar(j, n_ids_clinical, 'Matching clinical data');
        i_find = find(strcmp(pat_normalized, id_list_normalized{j}));
        if ~isempty(i_find)
            i_find = i_find(1);
            lf(j) = T.LocalOrRegionalFailure(i_find);
            % --- Competing risk classification ---
            % Clinical decision: a patient with documented local/regional
            % failure (lf==1) is always coded as an *event* regardless of
            % cause of death.  Rationale: imaging-confirmed LF represents
            % the endpoint of interest; subsequent non-cancer death does
            % not negate the failure observation.
            %
            % Competing risk (lf==2) is assigned ONLY when:
            %   (a) NO local/regional failure was documented (lf==0), AND
            %   (b) a non-cancer cause of death is recorded.
            % This means patients who had LF AND later died of non-cancer
            % causes are counted as events (lf==1), not competing risks.
            if ismember('CauseOfDeath', T.Properties.VariableNames)
                cod_raw = T.CauseOfDeath(i_find);
                if iscell(cod_raw), cod_raw = cod_raw{1}; end
                cod = char(string(cod_raw));  % handles string, categorical, numeric, missing
                cod_lower = strtrim(lower(cod));
                is_unknown = isempty(cod_lower) || strcmp(cod_lower, 'unknown') || ...
                             strcmp(cod_lower, 'pending') || strcmp(cod_lower, 'n/a');
                if lf(j) == 0 && ~isempty(cod) && ~is_unknown && isempty(strfind(cod_lower, 'cancer'))
                    lf(j) = 2; % Competing risk: non-cancer death without LF
                end
            end
            lf_date(j) = T.LocoregionalFailureDateOfLocalOrRegionalFailure(i_find);
            % Use the later of the two censor dates, but if one is NaT use the other.
            lf_cens = T.LocalFailureDateOfLocalFailureOrCensor(i_find);
            rf_cens = T.RegionalFailureDateOfRegionalFailureOrCensor(i_find);
            if isnat(lf_cens) && isnat(rf_cens)
                censor_date(j) = NaT;
            elseif isnat(lf_cens)
                censor_date(j) = rf_cens;
            elseif isnat(rf_cens)
                censor_date(j) = lf_cens;
            else
                censor_date(j) = max(lf_cens, rf_cens);
            end
            rtstartdate(j) = T.RTStartDate(i_find);
            rtenddate(j) = T.RTStopDate(i_find);
        end
    end
end

% Warn if competing risks could not be identified because the clinical
% spreadsheet lacks a CauseOfDeath column.  Without it, non-cancer deaths
% are treated as administrative censoring, which inflates the CSH estimate.
if ~exist('OCTAVE_VERSION', 'builtin') && ~ismember('CauseOfDeath', T.Properties.VariableNames)
    warning('metrics_baseline:noCauseOfDeath', ...
        'CauseOfDeath column not found in clinical spreadsheet. Competing risks cannot be identified; CSH survival analysis may be biased.');
end

if exist('OCTAVE_VERSION', 'builtin')
    total_time = m_total_time;
    total_follow_up_time = m_total_follow_up_time;
else
    total_time = days(lf_date - rtenddate);
    total_follow_up_time = days(censor_date - rtenddate);
end

fprintf('  --- DEEP LEARNING RIGOR AUDIT ---\n');
manifest_file = fullfile(config_struct.dataloc, 'dl_validation_manifest.mat');
dl_provenance = load_dl_provenance(manifest_file);

dwi_type_names = {'Standard', 'dnCNN', 'IVIMnet'};
if isfield(config_struct, 'output_folder')
    output_folder = config_struct.output_folder;
else
    output_folder = fullfile(config_struct.dataloc, 'saved_figures');
end
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
dtype = config_struct.dwi_types_to_run;
diary_file = fullfile(output_folder, ['metrics_baseline_output_' dwi_type_names{dtype} '.txt']);
if exist(diary_file, 'file'), delete(diary_file); end
diary(diary_file);
w_state_baseline = warning('off', 'stats:glmfit:IterationLimit');
prev_fig_vis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');  

if ~exist('nTp', 'var') || isempty(nTp)
    nTp = size(adc_mean, 2);
end

% Set up common variables for the next scripts
dtype_label = dwi_type_names{dtype};
exclude_missing_baseline = true;
valid_baseline = ~isnan(gtv_vol(:,1)) & ~isnan(adc_mean(:,1,dtype));
exclude_outliers = true;

% Outlier detection is outcome-BLINDED: thresholds are derived solely from
% baseline imaging metrics (IQR fences) without reference to the lf outcome
% variable.  However, differential removal across outcome groups can still
% introduce selection bias.  We log the outcome distribution of flagged
% outliers below so the researcher can verify no systematic imbalance.
baseline_metrics_oi   = {adc_mean(:,1,dtype), d_mean(:,1,dtype), f_mean(:,1,dtype), dstar_mean(:,1,dtype)};
baseline_metric_names = {'ADC', 'D', 'f', 'D*'};
is_outlier = false(size(lf));
n_baseline_metrics = numel(baseline_metrics_oi);
for metric_idx = 1:n_baseline_metrics
    text_progress_bar(metric_idx, n_baseline_metrics, 'Detecting outliers');
    col = baseline_metrics_oi{metric_idx};
    col_clean = col(~isnan(col));
    if numel(col_clean) < 3, continue; end
    med_val = median(col_clean);
    iqr_val = iqr(col_clean);
    if iqr_val == 0, continue; end
    lower_fence = med_val - 3 * iqr_val;
    upper_fence = med_val + 3 * iqr_val;
    outlier_flags = (col < lower_fence | col > upper_fence) & ~isnan(col);
    if any(outlier_flags)
        n_out = sum(outlier_flags);
        n_out_lf = sum(outlier_flags & (lf == 1));
        n_out_lc = sum(outlier_flags & (lf == 0));
        n_out_cr = sum(outlier_flags & (lf == 2));
        fprintf('  Outlier flag (%s): %d flagged (LF=%d, LC=%d, CR=%d)\n', ...
            baseline_metric_names{metric_idx}, n_out, n_out_lf, n_out_lc, n_out_cr);
    end
    is_outlier = is_outlier | outlier_flags;
end
n_total_outliers = sum(is_outlier);
if n_total_outliers > 0
    fprintf('  Total outliers removed: %d / %d (%.1f%%)\n', ...
        n_total_outliers, numel(lf), 100*n_total_outliers/numel(lf));
end
non_outlier = ~is_outlier;

m_lf                   = lf;
m_total_time           = total_time;
m_total_follow_up_time = total_follow_up_time;
m_gtv_vol              = gtv_vol;
m_adc_mean             = adc_mean;
m_d_mean               = d_mean;
m_f_mean               = f_mean;
m_dstar_mean           = dstar_mean;
m_id_list              = id_list;
m_mrn_list             = mrn_list;
m_d95_gtvp             = d95_gtvp;
m_v50gy_gtvp           = v50gy_gtvp;
m_data_vectors_gtvp    = data_vectors_gtvp;

if size(m_gtv_vol, 2) < nTp, m_gtv_vol = [m_gtv_vol, nan(size(m_gtv_vol, 1), nTp - size(m_gtv_vol, 2))]; end
if size(m_d95_gtvp, 2) < nTp, m_d95_gtvp = [m_d95_gtvp, nan(size(m_d95_gtvp, 1), nTp - size(m_d95_gtvp, 2))]; end
if size(m_v50gy_gtvp, 2) < nTp, m_v50gy_gtvp = [m_v50gy_gtvp, nan(size(m_v50gy_gtvp, 1), nTp - size(m_v50gy_gtvp, 2))]; end
if size(dmean_gtvp, 2) < nTp, dmean_gtvp = [dmean_gtvp, nan(size(dmean_gtvp, 1), nTp - size(dmean_gtvp, 2))]; end


if exclude_missing_baseline
    % Report imaging-based exclusion statistics for informative censoring assessment.
    % Patients are excluded based on imaging availability (valid_baseline), which
    % could be correlated with disease severity.  If event rates differ between
    % excluded and included patients, results may not generalise.
    n_excluded_baseline = sum(~valid_baseline);
    n_total_cohort = numel(valid_baseline);
    if n_excluded_baseline > 0
        fprintf('  ⚠️  Excluded %d/%d patients due to missing baseline imaging (GTV or ADC).\n', ...
            n_excluded_baseline, n_total_cohort);
        lf_excluded = lf(~valid_baseline);
        lf_included = lf(valid_baseline);
        lf_rate_excl = 100 * sum(lf_excluded == 1) / max(1, sum(isfinite(lf_excluded)));
        lf_rate_incl = 100 * sum(lf_included == 1) / max(1, sum(isfinite(lf_included)));
        fprintf('      LF rate: included=%.1f%%, excluded=%.1f%%\n', lf_rate_incl, lf_rate_excl);
        if abs(lf_rate_excl - lf_rate_incl) > 10
            fprintf('  ⚠️  >10pp difference in LF rates between imaging-excluded and included patients.\n');
            fprintf('      Potential informative censoring — interpret downstream survival results with caution.\n');
        end
    end

    m_lf                   = m_lf(valid_baseline);
    m_total_time           = m_total_time(valid_baseline);
    m_total_follow_up_time = m_total_follow_up_time(valid_baseline);
    m_gtv_vol              = m_gtv_vol(valid_baseline,:);
    m_adc_mean             = m_adc_mean(valid_baseline,:,:);
    m_d_mean               = m_d_mean(valid_baseline,:,:);
    m_f_mean               = m_f_mean(valid_baseline,:,:);
    m_dstar_mean           = m_dstar_mean(valid_baseline,:,:);
    m_id_list              = m_id_list(valid_baseline);
    m_mrn_list             = m_mrn_list(valid_baseline);
    m_d95_gtvp             = m_d95_gtvp(valid_baseline,:);
    m_v50gy_gtvp           = m_v50gy_gtvp(valid_baseline,:);
    m_data_vectors_gtvp    = m_data_vectors_gtvp(valid_baseline,:,:);
end

if exclude_outliers
    if exclude_missing_baseline
        outlier_current = ~non_outlier(valid_baseline);
    else
        outlier_current = ~non_outlier;
    end
    outlier_ids = m_id_list(outlier_current);
    if ~isempty(outlier_ids)
        fprintf('  ⚠️  Removed %d patients as outliers (IDs: %s)\n', numel(outlier_ids), strjoin(outlier_ids, ', '));
    end
    m_adc_mean(outlier_current,:,:)             = NaN;
    m_d_mean(outlier_current,:,:)               = NaN;
    m_f_mean(outlier_current,:,:)               = NaN;
    m_dstar_mean(outlier_current,:,:)           = NaN;
    m_gtv_vol(outlier_current,:)                = NaN;
    m_d95_gtvp(outlier_current,:)               = NaN;
    m_v50gy_gtvp(outlier_current,:)             = NaN;
end

ADC_abs = m_adc_mean(:,:,dtype);
D_abs   = m_d_mean(:,:,dtype);
f_abs   = m_f_mean(:,:,dtype);
Dstar_abs = m_dstar_mean(:,:,dtype);

% Use a data-adaptive epsilon (1% of baseline IQR) to prevent inflated
% percent changes when baseline values are near zero, while remaining
% proportional to the actual measurement scale.
% NOTE: Outlier patients were NaN-ified above (lines 318-324), so the
% isfinite() filter correctly excludes them from epsilon estimation.
adc_iqr_raw  = iqr(ADC_abs(isfinite(ADC_abs(:,1)), 1));
d_iqr_raw    = iqr(D_abs(isfinite(D_abs(:,1)), 1));
dstar_iqr_raw = iqr(Dstar_abs(isfinite(Dstar_abs(:,1)), 1));
% max(scalar, NaN) returns NaN in MATLAB's two-input form; guard explicitly.
if isnan(adc_iqr_raw),  adc_iqr_raw  = 0; end
if isnan(d_iqr_raw),    d_iqr_raw    = 0; end
if isnan(dstar_iqr_raw), dstar_iqr_raw = 0; end
adc_eps  = max(1e-8, 0.01 * adc_iqr_raw);
d_eps    = max(1e-8, 0.01 * d_iqr_raw);
dstar_eps = max(1e-8, 0.01 * dstar_iqr_raw);
% Exclude patients with near-zero or negative baselines from percent change
% computation to avoid sign-flipped or inflated ratios.  These patients get
% NaN percent change so they are excluded from downstream group comparisons.
% Known trade-off: this creates a discontinuity — baselines just below
% epsilon are excluded (NaN), while baselines just above may produce large
% ratios that are subsequently winsorized to ±500%.
adc_bl = ADC_abs(:,1);  adc_bl(adc_bl < adc_eps) = NaN;
d_bl   = D_abs(:,1);    d_bl(d_bl < d_eps) = NaN;
dstar_bl = Dstar_abs(:,1); dstar_bl(dstar_bl < dstar_eps) = NaN;
ADC_pct = ((ADC_abs - ADC_abs(:,1)) ./ adc_bl) * 100;
D_pct   = ((D_abs - D_abs(:,1)) ./ d_bl) * 100;
f_delta = (f_abs - f_abs(:,1));
Dstar_pct = ((Dstar_abs - Dstar_abs(:,1)) ./ dstar_bl) * 100;

% Winsorize percent changes at ±500% to limit influence of near-zero
% baselines that passed the epsilon filter but still produce extreme ratios.
pct_clip = 500;
ADC_pct(ADC_pct < -pct_clip) = -pct_clip;  ADC_pct(ADC_pct > pct_clip) = pct_clip;
D_pct(D_pct < -pct_clip) = -pct_clip;      D_pct(D_pct > pct_clip) = pct_clip;
Dstar_pct(Dstar_pct < -pct_clip) = -pct_clip;  Dstar_pct(Dstar_pct > pct_clip) = pct_clip;

valid_pts = isfinite(m_lf);
lf_group = m_lf(valid_pts);

% Re-organize metrics into distinct Sets:
% NOTE: f_delta is absolute change (f range [0,1]), not percent change.
% It is grouped with percent-change metrics for convenience in downstream
% univariate analysis (each metric tested independently).  In multivariate
% models (Cox, GLME), scale_td_panel z-scores all features to comparable
% scales, so the raw-scale difference does not affect model coefficients.
% The set_names label already distinguishes it as '\Delta f (abs)'.
metric_sets = {
    {ADC_abs, D_abs, f_abs, Dstar_abs}, ...          % Set 1
    {ADC_pct, D_pct, f_delta, Dstar_pct} ...          % Set 2
};

set_names = {
    {'ADC Absolute', 'D Absolute', 'f Absolute', 'D* Absolute'}, ...
    {'\Delta ADC (%)', '\Delta D (%)', '\Delta f (abs)', '\Delta D* (%)'} ...
};
time_labels = {'Fx1', 'Fx2', 'Fx3', 'Fx4', 'Fx5', 'Post'};

% Restore global state modified at the top of this function
diary off;
warning(w_state_baseline);
set(0, 'DefaultFigureVisible', prev_fig_vis);
end
