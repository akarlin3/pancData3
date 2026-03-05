function [m_lf, m_total_time, m_total_follow_up_time, m_gtv_vol, m_adc_mean, m_d_mean, m_f_mean, m_dstar_mean, m_id_list, m_mrn_list, m_d95_gtvp, m_v50gy_gtvp, m_data_vectors_gtvp, lf_group, valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_delta, Dstar_pct, nTp, metric_sets, set_names, time_labels, dtype_label, dl_provenance] = metrics_baseline(data_vectors_gtvp, data_vectors_gtvn, summary_metrics, config_struct)
% METRICS_BASELINE — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% Part 1/5 of the metrics step. Compiles baseline measures, cleans outliers,
% computes relative changes (percent delta), and groups metric sets for later steps.
%
% ANALYTICAL OVERVIEW:
%   This function serves as the data preparation gateway for all downstream
%   statistical and survival analyses.  It performs four critical tasks:
%
%   1. REPEATABILITY ANALYSIS — Quantifies measurement noise via within-patient
%      coefficient of variation (wCV) from same-day repeat scans.  This establishes
%      the minimum detectable change (Coefficient of Reproducibility) for each
%      biomarker, which is essential for distinguishing real treatment-induced
%      changes from random measurement fluctuation in pancreatic DWI (where
%      respiratory motion and susceptibility artifacts at the pancreas-air
%      interface are substantial noise sources).
%
%   2. CLINICAL OUTCOME LINKAGE — Matches imaging biomarkers to clinical
%      endpoints (locoregional failure, competing risks) from the clinical
%      spreadsheet.  The competing risk framework is critical because pancreatic
%      cancer patients frequently die of systemic disease or treatment
%      complications before local failure can be observed.  Ignoring competing
%      risks inflates the Cause-Specific Hazard estimate for local failure.
%
%   3. OUTLIER REMOVAL AND COHORT FILTERING — Applies outcome-blinded IQR
%      fencing to remove biophysically implausible measurements before they
%      can distort group comparisons.  Patients missing baseline imaging are
%      excluded because percent-change metrics require a reference value.
%
%   4. PERCENT DELTA COMPUTATION — Computes treatment-induced changes relative
%      to baseline.  ADC, D, and D* use percent change; f uses absolute change
%      because f values near zero (typical range 0.05-0.15) make percent change
%      numerically unstable and clinically uninterpretable.
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

% =========================================================================
% EXTRACT PRE-COMPUTED SUMMARY METRICS
% =========================================================================
% These arrays were aggregated during the 'load' step by compute_summary_metrics.m.
% Each metric is stored as [nPatients x nTimepoints x nDwiTypes], where DWI
% types are 1=Standard, 2=dnCNN, 3=IVIMnet.  The third dimension allows the
% same patient cohort to be analysed under different denoising pipelines
% without re-running the expensive DICOM loading and model fitting steps.
%
% Key diffusion parameters:
%   ADC  — Apparent Diffusion Coefficient from mono-exponential fit.
%          Reflects total water mobility (cellularity + perfusion).
%   D    — True tissue diffusion from IVIM bi-exponential fit.
%          Isolates cellularity by separating out the perfusion component.
%   f    — Perfusion fraction from IVIM.  Represents the volume fraction of
%          capillary blood flow.  Ranges [0,1] and is often very small (~0.05-0.15).
%   D*   — Pseudo-diffusion coefficient from IVIM.  Reflects the speed of
%          microcirculatory blood flow.  Physiologically noisy and highly variable.
%
% Dose metrics:
%   d95_gtvp   — Minimum dose (Gy) to 95% of the primary GTV volume.
%   v50gy_gtvp — Fraction of the primary GTV receiving >= 50 Gy.
%   dmean_gtvp — Mean dose to the primary GTV.
%
% Repeatability arrays (*_rpt):
%   Same-day repeat scans used to quantify measurement reproducibility.
%   Stored as [nPatients x nRepeats x nDwiTypes].
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
% =========================================================================
% REPEATABILITY ANALYSIS: Within-Patient Coefficient of Variation (wCV)
% =========================================================================
% Quantifies measurement reproducibility from same-day repeat DWI scans.
% wCV = SD(repeats) / mean(repeats) for each patient and DWI type.
%
% Clinical rationale: wCV establishes a "noise floor" for each biomarker.
% A treatment-induced change must exceed the Coefficient of Reproducibility
% (CoR ≈ 1.96 * sqrt(2) * wCV) to be considered a real biological signal
% rather than measurement noise.  This is critical for pancreatic DWI
% because respiratory motion and susceptibility artifacts at the
% pancreas-air interface inflate measurement variability.
%
% Patients with fewer than 2 repeat scans (n_rpt < 2) cannot have wCV
% computed and are set to NaN.  The denominator floor (1e-10) prevents
% division by zero for biomarkers with near-zero means (f and D* can be
% very small in poorly perfused pancreatic tumours).
%
% wCV is computed per-patient (not pooled) because inter-patient variation
% in tumor composition would inflate a pooled wCV, making the noise floor
% appear artificially high.
wcv_denom_floor = 1e-10;
if exist('OCTAVE_VERSION', 'builtin')
    % Use reshape instead of squeeze to preserve [nPatients x nDwiTypes]
    % shape.  squeeze() collapses a [1 x 1 x 3] result to [3 x 1] when
    % nPatients==1, causing dimension mismatches with the [nPatients x 1]
    % masking vectors (n_rpt, denom floors).
    nPat_rpt_oct = size(adc_mean_rpt, 1);
    nDwi_rpt_oct = size(adc_mean_rpt, 3);
    adc_denom = reshape(nanmean(adc_mean_rpt,2), [nPat_rpt_oct, nDwi_rpt_oct]);
    adc_denom(abs(adc_denom) < wcv_denom_floor) = nan;
    adc_wCV = reshape(nanstd(adc_mean_rpt,0,2), [nPat_rpt_oct, nDwi_rpt_oct]) ./ adc_denom;
    adc_wCV(n_rpt<2, :) = nan;
    adc_sub_denom = reshape(nanmean(adc_sub_rpt,2), [nPat_rpt_oct, nDwi_rpt_oct]);
    adc_sub_denom(abs(adc_sub_denom) < wcv_denom_floor) = nan;
    adc_wCV_sub = reshape(nanstd(adc_sub_rpt,0,2), [nPat_rpt_oct, nDwi_rpt_oct]) ./ adc_sub_denom;
    adc_wCV_sub(n_rpt<2, :) = nan;
    d_denom = reshape(nanmean(d_mean_rpt,2), [nPat_rpt_oct, nDwi_rpt_oct]);
    d_denom(abs(d_denom) < wcv_denom_floor) = nan;
    d_wCV = reshape(nanstd(d_mean_rpt,0,2), [nPat_rpt_oct, nDwi_rpt_oct]) ./ d_denom;
    d_wCV(n_rpt<2, :) = nan;
    f_denom = reshape(nanmean(f_mean_rpt,2), [nPat_rpt_oct, nDwi_rpt_oct]);
    f_denom(abs(f_denom) < wcv_denom_floor) = nan;
    f_wCV = reshape(nanstd(f_mean_rpt,0,2), [nPat_rpt_oct, nDwi_rpt_oct]) ./ f_denom;
    f_wCV(n_rpt<2, :) = nan;
    dstar_denom = reshape(nanmean(dstar_mean_rpt,2), [nPat_rpt_oct, nDwi_rpt_oct]);
    dstar_denom(abs(dstar_denom) < wcv_denom_floor) = nan;
    dstar_wCV = reshape(nanstd(dstar_mean_rpt,0,2), [nPat_rpt_oct, nDwi_rpt_oct]) ./ dstar_denom;
    dstar_wCV(n_rpt<2, :) = nan;
else
    % Use reshape instead of squeeze to preserve [nPatients x nDwiTypes]
    % shape.  squeeze() collapses a [1 x 1 x 3] result to [3 x 1] when
    % nPatients==1, causing dimension mismatches with the [nPatients x 1]
    % masking vectors (n_rpt, denom floors).
    nPat_rpt = size(adc_mean_rpt, 1);
    nDwi_rpt = size(adc_mean_rpt, 3);
    adc_denom = reshape(mean(adc_mean_rpt,2,'omitnan'), [nPat_rpt, nDwi_rpt]);
    adc_denom(abs(adc_denom) < wcv_denom_floor) = nan;
    adc_wCV = reshape(std(adc_mean_rpt,0,2,'omitnan'), [nPat_rpt, nDwi_rpt]) ./ adc_denom;
    adc_wCV(n_rpt<2, :) = nan;
    adc_sub_denom = reshape(mean(adc_sub_rpt,2,'omitnan'), [nPat_rpt, nDwi_rpt]);
    adc_sub_denom(abs(adc_sub_denom) < wcv_denom_floor) = nan;
    adc_wCV_sub = reshape(std(adc_sub_rpt,0,2,'omitnan'), [nPat_rpt, nDwi_rpt]) ./ adc_sub_denom;
    adc_wCV_sub(n_rpt<2, :) = nan;
    d_denom = reshape(mean(d_mean_rpt,2,'omitnan'), [nPat_rpt, nDwi_rpt]);
    d_denom(abs(d_denom) < wcv_denom_floor) = nan;
    d_wCV = reshape(std(d_mean_rpt,0,2,'omitnan'), [nPat_rpt, nDwi_rpt]) ./ d_denom;
    d_wCV(n_rpt<2, :) = nan;
    f_denom = reshape(mean(f_mean_rpt,2,'omitnan'), [nPat_rpt, nDwi_rpt]);
    f_denom(abs(f_denom) < wcv_denom_floor) = nan;
    f_wCV = reshape(std(f_mean_rpt,0,2,'omitnan'), [nPat_rpt, nDwi_rpt]) ./ f_denom;
    f_wCV(n_rpt<2, :) = nan;
    dstar_denom = reshape(mean(dstar_mean_rpt,2,'omitnan'), [nPat_rpt, nDwi_rpt]);
    dstar_denom(abs(dstar_denom) < wcv_denom_floor) = nan;
    dstar_wCV = reshape(std(dstar_mean_rpt,0,2,'omitnan'), [nPat_rpt, nDwi_rpt]) ./ dstar_denom;
    dstar_wCV(n_rpt<2, :) = nan;
end

fprintf('  --- SECTION 2: Load Clinical Outcome Data ---\n');
% =========================================================================
% LOAD CLINICAL OUTCOME DATA FROM SPREADSHEET
% =========================================================================
% Links imaging biomarkers to clinical endpoints.  The primary endpoint is
% locoregional failure (LF) — confirmed progression of disease at or near
% the treated pancreatic tumour site.  This is the most clinically relevant
% endpoint for evaluating whether DWI biomarkers can predict treatment
% response to radiation therapy.
%
% The spreadsheet provides:
%   - Local/regional failure status (binary: 0=local control, 1=failure)
%   - Dates for failure, censoring, RT start/stop (used to compute
%     time-to-event in days from end of RT)
%   - Cause of death (for competing risk classification)
clinical_data_sheet = fullfile(dataloc, config_struct.clinical_data_sheet);

% Determine cause-of-death column name from config (default: 'CauseOfDeath')
if isfield(config_struct, 'cause_of_death_column')
    cod_column = config_struct.cause_of_death_column;
else
    cod_column = 'CauseOfDeath';
end

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
    T.CauseOfDeath = repmat({''}, size(id_list));

    T_prop = struct();
    T_prop.VariableNames = {'Pat', 'MRN', 'LocalOrRegionalFailure', 'LocoregionalFailureDateOfLocalOrRegionalFailure', 'LocalFailureDateOfLocalFailureOrCensor', 'RegionalFailureDateOfRegionalFailureOrCensor', 'RTStartDate', 'RTStopDate', 'CauseOfDeath'};
    T.Properties = T_prop;
else
    if isfield(config_struct, 'clinical_sheet_name')
        T = readtable(clinical_data_sheet,'Sheet', config_struct.clinical_sheet_name);
    else
        T = readtable(clinical_data_sheet);
    end
    % Rename the configured cause-of-death column to the canonical name
    if ~strcmp(cod_column, 'CauseOfDeath') && ismember(cod_column, T.Properties.VariableNames)
        T.Properties.VariableNames{strcmp(T.Properties.VariableNames, cod_column)} = 'CauseOfDeath';
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

    % Strip embedded single quotes that Excel may add to text cells,
    % matching the same normalization applied in load_dwi_data.m (line 263).
    pat_raw = T.Pat;
    if iscategorical(pat_raw), pat_raw = cellstr(pat_raw); end
    pat_raw = strrep(pat_raw, '''', '');
    pat_normalized = strrep(pat_raw, '_', '-');
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
                % Classify as competing risk if death was not from
                % pancreatic/biliary cancer.  The previous check for
                % 'cancer' anywhere in the string incorrectly treated
                % deaths from unrelated cancers (e.g., lung cancer) as
                % pancreatic cancer deaths.
                is_panc_cancer = ~isempty(regexp(cod_lower, 'pancrea|biliar|disease\s+progression', 'once'));
                if lf(j) == 0 && ~isempty(cod) && ~is_unknown && ~is_panc_cancer
                    lf(j) = 2; % Competing risk: non-pancreatic-cancer death without LF
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
        '%s column not found in clinical spreadsheet. Competing risks cannot be identified; CSH survival analysis may be biased.', cod_column);
end

% =========================================================================
% COMPUTE TIME-TO-EVENT RELATIVE TO END OF RADIATION THERAPY
% =========================================================================
% Time origin is the RT end date (not start date) because:
%   1. Diffusion changes are expected to evolve after the last RT fraction
%   2. Using RT-end as time zero makes time-to-event comparable across
%      patients regardless of RT duration (e.g., 5-fraction SBRT vs 25-fraction)
%   3. Events before RT-end would have negative times, which are nonsensical
%
% total_time: days from RT-end to local/regional failure date (for LF patients)
% total_follow_up_time: days from RT-end to last known follow-up (for censored patients)
if exist('OCTAVE_VERSION', 'builtin')
    total_time = m_total_time;
    total_follow_up_time = m_total_follow_up_time;
else
    total_time = days(lf_date - rtenddate);
    total_follow_up_time = days(censor_date - rtenddate);
end

fprintf('  --- DEEP LEARNING RIGOR AUDIT ---\n');
% =========================================================================
% DEEP LEARNING PROVENANCE CHECK
% =========================================================================
% When using DL-denoised data (dnCNN or IVIMnet), patients used to TRAIN
% the DL model must be excluded from the ANALYSIS cohort.  If a patient's
% images were used to train the denoiser, the denoised output for that
% patient benefits from having "seen" itself — an insidious form of data
% leakage that inflates apparent biomarker precision and downstream AUC.
% The provenance manifest records which patient IDs were in each DL
% training set.  This is checked again at the LOOCV level in
% metrics_stats_predictive.m for defense-in-depth.
manifest_file = fullfile(config_struct.dataloc, 'dl_validation_manifest.mat');
dl_provenance = load_dl_provenance(manifest_file);

dwi_type_names = {'Standard', 'dnCNN', 'IVIMnet'};
if isfield(config_struct, 'output_folder')
    output_folder = config_struct.output_folder;
else
    timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
    output_folder = fullfile(config_struct.dataloc, sprintf('saved_files_%s', timestamp_str));
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

% =========================================================================
% DEFINE COHORT INCLUSION CRITERIA AND OUTLIER POLICY
% =========================================================================
% Set up common variables for the next scripts
dtype_label = dwi_type_names{dtype};

% Exclude patients missing baseline (Fx1) imaging.  A valid baseline
% requires both a segmented GTV volume and a measurable ADC value.
% Without baseline data, percent-change metrics are undefined and
% survival models cannot adjust for pre-treatment tumour characteristics.
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

% =========================================================================
% PREPARE OUTPUT ARRAYS (PREFIX m_ = "metrics-stage" working copies)
% =========================================================================
% Create mutable copies of all arrays.  These will be filtered by
% valid_baseline and outlier masks below.  The m_ prefix distinguishes
% these working copies from the raw summary_metrics inputs.
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

% Pad arrays with NaN columns if fewer timepoints were collected than the
% maximum expected.  This ensures consistent matrix dimensions across
% patients with differing scan schedules (e.g., patients who missed late
% fractions or whose post-RT scan was not acquired).
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
        % Exclude competing risks (lf==2) from denominator: these patients
        % died of non-cancer causes before LF could be observed, so they
        % should not dilute the LF rate used for informative censoring assessment.
        lf_rate_excl = 100 * sum(lf_excluded == 1) / max(1, sum(lf_excluded <= 1 & isfinite(lf_excluded)));
        lf_rate_incl = 100 * sum(lf_included == 1) / max(1, sum(lf_included <= 1 & isfinite(lf_included)));
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

% Use fixed, physiologically motivated epsilon values to prevent inflated
% percent changes when baseline values are near zero.  Fixed thresholds
% improve reproducibility across cohorts (previously used data-adaptive
% 1% of IQR, which varied with cohort composition).
% Values correspond to ~1% of typical physiological range:
%   ADC: 0.001-0.003 mm^2/s → eps = 1e-5
%   D:   0.001-0.003 mm^2/s → eps = 1e-5
%   D*:  0.005-0.050 mm^2/s → eps = 5e-5
adc_eps  = 1e-5;
d_eps    = 1e-5;
dstar_eps = 5e-5;
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
time_labels = [arrayfun(@(x) sprintf('Fx%d', x), 1:(nTp-1), 'UniformOutput', false), {'Post'}];

% Restore global state modified at the top of this function
diary off;
warning(w_state_baseline);
set(0, 'DefaultFigureVisible', prev_fig_vis);
end
