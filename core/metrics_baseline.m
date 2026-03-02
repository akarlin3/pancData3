function [m_lf, m_total_time, m_total_follow_up_time, m_gtv_vol, m_adc_mean, m_d_mean, m_f_mean, m_dstar_mean, m_id_list, m_mrn_list, m_d95_gtvp, m_v50gy_gtvp, m_data_vectors_gtvp, lf_group, valid_pts, ADC_abs, D_abs, f_abs, Dstar_abs, ADC_pct, D_pct, f_pct, Dstar_pct, nTp, metric_sets, set_names, time_labels, dtype_label, dl_provenance] = metrics_baseline(data_vectors_gtvp, data_vectors_gtvn, summary_metrics, config_struct)
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
if exist('OCTAVE_VERSION', 'builtin')
    adc_wCV = squeeze(nanstd(adc_mean_rpt,0,2))./squeeze(nanmean(adc_mean_rpt,2));
    adc_wCV(n_rpt<2, :) = nan;
    adc_wCV_sub = squeeze(nanstd(adc_sub_rpt,0,2))./squeeze(nanmean(adc_sub_rpt,2));
    adc_wCV_sub(n_rpt<2, :) = nan;
    d_wCV = squeeze(nanstd(d_mean_rpt,0,2))./squeeze(nanmean(d_mean_rpt,2));
    d_wCV(n_rpt<2, :) = nan;
    f_wCV = squeeze(nanstd(f_mean_rpt,0,2))./squeeze(nanmean(f_mean_rpt,2));
    f_wCV(n_rpt<2, :) = nan;
    dstar_wCV = squeeze(nanstd(dstar_mean_rpt,0,2))./squeeze(nanmean(dstar_mean_rpt,2));
    dstar_wCV(n_rpt<2, :) = nan;
else
    adc_wCV = squeeze(std(adc_mean_rpt,0,2,'omitnan'))./squeeze(mean(adc_mean_rpt,2,'omitnan'));
    adc_wCV(n_rpt<2, :) = nan;
    adc_wCV_sub = squeeze(std(adc_sub_rpt,0,2,'omitnan'))./squeeze(mean(adc_sub_rpt,2,'omitnan'));
    adc_wCV_sub(n_rpt<2, :) = nan;
    d_wCV = squeeze(std(d_mean_rpt,0,2,'omitnan'))./squeeze(mean(d_mean_rpt,2,'omitnan'));
    d_wCV(n_rpt<2, :) = nan;
    f_wCV = squeeze(std(f_mean_rpt,0,2,'omitnan'))./squeeze(mean(f_mean_rpt,2,'omitnan'));
    f_wCV(n_rpt<2, :) = nan;
    dstar_wCV = squeeze(std(dstar_mean_rpt,0,2,'omitnan'))./squeeze(mean(dstar_mean_rpt,2,'omitnan'));
    dstar_wCV(n_rpt<2, :) = nan;
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

    for j = 1:length(id_list)
        i_find = find(~cellfun(@isempty, strfind(pat_normalized, id_list_normalized{j})));
        if ~isempty(i_find)
            i_find = i_find(1);
            lf(j) = T.LocalOrRegionalFailure(i_find);
            if ismember('CauseOfDeath', T.Properties.VariableNames)
                cod = T.CauseOfDeath{i_find};
                if lf(j) == 0 && ~isempty(cod) && isempty(strfind(lower(cod), 'cancer'))
                    lf(j) = 2; % Competing risk
                end
            end
            lf_date(j) = T.LocoregionalFailureDateOfLocalOrRegionalFailure(i_find);
            censor_date(j) = max(T.LocalFailureDateOfLocalFailureOrCensor(i_find),T.RegionalFailureDateOfRegionalFailureOrCensor(i_find));
            rtstartdate(j) = T.RTStartDate(i_find);
            rtenddate(j) = T.RTStopDate(i_find);
        end
    end
end

if exist('OCTAVE_VERSION', 'builtin')
    total_time = m_total_time;
    total_follow_up_time = m_total_follow_up_time;
else
    total_time = days(lf_date - rtenddate);
    total_follow_up_time = days(censor_date - rtenddate);
end

fprintf('  --- DEEP LEARNING RIGOR AUDIT ---\n');
manifest_file = fullfile(pwd, 'dl_validation_manifest.mat');
dl_provenance = load_dl_provenance(manifest_file);

dwi_type_names = {'Standard', 'dnCNN', 'IVIMnet'};
if isfield(config_struct, 'output_folder')
    output_folder = config_struct.output_folder;
else
    output_folder = fullfile(pwd, 'saved_figures');
end
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
diary_file = fullfile(output_folder, 'metrics_output.txt');
if exist(diary_file, 'file'), delete(diary_file); end
diary(diary_file);
warning('off', 'stats:glmfit:IterationLimit');
set(0, 'DefaultFigureVisible', 'off');  

if ~exist('nTp', 'var') || isempty(nTp)
    nTp = size(adc_mean, 2);
end

fprintf('  --- Baseline Data Completeness Check ---\n');
for j = 1:length(id_list)
    baseline_adc = adc_mean(j, 1, 1);
    baseline_vol = gtv_vol(j, 1);
end

% Set up common variables for the next scripts
dtype = config_struct.dwi_types_to_run;
dtype_label = dwi_type_names{dtype};
exclude_missing_baseline = true;
valid_baseline = ~isnan(gtv_vol(:,1)) & ~isnan(adc_mean(:,1,dtype));
exclude_outliers = true;

baseline_metrics_oi   = {adc_mean(:,1,dtype), d_mean(:,1,dtype), f_mean(:,1,dtype), dstar_mean(:,1,dtype)};
is_outlier = false(size(lf));
for metric_idx = 1:numel(baseline_metrics_oi)
    col = baseline_metrics_oi{metric_idx};
    col_clean = col(~isnan(col));
    if numel(col_clean) < 3, continue; end
    med_val = median(col_clean);
    iqr_val = iqr(col_clean);
    if iqr_val == 0, continue; end
    lower_fence = med_val - 3 * iqr_val;
    upper_fence = med_val + 3 * iqr_val;
    outlier_flags = (col < lower_fence | col > upper_fence) & ~isnan(col);
    is_outlier = is_outlier | outlier_flags;
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

epsilon = 1e-4;
ADC_pct = ((ADC_abs - ADC_abs(:,1)) ./ (ADC_abs(:,1) + epsilon)) * 100;
D_pct   = ((D_abs - D_abs(:,1)) ./ (D_abs(:,1) + epsilon)) * 100;
f_pct   = (f_abs - f_abs(:,1));
Dstar_pct = ((Dstar_abs - Dstar_abs(:,1)) ./ (Dstar_abs(:,1) + epsilon)) * 100;

valid_pts = isfinite(m_lf);
lf_group = m_lf(valid_pts);

% Re-organize metrics into 4 distinct Sets:
metric_sets = {
    {ADC_abs, D_abs, f_abs, Dstar_abs}, ...          % Set 1
    {ADC_pct, D_pct, f_pct, Dstar_pct} ...          % Set 2
};

set_names = {
    {'ADC Absolute', 'D Absolute', 'f Absolute', 'D* Absolute'}, ...
    {'\Delta ADC (%)', '\Delta D (%)', '\Delta f (abs)', '\Delta D* (%)'} ...
};
time_labels = {'Fx1', 'Fx2', 'Fx3', 'Fx4', 'Fx5', 'Post'};
end
