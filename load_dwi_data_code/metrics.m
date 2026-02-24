%% metrics.m — Pancreatic Cancer DWI/IVIM Treatment Response Analysis
% =========================================================================
% This script computes quantitative imaging biomarkers from diffusion-
% weighted MRI (DWI) and intravoxel incoherent motion (IVIM) data acquired
% during fractionated radiotherapy for pancreatic cancer.
%
% Pipeline overview:
%   1. Repeatability metrics (within-subject coefficient of variation, wCV)
%   2. Clinical outcome data loading (local failure status, time-to-event)
%   3. Longitudinal metric visualization (absolute & percent-change plots)
%   4. Dose–volume histogram (DVH) metrics for DWI-defined sub-volumes
%   5. Univariate analysis (Wilcoxon rank-sum boxplots: LC vs LF at each fraction)
%   6. Multiple-comparison corrections (FDR / Benjamini-Hochberg,
%      Holm-Bonferroni)
%   7. Elastic Net–regularized feature selection (α=0.5, L1+L2, 5-fold CV)
%   8. ROC analysis with Youden's J optimal cutoff
%   9. Multivariable logistic regression (with Firth penalized fallback)
%  10. Leave-pair-out (LPOCV) cross-validation
%  11. Kaplan-Meier survival analysis stratified by imaging biomarker
%  12. Representative parametric ADC maps (responder vs non-responder)
%  13. Sanity-check panels (volume change, heterogeneity, noise floor)
%
% Inputs (expected in workspace from upstream load_dwi_data pipeline):
%   adc_mean_rpt, adc_sub_rpt, d_mean_rpt, f_mean_rpt, dstar_mean_rpt
%       — Repeat-scan measurements [n_patients x n_repeats x n_dwi_types]
%   n_rpt              — Number of repeat scans per patient
%   adc_mean, d_mean, f_mean, dstar_mean
%       — Longitudinal metric arrays [n_patients x n_timepoints x n_types]
%   gtv_vol            — GTV volumes [n_patients x n_timepoints]
%   id_list, mrn_list  — Patient identifiers
%   d95_gtvp, v50gy_gtvp, data_vectors_gtvp — Dose and voxel-level data
%   dataloc            — Root data directory path
%   nTp                — Number of timepoints (fractions + post-RT)
%
% Outputs:
%   Figures saved to ./saved_figures/
%   CSV tables saved to dataloc
%   Diary log saved to ./saved_figures/metrics_output.txt
%
% Dependencies: MATLAB Statistics and Machine Learning Toolbox,
%               load_untouch_nii (NIfTI I/O), Image Processing Toolbox
% =========================================================================

%% ========================================================================
%  WORKSPACE CLEANUP
%  test.m runs as a script and leaks variables into the base workspace that
%  shadow built-in MATLAB functions used here (e.g. lines(), test(),
%  training()). Clear them defensively before metrics begins.
% =========================================================================
cleanup_names = {'lines', 'test', 'training', 'text', 'labels', 'cvp', 'k', 'n_pass', 'n_fail', ...
                 'metrics_code', 'vis_code'};
for ci = 1:numel(cleanup_names)
    if exist(cleanup_names{ci}, 'var'), clear(cleanup_names{ci}); end
end
clear cleanup_names ci

%% ========================================================================
%  SECTION 1: Repeatability — Within-Subject Coefficient of Variation (wCV)
% =========================================================================
% wCV quantifies scan-rescan reproducibility. It is defined as:
%   wCV = SD_within / Mean_within
% Only patients with >=2 Fx1 repeat scans (n_rpt >= 2) contribute.
% Results are reported as mean +/- SD across the cohort for each DWI type
% (Standard, dnCNN-denoised, IVIMnet).

% --- ADC repeatability (whole-GTV mean) ---
adc_wCV = squeeze(nanstd(adc_mean_rpt,[],2))./squeeze(nanmean(adc_mean_rpt,2));
adc_wCV(n_rpt<2) = nan;

% --- ADC repeatability (low-ADC subregion mean) ---
adc_wCV_sub = squeeze(nanstd(adc_sub_rpt,[],2))./squeeze(nanmean(adc_sub_rpt,2));
adc_wCV_sub(n_rpt<2) = nan;

% Print ADC wCV: whole-GTV and subregion, for Standard and dnCNN pipelines
fprintf('ADC   wCV (whole gtv) = %2.2f%s%2.2f%%, dnCNN: %2.2f%s%2.2f%%\n',nanmean(adc_wCV(:,1))*100,char(177),nanstd(adc_wCV(:,1))*100,nanmean(adc_wCV(:,2))*100,char(177),nanstd(adc_wCV(:,2))*100)
fprintf('ADC   wCV (subregion) = %2.2f%s%2.2f%%, dnCNN: %2.2f%s%2.2f%%\n',nanmean(adc_wCV_sub(:,1))*100,char(177),nanstd(adc_wCV_sub(:,1))*100,nanmean(adc_wCV_sub(:,2))*100,char(177),nanstd(adc_wCV_sub(:,2))*100)

% --- IVIM parameter repeatability (D, f, D*) ---
% Each has 3 DWI types: Standard, dnCNN, IVIMnet
d_wCV = squeeze(nanstd(d_mean_rpt,[],2))./squeeze(nanmean(d_mean_rpt,2));
d_wCV(repmat(n_rpt,[1,3])<2) = nan;  % Mask patients with <2 repeats

f_wCV = squeeze(nanstd(f_mean_rpt,[],2))./squeeze(nanmean(f_mean_rpt,2));
f_wCV(repmat(n_rpt,[1,3])<2) = nan;

dstar_wCV = squeeze(nanstd(dstar_mean_rpt,[],2))./squeeze(nanmean(dstar_mean_rpt,2));
dstar_wCV(repmat(n_rpt,[1,3])<2) = nan;

% Print IVIM wCV for all three processing pipelines (char(177) = ± symbol)
fprintf('D     wCV = %2.2f%s%2.2f%%, dnCNN: %2.2f%s%2.2f%%, IVIMnet: %2.2f%s%2.2f%%\n',nanmean(d_wCV(:,1))*100,char(177),nanstd(d_wCV(:,1))*100,nanmean(d_wCV(:,2))*100,char(177),nanstd(d_wCV(:,2))*100,nanmean(d_wCV(:,3))*100,char(177),nanstd(d_wCV(:,3))*100)
fprintf('f     wCV = %2.2f%s%2.2f%%, dnCNN: %2.2f%s%2.2f%%, IVIMnet: %2.2f%s%2.2f%%\n',nanmean(f_wCV(:,1))*100,char(177),nanstd(f_wCV(:,1))*100,nanmean(f_wCV(:,2))*100,char(177),nanstd(f_wCV(:,2))*100,nanmean(f_wCV(:,3))*100,char(177),nanstd(f_wCV(:,3))*100)
fprintf('Dstar wCV = %2.2f%s%2.2f%%, dnCNN: %2.2f%s%2.2f%%, IVIMnet: %2.2f%s%2.2f%%\n',nanmean(dstar_wCV(:,1))*100,char(177),nanstd(dstar_wCV(:,1))*100,nanmean(dstar_wCV(:,2))*100,char(177),nanstd(dstar_wCV(:,2))*100,nanmean(dstar_wCV(:,3))*100,char(177),nanstd(dstar_wCV(:,3))*100)

%% ========================================================================
%  SECTION 2: Load Clinical Outcome Data
% =========================================================================
% Reads the master clinical spreadsheet to extract:
%   - Local/regional failure status (binary: 0=LC, 1=LF)
%   - Failure and censor dates for time-to-event analysis
%   - RT start/stop dates for computing post-treatment intervals
clinical_data_sheet = [dataloc 'MASTER_pancreas_DWIanalysis.xlsx'];
T = readtable(clinical_data_sheet,'Sheet','Clin List_MR');

% Pre-allocate outcome arrays matched to the patient list
lf = nan(size(mrn_list));                % Local failure flag (0 or 1)
lf_date = NaT(size(mrn_list));           % Date of locoregional failure
censor_date = NaT(size(mrn_list));       % Last known follow-up / censor date
rtstartdate = NaT(size(mrn_list));       % Radiation therapy start date
rtenddate = NaT(size(mrn_list));         % Radiation therapy end date

% Match each study patient (by MRN) to the clinical spreadsheet
for j=1:length(mrn_list)
    i_find = find(ismember(T.MRN,str2num(mrn_list{j})));
    if ~isempty(i_find)
        i_find = i_find(1);  % Use first match if duplicates exist
        lf(j) = T.LocalOrRegionalFailure(i_find);
        lf_date(j) = T.LocoregionalFailureDateOfLocalOrRegionalFailure(i_find);
        % Take the later of local or regional censor dates
        censor_date(j) = max(T.LocalFailureDateOfLocalFailureOrCensor(i_find),T.RegionalFailureDateOfRegionalFailureOrCensor(i_find));
        rtstartdate(j) = T.RTStartDate(i_find);
        rtenddate(j) = T.RTStopDate(i_find);
    end
end

% Compute time-to-event variables (days from RT completion)
total_time = days(lf_date - rtenddate);               % Time to failure
total_follow_up_time = days(censor_date - rtenddate);  % Follow-up duration

% Summary statistics for the cohort
fprintf('LR observed in          %d / %d patients (%2.2f%%)\n',numel(lf(lf==1)),numel(lf(isfinite(lf))),100*numel(lf(lf==1))/numel(lf(isfinite(lf))));
fprintf('Median follow-up time = %2d days (%d - %d)\n',nanmedian(total_follow_up_time),nanmin(total_follow_up_time),nanmax(total_follow_up_time));
fprintf('Median time to LF     = %2d days (%d - %d)\n',nanmedian(total_time(lf==1)),nanmin(total_time(lf==1)),nanmax(total_time(lf==1)));

%% ========================================================================
%  SECTION 2.5: Deep Learning Rigor & Isolation Audit
% =========================================================================
% To prevent data leakage, deep learning models (dnCNN, IVIMnet) must NOT
% have been trained or fine-tuned on any patients held out during cross-
% validation. This section checks for a rigor manifest.

fprintf('\n--- DEEP LEARNING RIGOR AUDIT ---\n');
dl_provenance = struct();
manifest_file = fullfile(pwd, 'dl_validation_manifest.mat');

if exist('dl_provenance_workspace', 'var')
    dl_provenance = dl_provenance_workspace;
    fprintf('  Using DL provenance manifest from workspace.\n');
elseif exist(manifest_file, 'file')
    load(manifest_file, 'dl_provenance');
    fprintf('  Loaded DL provenance manifest: %s\n', manifest_file);
else
    fprintf('  [WARNING] No DL provenance manifest found. Assuming Standard weights.\n');
    fprintf('  [CRITICAL] Ensure dnCNN/IVIMnet weights are from 100%% independent datasets.\n');
    % Initialize empty manifest (no patients flagged as leaky)
    dl_provenance.dncnn_train_ids = {};
    dl_provenance.ivimnet_train_ids = {};
end
%% ========================================================================
%  SECTION 3: Pipeline Setup — Output Folder, Diary, and Figure Defaults
% =========================================================================
% Configure DWI processing type labels, create output directory for saved
% figures, start a diary log, and suppress figure display for batch mode.
dwi_type_names = {'Standard', 'dnCNN', 'IVIMnet'};
output_folder = fullfile(pwd, 'saved_figures');
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
diary_file = fullfile(output_folder, 'metrics_output.txt');
if exist(diary_file, 'file'), delete(diary_file); end
diary(diary_file);
warning('off', 'stats:glmfit:IterationLimit');
set(0, 'DefaultFigureVisible', 'off');  % Hide figures (save only)

% Robust calculation of number of timepoints (usually 6: Fx1-Fx5 + Post)
if ~exist('nTp', 'var') || isempty(nTp)
    nTp = size(adc_mean, 2);
end

%% ========================================================================
%  SECTION 4: Baseline Data Completeness Check
% =========================================================================
% Report patients missing Fx1 (pretreatment) ADC or GTV volume.
% Missing baseline data disqualifies a patient from percent-change analyses.
for j = 1:length(id_list)
    baseline_adc = adc_mean(j, 1, 1);
    baseline_vol = gtv_vol(j, 1);

    if isnan(baseline_adc) || isnan(baseline_vol)
        fprintf('Patient %s (MRN %s): Baseline ADC=%.4g, Baseline Vol=%.4g\n', ...
            id_list{j}, mrn_list{j}, baseline_adc, baseline_vol);

        if isempty(dwi_locations{j,1,1})
            fprintf('  -> No DWI data at Fx1\n');
        end
        if isempty(gtv_locations{j,1,1})
            fprintf('  -> No GTV mask at Fx1\n');
        end
    end
end

%% ========================================================================
%  SECTION 5: Longitudinal Metric Plotting (per DWI processing type)
% =========================================================================
% For each DWI type (Standard, dnCNN, IVIMnet):
%   - Optionally exclude patients with missing baseline data
%   - Optionally exclude outlier patients (Fx1 values > 3 IQR from median)
%   - Compute absolute values and percent-change from Fx1 for ADC, D, f, D*
%   - Generate 2x4 spaghetti + population-mean plots:
%     Top row: absolute values; Bottom row: percent change from baseline
for dtype = 1:3
dtype_label = dwi_type_names{dtype};
fprintf('\n=== Processing DWI Type %d: %s ===\n', dtype, dtype_label);

% Toggle baseline exclusion: if true, patients without Fx1 data are removed
exclude_missing_baseline = true;

% Identify patients who have both a valid Fx1 GTV volume and Fx1 ADC
valid_baseline = ~isnan(gtv_vol(:,1)) & ~isnan(adc_mean(:,1,dtype));
fprintf('Found %d/%d patients with missing baseline data\n', ...
    sum(~valid_baseline), length(valid_baseline));

% Toggle outlier exclusion: if true, patients with extreme Fx1 biomarker
% values (> 3 IQR from the cohort median) are removed from the analysis.
% Set to false to include all patients regardless of outlier status.
exclude_outliers = true;

% Identify outlier patients using a 3 IQR fence on each Fx1 biomarker.
% A patient is flagged as an outlier if ANY metric falls outside the fence.
baseline_metrics_oi   = {adc_mean(:,1,dtype), d_mean(:,1,dtype), ...
                          f_mean(:,1,dtype),   dstar_mean(:,1,dtype)};
baseline_names_oi     = {'ADC', 'D', 'f', 'D*'};
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
    outlier_indices = find(outlier_flags);
    for patient_idx = outlier_indices'
        fprintf('  Outlier: Patient %s  Fx1  %s = %.4g  (median=%.4g, 3*IQR fence=[%.4g, %.4g])\n', ...
            id_list{patient_idx}, baseline_names_oi{metric_idx}, col(patient_idx), ...
            med_val, lower_fence, upper_fence);
    end
    is_outlier = is_outlier | outlier_flags;
end
fprintf('Found %d/%d outlier patients at Fx1 (>3 IQR from median)\n', ...
    sum(is_outlier), length(is_outlier));
non_outlier = ~is_outlier;

% Create local (m_) copies of all arrays so exclusions do not alter
% upstream workspace variables used by subsequent DWI-type iterations
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
% Defensive padding: ensure longitudinal arrays have nTp columns to prevent 
% crashes during Post-RT (column 6) analysis when dose was only calculated 
% for Fractions 1-5.
long_vars = {'m_gtv_vol', 'm_d95_gtvp', 'm_v50gy_gtvp', 'dmean_gtvp'};
for lv = 1:numel(long_vars)
    vname = long_vars{lv};
    if exist(vname, 'var')
        curr_v = eval(vname);
        % curr_v is [nPat x nTimepoints]
        if size(curr_v, 2) < nTp
            % Pad with NaNs
            padding = nan(size(curr_v, 1), nTp - size(curr_v, 2));
            curr_v = [curr_v, padding];
            assignin('caller', vname, curr_v);
        end
    end
end

if exclude_missing_baseline
    fprintf('Excluding patients with missing baseline data\n');
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
    % Index the non_outlier mask relative to the patients still present
    % after baseline exclusion (or all patients if no baseline exclusion).
    if exclude_missing_baseline
        non_outlier_current = non_outlier(valid_baseline);
    else
        non_outlier_current = non_outlier;
    end
    n_excluded_oi = sum(~non_outlier_current);
    if n_excluded_oi > 0
        fprintf('Excluding %d outlier patient(s) from analysis\n', n_excluded_oi);
    end
    m_lf                   = m_lf(non_outlier_current);
    m_total_time           = m_total_time(non_outlier_current);
    m_total_follow_up_time = m_total_follow_up_time(non_outlier_current);
    m_gtv_vol              = m_gtv_vol(non_outlier_current,:);
    m_adc_mean             = m_adc_mean(non_outlier_current,:,:);
    m_d_mean               = m_d_mean(non_outlier_current,:,:);
    m_f_mean               = m_f_mean(non_outlier_current,:,:);
    m_dstar_mean           = m_dstar_mean(non_outlier_current,:,:);
    m_id_list              = m_id_list(non_outlier_current);
    m_mrn_list             = m_mrn_list(non_outlier_current);
    m_d95_gtvp             = m_d95_gtvp(non_outlier_current,:);
    m_v50gy_gtvp           = m_v50gy_gtvp(non_outlier_current,:);
    m_data_vectors_gtvp    = m_data_vectors_gtvp(non_outlier_current,:,:);
end

% Extract the 2D matrices [n_patients x n_timepoints] for the selected type
ADC_abs = m_adc_mean(:,:,dtype);
D_abs   = m_d_mean(:,:,dtype);
f_abs   = m_f_mean(:,:,dtype);
Dstar_abs = m_dstar_mean(:,:,dtype);

% Calculate percent change relative to pretreatment (Fx1)
% For f (bounded 0-1), compute the absolute difference instead of percent change
% For D*, add a small smoothing constant to the denominator to prevent noise spikes
epsilon = 1e-4;

ADC_pct = ((ADC_abs - ADC_abs(:,1)) ./ (ADC_abs(:,1) + epsilon)) * 100;
D_pct   = ((D_abs - D_abs(:,1)) ./ (D_abs(:,1) + epsilon)) * 100;
f_pct   = (f_abs - f_abs(:,1));
Dstar_pct = ((Dstar_abs - Dstar_abs(:,1)) ./ (Dstar_abs(:,1) + epsilon)) * 100;

% Group data for easy iteration
metrics_abs = {ADC_abs, D_abs, f_abs, Dstar_abs};
metrics_pct = {ADC_pct, D_pct, f_pct, Dstar_pct};
metric_names = {'ADC', 'D', 'f', 'D*'};
metric_units = {'mm^2/s', 'mm^2/s', 'Fraction', 'mm^2/s'};
x_labels = {'Fx1', 'Fx2', 'Fx3', 'Fx4', 'Fx5', 'Post'};
x_vals = 1:nTp;

% Create figure
figure('Name', ['Longitudinal Mean Metrics — ' dtype_label], 'Position', [100, 100, 1400, 700]);

for i = 1:4
    % -------------------------------------------------------------------
    % TOP ROW: Absolute Mean Values
    % -------------------------------------------------------------------
    subplot(2, 4, i);
    hold on;
    
    dat = metrics_abs{i};
    % Calculate population mean and standard error of the mean (SEM)
    pop_mean = nanmean(dat, 1);
    pop_se   = nanstd(dat, 0, 1) ./ sqrt(sum(~isnan(dat), 1));
    
    % Plot individual patient trajectories (spaghetti plot)
    plot(x_vals, dat', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    
    % Plot population average with error bars
    errorbar(x_vals, pop_mean, pop_se, '-k', 'LineWidth', 2, ...
        'Marker', 'o', 'MarkerFaceColor', 'k');
    
    % Formatting
    set(gca, 'XTick', x_vals, 'XTickLabel', x_labels, 'FontSize', 10);
    title(['Mean ', metric_names{i}], 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(metric_units{i});
    xlim([0.5, nTp+0.5]);
    grid on; box on;
    
    % -------------------------------------------------------------------
    % BOTTOM ROW: Percent Change from Fx1
    % -------------------------------------------------------------------
    subplot(2, 4, i+4);
    hold on;
    
    dat_pct = metrics_pct{i};
    % Calculate population mean percent change and SEM
    pop_mean_pct = nanmean(dat_pct, 1);
    pop_se_pct   = nanstd(dat_pct, 0, 1) ./ sqrt(sum(~isnan(dat_pct), 1));
    
    % Plot individual patient trajectories
    plot(x_vals, dat_pct', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    
    % Plot population average with error bars
    errorbar(x_vals, pop_mean_pct, pop_se_pct, '-r', 'LineWidth', 2, ...
        'Marker', 's', 'MarkerFaceColor', 'r');
    
    % Add a baseline reference line at 0% change
    yline(0, 'k--', 'LineWidth', 1.5);
    
    % Formatting
    set(gca, 'XTick', x_vals, 'XTickLabel', x_labels, 'FontSize', 10);
    if strcmp(metric_names{i}, 'f')
        title(['\Delta ', metric_names{i}, ' (abs)'], 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Absolute Change');
    else
        title(['\Delta ', metric_names{i}, ' (%)'], 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('% Change from Fx1');
    end
    xlim([0.5, nTp+0.5]);
    grid on; box on;
end

sgtitle(['Longitudinal Evolution of DWI and IVIM Metrics (' dtype_label ')'], 'FontSize', 16, 'FontWeight', 'bold');
% Shift subplot axes down to create visual separation between super-title and subplot titles
subplot_scale = 0.92;
allAx = findall(gcf, 'Type', 'Axes');
for k = 1:numel(allAx)
    pos = allAx(k).Position;
    allAx(k).Position = [pos(1), pos(2) * subplot_scale, pos(3), pos(4) * subplot_scale];
end
saveas(gcf, fullfile(output_folder, ['Longitudinal_Mean_Metrics_' dtype_label '.png']));
close(gcf);

%% ========================================================================
%  SECTION 6: Target Coverage — Sub-Volume Dose Metrics (D95, V50)
% =========================================================================
% For each patient and timepoint, define treatment-resistant sub-volumes
% using DWI/IVIM thresholds (low ADC / low D / low f / low D* voxels),
% then compute DVH metrics within those sub-volumes:
%   D95  — minimum dose covering 95% of the sub-volume (5th percentile)
%   V50  — percentage of sub-volume receiving >= 50 Gy
%
% Thresholds chosen based on literature values for restricted diffusion
% in pancreatic adenocarcinoma.
adc_thresh = 0.00115;   % ADC threshold (mm^2/s) — restricted diffusion
d_thresh = 0.001;       % D (true diffusion) threshold (mm^2/s)
f_thresh = 0.1;         % f (perfusion fraction) threshold — low perfusion
dstar_thresh = 0.01;    % D* (pseudo-diffusion) threshold (mm^2/s)

% Initialize storage arrays [n_patients x n_timepoints]
d95_adc_sub = nan(length(m_id_list), nTp);
v50_adc_sub = nan(length(m_id_list), nTp);

d95_d_sub = nan(length(m_id_list), nTp);
v50_d_sub = nan(length(m_id_list), nTp);

d95_f_sub = nan(length(m_id_list), nTp);
v50_f_sub = nan(length(m_id_list), nTp);

d95_dstar_sub = nan(length(m_id_list), nTp);
v50_dstar_sub = nan(length(m_id_list), nTp);

for j = 1:length(m_id_list)
    % Find original patient index to correctly map to unfiltered gtv_locations
    j_orig = find(strcmp(id_list, m_id_list{j}));
    
    % Restrict sub-volume DVH calculations to the baseline (Fx1) scan only.
    % WARNING: Longitudinal sub-volume dose metrics are uncorrected for spatial 
    % shifts (tumor deformation/migration) and require Deformable Image 
    % Registration (DIR) for validity.
    for k = 1
        % Extract vectors for standard DWI processing (dwi_type = 1, rpi = 1)
        adc_vec = m_data_vectors_gtvp(j,k,1).adc_vector;
        d_vec   = m_data_vectors_gtvp(j,k,1).d_vector;
        f_vec   = m_data_vectors_gtvp(j,k,1).f_vector;
        dstar_vec = m_data_vectors_gtvp(j,k,1).dstar_vector;
        dose_vec  = m_data_vectors_gtvp(j,k,1).dose_vector;
        
        % Proceed only if both dose and DWI data exist for this timepoint
        if ~isempty(dose_vec) && ~isempty(adc_vec)
            
            % Load 3D GTV mask to apply spatial operations
            gtv_mat = gtv_locations{j_orig, k, 1};
            has_3d = false;
            if ~isempty(gtv_mat)
                % Convert mapped network drive paths if necessary
                gtv_mat = strrep(gtv_mat, '/Volumes/aliottae/pancreas_dwi/', dataloc);
                gtv_mat = strrep(gtv_mat, '/', filesep);
                if exist(gtv_mat, 'file')
                    tmp = load(gtv_mat, 'Stvol3d');
                    gtv_mask_3d = tmp.Stvol3d;
                    if sum(gtv_mask_3d(:) == 1) == length(adc_vec)
                        has_3d = true;
                    end
                end
            end
            
            % Parameters for 3D morphological cleanup and statistical stability
            se = strel('sphere', 1);
            min_cc_voxels = 10;
            min_subvol_voxels = 100;
            
            % 1. ADC Sub-volume Coverage
            adc_mask_1d = adc_vec < adc_thresh;
            if has_3d
                vol_3d = false(size(gtv_mask_3d));
                vol_3d(gtv_mask_3d == 1) = adc_mask_1d;
                vol_3d = imclose(imopen(vol_3d, se), se);
                vol_3d = bwareaopen(vol_3d, min_cc_voxels);
                adc_mask_1d = vol_3d(gtv_mask_3d == 1);
            end
            dose_adc_sub = dose_vec(adc_mask_1d);
            if ~isempty(dose_adc_sub) && sum(adc_mask_1d) >= min_subvol_voxels
                d95_adc_sub(j,k) = prctile(dose_adc_sub, 5); % 5th percentile = D95
                v50_adc_sub(j,k) = sum(dose_adc_sub >= 50) / length(dose_adc_sub) * 100;
            end
            
            % 2. D Sub-volume Coverage
            d_mask_1d = d_vec < d_thresh;
            if has_3d
                vol_3d = false(size(gtv_mask_3d));
                vol_3d(gtv_mask_3d == 1) = d_mask_1d;
                vol_3d = imclose(imopen(vol_3d, se), se);
                vol_3d = bwareaopen(vol_3d, min_cc_voxels);
                d_mask_1d = vol_3d(gtv_mask_3d == 1);
            end
            dose_d_sub = dose_vec(d_mask_1d);
            if ~isempty(dose_d_sub) && sum(d_mask_1d) >= min_subvol_voxels
                d95_d_sub(j,k) = prctile(dose_d_sub, 5);
                v50_d_sub(j,k) = sum(dose_d_sub >= 50) / length(dose_d_sub) * 100;
            end
            
            % 3. f Sub-volume Coverage
            f_mask_1d = f_vec < f_thresh;
            if has_3d
                vol_3d = false(size(gtv_mask_3d));
                vol_3d(gtv_mask_3d == 1) = f_mask_1d;
                vol_3d = imclose(imopen(vol_3d, se), se);
                vol_3d = bwareaopen(vol_3d, min_cc_voxels);
                f_mask_1d = vol_3d(gtv_mask_3d == 1);
            end
            dose_f_sub = dose_vec(f_mask_1d);
            if ~isempty(dose_f_sub) && sum(f_mask_1d) >= min_subvol_voxels
                d95_f_sub(j,k) = prctile(dose_f_sub, 5);
                v50_f_sub(j,k) = sum(dose_f_sub >= 50) / length(dose_f_sub) * 100;
            end
            
            % 4. D* Sub-volume Coverage
            dstar_mask_1d = dstar_vec < dstar_thresh;
            if has_3d
                vol_3d = false(size(gtv_mask_3d));
                vol_3d(gtv_mask_3d == 1) = dstar_mask_1d;
                vol_3d = imclose(imopen(vol_3d, se), se);
                vol_3d = bwareaopen(vol_3d, min_cc_voxels);
                dstar_mask_1d = vol_3d(gtv_mask_3d == 1);
            end
            dose_dstar_sub = dose_vec(dstar_mask_1d);
            if ~isempty(dose_dstar_sub) && sum(dstar_mask_1d) >= min_subvol_voxels
                d95_dstar_sub(j,k) = prctile(dose_dstar_sub, 5);
                v50_dstar_sub(j,k) = sum(dose_dstar_sub >= 50) / length(dose_dstar_sub) * 100;
            end
            
        end
    end
end

%% ========================================================================
%  SECTION 7: Univariate Analysis — Metric Sets vs Local Failure (Wilcoxon Rank-Sum)
% =========================================================================
% For each metric × timepoint combination, run a Wilcoxon rank-sum test comparing
% local control (LC=0) vs local failure (LF=1). Results are displayed as
% multi-panel boxplot figures with p-values highlighted in red when < 0.05.
%
% Four metric sets are tested:
%   Set 1: Absolute DWI/IVIM values (ADC, D, f, D*)
%   Set 2: Percent change from Fx1 (ΔADC, ΔD, Δf, ΔD*)
%   Set 3: D95 for whole GTV and each DWI-defined resistant sub-volume
%   Set 4: V50 for whole GTV and each DWI-defined resistant sub-volume

% Exclude patients without clinical outcome data
valid_pts = isfinite(m_lf);
lf_group = m_lf(valid_pts);  % Binary outcome vector (0=LC, 1=LF)

% Re-organize metrics into 4 distinct Sets:
% 1. Absolute DWI values
% 2. Percent Change DWI values
% 3. D95 Dosimetry (Whole GTV + All Sub-volumes)
% 4. V50 Dosimetry (Whole GTV + All Sub-volumes)

metric_sets = {
    {ADC_abs, D_abs, f_abs, Dstar_abs}, ...          % Set 1
    {ADC_pct, D_pct, f_pct, Dstar_pct}, ...          % Set 2
    {m_d95_gtvp, d95_adc_sub, d95_d_sub, d95_f_sub, d95_dstar_sub}, ... % Set 3
    {m_v50gy_gtvp, v50_adc_sub, v50_d_sub, v50_f_sub, v50_dstar_sub}    % Set 4
};

set_names = {
    {'ADC Absolute', 'D Absolute', 'f Absolute', 'D* Absolute'}, ...
    {'\Delta ADC (%)', '\Delta D (%)', '\Delta f (abs)', '\Delta D* (%)'}, ...
    {'D95 Whole GTV', 'D95 Sub(ADC)', 'D95 Sub(D)', 'D95 Sub(f)', 'D95 Sub(D*)'}, ...
    {'V50 Whole GTV', 'V50 Sub(ADC)', 'V50 Sub(D)', 'V50 Sub(f)', 'V50 Sub(D*)'}
};

figure_titles = {
    '1. Absolute DWI/IVIM Metrics vs Local Failure (Wilcoxon Rank-Sum)', ...
    '2. Percent Change Metrics vs Local Failure (Wilcoxon Rank-Sum)', ...
    '3. Target Coverage (D95): Whole GTV vs Resistant Sub-volumes (Wilcoxon Rank-Sum)', ...
    '4. Target Coverage (V50): Whole GTV vs Resistant Sub-volumes (Wilcoxon Rank-Sum)'
};

time_labels = {'Fx1', 'Fx2', 'Fx3', 'Fx4', 'Fx5', 'Post'};

% Iterate through each set of metrics
for s = 1:length(metric_sets)
    current_metrics = metric_sets{s};
    current_names = set_names{s};
    
    % Create a large figure for this set
    fig = figure('Name', [figure_titles{s} ' — ' dtype_label], 'Position', [50, 50, 1600, 1000]);
    sgtitle([figure_titles{s} ' (' dtype_label ')'], 'FontSize', 16, 'FontWeight', 'bold');
    
    num_rows = length(current_metrics);
    % Determine columns based on available timepoints (Dosimetry usually has fewer columns)
    max_cols = 0;
    for m=1:num_rows
        max_cols = max(max_cols, size(current_metrics{m}, 2));
    end
    cols_to_plot = min(max_cols, length(time_labels));
    
    plot_idx = 1;
    
    for m = 1:num_rows
        metric_data = current_metrics{m};
        
        for tp = 1:cols_to_plot
            subplot(num_rows, cols_to_plot, plot_idx);
            
            % Check if data exists for this timepoint
            if tp > size(metric_data, 2)
                 axis off;
                 plot_idx = plot_idx + 1;
                 continue;
            end

            % Extract data for this fraction and valid patients
            y_raw = metric_data(valid_pts, tp);
            
            % Filter out NaNs for the statistical test
            has_data = ~isnan(y_raw);
            y = y_raw(has_data);
            g = lf_group(has_data);
            
            % Check if we have enough data in both groups to run Wilcoxon rank-sum test
            unique_groups = unique(g);
            if length(unique_groups) > 1 && length(y) > 2
                % Compute Wilcoxon rank-sum quietly
                p = ranksum(y(g==0), y(g==1));
                
                % Plot boxplot
                boxplot(y, g, 'Labels', {'LC (0)', 'LF (1)'});
                
                % Highlight significant p-values in red
                title_str = sprintf('%s - %s\np = %.3f', current_names{m}, time_labels{tp}, p);
                if p < 0.05
                    title(title_str, 'Color', 'r', 'FontWeight', 'bold');
                else
                    title(title_str, 'Color', 'k', 'FontWeight', 'normal');
                end
            else
                % Not enough data to compare
                title(sprintf('%s - %s\n(Insufficient Data)', current_names{m}, time_labels{tp}), 'FontSize', 8);
                axis off; 
            end
            
            if m == num_rows
                xlabel('Outcome');
            end
            if tp == 1
                ylabel(current_names{m});
            end
            
            grid on;
            plot_idx = plot_idx + 1;
        end
    end
    % Shift subplot axes down to create visual separation between super-title and subplot titles
    allAx = findall(fig, 'Type', 'Axes');
    for k = 1:numel(allAx)
        pos = allAx(k).Position;
        allAx(k).Position = [pos(1), pos(2) * subplot_scale, pos(3), pos(4) * subplot_scale];
    end
    saveas(gcf, fullfile(output_folder, sprintf('Metric_Set_%d_%s.png', s, dtype_label)));
    close(gcf);
end
%% ========================================================================
%  SECTION 8: Compile and Export All Significant (p < 0.05) Wilcoxon Rank-Sum Results
% =========================================================================
% Re-iterates through all metric sets and timepoints to collect every
% comparison that reached nominal significance (uncorrected p < 0.05).
% Stores metric name, timepoint, p-value, and group means for export.
% Initialize arrays to store significant results
sig_metric = {};
sig_fraction = {};
sig_pval = [];
sig_mean_LC = [];
sig_mean_LF = [];

% Iterate through the predefined metric sets
for s = 1:length(metric_sets)
    current_metrics = metric_sets{s};
    current_names = set_names{s};
    num_metrics = length(current_metrics);
    
    for m = 1:num_metrics
        metric_data = current_metrics{m};
        cols_to_plot = min(size(metric_data, 2), length(time_labels));
        
        for tp = 1:cols_to_plot
            % Extract data for this fraction and valid patients
            y_raw = metric_data(valid_pts, tp);
            has_data = ~isnan(y_raw);
            y = y_raw(has_data);
            g = lf_group(has_data);
            
            % Check if we have enough data in both groups
            unique_groups = unique(g);
            if length(unique_groups) > 1 && length(y) > 2
                p = ranksum(y(g==0), y(g==1));
                
                % Record the variable if it reaches statistical significance
                if p < 0.05
                    % Calculate group means for context
                    mean_LC = nanmean(y(g == 0));
                    mean_LF = nanmean(y(g == 1));
                    
                    % Append to storage arrays
                    sig_metric{end+1, 1} = current_names{m};
                    sig_fraction{end+1, 1} = time_labels{tp};
                    sig_pval(end+1, 1) = p;
                    sig_mean_LC(end+1, 1) = mean_LC;
                    sig_mean_LF(end+1, 1) = mean_LF;
                end
            end
        end
    end
end

% Construct table and export if any significant results were found
if ~isempty(sig_pval)
    sig_results_table = table(sig_metric, sig_fraction, sig_pval, sig_mean_LC, sig_mean_LF, ...
        'VariableNames', {'Metric', 'Timepoint', 'P_Value', 'Mean_LC', 'Mean_LF'});
    
    % Sort by p-value (most significant first)
    sig_results_table = sortrows(sig_results_table, 'P_Value');
    
    % Display in command window
    disp('----- Significant Findings (p < 0.05) -----');
    disp(sig_results_table);
    

    % Export to CSV
    export_filename = fullfile(dataloc, 'Significant_LF_Metrics.csv');
    writetable(sig_results_table, export_filename);
    fprintf('Saved significant results to: %s\n', export_filename);
else
    disp('No significant differences (p < 0.05) found between LC and LF groups for these metrics.');
end

%% ========================================================================
%  SECTION 9: FDR Correction — Per-Timepoint Benjamini-Hochberg
% =========================================================================
% BH is applied INDEPENDENTLY for each timepoint. Pooling across timepoints
% violates the positive-regression-dependency assumption in BH because
% serial daily measurements of the same patient are positively correlated.
% Within each timepoint, the family consists of one test per metric/metric-set.
if ~isempty(sig_pval)
    fprintf('\n----- PER-TIMEPOINT FDR (Benjamini-Hochberg, Q < 0.05) -----\n');
    
    for tp = 1:length(time_labels)
        tp_pvals  = [];
        tp_labels = {};
        
        for s = 1:length(metric_sets)
            current_metrics = metric_sets{s};
            current_names = set_names{s};
            for mi = 1:length(current_metrics)
                metric_data = current_metrics{mi};
                if tp > size(metric_data, 2), continue; end
                y_raw = metric_data(valid_pts, tp);
                g = lf_group(~isnan(y_raw));
                y = y_raw(~isnan(y_raw));
                if length(unique(g)) > 1 && length(y) > 2
                    tp_pvals(end+1, 1)  = ranksum(y(g==0), y(g==1));
                    tp_labels{end+1, 1} = current_names{mi};
                end
            end
        end
        
        if isempty(tp_pvals), continue; end
        
        % BH within this timepoint only
        n_tp = length(tp_pvals);
        [p_sort, sort_id] = sort(tp_pvals);
        q_tp = zeros(n_tp, 1);
        q_tp(n_tp) = p_sort(n_tp);
        for ii = n_tp-1:-1:1
            q_tp(ii) = min(q_tp(ii+1), p_sort(ii) * (n_tp / ii));
        end
        q_tp = min(q_tp, 1);
        q_unsorted = zeros(n_tp, 1);
        q_unsorted(sort_id) = q_tp;
        
        tp_table = table(tp_labels, tp_pvals, q_unsorted, ...
            'VariableNames', {'Metric', 'Raw_P', 'FDR_Q'});
        sig_tp = tp_table(tp_table.FDR_Q < 0.05, :);
        
        fprintf('\n  Timepoint: %s (family size = %d)\n', time_labels{tp}, n_tp);
        if isempty(sig_tp)
            fprintf('    None survived FDR correction.\n');
        else
            disp(sig_tp);
            writetable(sig_tp, fullfile(dataloc, sprintf('FDR_Sig_%s.csv', strrep(time_labels{tp},' ','_'))));
        end
    end
end

%% ========================================================================
%  SECTION 10: Per-Timepoint Analysis Loop (Fx2 and Fx3)
% =========================================================================
% For each target fraction (Fx2, Fx3), this loop performs:
%   a) Elastic Net–regularized feature selection (α=0.5, binomial logistic, 5-fold CV)
%   b) ROC analysis with Youden's J optimal cutoff
%   c) Multivariable logistic regression (Firth penalized fallback)
%   d) FDR and Holm-Bonferroni multiple-comparison corrections
%   e) Leave-pair-out cross-validation (LPOCV) for unbiased AUC
%   f) 2D scatter plot with logistic decision boundary
%   g) Kaplan-Meier survival analysis
%   h) Correlation matrix of significant features
%   i) Representative parametric ADC maps (responder vs non-responder)
for target_fx = 2:nTp
fx_label = x_labels{target_fx};
fprintf('\n=== Analyzing %s ===\n', fx_label);

%% ---------- Elastic Net Feature Selection at This Timepoint ----------
% Determine which DWI/IVIM metrics are informative for predicting LF.
% A base variable (ADC, D, f, D*) is flagged if Elastic Net selects either
% its absolute or percent-change form.
base_metric_names_all = {'ADC', 'D', 'f', 'D*'};
all_abs_data = {ADC_abs, D_abs, f_abs, Dstar_abs};       % Absolute values
all_pct_data = {ADC_pct, D_pct, f_pct, Dstar_pct};       % Percent change

sig_flags = false(1, 4);   % Will be set true for Elastic Net–selected metrics
sig_p_best = ones(1, 4);   % Best univariate p-value (for downstream sorting)

%% --- Elastic Net Feature Selection ---
    % Construct feature matrix: 18 columns = [8 Imaging | 10 Dosimetry/Sub-volumes]
    % Columns:
    % 1-4: Absolute Imaging (ADC, D, f, D*)
    % 5-8: Percent-change Imaging (ADC, D, f, D*)
    % 9-10: Whole-GTV Dosimetry (D95, V50)
    % 11-12: ADC Sub-volume (D95, V50)
    % 13-14: D Sub-volume (D95, V50)
    % 15-16: f Sub-volume (D95, V50)
    % 17-18: D* Sub-volume (D95, V50)
    X_lasso_all = [ADC_abs(valid_pts, target_fx), D_abs(valid_pts, target_fx), ...
                   f_abs(valid_pts, target_fx),   Dstar_abs(valid_pts, target_fx), ...
                   ADC_pct(valid_pts, target_fx), D_pct(valid_pts, target_fx), ...
                   f_pct(valid_pts, target_fx),   Dstar_pct(valid_pts, target_fx), ...
                   m_d95_gtvp(valid_pts, target_fx), m_v50gy_gtvp(valid_pts, target_fx), ...
                   d95_adc_sub(valid_pts, target_fx), v50_adc_sub(valid_pts, target_fx), ...
                   d95_d_sub(valid_pts, target_fx),   v50_d_sub(valid_pts, target_fx), ...
                   d95_f_sub(valid_pts, target_fx),   v50_f_sub(valid_pts, target_fx), ...
                   d95_dstar_sub(valid_pts, target_fx), v50_dstar_sub(valid_pts, target_fx)];
               
    feat_names_lasso = {'ADC_Abs', 'D_Abs', 'f_Abs', 'Dstar_Abs', ...
                        'ADC_Pct', 'D_Pct', 'f_Pct', 'Dstar_Pct', ...
                        'D95_GTVp', 'V50_GTVp', ...
                        'D95_Sub_ADC', 'V50_Sub_ADC', ...
                        'D95_Sub_D', 'V50_Sub_D', ...
                        'D95_Sub_f', 'V50_Sub_f', ...
                        'D95_Sub_Dstar', 'V50_Sub_Dstar'};

    original_feature_indices = 1:18;

    % --- Feature Exclusion Policy ---
    % Sub-volume dosimetry features (indices 11-18) require Deformable Image
    % Registration (DIR) to correct for tumor deformation between fractions.
    % DIR is now implemented via apply_dir_mask_propagation.m (imregdemons);
    % all 18 features are now eligible for LASSO selection.
    
    % Exclude all dosimetry features for Post-RT scan (no meaningful dose reference)
    if target_fx == 6
        X_lasso_all = X_lasso_all(:, 1:8);
        feat_names_lasso = feat_names_lasso(1:8);
        original_feature_indices = original_feature_indices(1:8);
    end
    
    % Ensure knn_impute_train_test is never fed columns composed entirely of NaNs
    valid_cols = ~all(isnan(X_lasso_all), 1);
    X_lasso_all = X_lasso_all(:, valid_cols);
    feat_names_lasso = feat_names_lasso(valid_cols);
    original_feature_indices = original_feature_indices(valid_cols);

    y_lasso_all = lf_group;

    % --- Imputation strategy to prevent complete-case attrition ---
    % Step 1: Exclude patients missing ALL imaging data (entire row NaN)
    base_cols = min(8, size(X_lasso_all, 2));
    has_any_imaging = any(~isnan(X_lasso_all(:, 1:base_cols)), 2);
    % Step 2: Also require a valid clinical outcome
    impute_mask = has_any_imaging & ~isnan(y_lasso_all);
    X_impute = X_lasso_all(impute_mask, :);
    y_clean  = y_lasso_all(impute_mask);

    % --- Proper 5-fold CV for Lambda Selection to prevent Data Leakage ---
    rng(42); % Set seed for reproducibility
    cvp = cvpartition(y_clean, 'KFold', 5);
    n_lambdas = 25;
    
    % Manual CV to find optimal lambda without data leakage
    common_Lambda = [];
    cv_failed = false;
    
    warning('off', 'all');
    for cv_i = 1:cvp.NumTestSets
        tr_idx = training(cvp, cv_i);
        te_idx = test(cvp, cv_i);
        X_tr = X_impute(tr_idx, :); y_tr = y_clean(tr_idx);
        X_te = X_impute(te_idx, :); y_te = y_clean(te_idx);
        
        % 1. Impute strictly using training fold
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr, X_te, 5);
        
        % 2. Pre-filter features with Pearson |r| > 0.8 using strictly training data
        keep_fold = filter_collinear_features(X_tr_imp, y_tr);
        X_tr_kept = X_tr_imp(:, keep_fold);
        X_te_kept = X_te_imp(:, keep_fold);
        
        try
            if cv_i == 1
                [B_fold, FitInfo_fold] = lassoglm(X_tr_kept, y_tr, 'binomial', ...
                    'Alpha', 0.5, 'NumLambda', n_lambdas, 'Standardize', true, 'MaxIter', 1000000);
                common_Lambda = FitInfo_fold.Lambda;
                all_deviance = zeros(length(common_Lambda), cvp.NumTestSets);
            else
                [B_fold, FitInfo_fold] = lassoglm(X_tr_kept, y_tr, 'binomial', ...
                    'Alpha', 0.5, 'Lambda', common_Lambda, 'Standardize', true, 'MaxIter', 1000000);
            end
            
            % Evaluate binomial deviance on validation fold
            eta = X_te_kept * B_fold + FitInfo_fold.Intercept;
            p = 1 ./ (1 + exp(-eta));
            p = max(min(p, 1 - 1e-10), 1e-10); % Cap probabilities to avoid log(0)
            y_te_mat = repmat(y_te, 1, length(common_Lambda));
            
            % -2 * Log-likelihood for binomial distribution
            dev_val = -2 * sum(y_te_mat .* log(p) + (1 - y_te_mat) .* log(1 - p), 1);
            all_deviance(:, cv_i) = dev_val';
        catch
            cv_failed = true;
            break;
        end
    end
    warning('on', 'all');
    
    if ~cv_failed && ~isempty(common_Lambda)
        % Find optimal lambda (minimum average deviance)
        mean_deviance = mean(all_deviance, 2);
        [~, idx_min] = min(mean_deviance);
        opt_lambda = common_Lambda(idx_min);
        
        % Retrain final model on all data (imputed with full dataset context) using the optimal lambda.
        % We omit the global collinearity filter here to avoid data leakage in the final reporting
        % model; the Elastic Net penalty (Alpha=0.5) will handle multicollinearity during fit.
        X_clean_all = knn_impute_train_test(X_impute, [], 5);
        
        try
            [B_final, FitInfo_final] = lassoglm(X_clean_all, y_clean, 'binomial', ...
                'Alpha', 0.5, 'Lambda', opt_lambda, 'Standardize', true, 'MaxIter', 1000000);
            
            coefs_en = B_final;
            
            % Map selected indices back to original 1..18 mapping
            selected_indices = find(coefs_en ~= 0);
            selected_indices = original_feature_indices(selected_indices);
            
            fprintf('Elastic Net Selected Features for %s (Opt Lambda=%.4f): %s\n', ...
                fx_label, opt_lambda, strjoin(feat_names_lasso(selected_indices), ', '));
        catch
            cv_failed = true;
        end
    end
    
    if cv_failed || isempty(common_Lambda)
        fprintf('Elastic Net failed to converge for %s. Fallback to empty selection.\n', fx_label);
        selected_indices = [];
    end

    %% --- Rigorous Nested LOOCV for Unbiased Risk Scores ---
    % This loop computes out-of-fold risk scores for EVERY patient by re-running
    % the entire pipeline (imputation, corr filter, LASSO selection) on N-1 folds.
    % These unbiased scores are used for both Kaplan-Meier and Cox HRs.
    n_pts_impute = size(X_impute, 1);
    risk_scores_oof = nan(n_pts_impute, 1);
    is_high_risk_oof = false(n_pts_impute, 1);
    
    % Pre-map patient IDs for leakage auditing
    id_list_valid = id_list(valid_pts);
    id_list_impute = id_list_valid(impute_mask);

    fprintf('  Generating unbiased out-of-fold risk scores via nested LOOCV...\n');
    for loo_i = 1:n_pts_impute
        % --- Deep Learning Isolation Check ---
        pat_id_i = id_list_impute{loo_i};
        is_leaky = false;
        if dtype == 2 % dnCNN
            if any(strcmp(dl_provenance.dncnn_train_ids, pat_id_i))
                is_leaky = true;
            end
        elseif dtype == 3 % IVIMnet
            if any(strcmp(dl_provenance.ivimnet_train_ids, pat_id_i))
                is_leaky = true;
            end
        end
        
        if is_leaky
            error('DATA LEAKAGE DETECTED: Patient %s was used to train the %s model. Fundamental isolation broken.', ...
                pat_id_i, dtype_label);
        end

        % Create train fold (N-1)
        train_mask = true(n_pts_impute, 1);
        train_mask(loo_i) = false;
        X_tr_fold = X_impute(train_mask, :);
        y_tr_fold = y_clean(train_mask);
        X_te_fold = X_impute(loo_i, :);
        
        % 1. Impute using KNN fitted strictly on training data
        [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr_fold, X_te_fold, 5);
        
        % 2. Pre-filter: drop features with Pearson |r| > 0.8 using ONLY training data
        % (mirrors the unified collinearity filter)
        keep_fold = filter_collinear_features(X_tr_imp, y_tr_fold);
        X_tr_kept = X_tr_imp(:, keep_fold);
        X_te_kept = X_te_imp(:, keep_fold);
        
        % 3. Nested 5-fold CV for Lambda selection on the train fold
        warning('off', 'all');
        try
            [B_loo, FitInfo_loo] = lassoglm(X_tr_kept, y_tr_fold, 'binomial', ...
                'Alpha', 0.5, 'CV', 5, 'NumLambda', 25, 'Standardize', true, 'MaxIter', 1000000);
            best_idx = FitInfo_loo.IndexMinDeviance;
            coefs_loo = B_loo(:, best_idx);
            intercept_loo = FitInfo_loo.Intercept(best_idx);
        catch
            coefs_loo = zeros(size(X_tr_kept, 2), 1);
            intercept_loo = 0;
        end
        warning('on', 'all');
        
        % 4. Score the held-out patient
        risk_scores_oof(loo_i) = X_te_kept * coefs_loo + intercept_loo;
        
        train_median = median(X_tr_kept * coefs_loo + intercept_loo);
        is_high_risk_oof(loo_i) = risk_scores_oof(loo_i) > train_median;
    end
    
    % Map back to original valid_pts indices (fill NaNs for excluded patients)
    risk_scores_all = nan(sum(valid_pts), 1);
    risk_scores_all(impute_mask) = risk_scores_oof;
    
    % Define high-risk group using the unbiased scores (median split)
    is_high_risk = false(sum(valid_pts), 1);
    is_high_risk(impute_mask) = is_high_risk_oof;

    %% --- End Elastic Net Feature Selection ---
    

% -----------------------------------------------------------------------
% Post–Elastic Net feature list: all 18 candidates treated as fully independent.
% Indices 1-4   = Absolute Imaging forms
% Indices 5-8   = Percent-Change Imaging forms
% Indices 9-10  = Whole-GTV Dosimetry (D95, V50)
% Indices 11-18 = Sub-volume Dosimetry (D95, V50 pairs)
% -----------------------------------------------------------------------
all_feat_data  = {ADC_abs,       D_abs,       f_abs,       Dstar_abs, ...
                  ADC_pct,       D_pct,       f_pct,       Dstar_pct, ...
                  m_d95_gtvp,    m_v50gy_gtvp, ...
                  d95_adc_sub,   v50_adc_sub, ...
                  d95_d_sub,     v50_d_sub, ...
                  d95_f_sub,     v50_f_sub, ...
                  d95_dstar_sub, v50_dstar_sub};

all_feat_names = {'ADC',         'D',         'f',         'D*', ...
                  'ADC',         'D',         'f',         'D*', ...
                  'D95 GTVp',    'V50 GTVp', ...
                  'D95 Sub(ADC)','V50 Sub(ADC)', ...
                  'D95 Sub(D)',  'V50 Sub(D)', ...
                  'D95 Sub(f)',  'V50 Sub(f)', ...
                  'D95 Sub(D*)', 'V50 Sub(D*)'};

all_feat_is_abs = [true          true         true         true  ...
                   false         false        false        false ...
                   true          false        true         false ...
                   true          false        true         false ...
                   true          false];

all_feat_disp  = {'Abs ADC',     'Abs D',     'Abs f',     'Abs D*', ...
                  '\Delta ADC',  '\Delta D',  '\Delta f',  '\Delta D*', ...
                  'D95 GTVp',    'V50 GTVp', ...
                  'D95 Sub(ADC)','V50 Sub(ADC)', ...
                  'D95 Sub(D)',  'V50 Sub(D)', ...
                  'D95 Sub(f)',  'V50 Sub(f)', ...
                  'D95 Sub(D*)', 'V50 Sub(D*)'};

all_feat_units = {'mm^2/s',      'mm^2/s',    'Fraction',  'mm^2/s', ...
                  '%',           '%',         '%',         '%', ...
                  'Gy',          '%', ...
                  'Gy',          '%', ...
                  'Gy',          '%', ...
                  'Gy',          '%', ...
                  'Gy',          '%'};

% n_sig = number of features LASSO selected (up to 18)
n_sig = length(selected_indices);

% Build parallel arrays consumed by all downstream loops
sig_data_selected = cell(1, n_sig);
sig_abs_data      = cell(1, n_sig);
sig_pct_data      = cell(1, n_sig);
sig_names         = cell(1, n_sig);
sig_is_abs        = false(1, n_sig);
sig_is_pct_imaging = false(1, n_sig);
sig_disp_names    = cell(1, n_sig);
sig_units         = cell(1, n_sig);

for si = 1:n_sig
    fi = selected_indices(si);              % direct column index (1-18)
    sig_data_selected{si} = all_feat_data{fi};
    sig_names{si}         = all_feat_names{fi};
    sig_is_abs(si)        = all_feat_is_abs(fi);
    sig_is_pct_imaging(si) = (fi >= 5 && fi <= 8);
    sig_disp_names{si}    = all_feat_disp{fi};
    sig_units{si}         = all_feat_units{fi};
    
    % Carry both Abs and Pct forms for dose-response / Table-1 use (Imaging only)
    if fi <= 4
        sig_abs_data{si} = all_feat_data{fi};
        sig_pct_data{si} = all_feat_data{fi + 4};
    elseif fi >= 5 && fi <= 8
        sig_abs_data{si} = all_feat_data{fi - 4};
        sig_pct_data{si} = all_feat_data{fi};
    else
        % For dosimetry, we don't necessarily have abs/pct pairs in the same way
        sig_abs_data{si} = all_feat_data{fi};
        sig_pct_data{si} = all_feat_data{fi};
    end
end

fprintf('Significant variables at %s: ', fx_label);
if n_sig == 0
    fprintf('NONE. Skipping downstream analyses for %s.\n', fx_label);
    continue;
else
    fprintf('%s\n', strjoin(sig_disp_names, ', '));
end

% Build time-to-event vector for all valid patients
times_km = m_total_time;
times_km(m_lf==0) = m_total_follow_up_time(m_lf==0);
events_km = m_lf;

% Subset to valid patients for this specific fraction
times_km = times_km(valid_pts);
events_km = events_km(valid_pts);

%% ---------- PRIMARY ROC Analysis (LOOCV Out-of-Fold Risk Scores) ----------
% The full ROC curve and optimal decision threshold are computed exclusively
% from risk_scores_all — the single, non-averaged out-of-fold prediction
% produced by the nested LOOCV loop above. Each patient is scored exactly
% once by a model that never saw their data.
%
% LPOCV (below) is RESERVED for the paired concordance AUC scalar only.
% It does NOT determine the ROC curve or the optimum cutoff.
%
% Optimal cutoff is determined by Youden's J: max(Sensitivity - FPR).

labels = lf_group; % 0 = LC, 1 = LF
valid_roc = ~isnan(risk_scores_all) & ~isnan(labels);

if sum(valid_roc) > 0
    % Use the LOOCV out-of-fold risk scores directly — no GLM refitting.
    [roc_X, roc_Y, roc_T, roc_AUC] = perfcurve(labels(valid_roc), risk_scores_all(valid_roc), 1);
    
    % Optimal cutoff via Youden's J Statistic (maximise Sensitivity + Specificity - 1)
    [~, roc_opt_idx] = max(roc_Y - roc_X);
    roc_opt_thresh = roc_T(roc_opt_idx);
    
    % --- Plotting the ROC Curve ---
    figure('Name', ['ROC Analysis - ' fx_label ' — ' dtype_label], 'Position', [200, 200, 700, 600]);
    hold on;
    
    plot(roc_X, roc_Y, 'b-', 'LineWidth', 2.5);
    leg_entries = {sprintf('LOOCV OOF Risk Score (AUC = %.3f)', roc_AUC)};
    
    plot(roc_X(roc_opt_idx), roc_Y(roc_opt_idx), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    leg_entries{end+1} = sprintf('Youden Cutoff (score = %.3f)', roc_opt_thresh);
    
    plot([0 1], [0 1], 'k--', 'LineWidth', 1.5);
    leg_entries{end+1} = 'Random Guess';
    
    xlabel('False Positive Rate (1 - Specificity)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('True Positive Rate (Sensitivity)', 'FontSize', 12, 'FontWeight', 'bold');
    title(['PRIMARY ROC Curve: LOOCV Out-of-Fold Risk Score (' fx_label ', ' dtype_label ')'], 'FontSize', 14);
    
    legend(leg_entries, 'Location', 'SouthEast', 'FontSize', 11);
    grid on; box on;
    hold off;
    saveas(gcf, fullfile(output_folder, ['ROC_OOF_Risk_Score_' fx_label '_' dtype_label '.png']));
    close(gcf);
    
    % Print results to the command window
    fprintf('\n--- PRIMARY ROC ANALYSIS (LOOCV Out-of-Fold Risk Score) for %s ---\n', fx_label);
    fprintf('  AUC  = %.3f\n  Youden Optimal Score Cutoff = %.4f\n', roc_AUC, roc_opt_thresh);
    fprintf('  Sensitivity = %.1f%%  |  Specificity = %.1f%%\n\n', ...
        roc_Y(roc_opt_idx)*100, (1-roc_X(roc_opt_idx))*100);
else
    roc_AUC = NaN; roc_opt_thresh = NaN;
    fprintf('\n--- PRIMARY ROC ANALYSIS (LOOCV OOF) for %s ---\n', fx_label);
    fprintf('Insufficient data for out-of-fold ROC analysis.\n\n');
end



%% ---------- Harrell's C-index Evaluation (Out-of-Fold Risk Scores) ----------
% Evaluates the continuous out-of-fold risk scores directly. Harrell's C-index 
% (concordance index) explicitly drops mathematically "incomparable pairs" 
% where a patient's censoring time precedes the paired event time.

fprintf('\n--- UNBIASED EVALUATION: Harrell''s C-index (LOOCV OOF Scores) ---\n');

% 1. Extract valid pairs with out-of-fold scores
valid_idx_c = find(~isnan(risk_scores_all) & ~isnan(times_km) & ~isnan(events_km));
n_eval = length(valid_idx_c);

if n_eval < 2
    fprintf('  Insufficient data for C-index calculation.\n');
else
    concordant = 0;
    ties = 0;
    total_comparable = 0;
    
    for i = 1:n_eval
        ii = valid_idx_c(i);
        ti = times_km(ii);
        ei = events_km(ii);
        si = risk_scores_all(ii);
        
        for j = i+1:n_eval
            jj = valid_idx_c(j);
            tj = times_km(jj);
            ej = events_km(jj);
            sj = risk_scores_all(jj);
            
            % Harrell's Comparability Rules:
            % 1. If both failed, comparable.
            % 2. If one failed at T_i and the other is censored at T_j > T_i, comparable.
            % Otherwise, incomparable.
            
            comparable = false;
            better_survival_idx = 0;
            worse_survival_idx = 0;
            
            if ti < tj
                if ei == 1
                    comparable = true;
                    worse_survival_idx = ii;
                    better_survival_idx = jj;
                end
            elseif tj < ti
                if ej == 1
                    comparable = true;
                    worse_survival_idx = jj;
                    better_survival_idx = ii;
                end
            else % ti == tj
                if ei == 1 && ej == 1
                    % Both failed at same time: technically a tie in time.
                    % Typically dropped or handled as concordant if scores match.
                    % Following "explicilty drop incomparable pairs" for ties in time.
                end
            end
            
            if comparable
                total_comparable = total_comparable + 1;
                score_worse = risk_scores_all(worse_survival_idx);
                score_better = risk_scores_all(better_survival_idx);
                
                % C-index concordance: Higher risk score should have shorter survival (worse)
                if score_worse > score_better
                    concordant = concordant + 1;
                elseif score_worse == score_better
                    ties = ties + 1;
                end
            end
        end
    end
    
    if total_comparable > 0
        c_index = (concordant + 0.5 * ties) / total_comparable;
        fprintf('  Evaluated Comparable Pairs: %d\n', total_comparable);
        fprintf('  Harrell''s C-index: %.3f\n', c_index);
    else
        fprintf('  No comparable pairs found for C-index calculation.\n');
    end
end

fprintf('------------------------------------------------\n');

%% ---------- 2D Scatter Plot with Logistic Decision Boundary ----------
% Visualizes the two-feature space for all pairs of significant variables
% (e.g. %-change vars) with LC (blue) and LF (red) patient points.
% Overlays the linear decision boundary from the logistic
% regression model: intercept + β1*x + β2*y = 0.
if n_sig >= 2
    for fi = 1:(n_sig-1)
        for fj = (fi+1):n_sig
            figure('Name', sprintf('2D Feature Space %s vs %s %s — %s', sig_names{fi}, sig_names{fj}, fx_label, dtype_label), 'Position', [100, 100, 800, 600]);
            hold on;
            
            x_val = sig_data_selected{fi}(valid_pts, target_fx);
            y_val = sig_data_selected{fj}(valid_pts, target_fx);
            group = lf_group;

            scatter(x_val(group==0), y_val(group==0), 80, 'b', 'filled', 'MarkerEdgeColor', 'k');
            scatter(x_val(group==1), y_val(group==1), 80, 'r', 'filled', 'MarkerEdgeColor', 'k');

            % Draw Decision Boundary (from logistic regression model on these 2 vars)
            mdl = fitglm([x_val, y_val], group, 'Distribution', 'binomial', 'Options', statset('MaxIter', 1000000));
            coefs = mdl.Coefficients.Estimate;
            x_range = linspace(min(x_val), max(x_val), 100);
            y_boundary = -(coefs(1) + coefs(2)*x_range) / coefs(3);
            plot(x_range, y_boundary, 'k--', 'LineWidth', 2);

            if sig_is_abs(fi), xl = sprintf('%s at %s (%s)', sig_names{fi}, fx_label, sig_units{fi}); else, xl = sprintf('\\Delta %s at %s (%s)', sig_names{fi}, fx_label, sig_units{fi}); end
            if sig_is_abs(fj), yl = sprintf('%s at %s (%s)', sig_names{fj}, fx_label, sig_units{fj}); else, yl = sprintf('\\Delta %s at %s (%s)', sig_names{fj}, fx_label, sig_units{fj}); end
            xlabel(xl, 'FontSize', 12, 'FontWeight', 'bold');
            ylabel(yl, 'FontSize', 12, 'FontWeight', 'bold');
            title(sprintf('Biomarker Interaction: Separation of LC vs LF (%s, %s)', fx_label, dtype_label), 'FontSize', 14);
            legend({'Local Control', 'Local Failure', 'Logistic Decision Boundary'}, 'Location', 'NorthWest');
            
            grid on; box on;
            xline(0, 'k-', 'Alpha', 0.2); yline(0, 'k-', 'Alpha', 0.2);
            
            safe_name1 = strrep(sig_names{fi}, '*', 'star');
            safe_name2 = strrep(sig_names{fj}, '*', 'star');
            saveas(gcf, fullfile(output_folder, sprintf('2D_Space_%s_vs_%s_%s_%s.png', safe_name1, safe_name2, fx_label, dtype_label)));
            close(gcf);
        end
    end
else
    fprintf('Skipping 2D scatter plots: requires at least 2 significant variables.\n');
end

%% ---------- Kaplan-Meier Survival Analysis ----------
% Uses nested Leave-One-Out Cross-Validation (LOOCV) to prevent data leakage. 
% For each held-out patient, a multivariable linear risk score is computed 
% using the non-zero Elastic Net coefficients on the N-1 fold.
% The median risk score of the training fold determines the classification threshold.

% Build the full candidate feature matrix for all valid patients (18 columns)
km_X_all = [ADC_abs(valid_pts, target_fx), D_abs(valid_pts, target_fx), ...
            f_abs(valid_pts, target_fx),   Dstar_abs(valid_pts, target_fx), ...
            ADC_pct(valid_pts, target_fx), D_pct(valid_pts, target_fx), ...
            f_pct(valid_pts, target_fx),   Dstar_pct(valid_pts, target_fx), ...
            m_d95_gtvp(valid_pts, target_fx), m_v50gy_gtvp(valid_pts, target_fx), ...
            d95_adc_sub(valid_pts, target_fx), v50_adc_sub(valid_pts, target_fx), ...
            d95_d_sub(valid_pts, target_fx),   v50_d_sub(valid_pts, target_fx), ...
            d95_f_sub(valid_pts, target_fx),   v50_f_sub(valid_pts, target_fx), ...
            d95_dstar_sub(valid_pts, target_fx), v50_dstar_sub(valid_pts, target_fx)];
km_y_all = lf_group;
n_km = length(km_y_all);

% (LOOCV logic moved to Elastic Net section to provide unified unbiased risk scores)

% (Survival data extraction moved earlier to support C-index calculation)

% Compute and plot Kaplan-Meier survival curves for each risk group
figure('Name', ['Kaplan-Meier ' fx_label ' — ' dtype_label], 'Position', [100, 100, 700, 600]);

[f, x, flow, fup] = ecdf(times_km, 'Censoring', ~events_km, 'Function', 'survivor', 'Alpha', 0.05);
if sum(is_high_risk) > 0 && sum(~is_high_risk) > 0
    [f1, x1] = ecdf(times_km(is_high_risk), 'Censoring', ~events_km(is_high_risk), 'Function', 'survivor');
    [f2, x2] = ecdf(times_km(~is_high_risk), 'Censoring', ~events_km(~is_high_risk), 'Function', 'survivor');
    stairs(x1, f1, 'r-', 'LineWidth', 2.5); hold on;
    stairs(x2, f2, 'b-', 'LineWidth', 2.5);
    legend({'High Risk', 'Low Risk'}, 'Location', 'SouthWest');
    
    % Perform log-rank test approximation via univariate Cox regression
    [~, ~, ~, stats] = coxphfit(is_high_risk, times_km, 'Censoring', ~events_km);
    p_val_km = stats.p;
    
    title_str = sprintf('Kaplan-Meier (LOOCV): Stratified by Multivariable Risk Score (%s) (%s)\nLog-Rank p = %.4f', fx_label, dtype_label, p_val_km);
else
    fprintf('Warning: Risk groups are degenerate. Skipping KM stratification.\n');
    title_str = sprintf('Kaplan-Meier (LOOCV): Stratified by Multivariable Risk Score (%s) (%s)', fx_label, dtype_label);
end
xlabel('Time to Local Failure (Days)', 'FontSize', 12);
ylabel('Local Control Probability', 'FontSize', 12);
title(title_str, 'FontSize', 14);
grid on; axis([0 max(times)+50 0 1.05]);
saveas(gcf, fullfile(output_folder, ['Kaplan_Meier_' fx_label '_' dtype_label '.png']));
close(gcf);

fprintf('Kaplan-Meier (LOOCV) generated using multivariable risk score. Unbiased out-of-fold risk assignment.\n');

%% ---------- Correlation Matrix of Significant Features ----------
% Computes Pearson correlation coefficients between all pairs of
% Elastic Net–selected percent-change variables. Displayed as a jet-coloured
% heatmap with R values in [-1, 1]. High inter-feature correlation
% suggests redundancy and potential multicollinearity in the combined model.
% Build feature matrix and names from significant percent-change variables
figure('Name', ['Correlation Matrix ' fx_label ' — ' dtype_label], 'Position', [100, 100, 600, 500]);
feats = zeros(sum(valid_pts), n_sig);
feat_names = cell(1, n_sig);
for vi = 1:n_sig
    feats(:, vi) = sig_data_selected{vi}(valid_pts, target_fx);
    if sig_is_abs(vi), feat_names{vi} = ['Abs ' sig_names{vi}]; else, feat_names{vi} = ['\Delta ' sig_names{vi}]; end
end
[R, P] = corrcoef(feats, 'Rows', 'complete');
heatmap(feat_names, feat_names, R, 'ColorLimits', [-1 1], 'Colormap', jet);
if n_sig >= 2
    title(['Feature Correlation (' fx_label ') (' dtype_label '). R(' sig_disp_names{1} ', ' sig_disp_names{2} ') = ' num2str(R(1,2), '%.2f')]);
else
    title(['Feature Correlation (' fx_label ') (' dtype_label ')']);
end
saveas(gcf, fullfile(output_folder, ['Correlation_Matrix_' fx_label '_' dtype_label '.png']));
close(gcf);

%% ---------- Representative Parametric ADC Maps ----------
for vi = 1:n_sig
    if selected_indices(vi) > 8
        fprintf('  Skipping Representative Map for Dosimetry/Sub-volume metric: %s\n', sig_names{vi});
        continue;
    end

    base_idx = mod(selected_indices(vi)-1, 4) + 1;
    switch base_idx
        case 1
            map_name = 'ADC'; clims = [0 2.5e-3]; diff_clims = [-1e-3 1e-3];
        case 2
            map_name = 'D'; clims = [0 2.5e-3]; diff_clims = [-1e-3 1e-3];
        case 3
            map_name = 'f'; clims = [0 0.5]; diff_clims = [-0.2 0.2];
        case 4
            map_name = 'D*'; clims = [0 0.05]; diff_clims = [-0.02 0.02];
    end

    curr_sig_pct_full = sig_data_selected{vi};
    curr_sig_name = sig_names{vi};
    if sig_is_abs(vi)
        curr_sig_disp = ['Abs ' curr_sig_name];
        curr_sig_file = ['Abs_' curr_sig_name];
    else
        curr_sig_disp = ['\Delta ' curr_sig_name];
        curr_sig_file = ['Delta_' curr_sig_name];
    end

fprintf('Generating Representative Maps for %s...\n', curr_sig_disp);

valid_for_plot = isfinite(m_lf) & ~isnan(curr_sig_pct_full(:,target_fx));

% Identify Responders (Local Control)
lc_candidates = find(m_lf == 0 & valid_for_plot);
if isempty(lc_candidates), continue; end
[~, best_rel_idx] = max(curr_sig_pct_full(lc_candidates, target_fx));
idx_responder = lc_candidates(best_rel_idx);

% Identify Non-Responders (Local Failure)
lf_candidates = find(m_lf == 1 & valid_for_plot);
if isempty(lf_candidates), continue; end
[~, worst_rel_idx] = min(curr_sig_pct_full(lf_candidates, target_fx));
idx_nonresponder = lf_candidates(worst_rel_idx);

patients_to_plot = [idx_responder, idx_nonresponder];
titles = {'Responder (Local Control)', 'Non-Responder (Local Failure)'};

figure('Name', ['Representative ADC Evolution ' curr_sig_disp ' ' fx_label ' — ' dtype_label], 'Position', [50, 50, 1200, 800]);
colormap(jet);

for p = 1:2
    pt_idx = patients_to_plot(p);
    fprintf('Processing Patient %d (%s)...\n', pt_idx, titles{p});
    
    slice_fx1 = []; slice_fx3 = [];
    mask_slice_fx1 = []; mask_slice_fx3 = [];
    
    for t = 1:2
        tp_num = [1, target_fx];  % Compare Fx1 vs target fraction
        tp = tp_num(t);
        
        % Construct NIfTI file paths for this patient/timepoint
        basefolder = [dataloc m_id_list{pt_idx} '/'];
        nii_path = [basefolder '/nii/'];
        
        % Set scan and GTV mask filenames based on timepoint
        if tp == 1
            scanID = 'fx1_dwi1'; gtvname = 'fx1_gtv1';
        else
            scanID = sprintf('fx%d_dwi1', target_fx); gtvname = sprintf('fx%d_gtv1', target_fx);
        end
        
        dwi_file = [nii_path scanID '.nii.gz'];   % 4D DWI volume
        bval_file = [nii_path scanID '.bval'];     % b-value text file
        gtv_file = [nii_path gtvname '.nii.gz'];   % Binary GTV mask

        if ~exist(dwi_file, 'file') || ~exist(bval_file, 'file')
            fprintf('  Skipping: Missing files for Pt %d Tp %d\n', pt_idx, tp);
            continue;
        end

        % Load NIfTI data; rot90 for radiological display convention
        dwi_nii = load_untouch_nii(dwi_file);
        dwi_img = double(rot90(dwi_nii.img)); 
        gtv_nii = load_untouch_nii(gtv_file);
        gtv_img = double(rot90(gtv_nii.img));
        
        % Read b-values from text file (space-delimited single line)
        fid = fopen(bval_file);
        bvals = str2num(fgetl(fid)); 
        fclose(fid);
        bvals = bvals(:); % Column vector
        
        % Validate b-values against expected study protocol
        expected_bvals = [0; 30; 150; 550];
        if size(dwi_img, 4) ~= length(bvals) || ~isequal(sort(bvals), expected_bvals)
            fprintf('  Protocol deviation: Pt %d Tp %d has non-standard b-values %s — excluding from analysis.\n', ...
                pt_idx, tp, mat2str(bvals'));
            continue;
        end
        
        % Select Central Slice (Max Tumor Area)
        gtv_areas = squeeze(sum(sum(gtv_img,1),2));
        [~, z_slice] = max(gtv_areas);
        
        slice_gtv = squeeze(gtv_img(:,:,z_slice));
        
        % Compute required parametric map and extract slice
        switch base_idx
            case 1 % ADC
                map_vol = fit_adc_mono(dwi_img, bvals);
                slice_map = squeeze(map_vol(:,:,z_slice));
                slice_map(slice_map < 0) = 0; 
                slice_map(slice_map > 3e-3) = 3e-3;
            case 2 % D
                opts = struct(); opts.bthr = 100;
                slice_dwi = dwi_img(:,:,z_slice,:);
                mask_ivim = true(size(slice_dwi,1), size(slice_dwi,2), 1);
                ivim_fit = IVIMmodelfit(slice_dwi, bvals, "seg", mask_ivim, opts);
                slice_map = squeeze(ivim_fit(:,:,1,1));
                slice_map(slice_map < 0) = 0; 
                slice_map(slice_map > 3e-3) = 3e-3;
            case 3 % f
                opts = struct(); opts.bthr = 100;
                slice_dwi = dwi_img(:,:,z_slice,:);
                mask_ivim = true(size(slice_dwi,1), size(slice_dwi,2), 1);
                ivim_fit = IVIMmodelfit(slice_dwi, bvals, "seg", mask_ivim, opts);
                slice_map = squeeze(ivim_fit(:,:,1,3));
                slice_map(slice_map < 0) = 0;
                slice_map(slice_map > 1) = 1;
            case 4 % D*
                opts = struct(); opts.bthr = 100;
                slice_dwi = dwi_img(:,:,z_slice,:);
                mask_ivim = true(size(slice_dwi,1), size(slice_dwi,2), 1);
                ivim_fit = IVIMmodelfit(slice_dwi, bvals, "seg", mask_ivim, opts);
                slice_map = squeeze(ivim_fit(:,:,1,4));
                slice_map(slice_map < 0) = 0;
        end
        
        % Store Fx1 and target-Fx slices separately for later comparison
        if t == 1
            slice_fx1 = slice_map; mask_slice_fx1 = slice_gtv;
        else
            slice_fx3 = slice_map; mask_slice_fx3 = slice_gtv;
        end
    end
    
    % --- IMAGE PROCESSING & PLOTTING ---
    if ~isempty(slice_fx1) && ~isempty(slice_fx3)
        % Resize target-Fx slice to match Fx1 dimensions (if needed)
        sz1 = size(slice_fx1);
        if ~isequal(size(slice_fx3), sz1)
            slice_fx3 = imresize(slice_fx3, sz1, 'bilinear');
            mask_slice_fx3 = imresize(mask_slice_fx3, sz1, 'nearest');
        end
        
        % Crop to the union of Fx1 and target-Fx GTV masks + margin
        combined_mask = (mask_slice_fx1 > 0.5) | (mask_slice_fx3 > 0.5);
        [rows, cols] = find(combined_mask);
        
        if isempty(rows), rows=1:sz1(1); cols=1:sz1(2); end % Fallback if no mask
        
        margin = 15;  % Pixel padding around GTV bounding box
        r_min = max(1, min(rows)-margin); r_max = min(sz1(1), max(rows)+margin);
        c_min = max(1, min(cols)-margin); c_max = min(sz1(2), max(cols)+margin);
        
        crop_fx1 = slice_fx1(r_min:r_max, c_min:c_max);
        crop_fx3 = slice_fx3(r_min:r_max, c_min:c_max);
        
        % Pixel-wise difference map (positive = ADC increase = response)
        crop_diff = crop_fx3 - crop_fx1;
        
        % 2×3 subplot layout: top row = responder, bottom row = non-responder
        row_offset = (p-1)*3;
        
        subplot(2, 3, 1 + row_offset);
        imagesc(crop_fx1, clims); axis image; axis off;
        title([titles{p} ' - Fx1'], 'FontSize', 10, 'Interpreter', 'none');
        if p==1, ylabel([map_name ' Map']); end
        
        subplot(2, 3, 2 + row_offset);
        imagesc(crop_fx3, clims); axis image; axis off;
        title(fx_label, 'FontSize', 10);
        
        subplot(2, 3, 3 + row_offset);
        imagesc(crop_diff, diff_clims); axis image; axis off;
        title(['\Delta (' fx_label ' - Fx1)'], 'FontSize', 10);
        
        % Overlay Fx1 GTV contour on difference map (white outline)
        hold on; contour(mask_slice_fx1(r_min:r_max, c_min:c_max), [0.5 0.5], 'w', 'LineWidth', 1); hold off;
    end
end

% Add shared colorbars
c1 = colorbar('Position', [0.92 0.55 0.02 0.35]); ylabel(c1, map_name);
c2 = colorbar('Position', [0.92 0.11 0.02 0.35]); ylabel(c2, ['\Delta ' map_name]);
sgtitle(['Representative Longitudinal ' map_name ' Response for ' curr_sig_disp ' (' fx_label ', ' dtype_label ')'], 'FontSize', 14, 'FontWeight', 'bold');
safe_name = strrep(curr_sig_file, '*', 'star');
saveas(gcf, fullfile(output_folder, ['Representative_' strrep(map_name, '*', 'star') '_' safe_name '_' fx_label '_' dtype_label '.png']));
close(gcf);
end % loop over sig variables

%% 5. "Sanity Checks": Volume, Noise, and Heterogeneity
for vi = 1:n_sig
curr_sig_pct_full = sig_data_selected{vi};
curr_sig_name = sig_names{vi};
if sig_is_abs(vi)
    curr_sig_disp = ['Abs ' curr_sig_name];
    curr_sig_file = ['Abs_' curr_sig_name];
else
    curr_sig_disp = ['\Delta ' curr_sig_name];
    curr_sig_file = ['Delta_' curr_sig_name];
end

figure('Name', ['Sanity Checks ' curr_sig_disp ' ' fx_label ' — ' dtype_label], 'Position', [100, 100, 1200, 500]);
subplot(1, 3, 1);
vol_fx1 = m_gtv_vol(valid_pts, 1);          
vol_fx3 = m_gtv_vol(valid_pts, target_fx);  
vol_pct = (vol_fx3 - vol_fx1) ./ vol_fx1 * 100;  

boxplot(vol_pct, lf_group, 'Labels', {'LC (0)', 'LF (1)'});
ylabel(['% Change in GTV Volume (' fx_label ')']);
title('Confounder Check: Volume', 'FontSize', 12, 'FontWeight', 'bold');
p_vol = ranksum(vol_pct(lf_group==0), vol_pct(lf_group==1));

y_lim = ylim;
% Ensure axis has finite and distinct limits before placing text annotation
if numel(y_lim) >= 2 && all(isfinite(y_lim)) && y_lim(2) > y_lim(1)
    text(1.5, y_lim(1) + 0.9*(y_lim(2)-y_lim(1)), sprintf('p = %.3f', p_vol), ...
        'HorizontalAlignment', 'center', 'FontSize', 11);
end
grid on;
if p_vol > 0.05
    xlabel('Conclusion: No Volumetric Bias');
else
    xlabel('Warning: Volume is a Confounder');
end

% NOTE: Histogram kurtosis of trace-average ADC is NOT valid DKI.
% True DKI requires fitting the full kurtosis tensor W to raw directional
% diffusion gradients. ADC SD (intra-tumour voxel spread) is used instead
% as a model-free, valid heterogeneity proxy.
subplot(1, 3, 2);
sd_fx1  = adc_sd(valid_pts, 1, 1);
sd_fxN  = adc_sd(valid_pts, target_fx, 1);
sd_delta = sd_fxN - sd_fx1;

boxplot(sd_delta, lf_group, 'Labels', {'LC (0)', 'LF (1)'});
ylabel(['\Delta ADC SD (' fx_label ') [mm\u00b2/s]']);
title('Heterogeneity: ADC SD Change', 'FontSize', 12, 'FontWeight', 'bold');
p_sd = ranksum(sd_delta(lf_group==0), sd_delta(lf_group==1));

y_lim = ylim;
if numel(y_lim) >= 2 && all(isfinite(y_lim)) && y_lim(2) > y_lim(1)
    text(1.5, y_lim(1) + 0.9*(y_lim(2)-y_lim(1)), sprintf('p = %.3f', p_sd), ...
        'HorizontalAlignment', 'center', 'FontSize', 11);
end
grid on;

subplot(1, 3, 3);
hold on;
wcv_est = 2.8; 
cor_est = 1.96 * sqrt(2) * wcv_est;

x_scatter = ones(size(lf_group));
x_scatter(lf_group==1) = 2;
x_scatter = x_scatter + (rand(size(x_scatter))-0.5)*0.2;
scatter(x_scatter, curr_sig_pct_full(valid_pts, target_fx), 50, 'filled', 'MarkerEdgeColor', 'k');

% Check the base variable index (1=ADC, 2=D, 3=f, 4=D*, 9+=Dosimetry)
base_idx = mod(selected_indices(vi)-1, 4) + 1;

if sig_is_pct_imaging(vi) && base_idx == 1
    % Apply 7.8% noise floor bounds strictly for ADC percent change
    yfill = [-cor_est cor_est cor_est -cor_est];
    xfill = [0.5 0.5 2.5 2.5];
    fill(xfill, yfill, [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    yline(0, 'k-');
    yline(cor_est, 'k--', 'CoR (+7.8%)');
    yline(-cor_est, 'k--', 'CoR (-7.8%)');
elseif ~sig_is_abs(vi) && contains(sig_units{vi}, '%')
    % For other percentage features (D*, dose volume, etc.), just add a 0 line
    yline(0, 'k-', 'Alpha', 0.3);
end

xticks([1 2]); xticklabels({'LC', 'LF'});
ylbl = sprintf('%s (%s)', sig_disp_names{vi}, sig_units{vi});
ylabel(ylbl);
title('Signal vs. Noise Floor', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0.5 2.5]);

sgtitle(['Validation (' curr_sig_disp '): Volume, Texture, and Noise (' fx_label ', ' dtype_label ')'], 'FontSize', 14, 'FontWeight', 'bold');
allAx = findall(gcf, 'Type', 'Axes');
for k = 1:numel(allAx)
    pos = allAx(k).Position;
    allAx(k).Position = [pos(1), pos(2) * 0.92, pos(3), pos(4) * 0.92];
end
safe_name = strrep(curr_sig_file, '*', 'star');
saveas(gcf, fullfile(output_folder, ['Sanity_Checks_' safe_name '_' fx_label '_' dtype_label '.png']));
close(gcf);
end % loop over sig variables

%% 6. Final Data for Manuscript (Table 1, HR, and Sub-volumes)
fprintf('\n--- GENERATING TABLE 1 & MANUSCRIPT DATA ---\n');

% --- PART A: COHORT CHARACTERISTICS (Table 1) ---
% Check if LC and LF groups differ at baseline
% 1. Baseline Volume (cc)
vol_baseline = m_gtv_vol(valid_pts, 1); 
[p_vol_base, ~, stats_vol] = ranksum(vol_baseline(lf_group==0), vol_baseline(lf_group==1));
fprintf('Baseline Tumor Volume (cc):\n');
fprintf('  LC: %.1f ± %.1f\n', mean(vol_baseline(lf_group==0)), std(vol_baseline(lf_group==0)));
fprintf('  LF: %.1f ± %.1f\n', mean(vol_baseline(lf_group==1)), std(vol_baseline(lf_group==1)));
fprintf('  p-value (Wilcoxon): %.3f\n', p_vol_base);

% 2. Baseline values for each significant variable
for vi = 1:n_sig
    sig_baseline = sig_abs_data{vi}(valid_pts, 1);
    [p_base_vi, ~, ~] = ranksum(sig_baseline(lf_group==0), sig_baseline(lf_group==1));
    fprintf('\nBaseline %s:\n', sig_names{vi});
    fprintf('  LC: %.4f ± %.4f\n', mean(sig_baseline(lf_group==0), 'omitnan'), std(sig_baseline(lf_group==0), 'omitnan'));
    fprintf('  LF: %.4f ± %.4f\n', mean(sig_baseline(lf_group==1), 'omitnan'), std(sig_baseline(lf_group==1), 'omitnan'));
    fprintf('  p-value: %.3f\n', p_base_vi);
end

% --- PART B: EXPLORATORY HAZARD RATIOS (In-Sample Features) ---
fprintf('\n--- EXPLORATORY SURVIVAL ANALYSIS (In-Sample / Double-Dipped) ---\n');
for vi = 1:n_sig
    curr_sig_pct_full = sig_data_selected{vi}(valid_pts, target_fx);
    curr_sig_name = sig_names{vi};
    if sig_is_abs(vi), var_desc = 'Absolute'; else, var_desc = 'Continuous Delta'; end

    valid_cox = ~isnan(curr_sig_pct_full) & ~isnan(times) & ~isnan(events);
    cox_times = times(valid_cox);
    cox_events = events(valid_cox);
    cox_biomarker = curr_sig_pct_full(valid_cox);

    if sum(valid_cox) > 5
        warning('error', 'stats:coxphfit:FitWarning');
        warning('error', 'stats:coxphfit:IterationLimit');
        try
            [b, logl, H, stats] = coxphfit(cox_biomarker, cox_times, 'Censoring', ~cox_events, 'Options', statset('MaxIter', 1000000));
            is_firth = false;
        catch ME
            if contains(ME.identifier, 'FitWarning') || contains(ME.identifier, 'IterationLimit') || contains(lower(ME.message), 'perfect')
                mdl_firth = fitglm(cox_biomarker, cox_events, 'Distribution', 'binomial', 'LikelihoodPenalty', 'jeffreys-prior', 'Options', statset('MaxIter', 1000000));
                b = mdl_firth.Coefficients.Estimate(2:end);
                stats.se = mdl_firth.Coefficients.SE(2:end);
                stats.p = mdl_firth.Coefficients.pValue(2:end);
                stats.sschres = zeros(size(cox_biomarker,1), size(cox_biomarker,2));
                logl = NaN; H = [];
                is_firth = true;
            else
                rethrow(ME);
            end
        end
        warning('on', 'stats:coxphfit:FitWarning');
        warning('on', 'stats:coxphfit:IterationLimit');

        fprintf('Exploratory Biomarker: %s %s at %s\n', var_desc, curr_sig_name, fx_label);
        hr_unit = sig_units{vi};
        if is_firth
            fprintf('  [Cox fallback to Firth Logistic due to separation]\n');
            fprintf('Odds Ratio (OR per %s change): %.2f\n', hr_unit, exp(b));
        else
            fprintf('Hazard Ratio (HR per %s change): %.2f\n', hr_unit, exp(b));
        end
        fprintf('95%% Confidence Interval: %.2f - %.2f\n', exp(b - 1.96*stats.se), exp(b + 1.96*stats.se));
        fprintf('p-value: %.4f\n', stats.p);
        
        % --- Proportional Hazards Verification ---
        % Using scaled Schoenfeld residuals from stats.sschres
        % Need to extract individuals who had an event
        event_idx = (cox_events == 1);
        if sum(event_idx) >= 3
            sch_res = stats.sschres(event_idx);
            event_times = cox_times(event_idx);
            log_event_times = log(event_times);  % Log-transform to reduce leverage of late events
            
            % Correlation test between scaled Schoenfeld residuals and log(event times)
            % Using Spearman rank correlation for robustness to outliers.
            [r_ph, p_ph] = corr(log_event_times, sch_res, 'Type', 'Spearman');
            fprintf('Proportional Hazards Test (Spearman correlation, log time):\n');
            fprintf('  R = %.3f, p-value = %.4f\n', r_ph, p_ph);
            
            if p_ph < 0.05
                warning('Proportional hazards assumption may be violated for %s (p = %.3f)', curr_sig_name, p_ph);
            end
            fprintf('\n');

            % Plot Scaled Schoenfeld Residuals vs log(event time)
            figure('Name', ['Schoenfeld Residuals: ' var_desc ' ' curr_sig_name ' at ' fx_label], 'Position', [100, 100, 600, 500]);
            scatter(log_event_times, sch_res, 60, 'b', 'filled', 'MarkerEdgeColor', 'k'); hold on;
            
            % Add trend line
            p_fit = polyfit(log_event_times, sch_res, 1);
            x_trend = linspace(min(log_event_times), max(log_event_times), 100);
            y_trend = polyval(p_fit, x_trend);
            plot(x_trend, y_trend, 'r-', 'LineWidth', 2);
            
            yline(0, 'k--', 'LineWidth', 1.5);
            xlabel('log(Time to Event) [log(Days)]', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('Scaled Schoenfeld Residuals', 'FontSize', 12, 'FontWeight', 'bold');
            title({['Proportional Hazards Check: ' curr_sig_name], ['Spearman Corr (log time): r = ' num2str(r_ph, '%.3f') ', p = ' num2str(p_ph, '%.3f')]}, 'FontSize', 14);
            legend('Residuals', 'Trend', 'Zero Line', 'Location', 'Best');
            grid on;
            
            safe_name = strrep(curr_sig_name, '*', 'star');
            if sig_is_abs(vi), file_prefix = 'Abs_'; else, file_prefix = 'Delta_'; end
            saveas(gcf, fullfile(output_folder, ['Schoenfeld_' file_prefix safe_name '_' fx_label '_' dtype_label '.png']));
            close(gcf);
        else
            fprintf('Proportional Hazards Test: Not enough events to test assumption.\n\n');
        end
    else
        fprintf('Biomarker: %s %s at %s: Not enough complete data.\n\n', var_desc, curr_sig_name, fx_label);
    end
end

% --- Unbiased Out-of-Sample Hazard Ratio (Rigorous LOOCV) ---
fprintf('\n--- UNBIASED OUT-OF-SAMPLE HAZARD RATIO (LOOCV) ---\n');
valid_unbiased = ~isnan(risk_scores_all) & ~isnan(times_km) & ~isnan(events_km);
if sum(valid_unbiased) > 5
    unbiased_times = times_km;
    unbiased_events = events_km;
    warning('error', 'stats:coxphfit:FitWarning');
    warning('error', 'stats:coxphfit:IterationLimit');
    try
        [b_unbiased, logl_true, ~, stats_unbiased] = coxphfit(risk_scores_all(valid_unbiased), unbiased_times(valid_unbiased), ...
            'Censoring', ~unbiased_events(valid_unbiased), 'Options', statset('MaxIter', 1000000));
        is_firth_unbiased = false;
    catch ME
        if contains(ME.identifier, 'FitWarning') || contains(ME.identifier, 'IterationLimit') || contains(lower(ME.message), 'perfect')
            mdl_firth_unb = fitglm(risk_scores_all(valid_unbiased), unbiased_events(valid_unbiased), 'Distribution', 'binomial', 'LikelihoodPenalty', 'jeffreys-prior', 'Options', statset('MaxIter', 1000000));
            b_unbiased = mdl_firth_unb.Coefficients.Estimate(2:end);
            stats_unbiased.se = mdl_firth_unb.Coefficients.SE(2:end);
            stats_unbiased.p = mdl_firth_unb.Coefficients.pValue(2:end);
            logl_true = mdl_firth_unb.LogLikelihood;
            is_firth_unbiased = true;
        else
            rethrow(ME);
        end
    end
    warning('on', 'stats:coxphfit:FitWarning');
    warning('on', 'stats:coxphfit:IterationLimit');
    
    fprintf('Predictor: Unbiased Out-of-Fold Risk Score (%s)\n', fx_label);
    if is_firth_unbiased
        fprintf('  [Cox fallback to Firth Logistic due to separation]\n');
        fprintf('  Out-of-Sample Odds Ratio: %.2f\n', exp(b_unbiased));
    else
        fprintf('  Out-of-Sample Hazard Ratio: %.2f\n', exp(b_unbiased));
    end
    fprintf('  95%% Confidence Interval: %.2f - %.2f\n', ...
        exp(b_unbiased - 1.96*stats_unbiased.se), exp(b_unbiased + 1.96*stats_unbiased.se));
    fprintf('  Standard p-value (Violates Independence): %.4f\n', stats_unbiased.p);
    
    % --- Permutation Test ---
    fprintf('  Running Permutation Test (1,000 shuffles) to generate empirical p-value...\n');
    n_perms = 1000;
    logl_null = nan(n_perms, 1);
    
    valid_idx = find(valid_unbiased);
    n_valid = length(valid_idx);
    
    n_pts_impute = size(X_impute, 1);
    map_valid_to_impute = cumsum(impute_mask); 
    idx_y_clean = map_valid_to_impute(valid_idx); 
    
    parfor p_i = 1:n_perms
        % Shuffle survival times and event indicators for valid patients
        perm_idx = randperm(n_valid);
        t_perm = unbiased_times;
        e_perm = unbiased_events;
        t_perm(valid_idx) = unbiased_times(valid_idx(perm_idx));
        e_perm(valid_idx) = unbiased_events(valid_idx(perm_idx));
        
        y_perm_clean = y_clean;
        y_perm_clean(idx_y_clean) = e_perm(valid_idx);
        
        risk_scores_perm_oof = nan(n_pts_impute, 1);
        
        for loo_i = 1:n_pts_impute
            % --- Deep Learning Isolation Check ---
            pat_id_i = id_list_impute{loo_i};
            is_leaky = false;
            if dtype == 2 % dnCNN
                if any(strcmp(dl_provenance.dncnn_train_ids, pat_id_i))
                    is_leaky = true;
                end
            elseif dtype == 3 % IVIMnet
                if any(strcmp(dl_provenance.ivimnet_train_ids, pat_id_i))
                    is_leaky = true;
                end
            end
            
            if is_leaky
                % Avoid parfor crash: normally we'd error, but since it's a parallel loop
                % we'll just set a flag or let the inner logic fail. 
                % Actually, better to catch this before the parfor, but the LOOCV is nested.
                % We'll use a local error which parfor will catch.
                error('DATA LEAKAGE DETECTED (Permutation): Patient %s leaked into %s.', ...
                    pat_id_i, dtype_label);
            end

            train_mask = true(n_pts_impute, 1);
            train_mask(loo_i) = false;
            X_tr_fold = X_impute(train_mask, :);
            y_tr_fold = y_perm_clean(train_mask);
            X_te_fold = X_impute(loo_i, :);
            
            [X_tr_imp, X_te_imp] = knn_impute_train_test(X_tr_fold, X_te_fold, 5);
            keep_fold = filter_collinear_features(X_tr_imp, y_tr_fold);
            X_tr_kept = X_tr_imp(:, keep_fold);
            X_te_kept = X_te_imp(:, keep_fold);
            
            w_state = warning('off', 'all');
            try
                [B_loo, FitInfo_loo] = lassoglm(X_tr_kept, y_tr_fold, 'binomial', ...
                    'Alpha', 0.5, 'CV', 5, 'NumLambda', 25, 'Standardize', true, 'MaxIter', 1000000);
                best_idx = FitInfo_loo.IndexMinDeviance;
                coefs_loo = B_loo(:, best_idx);
                intercept_loo = FitInfo_loo.Intercept(best_idx);
            catch
                coefs_loo = zeros(size(X_tr_kept, 2), 1);
                intercept_loo = 0;
            end
            warning(w_state);
            
            risk_scores_perm_oof(loo_i) = X_te_kept * coefs_loo + intercept_loo;
        end
        
        risk_scores_perm_all = nan(sum(valid_pts), 1);
        risk_scores_perm_all(impute_mask) = risk_scores_perm_oof;
        
        w_state = warning('off', 'all');
        try
            if is_firth_unbiased
                mdl_null = fitglm(risk_scores_perm_all(valid_idx), e_perm(valid_idx), 'Distribution', 'binomial', 'LikelihoodPenalty', 'jeffreys-prior');
                logl_null(p_i) = mdl_null.LogLikelihood;
            else
                [~, logl_iter, ~, ~] = coxphfit(risk_scores_perm_all(valid_idx), t_perm(valid_idx), 'Censoring', ~e_perm(valid_idx), 'Options', statset('MaxIter', 100));
                logl_null(p_i) = logl_iter;
            end
        catch
            logl_null(p_i) = NaN;
        end
        warning(w_state);
    end
    
    emp_p = sum(logl_null(~isnan(logl_null)) >= logl_true) / sum(~isnan(logl_null));
    fprintf('  Empirical p-value: %.4f\n', emp_p);
    
    if emp_p < 0.05
        fprintf('  Conclusion: RIGOROUS. Model retains significance out-of-sample.\n');
    else
        fprintf('  Conclusion: NOT SIGNIFICANT. In-sample findings may be due to over-fitting.\n');
    end
else
    fprintf('Not enough data for unbiased HR calculation.\n');
end

% --- PART C: SUB-VOLUME STATISTICS ---
% How big were the "Resistant" (Low Diffusion) sub-volumes?
%% ---------- Time-Dependent Cox PH Model (Counting-Process / Start–Stop) ----------
% Each patient contributes one row per DWI scan interval, so the model
% tracks how the biomarker evolves up to the event rather than using only
% a static snapshot at target_fx.
%
% Panel notation:
%   t_start  = scan day (Fx1=0, Fx2≈5, Fx3≈10, Fx4≈15, Fx5≈20, Post≈90)
%   t_stop   = start of next interval or total_time if last valid scan
%   event    = 1 only on the terminal row of a local-failure patient
%
% Reference: Andersen & Gill (1982); Klein & Moeschberger §8.4.

fprintf('\n--- TIME-DEPENDENT COX PH MODEL (Counting Process) ---\n');

td_scan_days = [0, 5, 10, 15, 20, 90];   % update if exact scan dates are available

% Covariates: all four IVIM parameters (absolute, all fractions)
td_feat_arrays = { ADC_abs(valid_pts,:), D_abs(valid_pts,:), ...
                   f_abs(valid_pts,:),   Dstar_abs(valid_pts,:) };
td_feat_names  = {'ADC', 'D', 'f', 'D*'};
td_n_feat      = numel(td_feat_arrays);

td_lf       = m_lf(valid_pts);
td_tot_time = m_total_time(valid_pts);
% Censored patients use follow-up time; events use time-to-event
cens_mask_td = (td_lf == 0) & ~isnan(m_total_follow_up_time(valid_pts));
td_tot_time(cens_mask_td) = m_total_follow_up_time(valid_pts & (m_lf(:)==0) & ~isnan(m_total_follow_up_time(:)));

[X_td, t_start_td, t_stop_td, event_td, pat_id_td] = ...
    build_td_panel(td_feat_arrays, td_feat_names, td_lf, td_tot_time, nTp, td_scan_days);

td_ok = (sum(event_td) >= 3) && (size(X_td, 1) > td_n_feat + 1);
if ~td_ok
    fprintf('  Insufficient events (%d) or intervals for time-dependent Cox model.\n', sum(event_td));
else
    % ---- Fit the time-dependent Cox model --------------------------------
    % coxphfit accepts a two-column [t_start t_stop] matrix for start-stop data.
    warning('error', 'stats:coxphfit:FitWarning');
    warning('error', 'stats:coxphfit:IterationLimit');
    is_firth_td = false;
    try
        [b_td, logl_td, ~, stats_td] = coxphfit(X_td, [t_start_td, t_stop_td], ...
            'Censoring', ~event_td, 'Ties', 'breslow', ...
            'Options', statset('MaxIter', 1000000));
    catch ME_td
        if contains(ME_td.identifier, 'FitWarning') || ...
           contains(ME_td.identifier, 'IterationLimit') || ...
           contains(lower(ME_td.message), 'perfect')
            % Firth penalised logistic fallback: use last observed covariate per patient
            n_vp = sum(valid_pts);
            td_last_X = nan(n_vp, td_n_feat);
            for ii = 1:n_vp
                rows_ii = (pat_id_td == ii);
                if any(rows_ii)
                    td_last_X(ii, :) = X_td(find(rows_ii, 1, 'last'), :);
                end
            end
            valid_firth_td = ~any(isnan(td_last_X), 2) & ~isnan(td_lf);
            mdl_firth_td = fitglm(td_last_X(valid_firth_td,:), td_lf(valid_firth_td), ...
                'Distribution', 'binomial', 'LikelihoodPenalty', 'jeffreys-prior', ...
                'Options', statset('MaxIter', 1000000));
            b_td           = mdl_firth_td.Coefficients.Estimate(2:end);
            stats_td.se    = mdl_firth_td.Coefficients.SE(2:end);
            stats_td.p     = mdl_firth_td.Coefficients.pValue(2:end);
            stats_td.sschres = zeros(size(X_td,1), td_n_feat);
            logl_td        = mdl_firth_td.LogLikelihood;
            is_firth_td    = true;
        else
            warning('on', 'stats:coxphfit:FitWarning');
            warning('on', 'stats:coxphfit:IterationLimit');
            rethrow(ME_td);
        end
    end
    warning('on', 'stats:coxphfit:FitWarning');
    warning('on', 'stats:coxphfit:IterationLimit');

    % ---- Likelihood-ratio test vs. null model (no covariates) -----------
    try
        [~, logl_null_td] = coxphfit(zeros(size(X_td,1),1), [t_start_td, t_stop_td], ...
            'Censoring', ~event_td, 'Ties', 'breslow', ...
            'Options', statset('MaxIter', 100));
        LRT_stat = 2 * (logl_td - logl_null_td);
        LRT_p    = 1 - chi2cdf(LRT_stat, td_n_feat);
    catch
        LRT_stat = NaN; LRT_p = NaN;
    end

    % ---- Print results table -------------------------------------------
    if is_firth_td
        fprintf('  [Fallback to Firth logistic due to separation]\n');
        hr_label = 'OR ';
    else
        hr_label = 'HR ';
    end
    fprintf('  %-10s  %6s  %6s  %6s  %6s\n', 'Covariate', hr_label, 'CI_lo', 'CI_hi', 'p');
    fprintf('  %s\n', repmat('-', 1, 52));
    for fi = 1:td_n_feat
        hr_i  = exp(b_td(fi));
        ci_lo = exp(b_td(fi) - 1.96*stats_td.se(fi));
        ci_hi = exp(b_td(fi) + 1.96*stats_td.se(fi));
        fprintf('  %-10s  %6.3f  %6.3f  %6.3f  %6.4f\n', ...
            td_feat_names{fi}, hr_i, ci_lo, ci_hi, stats_td.p(fi));
    end
    if ~isnan(LRT_p)
        fprintf('  Global LRT: chi2(%d) = %.2f, p = %.4f\n', td_n_feat, LRT_stat, LRT_p);
    end

    % ---- Schoenfeld residuals: PH assumption check ----------------------
    if ~is_firth_td && isfield(stats_td, 'sschres')
        event_rows_td = find(event_td);
        if length(event_rows_td) >= 3
            figure('Name', ['TD Cox Schoenfeld — ' fx_label ' — ' dtype_label], ...
                'Position', [100, 100, 300*td_n_feat, 320]);
            for fi = 1:td_n_feat
                subplot(1, td_n_feat, fi);
                sch_fi   = stats_td.sschres(event_rows_td, fi);
                logt_fi  = log(t_stop_td(event_rows_td));
                scatter(logt_fi, sch_fi, 30, 'b', 'filled', 'MarkerEdgeColor', 'k');
                hold on;
                pf       = polyfit(logt_fi, sch_fi, 1);
                xt       = linspace(min(logt_fi), max(logt_fi), 50);
                plot(xt, polyval(pf, xt), 'r-', 'LineWidth', 1.5);
                yline(0, 'k--');
                [r_ph_fi, p_ph_fi] = corr(logt_fi, sch_fi, 'Type', 'Spearman');
                ph_col   = 'k'; if p_ph_fi < 0.05, ph_col = 'r'; end
                title({td_feat_names{fi}, ...
                    sprintf('r=%.2f p=%.3f', r_ph_fi, p_ph_fi)}, ...
                    'Color', ph_col, 'FontSize', 10);
                xlabel('log(t)'); ylabel('Sch. residual'); grid on;
            end
            sgtitle(['TD Cox PH Check (' fx_label ', ' dtype_label ')'], 'FontSize', 12);
            saveas(gcf, fullfile(output_folder, ...
                ['TDCox_Schoenfeld_' fx_label '_' dtype_label '.png']));
            close(gcf);
        end
    end

    fprintf('  Time-Dependent Cox model complete.\n\n');
end

% --- PART C: SUB-VOLUME STATISTICS ---
% We used: D < 0.001 and f < 0.1
% Let's look at Fx2 (where the Dose finding was significant)
% Note: 'ivim_sub_vol' contains volume in cc. 'm_gtv_vol' is total.
if target_fx <= size(ivim_sub_vol, 2)
    sub_vol_cc = ivim_sub_vol(valid_pts, target_fx, 1);
else
    sub_vol_cc = nan(sum(valid_pts), 1);
end
if target_fx <= size(m_gtv_vol, 2)
    total_vol_cc = m_gtv_vol(valid_pts, target_fx);
else
    total_vol_cc = nan(sum(valid_pts), 1);
end
sub_vol_pct = (sub_vol_cc ./ total_vol_cc) * 100;

fprintf('\n--- SUB-VOLUME CHARACTERISTICS (Fx2) ---\n');
fprintf('Resistant Sub-volume Size (%% of GTV):\n');
fprintf('  Mean: %.1f%%\n', nanmean(sub_vol_pct));
fprintf('  Range: %.1f%% - %.1f%%\n', nanmin(sub_vol_pct), nanmax(sub_vol_pct));
fprintf('  Absolute Size: %.1f ± %.1f cc\n', nanmean(sub_vol_cc), nanstd(sub_vol_cc));

if nanmean(sub_vol_pct) < 1
    fprintf('WARNING: Sub-volumes are very small (<1%%). Dose finding may be noise.\n');
else
    fprintf('VALIDATION: Sub-volumes are substantial (%.1f%%). Dose finding is physically meaningful.\n', nanmean(sub_vol_pct));
end

%% 7. Final "Defense" Analyses: ICC, Dose-Response, Multivariable Cox (Fully Fixed)
fprintf('\n--- GENERATING DEFENSE ANALYSES ---\n');

% --- B. DOSE-RESPONSE CORRELATION ---
try
    if target_fx <= size(dmean_gtvp, 2)
        dose_vals = dmean_gtvp(valid_pts, target_fx);
    else
        dose_vals = nan(sum(valid_pts), 1);
    end
    
    % Test dose-response for each significant variable
    for vi = 1:n_sig
        % Skip Dosimetry and Sub-volume metrics (indices 9-18)
        if selected_indices(vi) > 8
            continue;
        end
        
        delta_vi = sig_pct_data{vi}(valid_pts, target_fx);
        clean_mask = ~isnan(dose_vals) & ~isnan(delta_vi);
        x_dose = dose_vals(clean_mask);
        y_resp = delta_vi(clean_mask);
        
        if length(x_dose) > 2
            [r_dose, p_dose] = corr(x_dose, y_resp, 'Type', 'Spearman');
            fprintf('\nDose-Response Analysis (Delta %s):\n', sig_names{vi});
            fprintf('  Correlation: R = %.3f, p = %.3f\n', r_dose, p_dose);
            if p_dose > 0.05
                fprintf('  Conclusion: INDEPENDENT (Biological).\n');
            else
                fprintf('  Conclusion: CORRELATED (Dose-driven).\n');
            end
        end
    end
catch ME
    fprintf('Error in Dose-Response: %s\n', ME.message);
end

% --- C. PRIMARY MULTIVARIABLE COX REGRESSION (Robust Version) ---
try
    raw_lf = m_lf(valid_pts);
    raw_time_fail = m_total_time(valid_pts);
    raw_time_cens = m_total_follow_up_time(valid_pts);
    
    surv_time = raw_time_fail;
    surv_time(raw_lf == 0) = raw_time_cens(raw_lf == 0);
    surv_event = raw_lf;
    
    raw_baseline_vol = m_gtv_vol(valid_pts, 1);
    vol_vec = (raw_baseline_vol - mean(raw_baseline_vol, 'omitnan')) ./ std(raw_baseline_vol, 'omitnan');
    
    fprintf('\n--- PRIMARY MULTIVARIABLE COX REGRESSION (Unbiased Risk Score + Volume) ---\n');
    
    % risk_scores_all is aligned with sum(valid_pts) elements
    final_mask = ~isnan(surv_time) & ~isnan(surv_event) & ~isnan(risk_scores_all) & ~isnan(vol_vec);
    
    if sum(final_mask) > 5
        y_time = surv_time(final_mask);
        y_event = surv_event(final_mask);
        x_biomarker = risk_scores_all(final_mask); 
        x_vol = vol_vec(final_mask);
        
        % --- Collinearity Check ---
        [r_collin, ~] = corr(x_biomarker, x_vol, 'Type', 'Spearman');
        if abs(r_collin) > 0.70
            fprintf('  Predictor 1: Unbiased Out-of-Fold Risk Score\n');
            fprintf('  [WARNING] Critical multicollinearity with Baseline Volume (|r| = %.2f).\n', r_collin);
            fprintf('  Skipping multivariable fit to avoid mischaracterization as confounded.\n\n');
        else
            warning('error', 'stats:coxphfit:FitWarning');
            warning('error', 'stats:coxphfit:IterationLimit');
            try
                [b, logl, H, stats] = coxphfit([x_biomarker, x_vol], y_time, 'Censoring', ~y_event, 'Options', statset('MaxIter', 1000000));
                is_firth = false;
            catch ME
                if contains(ME.identifier, 'FitWarning') || contains(ME.identifier, 'IterationLimit') || contains(lower(ME.message), 'perfect')
                    mdl_firth = fitglm([x_biomarker, x_vol], y_event, 'Distribution', 'binomial', 'LikelihoodPenalty', 'jeffreys-prior', 'Options', statset('MaxIter', 1000000));
                    b = mdl_firth.Coefficients.Estimate(2:end);
                    stats.se = mdl_firth.Coefficients.SE(2:end);
                    stats.p = mdl_firth.Coefficients.pValue(2:end);
                    stats.sschres = zeros(size(x_biomarker,1), 2);
                    logl = NaN; H = [];
                    is_firth = true;
                else
                    rethrow(ME);
                end
            end
            warning('on', 'stats:coxphfit:FitWarning');
            warning('on', 'stats:coxphfit:IterationLimit');
            
            fprintf('  Predictor 1: Unbiased Out-of-Fold Risk Score\n');
            if is_firth
                fprintf('  [Cox fallback to Firth Logistic due to separation]\n');
                fprintf('     OR: %.4f (p=%.4f)\n', exp(b(1)), stats.p(1));
                fprintf('  Predictor 2: Baseline Tumor Volume (Z-score)\n');
                fprintf('     OR: %.4f (p=%.4f)\n', exp(b(2)), stats.p(2));
            else
                fprintf('     HR: %.4f (p=%.4f)\n', exp(b(1)), stats.p(1));
                fprintf('  Predictor 2: Baseline Tumor Volume (Z-score)\n');
                fprintf('     HR: %.4f (p=%.4f)\n', exp(b(2)), stats.p(2));
            end
            
            if stats.p(1) < 0.05
                fprintf('  Conclusion: PRIMARY CLAIM ROBUST. Predictive independent of volume.\n');
            else
                fprintf('  Conclusion: CONFOUNDED. Significance lost after adjustment for baseline volume.\n');
            end
            
            % --- Proportional Hazards Verification (Multivariable) ---
            event_idx = (y_event == 1);
            if sum(event_idx) >= 3
                sch_res_bio = stats.sschres(:, 1);
                sch_res_vol = stats.sschres(:, 2);
                event_times = y_time(event_idx);
                log_event_times = log(event_times);  % Log-transform to reduce leverage of late events
                
                if size(stats.sschres, 1) == length(log_event_times)
                    % Spearman rank correlation on log(event_times) for robustness
                    [r_ph_bio, p_ph_bio] = corr(log_event_times, sch_res_bio, 'Type', 'Spearman');
                    [r_ph_vol, p_ph_vol] = corr(log_event_times, sch_res_vol, 'Type', 'Spearman');
                    
                    fprintf('  Proportional Hazards Test (Spearman correlation, log time):\n');
                    fprintf('     Variable 1 (Risk Score): R = %.3f, p-value = %.4f\n', r_ph_bio, p_ph_bio);
                    fprintf('     Variable 2 (Volume): R = %.3f, p-value = %.4f\n', r_ph_vol, p_ph_vol);
                else
                    fprintf('  Proportional Hazards Test: Residual count (%d) mismatch with event count (%d). Skipping.\n', size(stats.sschres, 1), length(log_event_times));
                end
                
                if p_ph_bio < 0.05
                    warning('Proportional hazards assumption may be violated for multivariable Risk Score (p = %.3f)', p_ph_bio);
                end
                if p_ph_vol < 0.05
                    warning('Proportional hazards assumption may be violated for multivariable Volume (p = %.3f)', p_ph_vol);
                end
                
                % Plot for Risk Score vs log(event time)
                figure('Name', ['MV Schoenfeld: Risk_Score at ' fx_label], 'Position', [100, 100, 600, 500]);
                scatter(log_event_times, sch_res_bio, 60, 'b', 'filled', 'MarkerEdgeColor', 'k'); hold on;
                p_fit = polyfit(log_event_times, sch_res_bio, 1);
                x_trend = linspace(min(log_event_times), max(log_event_times), 100);
                y_trend = polyval(p_fit, x_trend);
                plot(x_trend, y_trend, 'r-', 'LineWidth', 2);
                yline(0, 'k--', 'LineWidth', 1.5);
                xlabel('log(Time to Event) [log(Days)]', 'FontSize', 12, 'FontWeight', 'bold');
                ylabel('Scaled Schoenfeld Residuals', 'FontSize', 12, 'FontWeight', 'bold');
                title({['MV Proportional Hazards Check: Unbiased Risk Score'], ['Spearman Corr (log time): r = ' num2str(r_ph_bio, '%.3f') ', p = ' num2str(p_ph_bio, '%.3f')]}, 'FontSize', 14);
                legend('Residuals', 'Trend', 'Zero Line', 'Location', 'Best');
                grid on;
                saveas(gcf, fullfile(output_folder, ['MV_Schoenfeld_RiskScore_' fx_label '_' dtype_label '.png']));
                close(gcf);
            else
                fprintf('  Proportional Hazards Test: Not enough events to test assumption.\n');
            end
            fprintf('\n');
    end
    else
        fprintf('  Not enough complete data for multivariable Cox.\n\n');
    end

catch ME
    fprintf('Error in Cox Regression: %s\n', ME.message);
end
end % for target_fx

%% ---------- Q-Value (FDR) Analysis — Per-Timepoint BH ----------
% BH is applied SEPARATELY per timepoint to avoid inflating the penalty
% by pooling serial correlated measurements of the same patients.
% Within each timepoint the family = one Wilcoxon test per metric/set.
fprintf('\n--- PER-TIMEPOINT FDR + HOLM-BONFERRONI ANALYSIS ---\n');

for tp = 1:length(time_labels)
    tp_pvals  = [];
    tp_labels = {};
    
    for s = 1:length(metric_sets)
        current_metrics = metric_sets{s};
        current_names   = set_names{s};
        for mi = 1:length(current_metrics)
            metric_data = current_metrics{mi};
            if tp > size(metric_data, 2), continue; end
            y_raw = metric_data(valid_pts, tp);
            g = lf_group(~isnan(y_raw));
            y = y_raw(~isnan(y_raw));
            if length(unique(g)) > 1 && length(y) > 2
                tp_pvals(end+1, 1)  = ranksum(y(g==0), y(g==1));
                tp_labels{end+1, 1} = [current_names{mi} ' (' time_labels{tp} ')'];
            end
        end
    end
    
    if isempty(tp_pvals), continue; end
    n_tp = length(tp_pvals);
    alpha = 0.05;
    
    % --- BH FDR ---
    [p_sort, sort_id] = sort(tp_pvals);
    q_tp = zeros(n_tp, 1);
    q_tp(n_tp) = p_sort(n_tp);
    for ii = n_tp-1:-1:1
        q_tp(ii) = min(q_tp(ii+1), p_sort(ii) * (n_tp / ii));
    end
    q_tp = min(q_tp, 1);
    q_unsorted = zeros(n_tp, 1);
    q_unsorted(sort_id) = q_tp;
    
    tp_table = table(tp_labels, tp_pvals, q_unsorted, ...
        'VariableNames', {'Metric', 'P_Value', 'Q_Value'});
    tp_table = sortrows(tp_table, 'P_Value');
    
    % --- Holm-Bonferroni (FWER) ---
    holm_thresh = alpha ./ (n_tp + 1 - (1:n_tp)');
    is_sig_holm = false(n_tp, 1);
    for k = 1:n_tp
        if p_sort(k) < holm_thresh(k), is_sig_holm(k) = true; else, break; end
    end
    
    fprintf('\n=== %s (family n=%d) ===\n', time_labels{tp}, n_tp);
    fprintf('  Top 5 by p-value:\n');
    disp(tp_table(1:min(5, height(tp_table)), :));
    if any(is_sig_holm)
        fprintf('  Holm-Bonferroni: %d test(s) survived FWER control.\n', sum(is_sig_holm));
    else
        fprintf('  Holm-Bonferroni: No findings survived FWER (best p=%.4f vs threshold=%.4f).\n', ...
            p_sort(1), holm_thresh(1));
    end
    
    % Plot P vs Q for this timepoint
    figure('Name', ['P vs Q — ' time_labels{tp} ' — ' dtype_label], 'Position', [300 300 600 500]);
    scatter(tp_table.P_Value, tp_table.Q_Value, 50, 'filled');
    xlabel('Raw P-Value'); ylabel('Q-Value (FDR)');
    title(['P vs Q: ' time_labels{tp} ' (' dtype_label ')']);
    yline(0.05, 'r--', 'Q = 0.05');
    grid on;
    safe_tp = strrep(time_labels{tp}, ' ', '_');
    saveas(gcf, fullfile(output_folder, ['P_vs_Q_' safe_tp '_' dtype_label '.png']));
    close(gcf);
end

%% ---------- 8. Longitudinal Mixed-Effects Model (GLME) ----------
fprintf('\n--- LONGITUDINAL MIXED-EFFECTS MODEL (GLME) ---\n');
% Build long-format table across all timepoints (1 to nTp) for patients with known outcomes
long_PatientID = [];
long_Timepoint = [];
long_ADC = [];
long_D = [];
long_f = [];
long_Dstar = [];
long_LF = [];

patient_indices = find(valid_pts);
for i = 1:length(patient_indices)
    p_idx = patient_indices(i);
    for t = 1:nTp
        % Only include if at least one metric is not NaN
        if ~isnan(ADC_abs(p_idx, t)) || ~isnan(D_abs(p_idx, t)) || ~isnan(f_abs(p_idx, t)) || ~isnan(Dstar_abs(p_idx, t))
            long_PatientID = [long_PatientID; i]; 
            long_Timepoint = [long_Timepoint; t];
            long_ADC = [long_ADC; ADC_abs(p_idx, t)];
            long_D = [long_D; D_abs(p_idx, t)];
            long_f = [long_f; f_abs(p_idx, t)];
            long_Dstar = [long_Dstar; Dstar_abs(p_idx, t)];
            long_LF = [long_LF; lf_group(i)];
        end
    end
end

glme_table = table(categorical(long_PatientID), long_Timepoint, ...
    long_ADC, long_D, long_f, long_Dstar, long_LF, ...
    'VariableNames', {'PatientID', 'Timepoint', 'ADC', 'D', 'f', 'Dstar', 'LF'});

clean_idx = ~isnan(glme_table.ADC) & ~isnan(glme_table.D) & ~isnan(glme_table.f) & ~isnan(glme_table.Dstar);
glme_table_clean = glme_table(clean_idx, :);
baseline_idx = glme_table_clean.Timepoint == 1;

mean_ADC_base = mean(glme_table_clean.ADC(baseline_idx));
std_ADC_base = std(glme_table_clean.ADC(baseline_idx));
mean_D_base = mean(glme_table_clean.D(baseline_idx));
std_D_base = std(glme_table_clean.D(baseline_idx));
mean_f_base = mean(glme_table_clean.f(baseline_idx));
std_f_base = std(glme_table_clean.f(baseline_idx));
mean_Dstar_base = mean(glme_table_clean.Dstar(baseline_idx));
std_Dstar_base = std(glme_table_clean.Dstar(baseline_idx));

glme_table_clean.ADC_z = (glme_table_clean.ADC - mean_ADC_base) / std_ADC_base;
glme_table_clean.D_z = (glme_table_clean.D - mean_D_base) / std_D_base;
glme_table_clean.f_z = (glme_table_clean.f - mean_f_base) / std_f_base;
glme_table_clean.Dstar_z = (glme_table_clean.Dstar - mean_Dstar_base) / std_Dstar_base;

% Ensure LF and Timepoint are categorical to properly evaluate interaction
glme_table_clean.LF = categorical(glme_table_clean.LF);
glme_table_clean.Timepoint = categorical(glme_table_clean.Timepoint);

biomarkers = {'ADC_z', 'D_z', 'f_z', 'Dstar_z'};
warning('off', 'all');
for b = 1:length(biomarkers)
    bm = biomarkers{b};
    formula = sprintf('%s ~ 1 + LF * Timepoint + (1|PatientID)', bm);
    try
        glme = fitglme(glme_table_clean, formula, 'OptimizerOptions', statset('MaxIter', 1000000));
        fprintf('\n--- %s ---\n', formula);
        
        anova_res = anova(glme);
        row_idx = find(strcmp(anova_res.Term, 'LF:Timepoint'));
        if ~isempty(row_idx)
            pval = anova_res.pValue(row_idx);
            fprintf('Interaction P-Value (LF * Timepoint): %.4f\n', pval);
            if pval < 0.05
                fprintf('  -> SIGNIFICANT DIFFERENCE in trajectory between LC and LF groups.\n');
            else
                fprintf('  -> No significant difference in trajectory between LC and LF groups.\n');
            end
        else
            disp(anova_res);
        end
    catch ME
        fprintf('GLME model for %s failed to converge: %s\n', bm, ME.message);
    end
end
warning('on', 'all');

fprintf('\n=== Completed DWI Type %d: %s ===\n', dtype, dtype_label);
end % for dtype
diary off

    % Filters collinear features using Pearson |r| > 0.8 threshold.
    % When two features are highly correlated, the one with the higher univariate
    % Wilcoxon rank-sum p-value (less significant) is dropped.
    R = corrcoef(X);
    drop_idx = false(1, size(X, 2));
    for fi = 1:size(X, 2)
        if drop_idx(fi), continue; end
        for fj = fi+1:size(X, 2)
            if abs(R(fi, fj)) > 0.8
                p_fi = ranksum(X(y==0, fi), X(y==1, fi));
                p_fj = ranksum(X(y==0, fj), X(y==1, fj));
                if p_fj >= p_fi
                    drop_idx(fj) = true;
                else
                    drop_idx(fi) = true;
                    break;
                end
            end
        end
    end
    keep_idx = find(~drop_idx);
end