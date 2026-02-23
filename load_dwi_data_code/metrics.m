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
%   7. LASSO-regularized feature selection (binomial logistic, 5-fold CV)
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
fprintf('Median time to LF     = %2d days (%d - %d)\n',nanmedian(total_time(lf==1)),nanmin(total_time(lf==1)),nanmax(total_time(lf==1)));%% ========================================================================
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
% (Debug prints removed)
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
% Formula: ((Fx_n - Fx_1) ./ (0.5 * (Fx_n + Fx_1))) * 100
ADC_pct = ((ADC_abs - ADC_abs(:,1)) ./ (0.5 * (ADC_abs + ADC_abs(:,1)))) * 100;
D_pct   = ((D_abs - D_abs(:,1)) ./ (0.5 * (D_abs + D_abs(:,1)))) * 100;
f_pct   = ((f_abs - f_abs(:,1)) ./ (0.5 * (f_abs + f_abs(:,1)))) * 100;
Dstar_pct = ((Dstar_abs - Dstar_abs(:,1)) ./ (0.5 * (Dstar_abs + Dstar_abs(:,1)))) * 100;

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
    title(['\Delta ', metric_names{i}, ' (%)'], 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('% Change from Fx1');
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
    for k = 1:nTp
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
            
            % Parameters for 3D morphological cleanup
            se = strel('sphere', 1);
            min_cc_voxels = 10;
            
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
            if ~isempty(dose_adc_sub)
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
            if ~isempty(dose_d_sub)
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
            if ~isempty(dose_f_sub)
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
            if ~isempty(dose_dstar_sub)
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
    {'\Delta ADC (%)', '\Delta D (%)', '\Delta f (%)', '\Delta D* (%)'}, ...
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
%  SECTION 9: FDR Correction — Benjamini-Hochberg Procedure
% =========================================================================
% Controls the false discovery rate (FDR) across ALL tested comparisons.
% Steps:
%   1. Re-collect every raw p-value (not just p < 0.05)
%   2. Rank p-values, compute BH q-values: q_k = p_k * (m / k)
%   3. Enforce monotonicity from bottom up
%   4. Report comparisons surviving FDR < 0.05
if ~isempty(sig_pval)
    % 1. Get all p-values (not just the significant ones)
    % We need to re-run the loop to capture ALL p-values, not just p<0.05
    all_pvals = [];
    all_labels = {};
    
    for s = 1:length(metric_sets)
        current_metrics = metric_sets{s};
        current_names = set_names{s};
        for m = 1:length(current_metrics)
            metric_data = current_metrics{m};
            cols = min(size(metric_data, 2), length(time_labels));
            for tp = 1:cols
                y_raw = metric_data(valid_pts, tp);
                g = lf_group(~isnan(y_raw));
                y = y_raw(~isnan(y_raw));
                if length(unique(g)) > 1 && length(y) > 2
                    p = ranksum(y(g==0), y(g==1));
                    all_pvals(end+1, 1) = p;
                    all_labels{end+1, 1} = [current_names{m} ' (' time_labels{tp} ')'];
                end
            end
        end
    end
    
    % 2. Calculate FDR adjusted p-values (q-values)
    m = length(all_pvals);             % Total number of tests
    [sorted_p, sort_ids] = sort(all_pvals);
    
    % Benjamini-Hochberg Critical Value calculation: (rank / total_tests) * alpha
    % We calculate the q-value (adjusted p) directly:
    q_values = zeros(size(sorted_p));
    q_values(m) = sorted_p(m);
    for i = m-1:-1:1
        q_values(i) = min(q_values(i+1), sorted_p(i) * (m/i));
    end
    
    % Map back to original order
    q_values_unsorted(sort_ids) = q_values;
    
    % 3. Display Robust Results (FDR < 0.05)
    results_table = table(all_labels, all_pvals, q_values_unsorted', ...
        'VariableNames', {'Comparison', 'Raw_P', 'FDR_Q'});
    
    significant_fdr = results_table(results_table.FDR_Q < 0.05, :);
    
    fprintf('\n----- FDR CORRECTED SIGNIFICANT RESULTS (Q < 0.05) -----\n');
    if isempty(significant_fdr)
        disp('None of the findings survived FDR correction.');
    else
        disp(significant_fdr);
        writetable(significant_fdr, [dataloc 'FDR_Significant_Results.csv']);
    end
end

%% ========================================================================
%  SECTION 10: Per-Timepoint Analysis Loop (Fx2 and Fx3)
% =========================================================================
% For each target fraction (Fx2, Fx3), this loop performs:
%   a) LASSO-regularized feature selection (binomial logistic, 5-fold CV)
%   b) ROC analysis with Youden's J optimal cutoff
%   c) Multivariable logistic regression (Firth penalized fallback)
%   d) FDR and Holm-Bonferroni multiple-comparison corrections
%   e) Leave-pair-out cross-validation (LPOCV) for unbiased AUC
%   f) 2D scatter plot with logistic decision boundary
%   g) Kaplan-Meier survival analysis
%   h) Correlation matrix of significant features
%   i) Representative parametric ADC maps (responder vs non-responder)
for target_fx = 2:nTp
fx_label = sprintf('Fx%d', target_fx);
fprintf('\n=== Analyzing %s ===\n', fx_label);

%% ---------- LASSO Feature Selection at This Timepoint ----------
% Determine which DWI/IVIM metrics are informative for predicting LF.
% A base variable (ADC, D, f, D*) is flagged if LASSO selects either
% its absolute or percent-change form.
base_metric_names_all = {'ADC', 'D', 'f', 'D*'};
all_abs_data = {ADC_abs, D_abs, f_abs, Dstar_abs};       % Absolute values
all_pct_data = {ADC_pct, D_pct, f_pct, Dstar_pct};       % Percent change

sig_flags = false(1, 4);   % Will be set true for LASSO-selected metrics
sig_p_best = ones(1, 4);   % Best univariate p-value (for downstream sorting)

%% --- LASSO Feature Selection ---
    % Construct feature matrix: 8 columns = [4 Absolute | 4 %-Change]
    % Rows = patients with valid clinical outcome data
    X_lasso = [ADC_abs(valid_pts, target_fx), D_abs(valid_pts, target_fx), ...
               f_abs(valid_pts, target_fx),   Dstar_abs(valid_pts, target_fx), ...
               ADC_pct(valid_pts, target_fx), D_pct(valid_pts, target_fx), ...
               f_pct(valid_pts, target_fx),   Dstar_pct(valid_pts, target_fx)];
           
    feat_names_lasso = {'ADC_Abs', 'D_Abs', 'f_Abs', 'Dstar_Abs', ...
                        'ADC_Pct', 'D_Pct', 'f_Pct', 'Dstar_Pct'};
    y_lasso = lf_group;

    % --- Imputation strategy to prevent complete-case attrition ---
    % Step 1: Exclude patients missing ALL imaging data (entire row NaN)
    has_any_imaging = any(~isnan(X_lasso), 2);
    % Step 2: Also require a valid clinical outcome
    impute_mask = has_any_imaging & ~isnan(y_lasso);
    X_impute = X_lasso(impute_mask, :);
    y_clean  = y_lasso(impute_mask);

    % --- Proper 5-fold CV for Lambda Selection to prevent Data Leakage ---
    rng(42); % Set seed for reproducibility
    cvp = cvpartition(y_clean, 'KFold', 5);
    n_lambdas = 25;

    % Generate common lambda sequence using dummy full-data pass
    dummy_clean = knn_impute_train_test(X_impute, [], 5);
    [~, FitInfo_dummy] = lassoglm(dummy_clean, y_clean, 'binomial', 'NumLambda', n_lambdas, 'Standardize', true, 'MaxIter', 1000000);
    if length(FitInfo_dummy.Lambda) < n_lambdas
        common_Lambda = [FitInfo_dummy.Lambda, zeros(1, n_lambdas - length(FitInfo_dummy.Lambda))];
    else
        common_Lambda = FitInfo_dummy.Lambda;
    end

    dev = zeros(n_lambdas, cvp.NumTestSets);
    for fold_i = 1:cvp.NumTestSets
        tr_idx = training(cvp, fold_i);
        te_idx = test(cvp, fold_i);
        
        X_tr = X_impute(tr_idx, :);
        X_te = X_impute(te_idx, :);
        y_tr = y_clean(tr_idx);
        y_te = y_clean(te_idx);
        
        % Impute using KNN fitted strictly on training data
        [X_tr, X_te] = knn_impute_train_test(X_tr, X_te, 5);
        
        % Pre-filter: drop features with Pearson |r| > 0.8 using ONLY training data
        R_fold = corrcoef(X_tr);
        drop_fold = false(1, size(X_tr, 2));
        for fi = 1:size(X_tr, 2)
            if drop_fold(fi), continue; end
            for fj = fi+1:size(X_tr, 2)
                if abs(R_fold(fi, fj)) > 0.8
                    drop_fold(fj) = true;
                end
            end
        end
        keep_fold = find(~drop_fold);
        X_tr_kept = X_tr(:, keep_fold);
        X_te_kept = X_te(:, keep_fold);
        
        try
            [B_fold, FitInfo_fold] = lassoglm(X_tr_kept, y_tr, 'binomial', ...
                'Lambda', common_Lambda, 'Standardize', true, 'MaxIter', 1000000);
            
            for lam_idx = 1:length(common_Lambda)
                if common_Lambda(lam_idx) == 0, continue; end
                coefs = B_fold(:, lam_idx);
                intercept = FitInfo_fold.Intercept(lam_idx);
                
                eta = X_te_kept * coefs + intercept;
                mu = 1 ./ (1 + exp(-eta));
                mu(mu == 0) = eps; mu(mu == 1) = 1-eps; 
                
                dev(lam_idx, fold_i) = 2 * sum(-y_te .* log(mu) - (1 - y_te) .* log(1 - mu));
            end
        catch
            dev(:, fold_i) = Inf; 
        end
    end

    % Find optimal lambda
    mean_dev = mean(dev, 2);
    [~, best_lam_idx] = min(mean_dev);
    best_lambda = common_Lambda(best_lam_idx);

    % --- Final Model Fit on Full Data ---
    X_clean = knn_impute_train_test(X_impute, [], 5);
    R_corr = corrcoef(X_clean);
    drop_flag = false(1, size(X_clean, 2));
    for fi = 1:size(X_clean, 2)
        if drop_flag(fi), continue; end
        for fj = fi+1:size(X_clean, 2)
            if abs(R_corr(fi, fj)) > 0.8
                drop_flag(fj) = true;
                fprintf('  Dropping %s (|r| = %.2f with %s)\n', ...
                    feat_names_lasso{fj}, abs(R_corr(fi, fj)), feat_names_lasso{fi});
            end
        end
    end
    keep_idx = find(~drop_flag);
    X_clean = X_clean(:, keep_idx);
    feat_names_lasso_kept = feat_names_lasso(keep_idx);
    fprintf('  Retained %d / %d features after corr filter\n', length(keep_idx), length(feat_names_lasso));

    try
        [B_lasso, FitInfo] = lassoglm(X_clean, y_clean, 'binomial', ...
            'Lambda', best_lambda, 'Standardize', true, 'MaxIter', 1000000);
        
        coefs_lasso = B_lasso(:, 1);
        selected_kept = find(coefs_lasso ~= 0);
        selected_indices = keep_idx(selected_kept);
        
        fprintf('Elastic Net Selected Features for %s: %s\n', ...
            fx_label, strjoin(feat_names_lasso(selected_indices), ', '));
    catch
        fprintf('Elastic Net failed to converge. Fallback to empty selection.\n');
        selected_indices = [];
    end

    %% --- End LASSO Feature Selection ---

% -----------------------------------------------------------------------
% Post-LASSO feature list: all 8 candidates treated as fully independent.
% Indices 1-4 = Absolute forms; indices 5-8 = Percent-Change forms.
% Both forms of ANY metric may be retained simultaneously.
% -----------------------------------------------------------------------
all_feat_data  = {ADC_abs,       D_abs,       f_abs,       Dstar_abs, ...
                  ADC_pct,       D_pct,       f_pct,       Dstar_pct};
all_feat_names = {'ADC',         'D',         'f',         'D*', ...
                  'ADC',         'D',         'f',         'D*'};
all_feat_is_abs = [true          true         true         true  ...
                   false         false        false        false ];
all_feat_disp  = {'Abs ADC',     'Abs D',     'Abs f',     'Abs D*', ...
                  '\Delta ADC',  '\Delta D',  '\Delta f',  '\Delta D*'};

% n_sig = number of features LASSO selected (up to 8, not capped at 4)
n_sig = length(selected_indices);

% Build parallel arrays consumed by all downstream loops
sig_data_selected = cell(1, n_sig);
sig_abs_data      = cell(1, n_sig);
sig_pct_data      = cell(1, n_sig);
sig_names         = cell(1, n_sig);
sig_is_abs        = false(1, n_sig);
sig_disp_names    = cell(1, n_sig);

for si = 1:n_sig
    fi = selected_indices(si);              % direct column index (1-8)
    sig_data_selected{si} = all_feat_data{fi};
    sig_names{si}         = all_feat_names{fi};
    sig_is_abs(si)        = all_feat_is_abs(fi);
    sig_disp_names{si}    = all_feat_disp{fi};
    % Carry both Abs and Pct forms for dose-response / Table-1 use
    if fi <= 4
        sig_abs_data{si} = all_feat_data{fi};
        sig_pct_data{si} = all_feat_data{fi + 4};
    else
        sig_abs_data{si} = all_feat_data{fi - 4};
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

%% ---------- ROC Analysis (Receiver Operating Characteristic) ----------
% Computes ROC curves for each significant percent-change variable at the
% target fraction. Fits a univariate logistic regression to handle the
% direction of association automatically, then passes predicted
% probabilities to perfcurve.
% Youden's J statistic (max[Sensitivity - FPR]) determines optimal cutoff.

% Extract the percent change data for significant variables at target fraction
labels = lf_group; % 0 = LC, 1 = LF

% Pre-allocate ROC storage: FPR (X), TPR (Y), thresholds (T), AUC, cutoff
% Note: clear 'lines' to prevent any workspace variable named 'lines' from
% shadowing the built-in lines() colormap function (test.m can leak this).
if exist('lines', 'var'), clear lines; end
roc_colors = lines(n_sig);       % Distinct colours per variable
roc_X = cell(1, n_sig);          % False positive rate vectors
roc_Y = cell(1, n_sig);          % True positive rate vectors
roc_T = cell(1, n_sig);          % Threshold vectors
roc_AUC = zeros(1, n_sig);       % Area under the curve
roc_opt_idx = zeros(1, n_sig);   % Index of Youden-optimal point
roc_opt_thresh = zeros(1, n_sig);% Optimal %-change cutoff (un-negated)
sig_fx_data = cell(1, n_sig);    % Cleaned predictor data at target Fx

for vi = 1:n_sig
    sig_fx_data{vi} = sig_data_selected{vi}(valid_pts, target_fx);
    valid_vi = ~isnan(sig_fx_data{vi}) & ~isnan(labels);
    data_clean = sig_fx_data{vi}(valid_vi);
    labels_clean = labels(valid_vi);

    % Fit univariate logistic regression to determine direction automatically
    mdl_roc = fitglm(data_clean, labels_clean, 'Distribution', 'binomial', 'Options', statset('MaxIter', 1000000));
    [roc_X{vi}, roc_Y{vi}, roc_T{vi}, roc_AUC(vi)] = perfcurve(labels_clean, mdl_roc.Fitted.Probability, 1);

    % Optimal cutoff via Youden's J Statistic
    [~, roc_opt_idx(vi)] = max(roc_Y{vi} - roc_X{vi});
    roc_opt_thresh(vi) = roc_T{vi}(roc_opt_idx(vi));
end

% --- Plotting the ROC Curves ---
figure('Name', ['ROC Analysis - ' fx_label ' — ' dtype_label], 'Position', [200, 200, 700, 600]);
hold on;

leg_entries = {};
feat_disp_names = sig_disp_names;   % pre-built during feature-mapping above
for vi = 1:n_sig
    plot(roc_X{vi}, roc_Y{vi}, '-', 'Color', roc_colors(vi,:), 'LineWidth', 2);
    leg_entries{end+1} = sprintf('%s (AUC = %.3f)', feat_disp_names{vi}, roc_AUC(vi));
end
for vi = 1:n_sig
    plot(roc_X{vi}(roc_opt_idx(vi)), roc_Y{vi}(roc_opt_idx(vi)), 'o', ...
        'Color', roc_colors(vi,:), 'MarkerSize', 8, 'MarkerFaceColor', roc_colors(vi,:));
    leg_entries{end+1} = sprintf('Optimal %s Cutoff', feat_disp_names{vi});
end

plot([0 1], [0 1], 'k--', 'LineWidth', 1.5);
leg_entries{end+1} = 'Random Guess';

xlabel('False Positive Rate (1 - Specificity)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('True Positive Rate (Sensitivity)', 'FontSize', 12, 'FontWeight', 'bold');
sig_names_str = strjoin(feat_disp_names, ' and ');
title(['ROC Curve: ' sig_names_str ' Predicting Local Failure (' fx_label ', ' dtype_label ')'], 'FontSize', 14);

legend(leg_entries, 'Location', 'SouthEast', 'FontSize', 11);
grid on; box on;
hold off;
saveas(gcf, fullfile(output_folder, ['ROC_' fx_label '_' dtype_label '.png']));
close(gcf);

% Print results to the command window
fprintf('\n--- ROC Analysis Results for %s ---\n', fx_label);
for vi = 1:n_sig
    if sig_is_abs(vi), unit_str = ''; else, unit_str = '%%'; end
    fprintf('%s: \n  AUC = %.3f \n  Optimal Threshold (Youden) = %.2f%s\n', ...
        feat_disp_names{vi}, roc_AUC(vi), roc_opt_thresh(vi), unit_str);
    fprintf('  Sensitivity at Threshold = %.2f%%\n  Specificity at Threshold = %.2f%%\n\n', ...
        roc_Y{vi}(roc_opt_idx(vi))*100, (1-roc_X{vi}(roc_opt_idx(vi)))*100);
end



%% ---------- Leave-Pair-Out Cross-Validation (LPOCV) ----------
% LPOCV provides an unbiased AUC estimate by testing every possible
% (LF, LC) pair. For each pair, a model is trained on all other patients
% and scores the held-out pair. A pair is concordant if the model assigns
% a higher P(LF) to the actual LF patient. The proportion of concordant
% pairs equals the AUC. A heatmap shows which specific pairs were
% correctly ranked.
% We test every possible pair of (1 LF vs 1 LC).

% 1. Setup Data (all candidate features for unbiased feature selection)
labels = lf_group;

% Build predictor matrix from ALL candidate features (8 columns):
% [4 Absolute | 4 %-Change] to allow LASSO to select within each fold
lpocv_cols = [ADC_abs(valid_pts, target_fx), D_abs(valid_pts, target_fx), ...
              f_abs(valid_pts, target_fx),   Dstar_abs(valid_pts, target_fx), ...
              ADC_pct(valid_pts, target_fx), D_pct(valid_pts, target_fx), ...
              f_pct(valid_pts, target_fx),   Dstar_pct(valid_pts, target_fx)];

% Filter valid rows: Require valid outcome and at least some imaging data
has_any_imaging = any(~isnan(lpocv_cols), 2);
valid_rows = ~isnan(labels) & has_any_imaging;
X_data = lpocv_cols(valid_rows, :);
Y_data = labels(valid_rows);

% 2. Identify Indices for Each Class
idx_lf = find(Y_data == 1); % Local Failures
idx_lc = find(Y_data == 0); % Local Controls
n_lf = length(idx_lf);
n_lc = length(idx_lc);

fprintf('\n--- Starting LPOCV (Leave-Pair-Out) ---\n');
fprintf('Testing %d pairs (%d LF x %d LC)...\n', n_lf * n_lc, n_lf, n_lc);

% 3. Iterate Through Every Pair
pair_scores = zeros(n_lf, n_lc); % Store 1 (correct), 0 (wrong), 0.5 (tie)

for i = 1:n_lf
    for j = 1:n_lc
        % The specific pair to leave out
        current_lf = idx_lf(i);
        current_lc = idx_lc(j);
        
        % Training set: All data EXCEPT these two
        train_mask = true(length(Y_data), 1);
        train_mask([current_lf, current_lc]) = false;
        
        X_train = X_data(train_mask, :);
        Y_train = Y_data(train_mask);
        
        % Impute train and test sets using KNN fitted strictly on training data
        X_test_pair = X_data([current_lf, current_lc], :);
        [X_train, X_test_pair] = knn_impute_train_test(X_train, X_test_pair, 5);
        
        % Pre-filter: drop features with Pearson |r| > 0.8 using ONLY training data
        % (mirrors the identical filter in the main LASSO CV fold loop)
        R_fold_lp = corrcoef(X_train);
        drop_fold_lp = false(1, size(X_train, 2));
        for fi_lp = 1:size(X_train, 2)
            if drop_fold_lp(fi_lp), continue; end
            for fj_lp = fi_lp+1:size(X_train, 2)
                if abs(R_fold_lp(fi_lp, fj_lp)) > 0.8
                    drop_fold_lp(fj_lp) = true;
                end
            end
        end
        keep_fold_lpocv = find(~drop_fold_lp);
        X_train_kept = X_train(:, keep_fold_lpocv);
        X_test_pair_kept = X_test_pair(:, keep_fold_lpocv);
        
        % Feature selection via Elastic Net on correlation-filtered training data only (no leakage)
        warning('off', 'all');
        try
            [B_cv, FitInfo_cv] = lassoglm(X_train_kept, Y_train, 'binomial', ...
                'Alpha', 0.5, 'CV', 5, 'NumLambda', 25, 'Standardize', true, 'MaxIter', 1000000);
            idx_min = FitInfo_cv.IndexMinDeviance;
            % sel_features indexes into the kept (filtered) columns of X_train_kept
            sel_features = find(B_cv(:, idx_min) ~= 0);
        catch
            sel_features = [];
        end
        warning('on', 'all');
        
        % If LASSO selects zero features, default to a tie (AUC = 0.5)
        if isempty(sel_features)
            pair_scores(i,j) = 0.5;
            continue;
        end
        
        % Train logistic regression on selected features only
        % Suppress iteration-limit and perfect-separation warnings that
        % commonly occur with small leave-pair-out training sets.
        glmfit_opts = statset('MaxIter', 10000);
        warning('off', 'stats:glmfit:PerfectSeparation');
        warning('off', 'stats:glmfit:IterationLimit');
        try
            mdl_cv = fitglm(X_train_kept(:, sel_features), Y_train, 'Distribution', 'binomial', 'Options', glmfit_opts);
        catch ME
            if strcmp(ME.identifier, 'stats:glmfit:PerfectSeparation') || ...
               contains(ME.message, 'perfectly separate')
                mdl_cv = fitglm(X_train_kept(:, sel_features), Y_train, 'Distribution', 'binomial', 'LikelihoodPenalty', 'jeffreys-prior', 'Options', statset('MaxIter', 1000000));
            else
                warning('on', 'stats:glmfit:PerfectSeparation');
                warning('on', 'stats:glmfit:IterationLimit');
                pair_scores(i,j) = 0.5;
                continue;
            end
        end
        warning('on', 'stats:glmfit:PerfectSeparation');
        warning('on', 'stats:glmfit:IterationLimit');
        
        % Predict Probabilities for the held-out pair
        % (X_test_pair_kept is already column-filtered; sel_features
        %  indexes into those filtered columns, matching X_train_kept)
        prob_lf = predict(mdl_cv, X_test_pair_kept(1, sel_features));
        prob_lc = predict(mdl_cv, X_test_pair_kept(2, sel_features));
        
        % Score the Pair (Concordance)
        if prob_lf > prob_lc
            pair_scores(i,j) = 1;   % Correct ranking
        elseif prob_lf == prob_lc
            pair_scores(i,j) = 0.5; % Tie
        else
            pair_scores(i,j) = 0;   % Incorrect ranking
        end
    end
end

% 4. Calculate LPOCV AUC
% The AUC is simply the average of these pairwise scores
auc_lpocv = mean(pair_scores(:));

fprintf('LPOCV AUC = %.3f\n', auc_lpocv);

% Final Validation Output
fprintf('------------------------------------------------\n');
fprintf('Unbiased Validation (LPOCV):\n');
fprintf('  LPOCV AUC: %.3f\n', auc_lpocv);
fprintf('  (Note: Optimistic training AUCs have been removed to prevent overfitting bias)\n');
fprintf('------------------------------------------------\n');

% 5. Visualization of the "Pairwise Matrix"
% This heatmap shows which specific pairs were hard to classify
figure('Name', ['LPOCV Pairwise Success ' fx_label ' — ' dtype_label], 'Position', [400, 400, 600, 500]);
imagesc(pair_scores);
colormap(gray); % White = Correct (1), Black = Wrong (0)
title(['LPOCV Pairwise Ranking ' fx_label ' (AUC = ' num2str(auc_lpocv, '%.3f') ') (' dtype_label ')']);
xlabel('Local Control Patients (Index)');
ylabel('Local Failure Patients (Index)');
colorbar;
axis square;
saveas(gcf, fullfile(output_folder, ['LPOCV_Pairwise_' fx_label '_' dtype_label '.png']));
close(gcf);

%% ---------- 2D Scatter Plot with Logistic Decision Boundary ----------
% Visualizes the two-feature space (first two significant %-change vars)
% with LC (blue) and LF (red) patient points. If exactly 2 features are
% significant, overlays the linear decision boundary from the logistic
% regression model: intercept + β1*x + β2*y = 0.
% Only applicable when at least 2 significant variables exist
if n_sig >= 2
    figure('Name', ['2D Feature Space ' fx_label ' — ' dtype_label], 'Position', [100, 100, 800, 600]);
    hold on;

    % Use the first two significant variables for the 2D plot
    x_val = sig_data_selected{1}(valid_pts, target_fx);
    y_val = sig_data_selected{2}(valid_pts, target_fx);
    group = lf_group;

    scatter(x_val(group==0), y_val(group==0), 80, 'b', 'filled', 'MarkerEdgeColor', 'k');
    scatter(x_val(group==1), y_val(group==1), 80, 'r', 'filled', 'MarkerEdgeColor', 'k');

    % Draw Decision Boundary (from logistic regression model)
    if n_sig == 2
        mdl = fitglm([x_val, y_val], group, 'Distribution', 'binomial', 'Options', statset('MaxIter', 1000000));
        coefs = mdl.Coefficients.Estimate;
        x_range = linspace(min(x_val), max(x_val), 100);
        y_boundary = -(coefs(1) + coefs(2)*x_range) / coefs(3);
        plot(x_range, y_boundary, 'k--', 'LineWidth', 2);
    end

    if sig_is_abs(1), xl = ['Abs ' sig_names{1} ' at ' fx_label]; else, xl = ['\Delta ' sig_names{1} ' (%) at ' fx_label]; end
    if sig_is_abs(2), yl = ['Abs ' sig_names{2} ' at ' fx_label]; else, yl = ['\Delta ' sig_names{2} ' (%) at ' fx_label]; end
    xlabel(xl, 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(yl, 'FontSize', 12, 'FontWeight', 'bold');
    title(['Biomarker Interaction: Separation of LC vs LF (' fx_label ', ' dtype_label ')'], 'FontSize', 14);
    if n_sig == 2
        legend({'Local Control', 'Local Failure', 'Logistic Decision Boundary'}, 'Location', 'NorthWest');
    else
        legend({'Local Control', 'Local Failure'}, 'Location', 'NorthWest');
    end
    grid on; box on;
    xline(0, 'k-', 'Alpha', 0.2); yline(0, 'k-', 'Alpha', 0.2);
    saveas(gcf, fullfile(output_folder, ['2D_Feature_Space_' fx_label '_' dtype_label '.png']));
    close(gcf);
else
    fprintf('Skipping 2D scatter plot: requires at least 2 significant variables.\n');
end

%% ---------- Kaplan-Meier Survival Analysis ----------
% Uses nested Leave-One-Out Cross-Validation (LOOCV) to prevent data leakage. 
% For each held-out patient, a multivariable linear risk score is computed 
% using the non-zero coefficients from lassoglm on the N-1 fold.
% The median risk score of the training fold determines the classification threshold.

% Build the full candidate feature matrix for all valid patients
km_X_all = [ADC_abs(valid_pts, target_fx), D_abs(valid_pts, target_fx), ...
            f_abs(valid_pts, target_fx),   Dstar_abs(valid_pts, target_fx), ...
            ADC_pct(valid_pts, target_fx), D_pct(valid_pts, target_fx), ...
            f_pct(valid_pts, target_fx),   Dstar_pct(valid_pts, target_fx)];
km_y_all = lf_group;
n_km = length(km_y_all);

% Build time-to-event vector
times = m_total_time;
times(m_lf==0) = m_total_follow_up_time(m_lf==0);
events = m_lf;  % Event indicator: 1 = failure, 0 = censored

times = times(valid_pts);
events = events(valid_pts);

% LOOCV
is_high_risk = false(n_km, 1);
risk_scores_all = nan(n_km, 1); 

for loo_i = 1:n_km
    train_mask = true(n_km, 1);
    train_mask(loo_i) = false;
    X_tr = km_X_all(train_mask, :);
    y_tr = km_y_all(train_mask);

    valid_tr = any(~isnan(X_tr), 2) & ~isnan(y_tr);
    X_tr = X_tr(valid_tr, :);
    y_tr = y_tr(valid_tr);
    
    % Imputation and Pre-filtering per train fold using KNN
    X_test_imp = km_X_all(loo_i, :);
    [X_tr_imp, X_test_imp] = knn_impute_train_test(X_tr, X_test_imp, 5);
    
    R_fold = corrcoef(X_tr_imp);
    drop_fold = false(1, size(X_tr_imp, 2));
    for fi = 1:size(X_tr_imp, 2)
        if drop_fold(fi), continue; end
        for fj = fi+1:size(X_tr_imp, 2)
            if abs(R_fold(fi, fj)) > 0.8
                drop_fold(fj) = true;
            end
        end
    end
    keep_fold = find(~drop_fold);
    X_tr_kept = X_tr_imp(:, keep_fold);
    X_te_kept = X_test_imp(:, keep_fold);

    warning('off', 'all');
    try
        [B_loo, FI_loo] = lassoglm(X_tr_kept, y_tr, 'binomial', ...
            'Alpha', 0.5, 'CV', 5, 'NumLambda', 25, 'Standardize', true, 'MaxIter', 1000000);
        coefs_loo = B_loo(:, FI_loo.IndexMinDeviance);
        intercept_loo = FI_loo.Intercept(FI_loo.IndexMinDeviance);
    catch
        coefs_loo = zeros(size(X_tr_kept, 2), 1);
        intercept_loo = 0;
    end
    warning('on', 'all');

    % Compute multivariable linear risk score on train fold
    risk_train = X_tr_kept * coefs_loo + intercept_loo;
    cutoff_loo = median(risk_train);

    % Compute risk score for held-out patient
    risk_test = X_te_kept * coefs_loo + intercept_loo;

    % High risk = higher probability of LF (score higher than median)
    is_high_risk(loo_i) = risk_test > cutoff_loo;
    risk_scores_all(loo_i) = risk_test;
end

% Compute and plot Kaplan-Meier survival curves for each risk group
figure('Name', ['Kaplan-Meier ' fx_label ' — ' dtype_label], 'Position', [100, 100, 700, 600]);
[f, x, flow, fup] = ecdf(times, 'Censoring', ~events, 'Function', 'survivor', 'Alpha', 0.05);
if sum(is_high_risk) > 0 && sum(~is_high_risk) > 0
    [f1, x1] = ecdf(times(is_high_risk), 'Censoring', ~events(is_high_risk), 'Function', 'survivor');
    [f2, x2] = ecdf(times(~is_high_risk), 'Censoring', ~events(~is_high_risk), 'Function', 'survivor');
    stairs(x1, f1, 'r-', 'LineWidth', 2.5); hold on;
    stairs(x2, f2, 'b-', 'LineWidth', 2.5);
    legend({'High Risk', 'Low Risk'}, 'Location', 'SouthWest');
else
    fprintf('Warning: Risk groups are degenerate. Skipping KM stratification.\n');
end
xlabel('Time to Local Failure (Days)', 'FontSize', 12);
ylabel('Local Control Probability', 'FontSize', 12);
title(['Kaplan-Meier (LOOCV): Stratified by Multivariable Risk Score (' fx_label ') (' dtype_label ')'], 'FontSize', 14);
grid on; axis([0 max(times)+50 0 1.05]);
saveas(gcf, fullfile(output_folder, ['Kaplan_Meier_' fx_label '_' dtype_label '.png']));
close(gcf);

fprintf('Kaplan-Meier (LOOCV) generated using multivariable risk score. Unbiased out-of-fold risk assignment.\n');

%% ---------- Correlation Matrix of Significant Features ----------
% Computes Pearson correlation coefficients between all pairs of
% LASSO-selected percent-change variables. Displayed as a jet-coloured
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
    title(['Feature Correlation (' fx_label ') (' dtype_label '). R(' feat_names{1} ', ' feat_names{2} ') = ' num2str(R(1,2), '%.2f')]);
else
    title(['Feature Correlation (' fx_label ') (' dtype_label ')']);
end
saveas(gcf, fullfile(output_folder, ['Correlation_Matrix_' fx_label '_' dtype_label '.png']));
close(gcf);

%% ---------- Representative Parametric ADC Maps ----------
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
        
        % Compute full ADC volume (log-linear OLS via fit_adc_mono) and extract slice
        adc_vol = fit_adc_mono(dwi_img, bvals);
        slice_adc = squeeze(adc_vol(:,:,z_slice));
        
        % Clamp ADC to physiological range [0, 3×10⁻³ mm²/s]
        slice_adc(slice_adc < 0) = 0; 
        slice_adc(slice_adc > 3e-3) = 3e-3;
        
        % Store Fx1 and target-Fx slices separately for later comparison
        if t == 1
            slice_fx1 = slice_adc; mask_slice_fx1 = slice_gtv;
        else
            slice_fx3 = slice_adc; mask_slice_fx3 = slice_gtv;
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
        imagesc(crop_fx1, [0 2.5e-3]); axis image; axis off;
        title([titles{p} ' - Fx1'], 'FontSize', 10, 'Interpreter', 'none');
        if p==1, ylabel('ADC Map'); end
        
        subplot(2, 3, 2 + row_offset);
        imagesc(crop_fx3, [0 2.5e-3]); axis image; axis off;
        title(fx_label, 'FontSize', 10);
        
        subplot(2, 3, 3 + row_offset);
        imagesc(crop_diff, [-1e-3 1e-3]); axis image; axis off;
        title(['\Delta (' fx_label ' - Fx1)'], 'FontSize', 10);
        
        % Overlay Fx1 GTV contour on difference map (white outline)
        hold on; contour(mask_slice_fx1(r_min:r_max, c_min:c_max), [0.5 0.5], 'w', 'LineWidth', 1); hold off;
    end
end

% Add shared colorbars: top = absolute ADC, bottom = difference (ΔADC)
c1 = colorbar('Position', [0.92 0.55 0.02 0.35]); ylabel(c1, 'ADC');
c2 = colorbar('Position', [0.92 0.11 0.02 0.35]); ylabel(c2, '\Delta ADC');
sgtitle(['Representative Longitudinal ADC Response for ' curr_sig_disp ' (' fx_label ', ' dtype_label ')'], 'FontSize', 14, 'FontWeight', 'bold');
safe_name = strrep(curr_sig_file, '*', 'star');
saveas(gcf, fullfile(output_folder, ['Representative_ADC_' safe_name '_' fx_label '_' dtype_label '.png']));
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

subplot(1, 3, 2);
kurt_fx1 = adc_kurt(valid_pts, 1, 1);            
kurt_fx3 = adc_kurt(valid_pts, target_fx, 1);     
kurt_delta = kurt_fx3 - kurt_fx1;                  

boxplot(kurt_delta, lf_group, 'Labels', {'LC (0)', 'LF (1)'});
ylabel(['\Delta ADC Kurtosis (' fx_label ')']);
title('Heterogeneity: Texture Change', 'FontSize', 12, 'FontWeight', 'bold');
p_kurt = ranksum(kurt_delta(lf_group==0), kurt_delta(lf_group==1));

y_lim = ylim;
text(1.5, y_lim(2)*0.9, sprintf('p = %.3f', p_kurt), 'HorizontalAlignment', 'center', 'FontSize', 11);
grid on;

subplot(1, 3, 3);
hold on;
wcv_est = 2.8; 
cor_est = 1.96 * sqrt(2) * wcv_est;

x_scatter = ones(size(lf_group));
x_scatter(lf_group==1) = 2;
x_scatter = x_scatter + (rand(size(x_scatter))-0.5)*0.2;
scatter(x_scatter, curr_sig_pct_full(valid_pts, target_fx), 50, 'filled', 'MarkerEdgeColor', 'k');

yfill = [-cor_est cor_est cor_est -cor_est];
xfill = [0.5 0.5 2.5 2.5];
fill(xfill, yfill, [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

yline(0, 'k-');
yline(cor_est, 'k--', 'CoR (+7.8%)');
yline(-cor_est, 'k--', 'CoR (-7.8%)');
xticks([1 2]); xticklabels({'LC', 'LF'});
if sig_is_abs(vi), ylbl = curr_sig_disp; else, ylbl = [curr_sig_disp ' (%)']; end
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

% --- PART B: HAZARD RATIOS (Cox Regression) ---
fprintf('\n--- SURVIVAL ANALYSIS (Prognostic Value) ---\n');
for vi = 1:n_sig
    curr_sig_pct_full = sig_data_selected{vi}(valid_pts, target_fx);
    curr_sig_name = sig_names{vi};
    if sig_is_abs(vi), var_desc = 'Absolute'; else, var_desc = 'Continuous Delta'; end

    valid_cox = ~isnan(curr_sig_pct_full) & ~isnan(times) & ~isnan(events);
    cox_times = times(valid_cox);
    cox_events = events(valid_cox);
    cox_biomarker = curr_sig_pct_full(valid_cox);

    if sum(valid_cox) > 5
        [b, logl, H, stats] = coxphfit(cox_biomarker, cox_times, 'Censoring', ~cox_events, 'Options', statset('MaxIter', 1000000));
        fprintf('Biomarker: %s %s at %s\n', var_desc, curr_sig_name, fx_label);
        if sig_is_abs(vi), hr_unit = '1 unit'; else, hr_unit = '1%%'; end
        fprintf('Hazard Ratio (HR per %s change): %.2f\n', hr_unit, exp(b));
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

% --- PART C: SUB-VOLUME STATISTICS ---
% How big were the "Resistant" (Low Diffusion) sub-volumes?
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

% --- C. MULTIVARIABLE COX REGRESSION (Robust Version) ---
try
    raw_lf = m_lf(valid_pts);
    raw_time_fail = m_total_time(valid_pts);
    raw_time_cens = m_total_follow_up_time(valid_pts);
    
    surv_time = raw_time_fail;
    surv_time(raw_lf == 0) = raw_time_cens(raw_lf == 0);
    surv_event = raw_lf;
    
    raw_baseline_vol = m_gtv_vol(valid_pts, 1);
    vol_vec = (raw_baseline_vol - mean(raw_baseline_vol, 'omitnan')) ./ std(raw_baseline_vol, 'omitnan');
    
    fprintf('\nMultivariable Cox Regression:\n');
    for vi = 1:n_sig
        curr_sig_pct_full = sig_data_selected{vi}(valid_pts, target_fx);
        curr_sig_name = sig_names{vi};
        if sig_is_abs(vi), var_desc = 'Absolute'; else, var_desc = 'Continuous Delta'; end
        
        final_mask = ~isnan(surv_time) & ~isnan(surv_event) & ~isnan(curr_sig_pct_full) & ~isnan(vol_vec);
        
        if sum(final_mask) > 5
            y_time = surv_time(final_mask);
            y_event = surv_event(final_mask);
            x_biomarker = curr_sig_pct_full(final_mask); 
            x_vol = vol_vec(final_mask);
            
            [b, logl, H, stats] = coxphfit([x_biomarker, x_vol], y_time, 'Censoring', ~y_event, 'Options', statset('MaxIter', 1000000));
            
            if sig_is_abs(vi), hr_unit = '1 unit'; else, hr_unit = '1%%'; end
            fprintf('  Variable 1: %s %s\n', var_desc, curr_sig_name);
            fprintf('     HR per %s change: %.4f (p=%.4f)\n', hr_unit, exp(b(1)), stats.p(1));
            fprintf('  Variable 2: Baseline Tumor Volume (Z-score)\n');
            fprintf('     HR: %.4f (p=%.4f)\n', exp(b(2)), stats.p(2));
            
            if stats.p(1) < 0.05
                fprintf('  Conclusion: ROBUST. Predictive independent of volume.\n');
            else
                fprintf('  Conclusion: CONFOUNDED. Significance lost after adjustment.\n');
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
                    fprintf('     Variable 1 (%s): R = %.3f, p-value = %.4f\n', curr_sig_name, r_ph_bio, p_ph_bio);
                    fprintf('     Variable 2 (Volume): R = %.3f, p-value = %.4f\n', r_ph_vol, p_ph_vol);
                else
                    fprintf('  Proportional Hazards Test: Residual count (%d) mismatch with event count (%d). Skipping.\n', size(stats.sschres, 1), length(log_event_times));
                end
                
                if p_ph_bio < 0.05
                    warning('Proportional hazards assumption may be violated for multivariable %s (p = %.3f)', curr_sig_name, p_ph_bio);
                end
                if p_ph_vol < 0.05
                    warning('Proportional hazards assumption may be violated for multivariable Volume (p = %.3f)', p_ph_vol);
                end
                
                % Plot for Biomarker vs log(event time)
                figure('Name', ['MV Schoenfeld: ' var_desc ' ' curr_sig_name ' at ' fx_label], 'Position', [100, 100, 600, 500]);
                scatter(log_event_times, sch_res_bio, 60, 'b', 'filled', 'MarkerEdgeColor', 'k'); hold on;
                p_fit = polyfit(log_event_times, sch_res_bio, 1);
                x_trend = linspace(min(log_event_times), max(log_event_times), 100);
                y_trend = polyval(p_fit, x_trend);
                plot(x_trend, y_trend, 'r-', 'LineWidth', 2);
                yline(0, 'k--', 'LineWidth', 1.5);
                xlabel('log(Time to Event) [log(Days)]', 'FontSize', 12, 'FontWeight', 'bold');
                ylabel('Scaled Schoenfeld Residuals', 'FontSize', 12, 'FontWeight', 'bold');
                title({['MV Proportional Hazards Check: ' curr_sig_name], ['Spearman Corr (log time): r = ' num2str(r_ph_bio, '%.3f') ', p = ' num2str(p_ph_bio, '%.3f')]}, 'FontSize', 14);
                legend('Residuals', 'Trend', 'Zero Line', 'Location', 'Best');
                grid on;
                safe_name = strrep(curr_sig_name, '*', 'star');
                if sig_is_abs(vi), file_prefix = 'Abs_'; else, file_prefix = 'Delta_'; end
                saveas(gcf, fullfile(output_folder, ['MV_Schoenfeld_' file_prefix safe_name '_' fx_label '_' dtype_label '.png']));
                close(gcf);
            else
                fprintf('  Proportional Hazards Test: Not enough events to test assumption.\n');
            end
            fprintf('\n');
        else
            fprintf('  Variable 1: %s %s: Not enough complete data.\n\n', var_desc, curr_sig_name);
        end
    end

catch ME
    fprintf('Error in Cox Regression: %s\n', ME.message);
end
end % for target_fx

%% ---------- Q-Value (FDR) Analysis — Full Metric Sweep ----------
% Re-gathers ALL raw p-values from the 4 metric sets × all timepoints,
% applies Benjamini-Hochberg to compute q-values, and generates a
% scatter plot of raw P vs adjusted Q to visualize inflation.
% Runs once globally (after all per-timepoint loops) to avoid redundant
% re-computation and duplication of figures.

% 1. Gather all raw p-values from the previous analysis
raw_p_vals = [];
metric_names = {};

for s = 1:length(metric_sets)
    current_metrics = metric_sets{s};
    current_names = set_names{s};
    for m = 1:length(current_metrics)
        metric_data = current_metrics{m};
        for tp = [2, 3]
            if tp > size(metric_data, 2)
                continue;
            end
            % Extract data
            y_raw = metric_data(valid_pts, tp);
            g = lf_group(~isnan(y_raw));
            y = y_raw(~isnan(y_raw));

            % Run ranksum if sufficient data exists
            if length(unique(g)) > 1 && length(y) > 2
                p = ranksum(y(g==0), y(g==1));
                raw_p_vals(end+1, 1) = p;
                metric_names{end+1, 1} = [current_names{m} ' (' time_labels{tp} ')'];
            end
        end
    end
end

% 2. Sort P-values to calculate Q-values (Benjamini-Hochberg)
[p_sorted, sort_idx] = sort(raw_p_vals);
m = length(p_sorted); % Total number of tests
q_vals = zeros(size(p_sorted));

% Calculate Q-values
% Formula: q_val = p_val * (Total_Tests / Rank)
for k = 1:m
    q_vals(k) = p_sorted(k) * (m / k);
end

% 3. Enforce monotonicity (a q-value cannot decrease as p-value increases)
for k = m-1:-1:1
    q_vals(k) = min(q_vals(k), q_vals(k+1));
end
% Cap q-values at 1.0
q_vals(q_vals > 1) = 1;

% 4. Create a Results Table
% Map back to original order to match names
q_vals_unsorted = zeros(size(raw_p_vals));
q_vals_unsorted(sort_idx) = q_vals;

results_table = table(metric_names, raw_p_vals, q_vals_unsorted, ...
    'VariableNames', {'Metric', 'P_Value', 'Q_Value'});

% Sort by P-value for readability
results_table = sortrows(results_table, 'P_Value');

% 5. Display the Top 10 Results
fprintf('\n--- Top 10 Results by Q-Value (FDR) ---\n');
num_to_disp = min(10, height(results_table));
disp(results_table(1:num_to_disp, :));

% 6. Plot P-values vs Q-values
figure('Name', ['P-value vs Q-value Analysis — ' dtype_label], 'Position', [300, 300, 600, 500]);
scatter(results_table.P_Value, results_table.Q_Value, 50, 'filled');
xlabel('Raw P-Value');
ylabel('Q-Value (FDR)');
title(['Significance Analysis: P vs Q (Global, ' dtype_label ')']);
yline(0.05, 'r--', 'Q = 0.05 (Significance Threshold)');
grid on;
saveas(gcf, fullfile(output_folder, ['P_vs_Q_Global_' dtype_label '.png']));
close(gcf);

%% ---------- Holm-Bonferroni Correction (FWER Control) ----------
% Controls the family-wise error rate (FWER) using the step-down Holm
% procedure. More powerful than classical Bonferroni but still strict.
% Algorithm: sort p-values ascending; compare the k-th smallest to
% α/(m+1−k). Reject until the first failure, then stop.
% Visualization: sorted p-values vs step-down thresholds.
% Runs once globally (after all per-timepoint loops).

% 1. Gather all raw p-values
raw_p_vals = [];
metric_labels = {};
for s = 1:length(metric_sets)
    current_metrics = metric_sets{s};
    current_names = set_names{s};
    for m = 1:length(current_metrics)
        metric_data = current_metrics{m};
        for tp = [2, 3]
            if tp > size(metric_data, 2)
                continue;
            end
            y_raw = metric_data(valid_pts, tp);
            g = lf_group(~isnan(y_raw));
            y = y_raw(~isnan(y_raw));
            if length(unique(g)) > 1 && length(y) > 2
                p = ranksum(y(g==0), y(g==1));
                raw_p_vals(end+1, 1) = p;
                metric_labels{end+1, 1} = [current_names{m} ' (' time_labels{tp} ')'];
            end
        end
    end
end

% 2. Sort p-values from smallest to largest
[p_sorted, sort_idx] = sort(raw_p_vals);
m = length(p_sorted);       % Total number of tests
alpha = 0.05;               % Desired FWER

% 3. Calculate Holm-Bonferroni Thresholds
% The critical value adjusts for the number of REMAINING tests
% Threshold_k = alpha / (m + 1 - k)
holm_thresholds = alpha ./ (m + 1 - (1:m)');

% 4. Determine Significance
% A test is significant if its p-value < holm_threshold
% AND all preceding tests were also significant.
is_sig_holm = false(size(p_sorted));
for k = 1:m
    if p_sorted(k) < holm_thresholds(k)
        is_sig_holm(k) = true;
    else
        % Stop at the first failure. All subsequent tests are non-significant.
        break;
    end
end

% 5. Map back to original labels and display
results_holm = table(metric_labels(sort_idx), p_sorted, holm_thresholds, is_sig_holm, ...
    'VariableNames', {'Metric', 'Raw_P', 'Holm_Threshold', 'Significant'});

fprintf('\n--- Holm-Bonferroni FWER Analysis ---\n');
fprintf('Family Size (m): %d tests\n', m);
fprintf('Uncorrected FWER Probability: %.2f%%\n', (1 - (1-alpha)^m)*100);
fprintf('Top 5 Results:\n');
num_to_disp_holm = min(5, height(results_holm));
disp(results_holm(1:num_to_disp_holm, :));

if any(is_sig_holm)
    fprintf('SUCCESS: One or more results survived FWER control.\n');
else
    fprintf('RESULT: No findings survived FWER control.\n');
    fprintf('The "Best" P-value (%.4f) did not beat the first hurdle (%.4f).\n', ...
        p_sorted(1), holm_thresholds(1));
end

% 6. Visualization
figure('Name', ['Holm-Bonferroni Step-Down — ' dtype_label], 'Position', [400, 400, 700, 500]);
plot(1:m, p_sorted, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
plot(1:m, holm_thresholds, 'r--', 'LineWidth', 2);
xlabel('Rank (k)');
ylabel('P-Value');
title(['Holm-Bonferroni Procedure (Global, ' dtype_label ')']);
legend({'Your Sorted P-values', 'Holm Step-Down Threshold'}, 'Location', 'NorthWest');
grid on;
saveas(gcf, fullfile(output_folder, ['Holm_Bonferroni_Global_' dtype_label '.png']));
close(gcf);

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

warning('off', 'all');
try
    glme = fitglme(glme_table_clean, ...
        'LF ~ 1 + ADC_z + D_z + f_z + Dstar_z + Timepoint + (1|PatientID)', ...
        'Distribution', 'Binomial', 'Link', 'logit', 'OptimizerOptions', statset('MaxIter', 1000000));
    disp(glme);
catch ME
    fprintf('GLME model failed to converge: %s\n', ME.message);
end
warning('on', 'all');

fprintf('\n=== Completed DWI Type %d: %s ===\n', dtype, dtype_label);
end % for dtype
diary off