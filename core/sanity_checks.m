function [is_valid, validation_msg, data_vectors_gtvp, data_vectors_gtvn] = sanity_checks(data_vectors_gtvp, data_vectors_gtvn, summary_metrics)
% SANITY_CHECKS — "Understand the Data" (Validation & Sanity Checking)
is_valid = true;
validation_msg = 'Passed';

id_list = summary_metrics.id_list;
mrn_list = summary_metrics.mrn_list;
adc_mean = summary_metrics.adc_mean;
d_mean = summary_metrics.d_mean;
f_mean = summary_metrics.f_mean;
dstar_mean = summary_metrics.dstar_mean;
d95_gtvp = summary_metrics.d95_gtvp;
% Run this script AFTER loading the saved workspace from
% load_dwi_data_forAvery.m (i.e., after the "reload saved result" and
% "pull in maps and compute longitudinal metrics" sections have executed).
%
% Required workspace variables:
%   data_vectors_gtvp  – struct array [nPat x nTp x nDwiType] holding
%                        per-voxel ADC, D, f, D*, and dose vectors
%   data_vectors_gtvn  – same structure for the GTV-normal region
%   id_list            – cell array of patient folder-name identifiers
%   mrn_list           – cell array of medical-record-number strings
%   adc_mean, d_mean, f_mean, dstar_mean
%                      – [nPat x nTp x nDwiType] summary arrays of mean
%                        diffusion biomarkers per patient/timepoint/type
%   lf                 – binary outcome vector (0 = LC, 1 = LF)
%   d95_gtvp, v50gy_gtvp
%                      – [nPat x nTp] RT-Dose coverage metrics sampled
%                        on the GTV mask
%   dwi_locations, rtdose_locations, gtv_locations
%                      – cell arrays of file paths (used for alignment check)
%
% This script performs three categories of sanity checks:
%   1. Verify Convergence — flags Inf, NaN, and negative diffusion values
%      in the per-voxel model-fit outputs to detect fitting failures.
%   2. Identify Missingness & Outliers — summarises NaN patterns across
%      patients/fractions and flags extreme biomarker values that deviate
%      more than 3 × IQR from the cohort median.
%   3. Data Alignment — confirms that the RT Dose and DWI volumes are
%      spatially registered by comparing voxel counts inside the GTV mask.

% Create output directory for saved figures and text logs
output_folder = fullfile(pwd, 'saved_figures');
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

% Write all console output to a diary file so results are reproducible.
% Delete any previous diary to avoid appending stale data.
diary_file = fullfile(output_folder, 'sanity_checks_output.txt');
if exist(diary_file, 'file'), delete(diary_file); end
diary(diary_file);

fprintf('\n======================================================\n');
fprintf('  SANITY CHECKS — Understand the Data\n');
fprintf('======================================================\n');

% Cohort dimensions
nPat = length(id_list);            % total number of patients
nTp  = size(adc_mean, 2);         % number of timepoints (typically 6: Fx1–Fx5 + Post)
fx_labels = {'Fx1','Fx2','Fx3','Fx4','Fx5','Post'};

%% -----------------------------------------------------------------------
%  1. VERIFY CONVERGENCE
%  Iterate over every patient × timepoint and inspect the per-voxel
%  vectors produced by the IVIM / ADC model fits.  Non-physical values
%  (Inf, NaN, or negative diffusion coefficients) indicate fitting
%  failures that could corrupt downstream summary statistics.
% -----------------------------------------------------------------------
fprintf('\n--- 1. Verify Convergence ---\n');

conv_issues = 0;   % running count of flagged patient-timepoint-metric tuples

for j = 1:nPat
    for k = 1:nTp
        % Pull the per-voxel result struct for the primary (rpi=1)
        % acquisition of this patient and timepoint.
        if j > size(data_vectors_gtvp, 1) || k > size(data_vectors_gtvp, 2)
            continue;
        end
        s = data_vectors_gtvp(j, k, 1);

        % Bundle the four voxel-level biomarker vectors into a struct
        % so we can loop over them programmatically.
        vectors = struct( ...
            'ADC',   s.adc_vector, ...
            'D',     s.d_vector,   ...
            'f',     s.f_vector,   ...
            'Dstar', s.dstar_vector);

        fields = fieldnames(vectors);
        for fi = 1:numel(fields)
            v = vectors.(fields{fi});
            if isempty(v), continue; end   % no data at this timepoint

            % Count non-physical voxel values
            n_inf    = sum(isinf(v));
            n_nan    = sum(isnan(v));
            n_neg    = sum(v < 0);
            n_total  = numel(v);

            if n_inf > 0 || n_nan > 0 || n_neg > 0
                fprintf('  Patient %s  %s  %s : ', ...
                    id_list{j}, fx_labels{min(k,numel(fx_labels))}, fields{fi});
                if n_inf > 0
                    fprintf('Inf=%d/%d (%.1f%%)  ', n_inf, n_total, 100*n_inf/n_total);
                end
                if n_nan > 0
                    fprintf('NaN=%d/%d (%.1f%%)  ', n_nan, n_total, 100*n_nan/n_total);
                end
                if n_neg > 0
                    fprintf('Neg=%d/%d (%.1f%%)  ', n_neg, n_total, 100*n_neg/n_total);
                end
                fprintf('\n');
                conv_issues = conv_issues + 1;
            end
        end
    end
end

if conv_issues == 0
    fprintf('  All voxel-level fit values are finite, non-NaN, and non-negative.\n');
else
    fprintf('  Total convergence flags raised: %d\n', conv_issues);
end

%% -----------------------------------------------------------------------
%  2. IDENTIFY MISSINGNESS & OUTLIERS
%  2a. Missingness — For each summary metric and fraction, count how many
%      patients have NaN values (i.e., missing or failed acquisitions).
%  2b. Outlier detection — Flag individual patient-fraction biomarker
%      values that lie more than 3 × IQR beyond the cohort median.
%      Such extreme values may indicate motion artefacts, mis-segmentation,
%      or fitting failures that were not caught in Step 1.
% -----------------------------------------------------------------------
fprintf('\n--- 2. Identify Missingness & Outliers ---\n');

% --- 2a. Fraction-level missingness summary ---
fprintf('\n  2a. Fraction-level missingness (NaN in summary arrays):\n');
fprintf('  %-12s  %s\n', 'Metric', strjoin(fx_labels,'   '));

% Evaluate missingness for each of the four summary biomarker arrays
% (using Standard DWI, dtype index = 1).
summary_metrics = {adc_mean(:,:,1), d_mean(:,:,1), f_mean(:,:,1), dstar_mean(:,:,1)};
summary_names   = {'ADC_mean', 'D_mean', 'f_mean', 'Dstar_mean'};

for mi = 1:numel(summary_metrics)
    dat = summary_metrics{mi};
    nCols = min(size(dat, 2), numel(fx_labels));
    counts = sum(isnan(dat(:, 1:nCols)), 1);
    parts = cell(1, nCols);
    for ci = 1:nCols
        parts{ci} = sprintf('%3d/%3d', counts(ci), nPat);
    end
    fprintf('  %-12s  %s\n', summary_names{mi}, strjoin(parts, '  '));
end

% RT Dose missingness — D95 is representative; if it is NaN the
% dose volume was not available for that fraction.
fprintf('\n  RT Dose NaN check (d95_gtvp):\n');
nDoseCols = size(d95_gtvp, 2);
for k = 1:min(nDoseCols, 5)
    n_nan_dose = sum(isnan(d95_gtvp(:, k)));
    fprintf('    %s : %d / %d patients missing dose data\n', ...
        fx_labels{k}, n_nan_dose, nPat);
end

% --- 2b. Outlier detection (values > 3 IQR from median) ---
% The 3 × IQR fence is more conservative than the standard 1.5 × IQR
% box-plot whiskers, so only truly extreme observations are flagged.
fprintf('\n  2b. Outlier detection (>3 IQR from cohort median):\n');

outlier_count = 0;
for mi = 1:numel(summary_metrics)
    dat = summary_metrics{mi};
    nCols = min(size(dat, 2), numel(fx_labels));
    for k = 1:nCols
        col = dat(:, k);
        col_clean = col(~isnan(col));
        if numel(col_clean) < 3, continue; end

        med_val = median(col_clean);
        iqr_val = iqr(col_clean);
        if iqr_val == 0, continue; end

        lower_fence = med_val - 3 * iqr_val;
        upper_fence = med_val + 3 * iqr_val;

        for j = 1:nPat
            val = dat(j, k);
            if isnan(val), continue; end
            if val < lower_fence || val > upper_fence
                fprintf('    OUTLIER: Patient %s  %s  %s = %.4g  (median=%.4g, 3*IQR fence=[%.4g, %.4g])\n', ...
                    id_list{j}, fx_labels{k}, summary_names{mi}, val, med_val, lower_fence, upper_fence);
                outlier_count = outlier_count + 1;
            end
        end
    end
end

if outlier_count == 0
    fprintf('    No extreme outliers detected.\n');
else
    fprintf('    Total outlier flags: %d  — review the original MRI for these cases.\n', outlier_count);
end

%% -----------------------------------------------------------------------
%  3. DATA ALIGNMENT
%  For each patient and fraction, compare the number of voxels in the
%  dose vector with the number in the ADC vector (both sampled on the
%  same GTV mask).  A mismatch indicates that the RT Dose and DWI
%  volumes were not resampled to the same grid — i.e., a registration
%  failure.  Additionally, if more than 10 % of the in-mask dose voxels
%  are NaN, the dose map may not fully cover the GTV, suggesting a
%  spatial mis-alignment between the treatment-planning CT and the MRI.
% -----------------------------------------------------------------------
fprintf('\n--- 3. Data Alignment (RT Dose ↔ DWI) ---\n');

align_issues = 0;   % running count of dimensional or NaN-fraction flags

for j = 1:nPat
    nFx = min(size(data_vectors_gtvp, 2), 5); % dose only for Fx1-5
    for k = 1:nFx
        if j > size(data_vectors_gtvp, 1) || k > size(data_vectors_gtvp, 2)
            continue;
        end
        s = data_vectors_gtvp(j, k, 1);

        dose_vec = s.dose_vector;
        adc_vec  = s.adc_vector;

        if isempty(dose_vec) || isempty(adc_vec), continue; end

        if numel(dose_vec) ~= numel(adc_vec)
            fprintf('  MISMATCH: Patient %s  %s : dose voxels=%d, ADC voxels=%d\n', ...
                id_list{j}, fx_labels{k}, numel(dose_vec), numel(adc_vec));
            align_issues = align_issues + 1;
        end

        % Check for NaN-only dose regions inside the mask
        frac_nan = sum(isnan(dose_vec)) / numel(dose_vec);
        if frac_nan > 0.1
            fprintf('  WARNING: Patient %s  %s : %.1f%% of in-mask dose voxels are NaN (possible mis-registration)\n', ...
                id_list{j}, fx_labels{k}, frac_nan * 100);
            align_issues = align_issues + 1;
        end
    end
end

if align_issues == 0
    fprintf('  All dose/DWI pairs are dimensionally consistent with <10%% NaN.\n');
else
    fprintf('  Total alignment flags: %d\n', align_issues);
end

fprintf('\n======================================================\n');
fprintf('  Sanity checks complete.\n');
fprintf('======================================================\n');
diary off

if align_issues == 0
    is_valid = true;
    validation_msg = sprintf('Passed all alignment checks. %d convergence warnings.', conv_issues);
else
    % We treat alignment mismatch as invalidating the run
    is_valid = false;
    validation_msg = sprintf('Failed %d alignment checks.', align_issues);
end
end