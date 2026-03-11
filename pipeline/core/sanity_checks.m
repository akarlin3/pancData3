function [is_valid, validation_msg, data_vectors_gtvp, data_vectors_gtvn] = sanity_checks(data_vectors_gtvp, data_vectors_gtvn, summary_metrics, config_struct)
% SANITY_CHECKS — "Understand the Data" (Validation & Sanity Checking)
% Author: Avery Karlin
%
% Inputs:
%   data_vectors_gtvp - Struct array holding primary GTV parameters
%   data_vectors_gtvn - Struct array holding nodal GTV parameters
%   summary_metrics   - Struct containing patient ID and metric statistics
%   config_struct     - Configuration struct defining paths and thresholds
%
% Outputs:
%   is_valid          - Boolean indicating if the data is safe to visualize/analyze
%   validation_msg    - Log message explaining issues
%   data_vectors_gtvp - Primary GTV parameter map (passed through un-modified or corrected)
%   data_vectors_gtvn - Nodal GTV parameter map (passed through un-modified or corrected)
%
% ANALYTICAL RATIONALE — DATA QUALITY ASSURANCE
%   Before any visualization or statistical analysis, data must be validated
%   to ensure that model fitting succeeded and that spatial registrations
%   are consistent. In diffusion MRI analysis, common failure modes include:
%
%   1. Model fitting failures: IVIM biexponential fitting is numerically
%      unstable, especially for D* estimation. Failed fits produce Inf, NaN,
%      or negative parameter values that would corrupt summary statistics.
%
%   2. Missing data: Pancreatic DWI acquisitions are frequently incomplete
%      (patient discomfort, scanner errors, motion rejection). Understanding
%      the pattern of missingness (random vs systematic) is essential for
%      choosing appropriate statistical methods (complete-case vs imputation).
%
%   3. Spatial misalignment: RT dose maps and DWI volumes must be spatially
%      co-registered for dose-response analysis. A dimensional mismatch
%      between dose and ADC vectors indicates that the dose resampling or
%      DIR registration failed, making voxel-level correlations invalid.
%
%   4. Excessive NaN: If >50% of voxels across the cohort are NaN, the
%      data is too sparse for reliable statistical analysis. This gate
%      prevents downstream modules from producing misleading results.
%

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
if nargin > 3 && isfield(config_struct, 'output_folder')
    output_folder = config_struct.output_folder;
else
    timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
    output_folder = fullfile(fileparts(mfilename('fullpath')), '..', sprintf('saved_files_%s', timestamp_str));
end
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
fx_labels = [arrayfun(@(x) sprintf('Fx%d', x), 1:(nTp-1), 'UniformOutput', false), {'Post'}];

%% -----------------------------------------------------------------------
fprintf('\n--- SECTION 1: Verify Convergence ---\n');
%  1. VERIFY CONVERGENCE
%  Iterate over every patient × timepoint and inspect the per-voxel
%  vectors produced by the IVIM / ADC model fits.  Non-physical values
%  (Inf, NaN, or negative diffusion coefficients) indicate fitting
%  failures that could corrupt downstream summary statistics.
% -----------------------------------------------------------------------
fprintf('\n--- 1. Verify Convergence ---\n');

conv_issues = 0;   % running count of flagged patient-timepoint-metric tuples

% Iterate over every patient x timepoint and check for non-physical values.
% In the IVIM model, all parameters must be physically non-negative:
%   D >= 0 (diffusion coefficient cannot be negative)
%   f in [0, 1] (perfusion fraction is a volume fraction)
%   D* >= 0 (pseudo-diffusion coefficient cannot be negative)
%   ADC >= 0 (apparent diffusion coefficient cannot be negative)
% Inf values indicate numerical overflow (typically from division by near-zero
% signal). NaN values indicate fitting did not converge or input data was
% invalid. Negative values indicate the fit converged to a non-physical
% solution (e.g., signal increasing with b-value due to motion artifacts).
for j = 1:nPat
    text_progress_bar(j, nPat, 'Checking convergence');
    for k = 1:nTp
        % Pull the per-voxel result struct for the primary (rpi=1)
        % acquisition of this patient and timepoint.
        if j > size(data_vectors_gtvp, 1) || k > size(data_vectors_gtvp, 2)
            continue;
        end
        s = data_vectors_gtvp(j, k, 1);

        % Use configuration target to check the corresponding fields
        if exist('config_struct', 'var') && isfield(config_struct, 'dwi_types_to_run') && isnumeric(config_struct.dwi_types_to_run)
            dtype_list = config_struct.dwi_types_to_run;
        else
            dtype_list = 1;
        end

        for dtype_idx = dtype_list(:)'
            switch dtype_idx
                case 1
                    vectors = struct('ADC', s.adc_vector, 'D', s.d_vector, 'f', s.f_vector, 'Dstar', s.dstar_vector);
                    dtype_name = '';
                case 2
                    vectors = struct('ADC', s.adc_vector_dncnn, 'D', s.d_vector_dncnn, 'f', s.f_vector_dncnn, 'Dstar', s.dstar_vector_dncnn);
                    dtype_name = ' (DnCNN)';
                case 3
                    vectors = struct('ADC', s.adc_vector, 'D', s.d_vector_ivimnet, 'f', s.f_vector_ivimnet, 'Dstar', s.dstar_vector_ivimnet);
                    dtype_name = ' (IVIM-NET)';
                otherwise
                    continue;
            end

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
                    fprintf('  Patient %s  %s  %s%s : ', ...
                        id_list{j}, fx_labels{min(k,numel(fx_labels))}, fields{fi}, dtype_name);
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
end

if conv_issues == 0
    fprintf('  All voxel-level fit values are finite, non-NaN, and non-negative.\n');
else
    fprintf('  Total convergence flags raised: %d\n', conv_issues);
end

%% -----------------------------------------------------------------------
fprintf('\n--- SECTION 2: Identify Missingness & Outliers ---\n');
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

% Evaluate missingness using the active DWI type from config, not hardcoded 1.
% This ensures that when running the dnCNN or IVIMnet pipeline, we check
% missingness in the corresponding parameter arrays rather than always
% looking at Standard pipeline data (which might be fully populated even
% when the DL pipeline has widespread failures).
if exist('config_struct', 'var') && isfield(config_struct, 'dwi_types_to_run') && isnumeric(config_struct.dwi_types_to_run)
    dtype_miss = config_struct.dwi_types_to_run(1);  % scalar: use first type
else
    dtype_miss = 1;
end
% Clamp to available dimension to avoid out-of-bounds indexing
dtype_miss = min(dtype_miss, size(adc_mean, 3));
% Extract the 2D (patients x timepoints) slice for each diffusion
% parameter from the 3D (patients x timepoints x DWI_type) arrays.
summary_arrs  = {adc_mean(:,:,dtype_miss), d_mean(:,:,dtype_miss), f_mean(:,:,dtype_miss), dstar_mean(:,:,dtype_miss)};
summary_names = {'ADC_mean', 'D_mean', 'f_mean', 'Dstar_mean'};

for mi = 1:numel(summary_arrs)
    dat = summary_arrs{mi};
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
% The 3 x IQR fence is deliberately more conservative than the standard
% 1.5 x IQR box-plot whiskers, so only truly extreme observations are
% flagged. In diffusion MRI of pancreatic tumors, moderate outliers are
% expected due to biological heterogeneity (e.g., cystic vs solid tumors,
% varying degrees of necrosis). Only values beyond 3 x IQR are likely
% to represent technical failures (mis-segmentation capturing normal
% bowel, severe motion artifacts, or fitting failures not caught in
% Step 1) rather than biological variation.
fprintf('\n  2b. Outlier detection (>3 IQR from cohort median):\n');

outlier_count = 0;
for mi = 1:numel(summary_arrs)
    dat = summary_arrs{mi};
    nCols = min(size(dat, 2), numel(fx_labels));
    for k = 1:nCols
        col = dat(:, k);
        col_clean = col(~isnan(col));
        if numel(col_clean) < 3, continue; end

        med_val = median(col_clean);
        iqr_val = iqr(col_clean);
        if iqr_val == 0 || isnan(iqr_val), continue; end

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
fprintf('\n--- SECTION 3: Data Alignment (RT Dose <-> DWI) ---\n');
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

dim_mismatches = 0;   % hard errors: dose/ADC vector length mismatch
nan_warnings   = 0;   % soft warnings: high NaN fraction in dose

% Check alignment only for Fx1-Fx5 (treatment fractions with dose maps).
% The post-RT scan (Fx6) does not have an associated RT dose because the
% treatment course is complete; dose-response analysis is only meaningful
% for on-treatment timepoints.
for j = 1:nPat
    text_progress_bar(j, nPat, 'Checking alignment');
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
            dim_mismatches = dim_mismatches + 1;
        end

        % Check for NaN-only dose regions inside the mask
        frac_nan = sum(isnan(dose_vec)) / numel(dose_vec);
        if frac_nan > 0.1
            nan_warnings = nan_warnings + 1;
        end
    end
end

if dim_mismatches == 0 && nan_warnings == 0
    fprintf('  All dose/DWI pairs are dimensionally consistent with <10%% NaN.\n');
else
    if dim_mismatches > 0
        fprintf('  ❌ Dimensional mismatches: %d (dose/ADC vector length differs)\n', dim_mismatches);
    end
    if nan_warnings > 0
        fprintf('  💡 NaN dose warnings: %d pairs have >10%% NaN in-mask dose voxels (partial RT dose coverage)\n', nan_warnings);
    end
end

fprintf('\n======================================================\n');
fprintf('  Sanity checks complete.\n');
fprintf('======================================================\n');

% --- 4. Excessive NaN check: fail if >50% of voxels are NaN across the cohort ---
% This is the final go/no-go gate before visualization and statistical
% analysis. If more than half of all voxels across the entire cohort are
% NaN, the data is too sparse for reliable analysis — summary statistics
% would be dominated by the few remaining non-NaN values, and group
% comparisons would lack statistical power. The 50% threshold is a
% pragmatic balance: some NaN is expected (missing fractions, failed fits),
% but >50% indicates a systematic pipeline failure that requires investigation.
% NOTE: diary stays open so the pass/fail result is captured in the log.
% Use the active DWI type's fields (not always Standard) so that dnCNN/IVIMnet
% runs validate the correct vectors.
max_nan_frac = 0.5;
excessive_nan = false;
if exist('config_struct', 'var') && isfield(config_struct, 'dwi_types_to_run') && isnumeric(config_struct.dwi_types_to_run)
    dtype_nan = config_struct.dwi_types_to_run;
    if ~isscalar(dtype_nan), dtype_nan = dtype_nan(1); end
else
    dtype_nan = 1;
end
switch dtype_nan
    case 2
        nan_check_fields = {'adc_vector_dncnn', 'd_vector_dncnn', 'f_vector_dncnn', 'dstar_vector_dncnn'};
        nan_check_names  = {'ADC (DnCNN)', 'D (DnCNN)', 'f (DnCNN)', 'D* (DnCNN)'};
    case 3
        nan_check_fields = {'adc_vector', 'd_vector_ivimnet', 'f_vector_ivimnet', 'dstar_vector_ivimnet'};
        nan_check_names  = {'ADC', 'D (IVIMnet)', 'f (IVIMnet)', 'D* (IVIMnet)'};
    otherwise
        nan_check_fields = {'adc_vector', 'd_vector', 'f_vector', 'dstar_vector'};
        nan_check_names  = {'ADC', 'D', 'f', 'D*'};
end
for fi = 1:numel(nan_check_fields)
    total_voxels = 0;
    total_nans   = 0;
    for j = 1:nPat
        for k = 1:min(size(data_vectors_gtvp, 2), nTp)
            if j > size(data_vectors_gtvp, 1) || k > size(data_vectors_gtvp, 2), continue; end
            s_nan = data_vectors_gtvp(j, k, 1);
            if ~isfield(s_nan, nan_check_fields{fi}), continue; end
            v = s_nan.(nan_check_fields{fi});
            if ~isempty(v)
                total_voxels = total_voxels + numel(v);
                total_nans   = total_nans + sum(isnan(v));
            end
        end
    end
    if total_voxels > 0
        frac = total_nans / total_voxels;
        if frac > max_nan_frac
            fprintf('  ❌ Excessive NaN fraction in %s: %.1f%% (threshold: %.0f%%)\n', ...
                nan_check_names{fi}, frac*100, max_nan_frac*100);
            excessive_nan = true;
            % Diagnostic: help identify root cause for dnCNN/IVIMnet failures
            if dtype_nan == 2 || dtype_nan == 3
                n_present = 0; n_allnan = 0;
                for jj = 1:nPat
                    for kk = 1:min(size(data_vectors_gtvp, 2), nTp)
                        if jj > size(data_vectors_gtvp, 1) || kk > size(data_vectors_gtvp, 2), continue; end
                        s_diag = data_vectors_gtvp(jj, kk, 1);
                        if isfield(s_diag, nan_check_fields{fi})
                            v_diag = s_diag.(nan_check_fields{fi});
                            if ~isempty(v_diag), n_present = n_present + 1; end
                            if ~isempty(v_diag) && all(isnan(v_diag)), n_allnan = n_allnan + 1; end
                        end
                    end
                end
                fprintf('      💡 %d scans have %s data; %d are all-NaN.\n', n_present, nan_check_names{fi}, n_allnan);
                if dtype_nan == 2
                    fprintf('      💡 Check: (1) dependencies/dncnn_model.mat exists, (2) GTV masks are valid, (3) pre-cached dncnn/*.nii.gz files.\n');
                else
                    fprintf('      💡 Check: (1) IVIMnet .mat parameter maps exist, (2) GTV masks are valid.\n');
                end
            end
        end
    end
end

% --- Final go/no-go decision ---
% is_valid controls whether downstream modules (visualize, metrics) proceed.
% Dimensional mismatches are hard failures (spatial registration broke).
% Excessive NaN is a hard failure (data too sparse for statistics).
% NaN dose warnings and convergence issues are soft warnings (common in
% clinical pancreatic DWI data due to motion artifacts and incomplete scans).
if dim_mismatches == 0 && ~excessive_nan
    is_valid = true;
    validation_msg = sprintf('Passed all alignment checks. %d convergence warnings, %d NaN dose warnings.', conv_issues, nan_warnings);
elseif excessive_nan
    is_valid = false;
    validation_msg = sprintf('Failed: >%.0f%% NaN voxels detected. %d dimensional mismatches.', max_nan_frac*100, dim_mismatches);
else
    % Only dimensional mismatches invalidate the run — NaN dose coverage
    % is common when the RT dose grid does not fully overlap the DWI FOV.
    is_valid = false;
    validation_msg = sprintf('Failed %d dimensional alignment checks.', dim_mismatches);
end

diary off
end