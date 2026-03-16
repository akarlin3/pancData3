function adc_out = compute_adc_metrics(config_struct, adc_vec, d_vec, f_vec, dstar_vec, ...
    adc_baseline, vox_vol, min_vox_hist, bin_edges, high_adc_thresh, adc_max, ...
    has_3d_iter, gtv_mask_3d, core_opts, k, j, data_vectors_gtvp, dwi_type)
% COMPUTE_ADC_METRICS — Computes ADC summary metrics for a single patient/timepoint/DWI-type
%
% Extracts GTV volume, whole-GTV ADC statistics, tumor core delineation via
% extract_tumor_core, restricted sub-volume, high-ADC sub-volume, fDM volume
% fractions, histogram, KS test vs baseline, and motion corruption flag.
%
% Inputs:
%   config_struct   - Configuration struct with thresholds and core_method
%   adc_vec         - ADC voxel vector for this patient/timepoint (column)
%   d_vec           - D voxel vector (may be empty)
%   f_vec           - f voxel vector (may be empty)
%   dstar_vec       - D* voxel vector (may be empty)
%   adc_baseline    - ADC voxel vector at Fx1 baseline (for KS test)
%   vox_vol         - Single voxel volume in cm^3
%   min_vox_hist    - Minimum voxels for kurtosis/skewness computation
%   bin_edges       - Histogram bin edges vector
%   high_adc_thresh - Threshold for high-ADC sub-volume
%   adc_max         - Maximum ADC for motion corruption detection
%   has_3d_iter     - Boolean: whether valid 3D mask is available for this iteration
%   gtv_mask_3d     - 3D logical mask (or empty)
%   core_opts       - Struct with timepoint_index and optional baseline vectors for fDM
%   k               - Timepoint index (for KS test skip at k==1)
%   j               - Patient index (for baseline vector lookup)
%   data_vectors_gtvp - Full data vectors struct (for baseline fDM lookup)
%   dwi_type        - DWI type index (1=Standard, 2=dnCNN, 3=IVIMnet)
%
% Outputs:
%   adc_out         - Struct with fields:
%     .gtv_vol_val        - GTV volume (numel(adc_vec)*vox_vol)
%     .adc_mean_val       - Mean ADC across whole GTV
%     .adc_kurt_val       - Kurtosis of ADC distribution
%     .adc_skew_val       - Skewness of ADC distribution
%     .adc_sd_val         - Standard deviation of ADC
%     .adc_sub_vol_val    - Restricted sub-volume in cm^3
%     .adc_sub_vol_pc_val - Restricted sub-volume as fraction of finite GTV
%     .adc_sub_mean_val   - Mean ADC within restricted sub-volume
%     .adc_sub_kurt_val   - Kurtosis within restricted sub-volume
%     .adc_sub_skew_val   - Skewness within restricted sub-volume
%     .high_adc_sub_vol_val    - High-ADC sub-volume in cm^3
%     .high_adc_sub_vol_pc_val - High-ADC sub-volume as fraction of finite GTV
%     .adc_histogram       - Laplace-smoothed histogram counts
%     .ks_stat_adc         - KS test statistic vs baseline (NaN if k==1)
%     .ks_pval_adc         - KS test p-value vs baseline (NaN if k==1)
%     .fx_corrupted_val    - Fraction of voxels above adc_max
%     .adc_vec_sub_mask    - Logical mask from extract_tumor_core (for IVIM reuse)
%     .fdm_responding_pc   - fDM responding fraction (NaN if not applicable)
%     .fdm_progressing_pc  - fDM progressing fraction (NaN if not applicable)
%     .fdm_stable_pc       - fDM stable fraction (NaN if not applicable)

% Initialize all output fields to NaN so that downstream code can safely
% aggregate across patients/timepoints without special-casing missing data.
adc_out = struct();
adc_out.gtv_vol_val = NaN;
adc_out.adc_mean_val = NaN;
adc_out.adc_kurt_val = NaN;
adc_out.adc_skew_val = NaN;
adc_out.adc_sd_val = NaN;
adc_out.adc_sub_vol_val = NaN;
adc_out.adc_sub_vol_pc_val = NaN;
adc_out.adc_sub_mean_val = NaN;
adc_out.adc_sub_kurt_val = NaN;
adc_out.adc_sub_skew_val = NaN;
adc_out.high_adc_sub_vol_val = NaN;
adc_out.high_adc_sub_vol_pc_val = NaN;
adc_out.adc_histogram = zeros(1, length(bin_edges)-1);
adc_out.ks_stat_adc = NaN;
adc_out.ks_pval_adc = NaN;
adc_out.fx_corrupted_val = NaN;
adc_out.adc_vec_sub_mask = false(size(adc_vec));
adc_out.fdm_responding_pc = NaN;
adc_out.fdm_progressing_pc = NaN;
adc_out.fdm_stable_pc = NaN;

if isempty(adc_vec)
    return;
end

n_finite_adc = sum(~isnan(adc_vec));

% GTV volume: use numel(adc_vec), not n_finite_adc, because the GTV
% mask defines the anatomical volume — NaN voxels (from failed ADC
% fits or artefacts) are still part of the tumour contour.
adc_out.gtv_vol_val = numel(adc_vec) * vox_vol;

% Whole-GTV ADC distribution statistics.
% Kurtosis and skewness characterize the shape of the ADC distribution:
%   - High kurtosis = heavy tails (heterogeneous tumor microenvironment)
%   - Negative skewness = tail toward low ADC (dominant restricted diffusion)
adc_out.adc_mean_val = nanmean_safe(adc_vec);
[adc_out.adc_kurt_val, adc_out.adc_skew_val] = compute_kurt_skew(adc_vec, min_vox_hist);
adc_out.adc_sd_val = nanstd_safe(adc_vec);

% CORE DELINEATION: Extract the restricted-diffusion tumor core using the
% configured method (threshold, Otsu, GMM, fDM, etc.). The returned mask
% is reused by compute_ivim_metrics for unified core methods.
adc_vec_sub_mask = extract_tumor_core(config_struct, adc_vec, d_vec, f_vec, dstar_vec, has_3d_iter, gtv_mask_3d, core_opts);
adc_vec_sub = adc_vec(adc_vec_sub_mask);           % ADC values within the core
adc_vec_high_sub = adc_vec(adc_vec > high_adc_thresh);  % High-ADC sub-volume (necrosis/edema)
adc_out.adc_vec_sub_mask = adc_vec_sub_mask;  % Export mask for IVIM reuse

% Compute Functional Diffusion Map (fDM) volume fractions when fDM core
% method is active and we are past baseline (k > 1). fDM classifies each
% voxel into responding/stable/progressing based on whether the diffusion
% parameter changed beyond the significance threshold from baseline.
if strcmpi(config_struct.core_method, 'fdm') && k > 1
    % Select which parameter to use for voxel-level change detection
    switch lower(config_struct.fdm_parameter)
        case 'adc'
            fdm_current = adc_vec;
            fdm_baseline = core_opts.baseline_adc_vec;
        case 'd'
            fdm_current = d_vec;
            fdm_baseline = core_opts.baseline_d_vec;
    end
    if ~isempty(fdm_baseline) && numel(fdm_baseline) == numel(fdm_current)
        % Use repeatability-derived Coefficient of Reproducibility (CoR)
        % when available; otherwise fall back to config threshold
        fdm_sig = config_struct.fdm_thresh;
        if isfield(core_opts, 'repeatability_cor') && ~isnan(core_opts.repeatability_cor)
            fdm_sig = core_opts.repeatability_cor;
        end
        delta_fdm = fdm_current - fdm_baseline;
        valid_delta = ~isnan(delta_fdm);
        n_valid_fdm = sum(valid_delta);
        if n_valid_fdm > 0
            % Responding: diffusivity increased beyond threshold (cell kill / edema)
            adc_out.fdm_responding_pc  = sum(delta_fdm(valid_delta) > fdm_sig) / n_valid_fdm;
            % Progressing: diffusivity decreased beyond threshold (increased cellularity)
            adc_out.fdm_progressing_pc = sum(delta_fdm(valid_delta) < -fdm_sig) / n_valid_fdm;
            % Stable: change within noise floor (no detectable treatment effect)
            adc_out.fdm_stable_pc      = sum(abs(delta_fdm(valid_delta)) <= fdm_sig) / n_valid_fdm;
        end
    end
end

% Restricted sub-volume: absolute volume (cm^3) and as a fraction of GTV.
adc_out.adc_sub_vol_val = sum(adc_vec_sub_mask) * vox_vol;
% Use count of finite (non-NaN) voxels as denominator so that voxels
% with failed ADC fits do not artificially deflate the percentage.
finite_vol = n_finite_adc * vox_vol;
if finite_vol > 0
    adc_out.adc_sub_vol_pc_val = adc_out.adc_sub_vol_val / finite_vol;
else
    adc_out.adc_sub_vol_pc_val = NaN;
end
if isempty(adc_vec_sub)
    adc_out.adc_sub_mean_val = NaN;
else
    adc_out.adc_sub_mean_val = nanmean_safe(adc_vec_sub);
end
[adc_out.adc_sub_kurt_val, adc_out.adc_sub_skew_val] = compute_kurt_skew(adc_vec_sub, min_vox_hist);

% Laplace-smoothed histogram: adds 1 to each bin count (additive smoothing)
% to avoid zero-probability bins in downstream KL-divergence or entropy calculations.
adc_out.adc_histogram = compute_histogram_laplace(adc_vec, bin_edges);
% NOTE: KS-test p-values are liberal because within-patient
% voxels are spatially autocorrelated (violates independence).
% Treat as descriptive, not inferential.
% Skip k==1: adc_vec and adc_baseline are the same data,
% so the KS test trivially returns stat=0, p=1.
if k > 1 && ~isempty(adc_baseline) && numel(adc_vec) >= min_vox_hist && numel(adc_baseline) >= min_vox_hist ...
        && any(~isnan(adc_vec)) && any(~isnan(adc_baseline))
    % Remove NaN before kstest2 (required for Octave compatibility)
    adc_vec_ks = adc_vec(~isnan(adc_vec));
    adc_bl_ks  = adc_baseline(~isnan(adc_baseline));
    [~, p, ks2stat] = kstest2(adc_vec_ks, adc_bl_ks);
    adc_out.ks_stat_adc = ks2stat;
    adc_out.ks_pval_adc = p;
end

% High-ADC sub-volume: identifies necrotic or edematous regions
% (unrestricted water diffusion above high_adc_thresh)
adc_out.high_adc_sub_vol_val = numel(adc_vec_high_sub) * vox_vol;
if finite_vol > 0
    adc_out.high_adc_sub_vol_pc_val = adc_out.high_adc_sub_vol_val / finite_vol;
else
    adc_out.high_adc_sub_vol_pc_val = NaN;
end
% Motion corruption flag: fraction of voxels with ADC above physical maximum.
% ADC values exceeding adc_max (typically 3e-3 mm^2/s) indicate motion
% artifact or failed fitting rather than genuine tissue properties.
if n_finite_adc > 0
    adc_out.fx_corrupted_val = sum(adc_vec > adc_max & ~isnan(adc_vec)) / n_finite_adc;
end

end

% Local helper functions (nanmean_safe, nanstd_safe, compute_kurt_skew,
% compute_histogram_laplace) have been extracted to pipeline/utils/ as
% shared utilities. These functions are on the MATLAB path and are
% accessible from within parfor loops as standalone .m files.
