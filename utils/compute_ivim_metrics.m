function ivim_out = compute_ivim_metrics(config_struct, d_vec, f_vec, dstar_vec, ...
    d_baseline, adc_vec_sub_mask, vox_vol, min_vox_hist, bin_edges, ...
    d_thresh, f_thresh, k)
% COMPUTE_IVIM_METRICS — Computes IVIM summary metrics for a single patient/timepoint/DWI-type
%
% Extracts failed-fit filtering, whole-GTV D/f/D* statistics, D sub-volume
% metrics, f sub-volume metrics, D/f histograms, and KS tests vs baseline.
%
% Inputs:
%   config_struct    - Configuration struct with core_method, fdm_parameter, etc.
%   d_vec            - D voxel vector for this patient/timepoint
%   f_vec            - f voxel vector (modified in-place for failed fits)
%   dstar_vec        - D* voxel vector (modified in-place for failed fits)
%   d_baseline       - D voxel vector at Fx1 baseline (for KS test)
%   adc_vec_sub_mask - Logical core mask from ADC computation (for unified methods)
%   vox_vol          - Single voxel volume in cm^3
%   min_vox_hist     - Minimum voxels for kurtosis/skewness computation
%   bin_edges        - Histogram bin edges vector
%   d_thresh         - Threshold for D sub-volume (restricted diffusion)
%   f_thresh         - Threshold for f sub-volume (low perfusion)
%   k                - Timepoint index (for KS test skip at k==1)
%
% Outputs:
%   ivim_out         - Struct with fields:
%     .f_vec           - Filtered f vector (failed fits set to NaN)
%     .dstar_vec       - Filtered D* vector (failed fits set to NaN)
%     .d_mean_val      - Mean D across whole GTV
%     .d_kurt_val      - Kurtosis of D distribution
%     .d_skew_val      - Skewness of D distribution
%     .d_sd_val        - Standard deviation of D
%     .d_histogram     - Laplace-smoothed D histogram
%     .ks_stat_d       - KS test statistic vs baseline (NaN if k==1)
%     .ks_pval_d       - KS test p-value vs baseline (NaN if k==1)
%     .d_sub_mean_val  - Mean D within restricted sub-volume
%     .d_sub_kurt_val  - Kurtosis of D sub-volume
%     .d_sub_skew_val  - Skewness of D sub-volume
%     .f_mean_val      - Mean f across whole GTV
%     .f_kurt_val      - Kurtosis of f distribution
%     .f_skew_val      - Skewness of f distribution
%     .f_sub_vol_val   - f sub-volume in cm^3
%     .dstar_mean_val  - Mean D* across whole GTV
%     .dstar_kurt_val  - Kurtosis of D* distribution
%     .dstar_skew_val  - Skewness of D* distribution

% Initialize all output fields to NaN so downstream aggregation handles
% missing IVIM data (e.g., Standard DWI type has no IVIM parameters)
% without special-casing.
ivim_out = struct();
ivim_out.f_vec = f_vec;
ivim_out.dstar_vec = dstar_vec;
ivim_out.d_mean_val = NaN;
ivim_out.d_kurt_val = NaN;
ivim_out.d_skew_val = NaN;
ivim_out.d_sd_val = NaN;
ivim_out.d_histogram = zeros(1, length(bin_edges)-1);
ivim_out.ks_stat_d = NaN;
ivim_out.ks_pval_d = NaN;
ivim_out.d_sub_mean_val = NaN;
ivim_out.d_sub_kurt_val = NaN;
ivim_out.d_sub_skew_val = NaN;
ivim_out.f_mean_val = NaN;
ivim_out.f_kurt_val = NaN;
ivim_out.f_skew_val = NaN;
ivim_out.f_sub_vol_val = NaN;
ivim_out.dstar_mean_val = NaN;
ivim_out.dstar_kurt_val = NaN;
ivim_out.dstar_skew_val = NaN;

if isempty(d_vec)
    return;
end

% Exclude exact-zero f values that co-occur with failed D* fits
% (D* == 0 or NaN), preserving genuine zero perfusion.
failed_fit = (f_vec == 0) & (isnan(dstar_vec) | dstar_vec == 0);
f_vec(failed_fit) = nan;
dstar_vec(failed_fit) = nan;
ivim_out.f_vec = f_vec;
ivim_out.dstar_vec = dstar_vec;

% For unified core methods (percentile, spectral, fdm),
% the core mask replaces individual parameter thresholds.
unified_methods = {'percentile', 'spectral', 'fdm'};
if any(strcmpi(config_struct.core_method, unified_methods))
    f_vec_sub = f_vec(adc_vec_sub_mask);
    ivim_out.f_sub_vol_val = sum(adc_vec_sub_mask) * vox_vol;
    d_vec_sub = d_vec(adc_vec_sub_mask);
else
    f_vec_sub = f_vec(f_vec < f_thresh);
    ivim_out.f_sub_vol_val = numel(f_vec_sub) * vox_vol;
    d_vec_sub = d_vec(d_vec < d_thresh);
end

ivim_out.d_mean_val = nanmean_safe(d_vec);
[ivim_out.d_kurt_val, ivim_out.d_skew_val] = compute_kurt_skew(d_vec, min_vox_hist);
ivim_out.d_sd_val = nanstd_safe(d_vec);

ivim_out.d_histogram = compute_histogram_laplace(d_vec, bin_edges);
% NOTE: KS-test p-values are liberal (see ADC comment in compute_adc_metrics).
% Skip k==1: d_vec and d_baseline are the same data.
if k > 1 && ~isempty(d_baseline) && numel(d_vec) >= min_vox_hist && numel(d_baseline) >= min_vox_hist ...
        && any(~isnan(d_vec)) && any(~isnan(d_baseline))
    % Remove NaN before kstest2 (consistent with ADC path)
    d_vec_ks = d_vec(~isnan(d_vec));
    d_bl_ks  = d_baseline(~isnan(d_baseline));
    [~, p, ks2stat] = kstest2(d_vec_ks, d_bl_ks);
    ivim_out.ks_stat_d = ks2stat;
    ivim_out.ks_pval_d = p;
end

if isempty(d_vec_sub)
    ivim_out.d_sub_mean_val = NaN;
else
    ivim_out.d_sub_mean_val = nanmean_safe(d_vec_sub);
end
[ivim_out.d_sub_kurt_val, ivim_out.d_sub_skew_val] = compute_kurt_skew(d_vec_sub, min_vox_hist);

ivim_out.f_mean_val = nanmean_safe(f_vec);
[ivim_out.f_kurt_val, ivim_out.f_skew_val] = compute_kurt_skew(f_vec, min_vox_hist);

ivim_out.dstar_mean_val = nanmean_safe(dstar_vec);
[ivim_out.dstar_kurt_val, ivim_out.dstar_skew_val] = compute_kurt_skew(dstar_vec, min_vox_hist);

end

%% --- Local helper functions (duplicated from compute_summary_metrics.m) ---

function result = nanmean_safe(v)
if exist('OCTAVE_VERSION', 'builtin')
    tmp = v(~isnan(v));
    if isempty(tmp)
        result = NaN;
    else
        result = mean(tmp);
    end
else
    result = nanmean(v);
end
end

function result = nanstd_safe(v)
if exist('OCTAVE_VERSION', 'builtin')
    tmp = v(~isnan(v));
    if isempty(tmp)
        result = NaN;
    else
        result = std(tmp);
    end
else
    result = nanstd(v);
end
end

function [kurt_val, skew_val] = compute_kurt_skew(v, min_vox_hist)
kurt_val = NaN;
skew_val = NaN;
if numel(v) >= min_vox_hist
    v_finite = v(~isnan(v));
    if numel(v_finite) >= min_vox_hist
        kurt_val = kurtosis(v_finite);
        skew_val = skewness(v_finite);
    end
end
end

function p1 = compute_histogram_laplace(vec, bin_edges)
if exist('OCTAVE_VERSION', 'builtin')
    vec_f = vec(~isnan(vec));
    c1 = histc(vec_f, bin_edges);
    c1(end-1) = c1(end-1) + c1(end);
    c1 = c1(1:end-1);
else
    [c1, ~] = histcounts(vec, bin_edges);
end
n_binned = sum(c1);
nbins = length(c1);
if n_binned > 0
    p1 = (c1 + 1) / (n_binned + nbins);
else
    p1 = zeros(size(c1));
end
end
