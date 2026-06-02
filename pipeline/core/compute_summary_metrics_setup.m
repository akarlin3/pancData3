function [adc_thresh, high_adc_thresh, d_thresh, f_thresh, min_vox_hist, ...
    nTp, nRpt, ...
    adc_sub_vol_pc, adc_sub_vol, adc_sub_mean, adc_sub_kurt, adc_sub_skew, ...
    high_adc_sub_vol, high_adc_sub_vol_pc, f_sub_vol, gtv_vol, ...
    fdm_responding_pc, fdm_progressing_pc, fdm_stable_pc, ...
    ALL_CORE_METHODS, n_all_methods, run_all_core, store_masks, per_method, ...
    adc_mean, adc_kurt, adc_skew, adc_sd, ...
    d_mean, d_kurt, d_skew, d_sd, d_sub_mean, d_sub_kurt, d_sub_skew, ...
    f_mean, f_kurt, f_skew, dstar_mean, dstar_kurt, dstar_skew, ...
    bin_edges, adc_histograms, d_histograms, ...
    ks_stats_adc, ks_pvals_adc, ks_stats_d, ks_pvals_d, ...
    fx_corrupted, adc_max, use_texture, texture_features, ...
    adc_mean_rpt, adc_sub_rpt, adc_sub_vol_rpt, adc_sub_vol_pc_rpt, ...
    fx_corrupted_rpt, d_mean_rpt, f_mean_rpt, dstar_mean_rpt, n_rpt, ...
    dice_rpt_adc, hd_max_rpt_adc, hd95_rpt_adc, ...
    dice_rpt_d, hd_max_rpt_d, hd95_rpt_d, ...
    dice_rpt_f, hd_max_rpt_f, hd95_rpt_f, ...
    dice_rpt_dstar, hd_max_rpt_dstar, hd95_rpt_dstar, ...
    morph_se, morph_min_cc, gtv_mask_cache, last_rpt_gtv_mat, last_rpt_gtv_mask] = ...
    compute_summary_metrics_setup(config_struct, data_vectors_gtvp, id_list)
% COMPUTE_SUMMARY_METRICS_SETUP — Threshold extraction and metric-array
% pre-allocation for compute_summary_metrics.
%
% Extracted verbatim from compute_summary_metrics.m (the threshold / array
% setup phase) to keep that file under the per-file line budget.  Reads the
% configuration thresholds and the cohort/timepoint dimensions, then
% pre-allocates every patient x timepoint x pipeline metric array (NaN-filled),
% the histogram/KS arrays, repeatability arrays, spatial-repeatability arrays,
% and the morphological / caching scratch state used by the main loop.
%
% per_method and texture_features are always returned (empty defaults when the
% corresponding feature is disabled) so the caller can capture them
% unconditionally; downstream code guards their use behind run_all_core /
% use_texture.

% ADC threshold for identifying "restricted diffusion" sub-volume.
% Voxels with ADC <= adc_thresh (typically ~1.0e-3 mm^2/s) represent
% regions of restricted diffusion, which in pancreatic tumors indicates
% high cellular density (viable tumor tissue). The restricted sub-volume
% is a candidate biomarker for treatment response: successful RT should
% kill tumor cells, reducing cellularity and increasing ADC, thereby
% shrinking the restricted sub-volume over time.
adc_thresh = config_struct.adc_thresh;

% Secondary ADC threshold for identifying "high ADC" sub-volume.
% Voxels with ADC > high_adc_thresh represent regions of high diffusivity
% (e.g., necrosis, edema, cystic change). Growth of this sub-volume
% during treatment may indicate tumor necrosis (positive response) or
% radiation-induced edema (confounding effect).
high_adc_thresh = config_struct.high_adc_thresh;

% IVIM-specific thresholds for sub-volume identification.
% d_thresh: true diffusion threshold (analogous to adc_thresh but for
%   the IVIM D parameter, which excludes perfusion contributions)
% f_thresh: perfusion fraction threshold — voxels with f < f_thresh
%   have low microvascularity, potentially indicating avascular necrotic
%   regions or hypoxic tumor zones resistant to RT
d_thresh = config_struct.d_thresh;
f_thresh = config_struct.f_thresh;

% Minimum number of finite voxels required to compute higher-order
% statistics (kurtosis, skewness). With too few voxels, these statistics
% are unreliable and highly sensitive to individual outlier voxels.
% kurtosis requires >= 4 data points by definition; the threshold is
% set higher for practical stability.
min_vox_hist = config_struct.min_vox_hist;

nTp = size(data_vectors_gtvp, 2);   % number of timepoints (Fx1–Fx5 + post)
nRpt = size(data_vectors_gtvp, 3);  % max number of repeat scans at Fx1

% --- Pre-allocate sub-volume metric arrays (patient x timepoint x pipeline) ---
% The 3rd dimension indexes the DWI processing pipeline:
%   1 = Standard (raw DWI + conventional fitting)
%   2 = DnCNN (denoised DWI + conventional fitting)
%   3 = IVIMnet (raw DWI + neural network fitting)
% NaN initialization ensures missing data propagates correctly through
% nanmean/nanstd without corrupting calculations with zeros.
adc_sub_vol_pc = nan(length(id_list),nTp,3);   % restricted sub-volume as % of GTV
adc_sub_vol = nan(length(id_list),nTp,3);       % restricted sub-volume in cm^3
adc_sub_mean = nan(length(id_list),nTp,3);       % mean ADC within restricted sub-volume
adc_sub_kurt = nan(length(id_list),nTp,3);       % kurtosis of restricted sub-volume ADC
adc_sub_skew = nan(length(id_list),nTp,3);       % skewness of restricted sub-volume ADC
high_adc_sub_vol = nan(length(id_list),nTp,3);   % high-ADC sub-volume in cm^3
high_adc_sub_vol_pc = nan(length(id_list),nTp,3); % high-ADC sub-volume as % of GTV
f_sub_vol = nan(length(id_list),nTp,3);           % low-perfusion sub-volume in cm^3
gtv_vol = nan(length(id_list),nTp);               % total GTV volume in cm^3 (pipeline-independent)

% --- fDM (Functional Diffusion Map) volume fractions ---
% Only populated when core_method = 'fdm'.  Each fraction represents the
% proportion of GTV voxels in each fDM class (responding/stable/progressing).
fdm_responding_pc = nan(length(id_list),nTp,3);
fdm_progressing_pc = nan(length(id_list),nTp,3);
fdm_stable_pc = nan(length(id_list),nTp,3);

% --- Multi-method core computation setup ---
% When run_all_core_methods is true, the pipeline computes sub-volume
% metrics for all 11 core delineation methods per patient/timepoint,
% storing them in a nested struct keyed by method name.
if isfield(config_struct, 'active_core_methods') && ~isempty(config_struct.active_core_methods)
    ALL_CORE_METHODS = config_struct.active_core_methods;
else
    ALL_CORE_METHODS = {'adc_threshold', 'd_threshold', 'df_intersection', ...
        'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', ...
        'percentile', 'spectral', 'fdm'};
end
n_all_methods = numel(ALL_CORE_METHODS);
run_all_core = isfield(config_struct, 'run_all_core_methods') && config_struct.run_all_core_methods;
store_masks = run_all_core && isfield(config_struct, 'store_core_masks') && config_struct.store_core_masks;

per_method = struct();
if run_all_core
    for m_init = 1:n_all_methods
        mname = ALL_CORE_METHODS{m_init};
        per_method.(mname).adc_sub_vol = nan(length(id_list), nTp, 3);
        per_method.(mname).adc_sub_vol_pc = nan(length(id_list), nTp, 3);
        per_method.(mname).adc_sub_mean = nan(length(id_list), nTp, 3);
        per_method.(mname).adc_sub_kurt = nan(length(id_list), nTp, 3);
        per_method.(mname).adc_sub_skew = nan(length(id_list), nTp, 3);
        per_method.(mname).f_sub_vol = nan(length(id_list), nTp, 3);
        per_method.(mname).d_sub_mean = nan(length(id_list), nTp, 3);
        per_method.(mname).d_sub_kurt = nan(length(id_list), nTp, 3);
        per_method.(mname).d_sub_skew = nan(length(id_list), nTp, 3);
        per_method.(mname).fdm_responding_pc = nan(length(id_list), nTp, 3);
        per_method.(mname).fdm_progressing_pc = nan(length(id_list), nTp, 3);
        per_method.(mname).fdm_stable_pc = nan(length(id_list), nTp, 3);
        if store_masks
            per_method.(mname).core_masks = cell(length(id_list), nTp);
        end
    end
end

% --- Whole-GTV summary statistics for ADC ---
adc_mean = nan(length(id_list),nTp,3);
adc_kurt = nan(length(id_list),nTp,3);
adc_skew = nan(length(id_list),nTp,3);
adc_sd = nan(length(id_list),nTp,3);

% --- Whole-GTV summary statistics for IVIM parameters (D, f, D*) ---
d_mean = nan(length(id_list),nTp,3);
d_kurt = nan(length(id_list),nTp,3);
d_skew = nan(length(id_list),nTp,3);
d_sd = nan(length(id_list),nTp,3);

% D sub-volume statistics
d_sub_mean = nan(length(id_list),nTp,3);
d_sub_kurt = nan(length(id_list),nTp,3);
d_sub_skew = nan(length(id_list),nTp,3);

% Perfusion fraction (f) summary statistics
f_mean = nan(length(id_list),nTp,3);
f_kurt = nan(length(id_list),nTp,3);
f_skew = nan(length(id_list),nTp,3);

% Pseudo-diffusion coefficient (D*) summary statistics
dstar_mean = nan(length(id_list),nTp,3);
dstar_kurt = nan(length(id_list),nTp,3);
dstar_skew = nan(length(id_list),nTp,3);

% --- Histogram and KS-test arrays for longitudinal distribution comparison ---
% Histogram bin edges span the physiological range of diffusion coefficients
% in soft tissue: 0 to 3.0e-3 mm^2/s in steps of 0.5e-4 mm^2/s (60 bins).
% This range covers both highly restricted tumor tissue (~0.5e-3) and
% free water/necrosis (~2.5e-3). The bin width of 0.5e-4 provides
% sufficient resolution to detect subtle distributional shifts during RT.
bin_edges = 0:0.5e-4:3e-3;
adc_histograms = nan(length(id_list),nTp,length(bin_edges)-1,3);  % smoothed probability distributions
d_histograms = nan(length(id_list),nTp,length(bin_edges)-1,3);

% KS (Kolmogorov-Smirnov) test statistics compare each timepoint's
% voxel distribution to the Fx1 baseline. The KS statistic (0-1) measures
% the maximum difference between cumulative distribution functions.
% A large KS statistic at Fx3 vs Fx1 indicates the tumor's diffusion
% profile has shifted substantially during treatment — potentially a
% more sensitive response indicator than comparing means alone, because
% it captures changes in distribution shape (not just location).
ks_stats_adc = nan(length(id_list),nTp,3);
ks_pvals_adc = nan(length(id_list),nTp,3);
ks_stats_d = nan(length(id_list),nTp,3);
ks_pvals_d = nan(length(id_list),nTp,3);

% --- Motion corruption flag ---
% Fraction of voxels with ADC > adc_max (an unrealistically high value,
% typically > 3.0e-3 mm^2/s, approaching free water diffusivity).
% High ADC values at many voxels indicate bulk patient motion during the
% DWI acquisition, which causes signal averaging across tissue boundaries
% and artificially inflated apparent diffusion. Pancreatic DWI is
% particularly susceptible to respiratory motion artifacts.
fx_corrupted = nan(length(id_list),nTp,3);
adc_max = config_struct.adc_max;

% --- Texture features (optional, when use_texture_features is enabled) ---
% GLCM and first-order texture features capture spatial heterogeneity
% patterns in parameter maps that summary statistics (mean, kurtosis)
% may miss. These include contrast, correlation, energy, homogeneity,
% entropy, and other radiomics descriptors.
use_texture = isfield(config_struct, 'use_texture_features') && config_struct.use_texture_features;
if use_texture
    texture_features = cell(length(id_list), nTp, 3);
else
    texture_features = {};
end

% --- Repeatability arrays (Fx1 repeat scans only, for wCV calculation) ---
% Within-session repeat scans at Fx1 allow computation of within-subject
% coefficient of variation (wCV = SD_within / mean), which quantifies
% the inherent measurement variability of each diffusion parameter.
% This is critical for interpreting longitudinal changes: only changes
% exceeding wCV can be attributed to treatment effects with confidence.
% For example, if ADC wCV is 5%, a 3% change between Fx1 and Fx2 is
% within noise, but a 15% change is likely a true biological effect.
adc_mean_rpt = nan(length(id_list),nRpt,3);     % mean ADC per repeat scan
adc_sub_rpt = nan(length(id_list),nRpt,3);       % mean ADC in restricted sub-volume per repeat
adc_sub_vol_rpt = nan(length(id_list),nRpt,3);   % restricted sub-volume (cm^3) per repeat
adc_sub_vol_pc_rpt = nan(length(id_list),nRpt,3); % restricted sub-volume (fraction of GTV) per repeat
fx_corrupted_rpt = nan(length(id_list),nRpt,3);   % motion corruption flag per repeat
d_mean_rpt = nan(length(id_list),nRpt,3);         % mean D per repeat scan
f_mean_rpt = nan(length(id_list),nRpt,3);         % mean f per repeat scan
dstar_mean_rpt = nan(length(id_list),nRpt,3);     % mean D* per repeat scan

% Count of available repeat scans per patient (for wCV denominator)
n_rpt = nan(length(id_list),1);

% --- Spatial repeatability: Dice and Hausdorff between repeat sub-volumes ---
% These metrics quantify whether the spatial definition of the sensitive
% sub-volume is reproducible across same-session repeat DWI acquisitions.
% High Dice and low Hausdorff indicate that the sub-volume boundary is
% stable despite measurement noise; low Dice or high Hausdorff indicates
% that parameter noise causes the threshold-defined sub-region to shift
% spatially between acquisitions.  This complements wCV (which quantifies
% scalar mean repeatability) by capturing spatial agreement.
dice_rpt_adc = nan(length(id_list), 3);
hd_max_rpt_adc = nan(length(id_list), 3);
hd95_rpt_adc = nan(length(id_list), 3);
dice_rpt_d = nan(length(id_list), 3);
hd_max_rpt_d = nan(length(id_list), 3);
hd95_rpt_d = nan(length(id_list), 3);
dice_rpt_f = nan(length(id_list), 3);
hd_max_rpt_f = nan(length(id_list), 3);
hd95_rpt_f = nan(length(id_list), 3);
dice_rpt_dstar = nan(length(id_list), 3);
hd_max_rpt_dstar = nan(length(id_list), 3);
hd95_rpt_dstar = nan(length(id_list), 3);

% Morphological structuring element for sub-volume cleanup (same as
% calculate_subvolume_metrics.m).  A 3D sphere of radius 1 voxel is used
% for morphological opening to remove single-voxel noise from threshold-
% defined sub-volumes.  Without this cleanup, isolated voxels that just
% happen to fall below the ADC threshold (due to noise, not biology)
% would inflate the sub-volume size and corrupt spatial repeatability.
if exist('OCTAVE_VERSION', 'builtin')
    % Octave lacks strel('sphere'), so we build a 6-connected cross kernel
    sphere_kernel = zeros(3, 3, 3);
    sphere_kernel(2,2,:) = 1; sphere_kernel(2,:,2) = 1; sphere_kernel(:,2,2) = 1;
    morph_se = strel('arbitrary', sphere_kernel);
else
    morph_se = strel('sphere', 1);
end
morph_min_cc = 10;  % minimum connected component size (voxels) — clusters
                    % smaller than this are discarded as noise artifacts

% Cache for GTV mask loading (avoids repeated disk I/O for same file).
% Using containers.Map to cache ALL previously loaded masks by filename,
% not just the last one.  Many patient-timepoint pairs reuse the Fx1
% reference mask, so a full cache eliminates redundant safe_load_mask calls.
gtv_mask_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
% Legacy single-entry cache passed to compute_spatial_repeatability
last_rpt_gtv_mat = '';
last_rpt_gtv_mask = [];

end
