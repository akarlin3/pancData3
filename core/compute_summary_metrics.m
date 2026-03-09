function summary_metrics = compute_summary_metrics(config_struct, data_vectors_gtvp, id_list, mrn_list, lf, immuno, gtv_locations, dwi_locations, dmean_gtvp, d95_gtvp, v50gy_gtvp, fx_dates)
% COMPUTE_SUMMARY_METRICS — Computes longitudinal summary metrics from data vectors
% Part of the load_dwi_data.m refactoring.
%
% Inputs:
%   config_struct     - Configuration struct with thresholds (adc_thresh, etc.)
%   data_vectors_gtvp - Struct array of voxel-level parameter vectors
%   id_list           - Cell array of patient folder IDs
%   mrn_list          - Cell array of patient MRNs
%   lf                - Array of local failure status
%   immuno            - Array of immunotherapy status
%   gtv_locations     - Cell array of GTV path locations
%   dwi_locations     - Cell array of DWI DICOM path locations
%   dmean_gtvp        - Array of mean dose inside GTV
%   d95_gtvp          - Array of D95 dose metric inside GTV
%   v50gy_gtvp        - Array of V50Gy dose metric inside GTV
%   fx_dates          - (Optional) Cell matrix of DICOM StudyDate strings
%                       (patients x fractions) from discover_patient_files
%
% Outputs:
%   summary_metrics   - Struct containing computed mean, kurtosis, skewness,
%                       SD, subsets, histogram features, predictability, etc.
%
% ANALYTICAL RATIONALE — VOXEL-TO-SUMMARY AGGREGATION
%   Voxel-level parameter maps contain thousands of values per patient per
%   timepoint, which are too granular for patient-level statistical modeling
%   (e.g., survival analysis, group comparisons). This function aggregates
%   voxel distributions into summary statistics that capture different
%   aspects of the intra-tumoral parameter distribution:
%
%   - Mean: Central tendency of the distribution. ADC_mean reflects overall
%     tumor cellularity; D_mean isolates true tissue diffusion from perfusion.
%
%   - Kurtosis: Tail heaviness of the distribution. High kurtosis indicates
%     a heterogeneous tumor with a mix of very high and very low diffusivity
%     regions — potentially indicating regions of necrosis (high D) coexisting
%     with dense viable tumor (low D).
%
%   - Skewness: Asymmetry of the distribution. Positive skew (tail toward
%     high values) may indicate emerging necrotic regions during treatment.
%     Negative skew (tail toward low values) may indicate therapy-resistant
%     dense tumor foci.
%
%   - Standard deviation: Spread of the distribution, directly reflecting
%     intra-tumoral heterogeneity. Heterogeneity is an independent prognostic
%     factor in many cancers.
%
%   - Sub-volume metrics: Volume of tumor below an ADC/D threshold
%     (restricted diffusion sub-volume) or above a threshold (high-ADC
%     sub-volume). These capture the proportion of the tumor that is
%     highly cellular vs necrotic/edematous.
%
%   - KS test statistics: Kolmogorov-Smirnov two-sample test comparing
%     each timepoint's voxel distribution to the baseline (Fx1) distribution.
%     The KS statistic quantifies the magnitude of distributional shift
%     during treatment, which may be more sensitive to treatment response
%     than mean changes alone.
%
%   - Repeatability metrics: Within-session repeat scans at Fx1 enable
%     calculation of within-subject coefficient of variation (wCV), which
%     defines the measurement noise floor. Only longitudinal changes
%     exceeding wCV can be confidently attributed to treatment effects.
%

if nargin < 12, fx_dates = {}; end

% Build a DWI-type-specific filename suffix so that Standard, dnCNN, and
% IVIMnet runs do not overwrite each other's cached summary metrics.
if isfield(config_struct, 'dwi_type_name')
    file_prefix = ['_' config_struct.dwi_type_name];
else
    file_prefix = '';
end
% Checkpoint file for summary metrics — allows skipping expensive
% recomputation on subsequent pipeline runs with the same cohort/config.
summary_metrics_file = fullfile(config_struct.dataloc, ['summary_metrics' file_prefix '.mat']);
if isfield(config_struct, 'use_checkpoints') && config_struct.use_checkpoints
    checkpoint_loaded = false;
    if exist(summary_metrics_file, 'file')
        fprintf('  [CHECKPOINT] Found existing %s. Loading and skipping metrics computation...\n', ['summary_metrics' file_prefix '.mat']);
        tmp_metrics = load(summary_metrics_file, 'summary_metrics');
        checkpoint_loaded = true;
    else
        fallback_metrics_file = fullfile(config_struct.dataloc, 'summary_metrics.mat');
        if exist(fallback_metrics_file, 'file')
            fprintf('  [CHECKPOINT] Specific %s not found but fallback %s exists. Loading and skipping metrics computation...\n', ['summary_metrics' file_prefix '.mat'], 'summary_metrics.mat');
            tmp_metrics = load(fallback_metrics_file, 'summary_metrics');
            checkpoint_loaded = true;
        end
    end
    if checkpoint_loaded
        % Validate checkpoint dimensions match current cohort before using
        nPat_expected = length(id_list);
        nRpt_expected = size(data_vectors_gtvp, 3);
        sm = tmp_metrics.summary_metrics;
        dims_ok = isfield(sm, 'adc_mean_rpt') && ...
                  isfield(sm, 'adc_sub_vol_rpt') && ...
                  size(sm.adc_mean_rpt, 1) == nPat_expected && ...
                  size(sm.adc_mean_rpt, 2) == nRpt_expected;
        if dims_ok
            summary_metrics = sm;
            % Ensure fx_dates survives stale checkpoints that pre-date its
            % addition.  The caller always passes the current fx_dates, so
            % graft it onto the checkpoint if missing.
            if ~isfield(summary_metrics, 'fx_dates')
                summary_metrics.fx_dates = fx_dates;
            end
            return;
        else
            fprintf('  [CHECKPOINT] Stale checkpoint (dimension mismatch). Recomputing...\n');
        end
    end
end

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
ALL_CORE_METHODS = {'adc_threshold', 'd_threshold', 'df_intersection', ...
    'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', ...
    'percentile', 'spectral', 'fdm'};
n_all_methods = numel(ALL_CORE_METHODS);
run_all_core = isfield(config_struct, 'run_all_core_methods') && config_struct.run_all_core_methods;
store_masks = run_all_core && isfield(config_struct, 'store_core_masks') && config_struct.store_core_masks;

if run_all_core
    per_method = struct();
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

% Cache for GTV mask loading (avoids repeated disk I/O for same file)
last_rpt_gtv_mat = '';
last_rpt_gtv_mask = [];

% --- Main analysis loop: patient × timepoint × DWI pipeline ---
% This triple-nested loop (patient × timepoint × DWI type) is the
% computational core that transforms voxel-level parameter vectors into
% the patient-level summary statistics used by all downstream modules.
n_patients_metrics = length(id_list);
for j=1:n_patients_metrics
    text_progress_bar(j, n_patients_metrics, 'Computing summary metrics');
    for k=1:nTp
        % Extract per-voxel volume (cm^3) for converting voxel counts to
        % physical volumes.  A scalar vox_vol indicates valid GTV data;
        % an empty or multi-element value indicates missing/corrupt data.
        if length(data_vectors_gtvp(j,k,1).vox_vol) == 1
            vox_vol = data_vectors_gtvp(j,k,1).vox_vol;
        else
            vox_vol = NaN;
        end

        % Load 3D GTV mask for this patient/timepoint (needed by
        % extract_tumor_core methods that require spatial context:
        % region_growing, active_contours, spectral with 3D cleanup).
        has_3d = false;
        gtv_mask_3d = [];
        if nargin >= 7 && ~isempty(gtv_locations) && ...
                size(gtv_locations, 1) >= j && size(gtv_locations, 2) >= k
            gtv_mat = gtv_locations{j, k, 1};
            if ~isempty(gtv_mat) && exist(gtv_mat, 'file')
                gtv_mask_3d = safe_load_mask(gtv_mat, 'Stvol3d');
                if ~isempty(gtv_mask_3d)
                    has_3d = true;
                end
            end
            % Fx1 mask fallback for DIR-warped timepoints:
            % At Fx2+, process_single_scan warps all parameter maps to
            % the Fx1 (baseline) coordinate frame via imregdemons, then
            % extracts voxel vectors using the Fx1 reference mask.  The
            % vectors therefore have length == sum(Fx1_mask), but
            % gtv_locations{j,k,1} still points to the native fraction's
            % mask file which has a different voxel count (the physician
            % re-contoured on the Fx2+ anatomy).  Substituting the Fx1
            % mask restores the correct 3D geometry for spatial methods
            % (active_contours, region_growing) that must map 1D vectors
            % back into 3D voxel positions.
            if has_3d && k > 1
                ref_vec = data_vectors_gtvp(j,k,1).adc_vector;
                if ~isempty(ref_vec) && sum(gtv_mask_3d(:) == 1) ~= numel(ref_vec)
                    fx1_mat = gtv_locations{j, 1, 1};
                    if ~isempty(fx1_mat) && exist(fx1_mat, 'file')
                        fx1_mask_3d = safe_load_mask(fx1_mat, 'Stvol3d');
                        if ~isempty(fx1_mask_3d) && sum(fx1_mask_3d(:) == 1) == numel(ref_vec)
                            gtv_mask_3d = fx1_mask_3d;
                        else
                            has_3d = false;
                            gtv_mask_3d = [];
                        end
                    else
                        has_3d = false;
                        gtv_mask_3d = [];
                    end
                end
            end
            % Diagnostic: log when 3D mask is unavailable so fallback
            % sources can be traced to specific patients/timepoints.
            if ~has_3d
                ref_vec_diag = data_vectors_gtvp(j,k,1).adc_vector;
                if ~isempty(ref_vec_diag)
                    fprintf('  [3D mask] %s Fx%d: no valid 3D mask (vec=%d voxels)\n', ...
                        id_list{j}, k, numel(ref_vec_diag));
                end
            end
        end

        for dwi_type = config_struct.dwi_types_to_run

            % Select the appropriate voxel vectors depending on pipeline.
            % Each DWI type stores its parameter vectors in different struct
            % fields (e.g., d_vector vs d_vector_dncnn vs d_vector_ivimnet).
            % The helper function abstracts this field selection.
            [adc_vec, d_vec, f_vec, dstar_vec] = select_dwi_vectors(data_vectors_gtvp, j, k, 1, dwi_type);
            % Baseline (Fx1) vectors are needed for KS tests (distributional
            % shift from baseline) and fDM (functional diffusion map) computation.
            [adc_baseline, d_baseline, ~, ~]   = select_dwi_vectors(data_vectors_gtvp, j, 1, 1, dwi_type);

            % Safety-net: confirm 3D mask voxel count matches this DWI
            % type's vector length.  The Fx1 fallback above resolves most
            % mismatches, but this guard catches residual edge cases (e.g.,
            % Fx1 mask also unavailable, or vector length differs across
            % DWI types due to pipeline-specific NaN pruning).
            has_3d_iter = has_3d;
            if has_3d && ~isempty(adc_vec) && sum(gtv_mask_3d(:) == 1) ~= numel(adc_vec)
                has_3d_iter = false;
            end

            % --- Compute ADC summary metrics for this patient/timepoint ---
            % Build core_opts for extract_tumor_core (baseline vectors for fDM)
            core_opts = struct();
            core_opts.timepoint_index = k;
            if k > 1
                switch dwi_type
                    case 1
                        core_opts.baseline_adc_vec = data_vectors_gtvp(j,1,1).adc_vector;
                        core_opts.baseline_d_vec = data_vectors_gtvp(j,1,1).d_vector;
                    case 2
                        core_opts.baseline_adc_vec = data_vectors_gtvp(j,1,1).adc_vector_dncnn;
                        core_opts.baseline_d_vec = data_vectors_gtvp(j,1,1).d_vector_dncnn;
                    case 3
                        core_opts.baseline_adc_vec = data_vectors_gtvp(j,1,1).adc_vector;
                        core_opts.baseline_d_vec = data_vectors_gtvp(j,1,1).d_vector_ivimnet;
                end
            end

            adc_out = compute_adc_metrics(config_struct, adc_vec, d_vec, f_vec, dstar_vec, ...
                adc_baseline, vox_vol, min_vox_hist, bin_edges, high_adc_thresh, adc_max, ...
                has_3d_iter, gtv_mask_3d, core_opts, k, j, data_vectors_gtvp, dwi_type);

            if ~isempty(adc_vec)
                n_finite_adc = sum(~isnan(adc_vec));
                if isnan(gtv_vol(j,k))
                    gtv_vol(j,k) = adc_out.gtv_vol_val;
                end
                adc_mean(j,k,dwi_type) = adc_out.adc_mean_val;
                adc_kurt(j,k,dwi_type) = adc_out.adc_kurt_val;
                adc_skew(j,k,dwi_type) = adc_out.adc_skew_val;
                adc_sd(j,k,dwi_type) = adc_out.adc_sd_val;
                adc_sub_vol(j,k,dwi_type) = adc_out.adc_sub_vol_val;
                adc_sub_vol_pc(j,k,dwi_type) = adc_out.adc_sub_vol_pc_val;
                adc_sub_mean(j,k,dwi_type) = adc_out.adc_sub_mean_val;
                adc_sub_kurt(j,k,dwi_type) = adc_out.adc_sub_kurt_val;
                adc_sub_skew(j,k,dwi_type) = adc_out.adc_sub_skew_val;
                high_adc_sub_vol(j,k,dwi_type) = adc_out.high_adc_sub_vol_val;
                high_adc_sub_vol_pc(j,k,dwi_type) = adc_out.high_adc_sub_vol_pc_val;
                adc_histograms(j,k,:,dwi_type) = adc_out.adc_histogram;
                ks_stats_adc(j,k,dwi_type) = adc_out.ks_stat_adc;
                ks_pvals_adc(j,k,dwi_type) = adc_out.ks_pval_adc;
                fx_corrupted(j,k,dwi_type) = adc_out.fx_corrupted_val;
                fdm_responding_pc(j,k,dwi_type) = adc_out.fdm_responding_pc;
                fdm_progressing_pc(j,k,dwi_type) = adc_out.fdm_progressing_pc;
                fdm_stable_pc(j,k,dwi_type) = adc_out.fdm_stable_pc;
                adc_vec_sub_mask = adc_out.adc_vec_sub_mask;
                finite_vol = n_finite_adc * vox_vol;
            end

            % --- Compute IVIM summary metrics (D, f, D*) ---
            ivim_out = compute_ivim_metrics(config_struct, d_vec, f_vec, dstar_vec, ...
                d_baseline, adc_out.adc_vec_sub_mask, vox_vol, min_vox_hist, bin_edges, ...
                d_thresh, f_thresh, k);

            if ~isempty(d_vec)
                d_mean(j,k,dwi_type) = ivim_out.d_mean_val;
                d_kurt(j,k,dwi_type) = ivim_out.d_kurt_val;
                d_skew(j,k,dwi_type) = ivim_out.d_skew_val;
                d_sd(j,k,dwi_type) = ivim_out.d_sd_val;
                d_histograms(j,k,:,dwi_type) = ivim_out.d_histogram;
                ks_stats_d(j,k,dwi_type) = ivim_out.ks_stat_d;
                ks_pvals_d(j,k,dwi_type) = ivim_out.ks_pval_d;
                d_sub_mean(j,k,dwi_type) = ivim_out.d_sub_mean_val;
                d_sub_kurt(j,k,dwi_type) = ivim_out.d_sub_kurt_val;
                d_sub_skew(j,k,dwi_type) = ivim_out.d_sub_skew_val;
                f_mean(j,k,dwi_type) = ivim_out.f_mean_val;
                f_kurt(j,k,dwi_type) = ivim_out.f_kurt_val;
                f_skew(j,k,dwi_type) = ivim_out.f_skew_val;
                f_sub_vol(j,k,dwi_type) = ivim_out.f_sub_vol_val;
                dstar_mean(j,k,dwi_type) = ivim_out.dstar_mean_val;
                dstar_kurt(j,k,dwi_type) = ivim_out.dstar_kurt_val;
                dstar_skew(j,k,dwi_type) = ivim_out.dstar_skew_val;
            end

            % --- Multi-method core metrics (when enabled) ---
            if run_all_core && ~isempty(adc_vec)
                per_method = compute_multi_core_metrics(per_method, config_struct, ...
                    ALL_CORE_METHODS, adc_vec, d_vec, f_vec, dstar_vec, ...
                    has_3d_iter, gtv_mask_3d, core_opts, ...
                    j, k, dwi_type, vox_vol, min_vox_hist, ...
                    f_thresh, d_thresh, finite_vol, store_masks);
            end

            % --- Repeatability analysis: extract metrics from Fx1 repeat scans ---
            % Only Fx1 (k==1) has multiple repeat acquisitions. These
            % back-to-back scans (same session, same setup) measure
            % inherent scan-to-scan variability. The resulting wCV
            % (within-subject coefficient of variation) determines the
            % minimum detectable change (MDC) for each parameter:
            %   MDC = wCV * 1.96 * sqrt(2)
            % Longitudinal changes smaller than MDC cannot be distinguished
            % from measurement noise at 95% confidence.
            if k==1
                rp_count = 0;
                for rpi=1:size(data_vectors_gtvp, 3)
                    [adc_vec, d_vec, f_vec, dstar_vec] = select_dwi_vectors(data_vectors_gtvp, j, k, rpi, dwi_type);

                    % Apply the same failed-fit filter used in the main
                    % metrics path (lines 297-299) so that repeatability
                    % wCV for f and D* is not biased by spurious zeros.
                    failed_fit_rpt = (f_vec == 0) & (isnan(dstar_vec) | dstar_vec == 0);
                    f_vec(failed_fit_rpt) = nan;
                    dstar_vec(failed_fit_rpt) = nan;

                    if ~isempty(adc_vec)
                        rp_count = rp_count+1;
                        adc_mean_rpt(j,rpi,dwi_type) = nanmean_safe(adc_vec);
                        n_finite_rpt = sum(~isnan(adc_vec));
                        if n_finite_rpt > 0
                            fx_corrupted_rpt(j,rpi,dwi_type) = sum(adc_vec > adc_max & ~isnan(adc_vec)) / n_finite_rpt;
                        end
                        % CORE DELINEATION METHOD ABSTRACTION
                        rpt_core_opts = struct('timepoint_index', k);
                        adc_vec_sub_mask_rpt = extract_tumor_core(config_struct, adc_vec, d_vec, f_vec, dstar_vec, has_3d_iter, gtv_mask_3d, rpt_core_opts);
                        adc_vec_sub = adc_vec(adc_vec_sub_mask_rpt);
                        if isempty(adc_vec_sub)
                            adc_sub_rpt(j,rpi,dwi_type) = NaN;
                        else
                            adc_sub_rpt(j,rpi,dwi_type) = nanmean_safe(adc_vec_sub);
                        end
                        % Restricted sub-volume size per repeat scan.
                        % Complements adc_sub_rpt (mean ADC within sub-volume)
                        % by tracking whether the sub-volume SIZE is reproducible.
                        adc_sub_vol_rpt(j,rpi,dwi_type) = numel(adc_vec_sub) * vox_vol;
                        if n_finite_rpt > 0
                            finite_vol_rpt = n_finite_rpt * vox_vol;
                            adc_sub_vol_pc_rpt(j,rpi,dwi_type) = adc_sub_vol_rpt(j,rpi,dwi_type) / finite_vol_rpt;
                        end
                    end

                    if ~isempty(d_vec)
                        % Count D/f/D* repeats when ADC is absent (e.g.,
                        % IVIMnet which reuses standard-pipeline ADC).
                        if isempty(adc_vec)
                            rp_count = rp_count + 1;
                        end
                        d_mean_rpt(j,rpi,dwi_type) = nanmean_safe(d_vec);
                        f_mean_rpt(j,rpi,dwi_type) = nanmean_safe(f_vec);
                        dstar_mean_rpt(j,rpi,dwi_type) = nanmean_safe(dstar_vec);
                    end
                end
                % Record repeat count from any DWI type (first non-zero wins).
                % Previously only dwi_type==1 populated n_rpt, leaving it
                % all-NaN for dnCNN/IVIMnet-only runs and breaking wCV.
                if isnan(n_rpt(j)) || n_rpt(j) == 0
                    n_rpt(j) = rp_count;
                end

                % --- Spatial repeatability: Dice & Hausdorff between repeat sub-volumes ---
                % Compare threshold-defined sensitive sub-volumes across all
                % pairs of Fx1 repeat scans to assess spatial reproducibility.
                if rp_count >= 2
                    [dice_rpt_adc(j, dwi_type), hd_max_rpt_adc(j, dwi_type), hd95_rpt_adc(j, dwi_type), ...
                     dice_rpt_d(j, dwi_type), hd_max_rpt_d(j, dwi_type), hd95_rpt_d(j, dwi_type), ...
                     dice_rpt_f(j, dwi_type), hd_max_rpt_f(j, dwi_type), hd95_rpt_f(j, dwi_type), ...
                     dice_rpt_dstar(j, dwi_type), hd_max_rpt_dstar(j, dwi_type), hd95_rpt_dstar(j, dwi_type)] = ...
                        compute_spatial_repeatability(data_vectors_gtvp, j, dwi_type, ...
                            gtv_locations, adc_thresh, d_thresh, f_thresh, config_struct.dstar_thresh, ...
                            morph_se, morph_min_cc, last_rpt_gtv_mat, last_rpt_gtv_mask);
                end
            end
        end
    end
end

% --- Package all computed metrics into a single output struct ---
% This struct is the primary interface between the 'load' step and all
% downstream analysis modules (sanity_checks, visualize_results,
% metrics_baseline, metrics_longitudinal, metrics_survival, etc.).
summary_metrics = struct();
summary_metrics.adc_mean = adc_mean;
summary_metrics.adc_kurt = adc_kurt;
summary_metrics.adc_skew = adc_skew;
summary_metrics.adc_sd = adc_sd;
summary_metrics.d_mean = d_mean;
summary_metrics.f_mean = f_mean;
summary_metrics.dstar_mean = dstar_mean;
summary_metrics.f_sub_vol = f_sub_vol;
summary_metrics.adc_sub_vol = adc_sub_vol;
summary_metrics.adc_sub_vol_pc = adc_sub_vol_pc;
summary_metrics.high_adc_sub_vol = high_adc_sub_vol;
summary_metrics.high_adc_sub_vol_pc = high_adc_sub_vol_pc;
summary_metrics.d_kurt = d_kurt;
summary_metrics.d_skew = d_skew;
summary_metrics.d_sd = d_sd;
summary_metrics.f_kurt = f_kurt;
summary_metrics.f_skew = f_skew;
summary_metrics.dstar_kurt = dstar_kurt;
summary_metrics.dstar_skew = dstar_skew;
summary_metrics.d_sub_mean = d_sub_mean;
summary_metrics.d_sub_kurt = d_sub_kurt;
summary_metrics.d_sub_skew = d_sub_skew;
summary_metrics.adc_histograms = adc_histograms;
summary_metrics.d_histograms = d_histograms;
summary_metrics.ks_stats_adc = ks_stats_adc;
summary_metrics.ks_pvals_adc = ks_pvals_adc;
summary_metrics.ks_stats_d = ks_stats_d;
summary_metrics.ks_pvals_d = ks_pvals_d;
summary_metrics.adc_sub_mean = adc_sub_mean;
summary_metrics.adc_sub_kurt = adc_sub_kurt;
summary_metrics.adc_sub_skew = adc_sub_skew;
summary_metrics.fx_corrupted = fx_corrupted;
summary_metrics.fx_corrupted_rpt = fx_corrupted_rpt;
summary_metrics.gtv_vol = gtv_vol;
summary_metrics.id_list = id_list;
summary_metrics.mrn_list = mrn_list;
summary_metrics.d95_gtvp = d95_gtvp;
summary_metrics.v50gy_gtvp = v50gy_gtvp;
summary_metrics.lf = lf;
summary_metrics.immuno = immuno;
summary_metrics.adc_mean_rpt = adc_mean_rpt;
summary_metrics.adc_sub_rpt = adc_sub_rpt;
summary_metrics.adc_sub_vol_rpt = adc_sub_vol_rpt;
summary_metrics.adc_sub_vol_pc_rpt = adc_sub_vol_pc_rpt;
summary_metrics.d_mean_rpt = d_mean_rpt;
summary_metrics.f_mean_rpt = f_mean_rpt;
summary_metrics.dstar_mean_rpt = dstar_mean_rpt;
summary_metrics.n_rpt = n_rpt;
summary_metrics.dice_rpt_adc = dice_rpt_adc;
summary_metrics.hd_max_rpt_adc = hd_max_rpt_adc;
summary_metrics.hd95_rpt_adc = hd95_rpt_adc;
summary_metrics.dice_rpt_d = dice_rpt_d;
summary_metrics.hd_max_rpt_d = hd_max_rpt_d;
summary_metrics.hd95_rpt_d = hd95_rpt_d;
summary_metrics.dice_rpt_f = dice_rpt_f;
summary_metrics.hd_max_rpt_f = hd_max_rpt_f;
summary_metrics.hd95_rpt_f = hd95_rpt_f;
summary_metrics.dice_rpt_dstar = dice_rpt_dstar;
summary_metrics.hd_max_rpt_dstar = hd_max_rpt_dstar;
summary_metrics.hd95_rpt_dstar = hd95_rpt_dstar;
summary_metrics.dmean_gtvp = dmean_gtvp;
summary_metrics.gtv_locations = gtv_locations;
summary_metrics.dwi_locations = dwi_locations;
summary_metrics.fx_dates = fx_dates;
summary_metrics.core_method = config_struct.core_method;
summary_metrics.fdm_responding_pc = fdm_responding_pc;
summary_metrics.fdm_progressing_pc = fdm_progressing_pc;
summary_metrics.fdm_stable_pc = fdm_stable_pc;

% Package per-method core metrics when multi-method computation was enabled
if run_all_core
    summary_metrics.all_core_metrics = per_method;
end

if isfield(config_struct, 'use_checkpoints') && config_struct.use_checkpoints
    fprintf('  [CHECKPOINT] Saving summary_metrics to %s...\n', summary_metrics_file);
    save(summary_metrics_file, 'summary_metrics');
end

end

function result = nanmean_safe(v)
% NANMEAN_SAFE — Octave-compatible NaN-ignoring mean.
% MATLAB's nanmean is in the Statistics Toolbox; Octave may lack it.
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
% NANSTD_SAFE — Octave-compatible NaN-ignoring standard deviation.
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
% COMPUTE_KURT_SKEW — Compute kurtosis and skewness with minimum sample guard.
% Returns NaN if the number of finite voxels is below min_vox_hist, because
% higher-order moments are unreliable with too few data points (kurtosis
% formally requires n >= 4, but practical stability needs more).
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
% COMPUTE_HISTOGRAM_LAPLACE — Laplace-smoothed probability distribution.
% Adds 1 pseudo-count per bin (Laplace smoothing / additive smoothing) to
% prevent zero-probability bins that would cause log(0) = -Inf in KL
% divergence or other information-theoretic comparisons.  The smoothing
% has negligible effect when n_binned >> nbins (typical for >100 voxels).
if exist('OCTAVE_VERSION', 'builtin')
    vec_f = vec(~isnan(vec));
    c1 = histc(vec_f, bin_edges);
    % histc includes a count for exact matches of the last edge; merge it
    % into the penultimate bin to match histcounts behavior.
    c1(end-1) = c1(end-1) + c1(end);
    c1 = c1(1:end-1);
else
    [c1, ~] = histcounts(vec, bin_edges);
end
n_binned = sum(c1);
nbins = length(c1);
if n_binned > 0
    % Laplace smoothing: P(bin) = (count + 1) / (total + n_bins)
    p1 = (c1 + 1) / (n_binned + nbins);
else
    p1 = zeros(size(c1));
end
end

function [adc_vec, d_vec, f_vec, dstar_vec] = select_dwi_vectors(data_vectors_gtvp, j, k, rpi, dwi_type)
% SELECT_DWI_VECTORS — Extract parameter vectors for the specified DWI pipeline.
%   dwi_type 1 (Standard): raw DWI with conventional mono/bi-exponential fits
%   dwi_type 2 (dnCNN): DnCNN-denoised DWI signal with conventional fits
%   dwi_type 3 (IVIMnet): raw DWI with neural-network IVIM parameter estimation
%   Note: IVIMnet reuses standard ADC (mono-exponential fit is model-free and
%   unaffected by the IVIM fitting method choice).
switch dwi_type
    case 1
        adc_vec   = data_vectors_gtvp(j,k,rpi).adc_vector;
        d_vec     = data_vectors_gtvp(j,k,rpi).d_vector;
        f_vec     = data_vectors_gtvp(j,k,rpi).f_vector;
        dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector;
    case 2
        adc_vec   = data_vectors_gtvp(j,k,rpi).adc_vector_dncnn;
        d_vec     = data_vectors_gtvp(j,k,rpi).d_vector_dncnn;
        f_vec     = data_vectors_gtvp(j,k,rpi).f_vector_dncnn;
        dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector_dncnn;
    case 3
        adc_vec   = data_vectors_gtvp(j,k,rpi).adc_vector;
        d_vec     = data_vectors_gtvp(j,k,rpi).d_vector_ivimnet;
        f_vec     = data_vectors_gtvp(j,k,rpi).f_vector_ivimnet;
        dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector_ivimnet;
end
end
