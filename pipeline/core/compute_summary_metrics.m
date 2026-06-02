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

% ============================================================
% SECTION INDEX (line numbers as of v2.4.0)
% ------------------------------------------------------------
%  No inline local functions — helpers in pipeline/utils/.
%    L  84 — DWI-type file prefix + checkpoint load (with staleness check)
%    L 189 — ADC / IVIM sub-volume thresholds
%    L 224 — Pre-allocate metric arrays (patient × timepoint × pipeline)
%    L 285 — Whole-GTV summary stats (ADC, D, f, D*)
%    L 312 — Histogram + KS-test arrays (longitudinal distribution shifts)
%    L 334 — Motion corruption flag
%    L 344 — Texture-feature setup (optional)
%    L 354 — Repeatability arrays (Fx1 repeats, wCV)
%    L 374 — Spatial repeatability (Dice/Hausdorff between repeat sub-volumes)
%    L 429 — Main analysis loop (patient × timepoint × DWI pipeline)
%    L 726 — Package output struct
% ============================================================

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
    % Invalidate the summary_metrics checkpoint when its upstream
    % dependency (pipeline_voxels_*.mat) is newer.  Without this guard, a
    % fresh voxel-cache run that re-discovers file paths (e.g. newly
    % contoured Fx1 GTV masks) is silently paired with a stale summary
    % that still reflects the old discovery state, producing all-NaN
    % downstream repeatability fields.
    voxel_cache_file = fullfile(config_struct.dataloc, ['pipeline_voxels' file_prefix '.mat']);
    is_stale_vs_vectors = false;
    if exist(summary_metrics_file, 'file') && exist(voxel_cache_file, 'file')
        sm_info = dir(summary_metrics_file);
        dv_info = dir(voxel_cache_file);
        if ~isempty(sm_info) && ~isempty(dv_info) && sm_info.datenum < dv_info.datenum
            is_stale_vs_vectors = true;
            fprintf('  [CHECKPOINT] %s is older than %s — ignoring stale checkpoint and recomputing.\n', ...
                ['summary_metrics' file_prefix '.mat'], ['pipeline_voxels' file_prefix '.mat']);
        end
    end
    if exist(summary_metrics_file, 'file') && ~is_stale_vs_vectors
        fprintf('  [CHECKPOINT] Found existing %s. Loading and skipping metrics computation...\n', ['summary_metrics' file_prefix '.mat']);
        tmp_metrics = load(summary_metrics_file, 'summary_metrics');
        checkpoint_loaded = true;
    elseif ~is_stale_vs_vectors
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

        % Invalidate checkpoints that pre-date the shared-Fx1 fallback
        % used by compute_spatial_repeatability / optimize_adc_threshold.
        % Trigger when: (a) dice_rpt_adc field is missing entirely (pre-v2.4
        % checkpoint), or (b) the column for the current DWI type is all-NaN.
        % In either case, only invalidate if the cohort actually has patients
        % with >=2 non-empty Fx1 repeat ADC vectors.
        rpt_stale = false;
        if dims_ok
            dwi_col = 1;
            if isfield(config_struct, 'dwi_types_to_run') && ...
                    ~isempty(config_struct.dwi_types_to_run)
                dwi_col = config_struct.dwi_types_to_run(1);
            end

            dice_all_nan = true;  % assume stale until proven otherwise
            if isfield(sm, 'dice_rpt_adc')
                dra = sm.dice_rpt_adc;
                if size(dra, 2) >= dwi_col && ~all(isnan(dra(:, dwi_col)))
                    dice_all_nan = false;
                end
            end

            if dice_all_nan
                for jj = 1:size(data_vectors_gtvp, 1)
                    n_nonempty = 0;
                    for rr = 1:size(data_vectors_gtvp, 3)
                        if ~isempty(data_vectors_gtvp(jj, 1, rr).adc_vector)
                            n_nonempty = n_nonempty + 1;
                        end
                    end
                    if n_nonempty >= 2
                        rpt_stale = true; break;
                    end
                end
            end
        end

        if dims_ok && ~rpt_stale
            summary_metrics = sm;
            % Ensure fx_dates survives stale checkpoints that pre-date its
            % addition.  The caller always passes the current fx_dates, so
            % graft it onto the checkpoint if missing.
            if ~isfield(summary_metrics, 'fx_dates')
                summary_metrics.fx_dates = fx_dates;
            end
            return;
        elseif rpt_stale
            fprintf('  [CHECKPOINT] Stale checkpoint (dice_rpt all-NaN but cohort has Fx1 repeats). Recomputing...\n');
        else
            fprintf('  [CHECKPOINT] Stale checkpoint (dimension mismatch). Recomputing...\n');
        end
    end
end

% Threshold extraction + metric-array pre-allocation are delegated to
% compute_summary_metrics_setup so this orchestrator stays under the
% per-file line budget.  The helper reads the configuration thresholds
% and cohort/timepoint dimensions, then pre-allocates every metric,
% histogram/KS, repeatability, and spatial-repeatability array plus the
% morphological / caching scratch state used by the main loop.  per_method
% and texture_features are always returned (empty defaults when disabled).
[adc_thresh, high_adc_thresh, d_thresh, f_thresh, min_vox_hist, ...
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
    compute_summary_metrics_setup(config_struct, data_vectors_gtvp, id_list);

% Reset the per-call trace file and persistent counter so each pipeline
% run starts with a clean repeat_dice_trace.txt.
trace_reset = fullfile(tempdir, 'repeat_dice_trace.txt');
if exist(trace_reset, 'file') == 2
    try; delete(trace_reset); catch; end
end
try; clear compute_spatial_repeatability; catch; end

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
                if gtv_mask_cache.isKey(gtv_mat)
                    gtv_mask_3d = gtv_mask_cache(gtv_mat);
                else
                    gtv_mask_3d = safe_load_mask(gtv_mat, 'Stvol3d');
                    gtv_mask_cache(gtv_mat) = gtv_mask_3d;
                end
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
                        if gtv_mask_cache.isKey(fx1_mat)
                            fx1_mask_3d = gtv_mask_cache(fx1_mat);
                        else
                            fx1_mask_3d = safe_load_mask(fx1_mat, 'Stvol3d');
                            gtv_mask_cache(fx1_mat) = fx1_mask_3d;
                        end
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

            % --- Texture features (when enabled) ---
            if use_texture && ~isempty(adc_vec) && has_3d_iter
                try
                    % Derive isotropic voxel spacing from voxel volume for shape features
                    if isfinite(vox_vol) && vox_vol > 0
                        iso_spacing = vox_vol^(1/3) * [1 1 1];
                    else
                        iso_spacing = [1 1 1];
                    end
                    % Pass quantization method from config (IBSI compliance)
                    if isfield(config_struct, 'texture_quantization_method')
                        tex_quant_method = config_struct.texture_quantization_method;
                    else
                        tex_quant_method = 'fixed_bin_number';
                    end
                    tex_3d = true;
                    if isfield(config_struct, 'texture_3d')
                        tex_3d = config_struct.texture_3d;
                    end
                    % Reconstruct the 3D parameter map from the 1D voxel
                    % vector and the 3D mask so that compute_texture_features
                    % receives spatially coherent data for GLCM/GLRLM.
                    adc_map_3d = NaN(size(gtv_mask_3d));
                    adc_map_3d(gtv_mask_3d) = adc_vec;
                    texture_features{j, k, dwi_type} = compute_texture_features(adc_map_3d, gtv_mask_3d, 32, iso_spacing, tex_quant_method, tex_3d);
                catch
                    % Texture extraction is non-fatal
                end
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
if isfield(config_struct, 'core_method')
    summary_metrics.core_method = config_struct.core_method;
else
    summary_metrics.core_method = 'adc_threshold';
end
summary_metrics.fdm_responding_pc = fdm_responding_pc;
summary_metrics.fdm_progressing_pc = fdm_progressing_pc;
summary_metrics.fdm_stable_pc = fdm_stable_pc;

% Package per-method core metrics when multi-method computation was enabled
if run_all_core
    summary_metrics.all_core_metrics = per_method;
end

% Package texture features when enabled
if use_texture
    summary_metrics.texture_features = texture_features;
end

% Repeat-Dice sanity summary so the log makes it obvious whether the
% shared-Fx1 fallback / path normalization actually succeeded this run.
try
    n_pat = size(summary_metrics.dice_rpt_adc, 1);
    n_dwi = size(summary_metrics.dice_rpt_adc, 2);
    diag_file = fullfile(config_struct.dataloc, 'repeat_dice_diagnostic.txt');
    diag_fid = fopen(diag_file, 'w');
    for dc = 1:n_dwi
        col = summary_metrics.dice_rpt_adc(:, dc);
        msg = sprintf('[REPEAT-DICE] dwi_col=%d: %d/%d patients with finite dice_rpt_adc', ...
            dc, sum(~isnan(col)), n_pat);
        fprintf('  %s\n', msg);
        if diag_fid > 0; fprintf(diag_fid, '%s\n', msg); end
    end
    if diag_fid > 0
        fprintf(diag_fid, 'Timestamp: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fclose(diag_fid);
        fprintf('  [REPEAT-DICE] Diagnostic written to %s\n', diag_file);
    end
    % Also copy the per-call trace file from tempdir (written by
    % compute_spatial_repeatability) into dataloc so the user can inspect
    % per-patient failure points (empty path / exist=false / voxel mismatch).
    trace_src = fullfile(tempdir, 'repeat_dice_trace.txt');
    if exist(trace_src, 'file') == 2
        trace_dst = fullfile(config_struct.dataloc, 'repeat_dice_trace.txt');
        try
            copyfile(trace_src, trace_dst);
            fprintf('  [REPEAT-DICE] Per-patient trace copied to %s\n', trace_dst);
        catch
        end
    end
catch
end

if isfield(config_struct, 'use_checkpoints') && config_struct.use_checkpoints
    fprintf('  [CHECKPOINT] Saving summary_metrics to %s...\n', summary_metrics_file);
    save(summary_metrics_file, 'summary_metrics');
end

end

% Local helper functions (nanmean_safe, nanstd_safe, compute_kurt_skew,
% compute_histogram_laplace) have been extracted to pipeline/utils/ as
% shared utilities used by compute_adc_metrics.m and compute_ivim_metrics.m.
% select_dwi_vectors is also a shared utility in utils/select_dwi_vectors.m.
