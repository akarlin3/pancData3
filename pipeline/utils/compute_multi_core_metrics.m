function per_method = compute_multi_core_metrics(per_method, config_struct, ...
    ALL_CORE_METHODS, adc_vec, d_vec, f_vec, dstar_vec, ...
    has_3d_iter, gtv_mask_3d, core_opts, ...
    j, k, dwi_type, vox_vol, min_vox_hist, ...
    f_thresh, d_thresh, finite_vol, store_masks)
% COMPUTE_MULTI_CORE_METRICS — Computes sub-volume metrics for configurable core methods
%
% Loops over all core delineation methods specified in config, calls extract_tumor_core for
% each, and computes per-method ADC/D/f sub-volume metrics plus optional
% fDM volume fractions. Optionally stores core masks.
%
% Inputs:
%   per_method       - Struct of pre-allocated per-method metric arrays (modified in-place)
%   config_struct    - Configuration struct with thresholds and settings
%   ALL_CORE_METHODS - Cell array of all method name strings (kept for backward compatibility)
%   adc_vec          - ADC voxel vector for this patient/timepoint
%   d_vec            - D voxel vector (may be empty)
%   f_vec            - f voxel vector (may be empty)
%   dstar_vec        - D* voxel vector (may be empty)
%   has_3d_iter      - Boolean: whether valid 3D mask is available
%   gtv_mask_3d      - 3D GTV mask (may be empty)
%   core_opts        - Struct with timepoint_index and optional baseline vectors
%   j                - Patient index
%   k                - Timepoint index
%   dwi_type         - DWI pipeline type index (1=Standard, 2=DnCNN, 3=IVIMnet)
%   vox_vol          - Single voxel volume in cm^3
%   min_vox_hist     - Minimum voxels for kurtosis/skewness computation
%   f_thresh         - f threshold for low-perfusion sub-volume
%   d_thresh         - D threshold for restricted diffusion sub-volume
%   finite_vol       - Finite (non-NaN) ADC volume in cm^3
%   store_masks      - Boolean: whether to store core masks in per_method
%
% Outputs:
%   per_method       - Updated struct with per-method metrics filled in for (j,k,dwi_type)

% Get core methods from config, with fallback to ALL_CORE_METHODS for backward compatibility
if isfield(config_struct, 'active_core_methods') && ~isempty(config_struct.active_core_methods)
    core_methods = config_struct.active_core_methods;
elseif isfield(config_struct, 'core_methods') && ~isempty(config_struct.core_methods)
    core_methods = config_struct.core_methods;
else
    core_methods = ALL_CORE_METHODS;
end

n_all_methods = numel(core_methods);
% "Unified" methods produce a single core mask applied to all parameters,
% as opposed to parameter-specific thresholding (where D, f, D* each get
% independent sub-volumes). This matters for computing D/f sub-volume
% metrics below: unified methods use the shared core mask, while others
% apply per-parameter thresholds independently.
unified_set = {'percentile', 'spectral', 'fdm'};

rng(42);  % reproducible clustering (GMM, k-means, spectral)
% Suppress expected warnings from extract_tumor_core when methods cannot
% run (e.g., too few voxels for spectral clustering, no 3D mask for
% active contours). These are expected in multi-method comparison runs
% where not every method can succeed for every patient.
prev_warn_csm = warning('query');
warning('off', 'extract_tumor_core:tooFewForSpectral');
warning('off', 'extract_tumor_core:no3DForActiveContours');
warning('off', 'extract_tumor_core:no3DForRegionGrowing');
warning('off', 'extract_tumor_core:fdmBaseline');
warning('off', 'extract_tumor_core:fdmNoBaseline');
warning('off', 'extract_tumor_core:noDValues');
warning('off', 'extract_tumor_core:noIVIMValues');
warning('off', 'extract_tumor_core:noSpectralCluster');

for m_idx = 1:n_all_methods
    mname = core_methods{m_idx};
    % Create a temporary config with this method name so extract_tumor_core
    % uses the correct delineation algorithm while preserving all other settings
    temp_cfg = config_struct;
    temp_cfg.core_method = mname;

    try
        % extract_tumor_core returns a logical mask (same length as adc_vec)
        % identifying voxels belonging to the tumor core sub-volume
        m_mask = extract_tumor_core(temp_cfg, adc_vec, d_vec, f_vec, dstar_vec, has_3d_iter, gtv_mask_3d, core_opts);
    catch
        % If the method fails (e.g., insufficient voxels), default to
        % an all-false mask so downstream metrics are NaN rather than erroring
        m_mask = false(size(adc_vec));
    end

    % Optionally store the binary core mask for downstream analysis
    % (e.g., Dice/Hausdorff comparison between methods via compare_core_methods)
    if store_masks
        per_method.(mname).core_masks{j, k} = m_mask;
    end

    % --- ADC sub-volume metrics ---
    % Extract ADC values within the core mask and compute summary statistics.
    % These capture the restricted-diffusion characteristics of the tumor core,
    % which correlates with cellularity and treatment resistance in pancreatic cancer.
    m_adc_sub = adc_vec(m_mask);
    per_method.(mname).adc_sub_vol(j,k,dwi_type) = sum(m_mask) * vox_vol;  % volume in cm^3
    if finite_vol > 0
        % Sub-volume as a fraction of total finite GTV volume
        per_method.(mname).adc_sub_vol_pc(j,k,dwi_type) = per_method.(mname).adc_sub_vol(j,k,dwi_type) / finite_vol;
    end
    if ~isempty(m_adc_sub)
        per_method.(mname).adc_sub_mean(j,k,dwi_type) = nanmean(m_adc_sub);
    end
    % Higher-order statistics (kurtosis, skewness) require a minimum number
    % of voxels to be meaningful; below min_vox_hist they are too noisy.
    % Kurtosis captures the "peakedness" of the ADC distribution (high
    % kurtosis = more extreme values), while skewness indicates asymmetry.
    if numel(m_adc_sub) >= min_vox_hist
        m_adc_finite = m_adc_sub(~isnan(m_adc_sub));
        if numel(m_adc_finite) >= min_vox_hist
            per_method.(mname).adc_sub_kurt(j,k,dwi_type) = kurtosis(m_adc_finite);
            per_method.(mname).adc_sub_skew(j,k,dwi_type) = skewness(m_adc_finite);
        end
    end

    % --- D/f sub-volume metrics (IVIM parameters) ---
    % D (true diffusion coefficient) and f (perfusion fraction) are IVIM
    % bi-exponential model parameters that separate tissue cellularity (D)
    % from microvascular perfusion (f). Only computed when IVIM data exists.
    if ~isempty(d_vec)
        if any(strcmpi(mname, unified_set))
            % Unified methods: use the same core mask for D and f sub-volumes.
            % This ensures spatial consistency — the "core" is the same region
            % regardless of which parameter is being summarized.
            m_f_sub = f_vec(m_mask); %#ok<NASGU>
            per_method.(mname).f_sub_vol(j,k,dwi_type) = sum(m_mask) * vox_vol;
            m_d_sub = d_vec(m_mask);
        else
            % Non-unified methods: apply independent per-parameter thresholds.
            % f < f_thresh identifies low-perfusion voxels (ischemic core),
            % D < d_thresh identifies restricted-diffusion voxels (high cellularity).
            m_f_sub = f_vec(f_vec < f_thresh); %#ok<NASGU>
            per_method.(mname).f_sub_vol(j,k,dwi_type) = numel(f_vec(f_vec < f_thresh)) * vox_vol;
            m_d_sub = d_vec(d_vec < d_thresh);
        end
        if ~isempty(m_d_sub)
            per_method.(mname).d_sub_mean(j,k,dwi_type) = nanmean(m_d_sub);
        end
        if numel(m_d_sub) >= min_vox_hist
            m_d_finite = m_d_sub(~isnan(m_d_sub));
            if numel(m_d_finite) >= min_vox_hist
                per_method.(mname).d_sub_kurt(j,k,dwi_type) = kurtosis(m_d_finite);
                per_method.(mname).d_sub_skew(j,k,dwi_type) = skewness(m_d_finite);
            end
        end
    end

    % --- fDM (functional Diffusion Map) volume fractions ---
    % fDM classifies each voxel as responding, stable, or progressing based
    % on the change in a diffusion parameter relative to baseline. This is
    % a per-voxel longitudinal analysis technique (Moffat et al., PNAS 2005)
    % that captures spatial heterogeneity in treatment response.
    % Only computed for the 'fdm' method and post-baseline timepoints (k>1).
    if strcmpi(mname, 'fdm') && k > 1 && isfield(core_opts, 'baseline_adc_vec')
        % Select the diffusion parameter for fDM analysis (ADC or D)
        switch lower(config_struct.fdm_parameter)
            case 'adc'
                fdm_cur = adc_vec;
                fdm_bl = core_opts.baseline_adc_vec;
            case 'd'
                fdm_cur = d_vec;
                fdm_bl = core_opts.baseline_d_vec;
            otherwise
                fdm_cur = adc_vec;
                fdm_bl = core_opts.baseline_adc_vec;
        end
        % Require matched voxel counts between baseline and current scan
        if ~isempty(fdm_bl) && numel(fdm_bl) == numel(fdm_cur)
            fdm_sig_pm = config_struct.fdm_thresh;  % significance threshold for voxel change
            delta_pm = fdm_cur - fdm_bl;  % per-voxel change from baseline
            valid_pm = ~isnan(delta_pm);
            n_valid_pm = sum(valid_pm);
            if n_valid_pm > 0
                % Classify voxels: responding (increased diffusion = cell death),
                % progressing (decreased diffusion = increased cellularity),
                % stable (change within noise threshold)
                per_method.(mname).fdm_responding_pc(j,k,dwi_type)  = sum(delta_pm(valid_pm) > fdm_sig_pm) / n_valid_pm;
                per_method.(mname).fdm_progressing_pc(j,k,dwi_type) = sum(delta_pm(valid_pm) < -fdm_sig_pm) / n_valid_pm;
                per_method.(mname).fdm_stable_pc(j,k,dwi_type)      = sum(abs(delta_pm(valid_pm)) <= fdm_sig_pm) / n_valid_pm;
            end
        end
    end
end

% Restore original warning state and clear stale warnings so the
% orchestrator doesn't re-log suppressed ones from extract_tumor_core
warning(prev_warn_csm);
lastwarn('');

end