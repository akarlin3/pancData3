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

if nargin < 12, fx_dates = {}; end

if isfield(config_struct, 'dwi_type_name')
    file_prefix = ['_' config_struct.dwi_type_name];
else
    file_prefix = '';
end
summary_metrics_file = fullfile(config_struct.dataloc, ['summary_metrics' file_prefix '.mat']);
if isfield(config_struct, 'use_checkpoints') && config_struct.use_checkpoints
    if exist(summary_metrics_file, 'file')
        fprintf('  [CHECKPOINT] Found existing %s. Loading and skipping metrics computation...\n', ['summary_metrics' file_prefix '.mat']);
        tmp_metrics = load(summary_metrics_file, 'summary_metrics');
        summary_metrics = tmp_metrics.summary_metrics;
        return;
    else
        fallback_metrics_file = fullfile(config_struct.dataloc, 'summary_metrics.mat');
        if exist(fallback_metrics_file, 'file')
            fprintf('  [CHECKPOINT] Specific %s not found but fallback %s exists. Loading and skipping metrics computation...\n', ['summary_metrics' file_prefix '.mat'], 'summary_metrics.mat');
            tmp_metrics = load(fallback_metrics_file, 'summary_metrics');
            summary_metrics = tmp_metrics.summary_metrics;
            return;
        end
    end
end

% ADC threshold for identifying "restricted diffusion" sub-volume
adc_thresh = config_struct.adc_thresh;

% Secondary ADC threshold for identifying "high ADC" sub-volume
high_adc_thresh = config_struct.high_adc_thresh;

% IVIM thresholds for sub-volume identification
d_thresh = config_struct.d_thresh;
f_thresh = config_struct.f_thresh;

% Minimum voxel threshold for higher-order histogram metrics
min_vox_hist = config_struct.min_vox_hist;

nTp = size(data_vectors_gtvp, 2);   % number of timepoints (Fx1–Fx5 + post)
nRpt = size(data_vectors_gtvp, 3);  % max number of repeat scans at Fx1

% --- Pre-allocate sub-volume metric arrays (patient × timepoint × pipeline) ---
adc_sub_vol_pc = nan(length(id_list),nTp,3);
adc_sub_vol = nan(length(id_list),nTp,3);
adc_sub_mean = nan(length(id_list),nTp,3);
adc_sub_kurt = nan(length(id_list),nTp,3);
adc_sub_skew = nan(length(id_list),nTp,3);
high_adc_sub_vol = nan(length(id_list),nTp,3);
high_adc_sub_vol_pc = nan(length(id_list),nTp,3);
f_sub_vol = nan(length(id_list),nTp,3);
gtv_vol = nan(length(id_list),nTp);

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
bin_edges = 0:0.5e-4:3e-3;
adc_histograms = nan(length(id_list),nTp,length(bin_edges)-1,3);
d_histograms = nan(length(id_list),nTp,length(bin_edges)-1,3);
ks_stats_adc = nan(length(id_list),nTp,3);
ks_pvals_adc = nan(length(id_list),nTp,3);
ks_stats_d = nan(length(id_list),nTp,3);
ks_pvals_d = nan(length(id_list),nTp,3);

% --- Motion corruption flag ---
fx_corrupted = nan(length(id_list),nTp,3);
adc_max = config_struct.adc_max;

% --- Repeatability arrays (Fx1 repeat scans only, for wCV calculation) ---
adc_mean_rpt = nan(length(id_list),nRpt,3);
adc_sub_rpt = nan(length(id_list),nRpt,3);
fx_corrupted_rpt = nan(length(id_list),nRpt,3);
d_mean_rpt = nan(length(id_list),nRpt,3);
f_mean_rpt = nan(length(id_list),nRpt,3);
dstar_mean_rpt = nan(length(id_list),nRpt,3);

% --- Pooled voxel vectors across all patients (for population-level analysis) ---
n_rpt = nan(length(id_list),1);

% --- Main analysis loop: patient × timepoint × DWI pipeline ---
n_patients_metrics = length(id_list);
for j=1:n_patients_metrics
    text_progress_bar(j, n_patients_metrics, 'Computing summary metrics');
    for k=1:nTp
        if length(data_vectors_gtvp(j,k,1).vox_vol) == 1
            vox_vol = data_vectors_gtvp(j,k,1).vox_vol;
        else
            vox_vol = NaN;
        end
        for dwi_type = config_struct.dwi_types_to_run

            % Select the appropriate voxel vectors depending on pipeline
            switch dwi_type
                case 1  % Standard (raw DWI)
                    adc_vec = data_vectors_gtvp(j,k,1).adc_vector;
                    d_vec = data_vectors_gtvp(j,k,1).d_vector;
                    f_vec = data_vectors_gtvp(j,k,1).f_vector;
                    dstar_vec = data_vectors_gtvp(j,k,1).dstar_vector;

                    adc_baseline = data_vectors_gtvp(j,1,1).adc_vector;
                    d_baseline = data_vectors_gtvp(j,1,1).d_vector;
                case 2  % DnCNN-denoised + conventional IVIM fit
                    adc_vec = data_vectors_gtvp(j,k,1).adc_vector_dncnn;
                    d_vec = data_vectors_gtvp(j,k,1).d_vector_dncnn;
                    f_vec = data_vectors_gtvp(j,k,1).f_vector_dncnn;
                    dstar_vec = data_vectors_gtvp(j,k,1).dstar_vector_dncnn;

                    adc_baseline = data_vectors_gtvp(j,1,1).adc_vector_dncnn;
                    d_baseline = data_vectors_gtvp(j,1,1).d_vector_dncnn;
                case 3  % IVIMnet deep-learning IVIM fit (ADC uses standard pipeline)
                    adc_vec = data_vectors_gtvp(j,k,1).adc_vector;
                    d_vec = data_vectors_gtvp(j,k,1).d_vector_ivimnet;
                    f_vec = data_vectors_gtvp(j,k,1).f_vector_ivimnet;
                    dstar_vec = data_vectors_gtvp(j,k,1).dstar_vector_ivimnet;

                    adc_baseline = data_vectors_gtvp(j,1,1).adc_vector;
                    d_baseline = data_vectors_gtvp(j,1,1).d_vector_ivimnet;
            end

            % --- Compute ADC summary metrics for this patient/timepoint ---
            if ~isempty(adc_vec)
                n_finite_adc = sum(~isnan(adc_vec));
                % Only set gtv_vol from the first dwi_type to avoid
                % overwriting with pipeline-specific voxel counts.  GTV
                % volume is a geometric property independent of DWI type.
                % Use numel(adc_vec), not n_finite_adc, because the GTV
                % mask defines the anatomical volume — NaN voxels (from
                % failed ADC fits or artefacts) are still part of the
                % tumour contour and should contribute to volume.
                if isnan(gtv_vol(j,k))
                    gtv_vol(j,k) = numel(adc_vec)*vox_vol;
                end
                if exist('OCTAVE_VERSION', 'builtin')
                    % Octave nanmean might fail if it's a vector of all nans or something else
                    tmp = adc_vec(~isnan(adc_vec));
                    if isempty(tmp)
                        adc_mean(j,k,dwi_type) = NaN;
                    else
                        adc_mean(j,k,dwi_type) = mean(tmp);
                    end
                else
                    adc_mean(j,k,dwi_type) = nanmean(adc_vec);
                end
                if numel(adc_vec) >= min_vox_hist
                    adc_vec_finite = adc_vec(~isnan(adc_vec));
                    if numel(adc_vec_finite) >= min_vox_hist
                        adc_kurt(j,k,dwi_type) = kurtosis(adc_vec_finite);
                        adc_skew(j,k,dwi_type) = skewness(adc_vec_finite);
                    end
                end

                if exist('OCTAVE_VERSION', 'builtin')
                    tmp = adc_vec(~isnan(adc_vec));
                    if isempty(tmp)
                        adc_sd(j,k,dwi_type) = NaN;
                    else
                        adc_sd(j,k,dwi_type) = std(tmp);
                    end
                else
                    adc_sd(j,k,dwi_type) = nanstd(adc_vec);
                end
                
                adc_vec_sub = adc_vec(adc_vec<adc_thresh);
                adc_vec_high_sub = adc_vec(adc_vec>high_adc_thresh);

                adc_sub_vol(j,k,dwi_type) = numel(adc_vec_sub)*vox_vol;
                % Use count of finite (non-NaN) voxels as denominator so
                % NaN voxels do not artificially deflate the percentage.
                finite_vol = n_finite_adc * vox_vol;
                if finite_vol > 0
                    adc_sub_vol_pc(j,k,dwi_type) = adc_sub_vol(j,k,dwi_type)/finite_vol;
                else
                    adc_sub_vol_pc(j,k,dwi_type) = NaN;
                end
                if isempty(adc_vec_sub)
                    adc_sub_mean(j,k,dwi_type) = NaN;
                else
                    if exist('OCTAVE_VERSION', 'builtin')
                        tmp = adc_vec_sub(~isnan(adc_vec_sub));
                        if isempty(tmp)
                            adc_sub_mean(j,k,dwi_type) = NaN;
                        else
                            adc_sub_mean(j,k,dwi_type) = mean(tmp);
                        end
                    else
                        adc_sub_mean(j,k,dwi_type) = nanmean(adc_vec_sub);
                    end
                end
                if numel(adc_vec_sub) >= min_vox_hist
                    adc_sub_finite = adc_vec_sub(~isnan(adc_vec_sub));
                    if numel(adc_sub_finite) >= min_vox_hist
                        adc_sub_kurt(j,k,dwi_type) = kurtosis(adc_sub_finite);
                        adc_sub_skew(j,k,dwi_type) = skewness(adc_sub_finite);
                    end
                end

                if exist('OCTAVE_VERSION', 'builtin')
                    adc_vec_hist = adc_vec(~isnan(adc_vec));
                    c1 = histc(adc_vec_hist, bin_edges);
                    c1(end-1) = c1(end-1) + c1(end);  % merge last-edge count into final bin (match histcounts)
                    c1 = c1(1:end-1);
                else
                    [c1, ~] = histcounts(adc_vec, bin_edges);
                end
                n_binned_adc = sum(c1);
                nbins_adc = length(c1);
                if n_binned_adc > 0
                    % Laplace (add-one) smoothing avoids zero-probability
                    % bins that would cause log(0) in KL divergence or
                    % other distribution distance metrics.  Machine epsilon
                    % (eps) is too small and distorts such measures.
                    p1 = (c1 + 1) / (n_binned_adc + nbins_adc);
                else
                    p1 = zeros(size(c1));
                end
                adc_histograms(j,k,:,dwi_type) = p1;
                % NOTE: KS-test p-values are liberal because within-patient
                % voxels are spatially autocorrelated (violates independence).
                % Treat as descriptive, not inferential.
                if ~isempty(adc_baseline) && numel(adc_vec) >= min_vox_hist && numel(adc_baseline) >= min_vox_hist ...
                        && any(~isnan(adc_vec)) && any(~isnan(adc_baseline))
                    % Remove NaN before kstest2 (required for Octave compatibility)
                    adc_vec_ks = adc_vec(~isnan(adc_vec));
                    adc_bl_ks  = adc_baseline(~isnan(adc_baseline));
                    [~,p,ks2stat] = kstest2(adc_vec_ks, adc_bl_ks);
                    ks_stats_adc(j,k,dwi_type) = ks2stat;
                    ks_pvals_adc(j,k,dwi_type) = p;
                end

                high_adc_sub_vol(j,k,dwi_type) = numel(adc_vec_high_sub)*vox_vol;
                if finite_vol > 0
                    high_adc_sub_vol_pc(j,k,dwi_type) = high_adc_sub_vol(j,k,dwi_type)/finite_vol;
                else
                    high_adc_sub_vol_pc(j,k,dwi_type) = NaN;
                end
                if n_finite_adc > 0
                    fx_corrupted(j,k,dwi_type) = sum(adc_vec > adc_max & ~isnan(adc_vec)) / n_finite_adc;
                end
            end

            % --- Compute IVIM summary metrics (D, f, D*) ---
            if ~isempty(d_vec)
                % Exclude only exact-zero f values that co-occur with failed D*
                % fits (D* == 0 or NaN), preserving genuine zero perfusion.
                failed_fit = (f_vec == 0) & (isnan(dstar_vec) | dstar_vec == 0);
                f_vec(failed_fit) = nan;

                f_vec_sub = f_vec(f_vec<f_thresh);
                f_sub_vol(j,k,dwi_type) = numel(f_vec_sub)*vox_vol;
                d_vec_sub = d_vec(d_vec<d_thresh);

                if exist('OCTAVE_VERSION', 'builtin')
                    tmp = d_vec(~isnan(d_vec));
                    if isempty(tmp)
                        d_mean(j,k,dwi_type) = NaN;
                    else
                        d_mean(j,k,dwi_type) = mean(tmp);
                    end
                else
                    d_mean(j,k,dwi_type) = nanmean(d_vec);
                end
                if numel(d_vec) >= min_vox_hist
                    d_vec_finite = d_vec(~isnan(d_vec));
                    if numel(d_vec_finite) >= min_vox_hist
                        d_kurt(j,k,dwi_type) = kurtosis(d_vec_finite);
                        d_skew(j,k,dwi_type) = skewness(d_vec_finite);
                    end
                end

                if exist('OCTAVE_VERSION', 'builtin')
                    tmp = d_vec(~isnan(d_vec));
                    if isempty(tmp)
                        d_sd(j,k,dwi_type) = NaN;
                    else
                        d_sd(j,k,dwi_type) = std(tmp);
                    end
                else
                    d_sd(j,k,dwi_type) = nanstd(d_vec);
                end

                if exist('OCTAVE_VERSION', 'builtin')
                    d_vec_hist = d_vec(~isnan(d_vec));
                    c1 = histc(d_vec_hist, bin_edges);
                    c1(end-1) = c1(end-1) + c1(end);  % merge last-edge count into final bin (match histcounts)
                    c1 = c1(1:end-1);
                else
                    [c1, ~] = histcounts(d_vec, bin_edges);
                end
                n_binned_d = sum(c1);
                nbins_d = length(c1);
                if n_binned_d > 0
                    p1 = (c1 + 1) / (n_binned_d + nbins_d);
                else
                    p1 = zeros(size(c1));
                end
                d_histograms(j,k,:,dwi_type) = p1;
                % NOTE: KS-test p-values are liberal (see ADC comment above).
                if ~isempty(d_baseline) && numel(d_vec) >= min_vox_hist && numel(d_baseline) >= min_vox_hist ...
                        && any(~isnan(d_vec)) && any(~isnan(d_baseline))
                    % Remove NaN before kstest2 (consistent with ADC path)
                    d_vec_ks = d_vec(~isnan(d_vec));
                    d_bl_ks  = d_baseline(~isnan(d_baseline));
                    [~,p,ks2stat] = kstest2(d_vec_ks, d_bl_ks);
                    ks_stats_d(j,k,dwi_type) = ks2stat;
                    ks_pvals_d(j,k,dwi_type) = p;
                end

                if isempty(d_vec_sub)
                    d_sub_mean(j,k,dwi_type) = NaN;
                else
                    if exist('OCTAVE_VERSION', 'builtin')
                        tmp = d_vec_sub(~isnan(d_vec_sub));
                        if isempty(tmp)
                            d_sub_mean(j,k,dwi_type) = NaN;
                        else
                            d_sub_mean(j,k,dwi_type) = mean(tmp);
                        end
                    else
                        d_sub_mean(j,k,dwi_type) = nanmean(d_vec_sub);
                    end
                end
                if numel(d_vec_sub) >= min_vox_hist
                    d_sub_finite = d_vec_sub(~isnan(d_vec_sub));
                    if numel(d_sub_finite) >= min_vox_hist
                        d_sub_kurt(j,k,dwi_type) = kurtosis(d_sub_finite);
                        d_sub_skew(j,k,dwi_type) = skewness(d_sub_finite);
                    end
                end

                if exist('OCTAVE_VERSION', 'builtin')
                    tmp = f_vec(~isnan(f_vec));
                    if isempty(tmp)
                        f_mean(j,k,dwi_type) = NaN;
                    else
                        f_mean(j,k,dwi_type) = mean(tmp);
                    end
                else
                    f_mean(j,k,dwi_type) = nanmean(f_vec);
                end
                if numel(f_vec) >= min_vox_hist
                    f_vec_finite = f_vec(~isnan(f_vec));
                    if numel(f_vec_finite) >= min_vox_hist
                        f_kurt(j,k,dwi_type) = kurtosis(f_vec_finite);
                        f_skew(j,k,dwi_type) = skewness(f_vec_finite);
                    end
                end

                if exist('OCTAVE_VERSION', 'builtin')
                    tmp = dstar_vec(~isnan(dstar_vec));
                    if isempty(tmp)
                        dstar_mean(j,k,dwi_type) = NaN;
                    else
                        dstar_mean(j,k,dwi_type) = mean(tmp);
                    end
                else
                    dstar_mean(j,k,dwi_type) = nanmean(dstar_vec);
                end
                if numel(dstar_vec) >= min_vox_hist
                    dstar_vec_finite = dstar_vec(~isnan(dstar_vec));
                    if numel(dstar_vec_finite) >= min_vox_hist
                        dstar_kurt(j,k,dwi_type) = kurtosis(dstar_vec_finite);
                        dstar_skew(j,k,dwi_type) = skewness(dstar_vec_finite);
                    end
                end
            end

            % --- Repeatability analysis: extract metrics from Fx1 repeat scans ---
            if k==1
                rp_count = 0;
                for rpi=1:size(data_vectors_gtvp, 3)
                    switch dwi_type
                        case 1
                            adc_vec = data_vectors_gtvp(j,k,rpi).adc_vector;
                            d_vec = data_vectors_gtvp(j,k,rpi).d_vector;
                            f_vec = data_vectors_gtvp(j,k,rpi).f_vector;
                            dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector;
                        case 2
                            adc_vec = data_vectors_gtvp(j,k,rpi).adc_vector_dncnn;
                            d_vec = data_vectors_gtvp(j,k,rpi).d_vector_dncnn;
                            f_vec = data_vectors_gtvp(j,k,rpi).f_vector_dncnn;
                            dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector_dncnn;
                        case 3
                            adc_vec = data_vectors_gtvp(j,k,rpi).adc_vector;
                            d_vec = data_vectors_gtvp(j,k,rpi).d_vector_ivimnet;
                            f_vec = data_vectors_gtvp(j,k,rpi).f_vector_ivimnet;
                            dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector_ivimnet;
                    end

                    if ~isempty(adc_vec)
                        rp_count = rp_count+1;
                        if exist('OCTAVE_VERSION', 'builtin')
                            tmp = adc_vec(~isnan(adc_vec));
                            if isempty(tmp)
                                adc_mean_rpt(j,rpi,dwi_type) = NaN;
                            else
                                adc_mean_rpt(j,rpi,dwi_type) = mean(tmp);
                            end
                        else
                            adc_mean_rpt(j,rpi,dwi_type) = nanmean(adc_vec);
                        end
                        n_finite_rpt = sum(~isnan(adc_vec));
                        if n_finite_rpt > 0
                            fx_corrupted_rpt(j,rpi,dwi_type) = sum(adc_vec > adc_max & ~isnan(adc_vec)) / n_finite_rpt;
                        end
                        adc_vec_sub = adc_vec(adc_vec<adc_thresh);
                        if isempty(adc_vec_sub)
                            adc_sub_rpt(j,rpi,dwi_type) = NaN;
                        else
                            if exist('OCTAVE_VERSION', 'builtin')
                                tmp = adc_vec_sub(~isnan(adc_vec_sub));
                                if isempty(tmp)
                                    adc_sub_rpt(j,rpi,dwi_type) = NaN;
                                else
                                    adc_sub_rpt(j,rpi,dwi_type) = mean(tmp);
                                end
                            else
                                adc_sub_rpt(j,rpi,dwi_type) = nanmean(adc_vec_sub);
                            end
                        end
                    end

                    if ~isempty(d_vec)
                        % Count D/f/D* repeats when ADC is absent (e.g.,
                        % IVIMnet which reuses standard-pipeline ADC).
                        if isempty(adc_vec)
                            rp_count = rp_count + 1;
                        end
                        if exist('OCTAVE_VERSION', 'builtin')
                            tmp = d_vec(~isnan(d_vec));
                            if isempty(tmp)
                                d_mean_rpt(j,rpi,dwi_type) = NaN;
                            else
                                d_mean_rpt(j,rpi,dwi_type) = mean(tmp);
                            end
                            tmp = f_vec(~isnan(f_vec));
                            if isempty(tmp)
                                f_mean_rpt(j,rpi,dwi_type) = NaN;
                            else
                                f_mean_rpt(j,rpi,dwi_type) = mean(tmp);
                            end
                            tmp = dstar_vec(~isnan(dstar_vec));
                            if isempty(tmp)
                                dstar_mean_rpt(j,rpi,dwi_type) = NaN;
                            else
                                dstar_mean_rpt(j,rpi,dwi_type) = mean(tmp);
                            end
                        else
                            d_mean_rpt(j,rpi,dwi_type) = nanmean(d_vec);
                            f_mean_rpt(j,rpi,dwi_type) = nanmean(f_vec);
                            dstar_mean_rpt(j,rpi,dwi_type) = nanmean(dstar_vec);
                        end
                    end
                end
                % Record repeat count from any DWI type (first non-zero wins).
                % Previously only dwi_type==1 populated n_rpt, leaving it
                % all-NaN for dnCNN/IVIMnet-only runs and breaking wCV.
                if isnan(n_rpt(j)) || n_rpt(j) == 0
                    n_rpt(j) = rp_count;
                end
            end
        end
    end
end

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
summary_metrics.d_mean_rpt = d_mean_rpt;
summary_metrics.f_mean_rpt = f_mean_rpt;
summary_metrics.dstar_mean_rpt = dstar_mean_rpt;
summary_metrics.n_rpt = n_rpt;
summary_metrics.dmean_gtvp = dmean_gtvp;
summary_metrics.gtv_locations = gtv_locations;
summary_metrics.dwi_locations = dwi_locations;
summary_metrics.fx_dates = fx_dates;

if isfield(config_struct, 'use_checkpoints') && config_struct.use_checkpoints
    fprintf('  [CHECKPOINT] Saving summary_metrics to %s...\n', summary_metrics_file);
    save(summary_metrics_file, 'summary_metrics');
end

end
