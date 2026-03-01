function summary_metrics = compute_summary_metrics(config_struct, data_vectors_gtvp, id_list, mrn_list, lf, immuno, gtv_locations, dwi_locations, dmean_gtvp, d95_gtvp, v50gy_gtvp)
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
%
% Outputs:
%   summary_metrics   - Struct containing computed mean, kurtosis, skewness,
%                       SD, subsets, histogram features, predictability, etc.
%

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

nTp = 6;   % number of timepoints (Fx1–Fx5 + post)
nRpt = 6;  % max number of repeat scans at Fx1

% --- Pre-allocate sub-volume metric arrays (patient × timepoint × pipeline) ---
adc_sub_vol_pc = nan(length(id_list),nTp,3);
adc_sub_vol = nan(length(id_list),nTp,3);
adc_sub_mean = nan(length(id_list),nTp,3);
adc_sub_kurt = nan(length(id_list),nTp,3);
adc_sub_skew = nan(length(id_list),nTp,3);
high_adc_sub_vol = nan(length(id_list),nTp,3);
ivim_sub_vol = nan(length(id_list),nTp,3);
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
for j=1:length(id_list)
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
                gtv_vol(j,k) = numel(adc_vec)*vox_vol;
                adc_mean(j,k,dwi_type) = nanmean(adc_vec);
                if numel(adc_vec) >= min_vox_hist
                    adc_kurt(j,k,dwi_type) = kurtosis(adc_vec);
                    adc_skew(j,k,dwi_type) = skewness(adc_vec);
                end
                adc_sd(j,k,dwi_type) = nanstd(adc_vec);
                
                adc_vec_sub = adc_vec(adc_vec<adc_thresh);
                adc_vec_high_sub = adc_vec(adc_vec>high_adc_thresh);

                adc_sub_vol(j,k,dwi_type) = numel(adc_vec_sub)*vox_vol;
                adc_sub_vol_pc(j,k,dwi_type) = adc_sub_vol(j,k,dwi_type)/gtv_vol(j,k);
                adc_sub_mean(j,k,dwi_type) = nanmean(adc_vec_sub);
                if numel(adc_vec_sub) >= min_vox_hist
                    adc_sub_kurt(j,k,dwi_type) = kurtosis(adc_vec_sub);
                    adc_sub_skew(j,k,dwi_type) = skewness(adc_vec_sub);
                end

                [c1, ~] = histcounts(adc_vec, bin_edges);
                p1 = c1 / numel(adc_vec); p1(p1==0)=eps;
                adc_histograms(j,k,:,dwi_type) = p1;
                if ~isempty(adc_baseline) && numel(adc_vec) >= min_vox_hist && numel(adc_baseline) >= min_vox_hist
                    [~,p,ks2stat] = kstest2(adc_vec,adc_baseline);
                    ks_stats_adc(j,k,dwi_type) = ks2stat;
                    ks_pvals_adc(j,k,dwi_type) = p;
                end

                high_adc_sub_vol(j,k,dwi_type) = numel(adc_vec_high_sub)*vox_vol;
                fx_corrupted(j,k,dwi_type) = numel(adc_vec(adc_vec>adc_max))/numel(adc_vec);
            end

            % --- Compute IVIM summary metrics (D, f, D*) ---
            if ~isempty(d_vec)
                f_vec(f_vec==0) = nan;  % zero perfusion fraction = failed fit

                ivim_vec_sub = d_vec(d_vec<d_thresh & f_vec<f_thresh);
                ivim_sub_vol(j,k,dwi_type) = numel(ivim_vec_sub)*vox_vol;
                d_vec_sub = d_vec(adc_vec<adc_thresh);

                d_mean(j,k,dwi_type) = nanmean(d_vec);
                if numel(d_vec) >= min_vox_hist
                    d_kurt(j,k,dwi_type) = kurtosis(d_vec);
                    d_skew(j,k,dwi_type) = skewness(d_vec);
                end
                d_sd(j,k,dwi_type) = nanstd(d_vec);

                [c1, ~] = histcounts(d_vec, bin_edges);
                p1 = c1 / numel(d_vec); p1(p1==0)=eps;
                d_histograms(j,k,:,dwi_type) = p1;
                if ~isempty(d_baseline) && numel(d_vec) >= min_vox_hist && numel(d_baseline) >= min_vox_hist
                    [~,p,ks2stat] = kstest2(d_vec,d_baseline);
                    ks_stats_d(j,k,dwi_type) = ks2stat;
                    ks_pvals_d(j,k,dwi_type) = p;
                end

                d_sub_mean(j,k,dwi_type) = nanmean(d_vec_sub);
                if numel(d_vec_sub) >= min_vox_hist
                    d_sub_kurt(j,k,dwi_type) = kurtosis(d_vec_sub);
                    d_sub_skew(j,k,dwi_type) = skewness(d_vec_sub);
                end

                f_mean(j,k,dwi_type) = nanmean(f_vec);
                if numel(f_vec) >= min_vox_hist
                    f_kurt(j,k,dwi_type) = kurtosis(f_vec);
                    f_skew(j,k,dwi_type) = skewness(f_vec);
                end

                dstar_mean(j,k,dwi_type) = nanmean(dstar_vec);
                if numel(dstar_vec) >= min_vox_hist
                    dstar_kurt(j,k,dwi_type) = kurtosis(dstar_vec);
                    dstar_skew(j,k,dwi_type) = skewness(dstar_vec);
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
                            adc_vec = [];
                            d_vec = data_vectors_gtvp(j,k,rpi).d_vector_ivimnet;
                            f_vec = data_vectors_gtvp(j,k,rpi).f_vector_ivimnet;
                            dstar_vec = data_vectors_gtvp(j,k,rpi).dstar_vector_ivimnet;
                    end

                    if ~isempty(adc_vec)
                        rp_count = rp_count+1;
                        adc_mean_rpt(j,rpi,dwi_type) = nanmean(adc_vec);
                        fx_corrupted_rpt(j,rpi,dwi_type) = numel(adc_vec(adc_vec>adc_max))/numel(adc_vec);
                        adc_vec_sub = adc_vec(adc_vec<adc_thresh);
                        adc_sub_rpt(j,rpi,dwi_type) = nanmean(adc_vec_sub);
                    end

                    if ~isempty(d_vec)
                        d_mean_rpt(j,rpi,dwi_type) = nanmean(d_vec);
                        f_mean_rpt(j,rpi,dwi_type) = nanmean(f_vec);
                        dstar_mean_rpt(j,rpi,dwi_type) = nanmean(dstar_vec);
                    end
                end
                if dwi_type==1
                    n_rpt(j) = rp_count;
                end
            end
        end
    end
end

summary_metrics = struct('adc_mean', adc_mean, 'adc_kurt', adc_kurt, 'adc_skew', adc_skew, 'adc_sd', adc_sd, 'd_mean', d_mean, 'f_mean', f_mean, 'dstar_mean', dstar_mean, ...
    'ivim_sub_vol', ivim_sub_vol, ...
    'adc_sub_vol', adc_sub_vol, 'adc_sub_vol_pc', adc_sub_vol_pc, 'high_adc_sub_vol', high_adc_sub_vol, 'd_kurt', d_kurt, 'd_skew', d_skew, 'd_sd', d_sd, ...
    'f_kurt', f_kurt, 'f_skew', f_skew, 'dstar_kurt', dstar_kurt, 'dstar_skew', dstar_skew, 'd_sub_mean', d_sub_mean, 'd_sub_kurt', d_sub_kurt, ...
    'd_sub_skew', d_sub_skew, 'adc_histograms', adc_histograms, 'd_histograms', d_histograms, 'ks_stats_adc', ks_stats_adc, 'ks_pvals_adc', ks_pvals_adc, ...
    'ks_stats_d', ks_stats_d, 'ks_pvals_d', ks_pvals_d, 'fx_corrupted', fx_corrupted, 'gtv_vol', gtv_vol, ...
    'id_list', {id_list}, 'mrn_list', {mrn_list}, 'd95_gtvp', d95_gtvp, 'v50gy_gtvp', v50gy_gtvp, 'lf', lf, 'immuno', immuno, ...
    'adc_mean_rpt', adc_mean_rpt, 'adc_sub_rpt', adc_sub_rpt, 'd_mean_rpt', d_mean_rpt, 'f_mean_rpt', f_mean_rpt, 'dstar_mean_rpt', dstar_mean_rpt, ...
    'n_rpt', n_rpt, 'dmean_gtvp', dmean_gtvp, 'gtv_locations', {gtv_locations}, 'dwi_locations', {dwi_locations});

if isfield(config_struct, 'use_checkpoints') && config_struct.use_checkpoints
    fprintf('  [CHECKPOINT] Saving summary_metrics to %s...\n', summary_metrics_file);
    save(summary_metrics_file, 'summary_metrics');
end

end
