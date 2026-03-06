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

if isfield(config_struct, 'dwi_type_name')
    file_prefix = ['_' config_struct.dwi_type_name];
else
    file_prefix = '';
end
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
% calculate_subvolume_metrics.m).
if exist('OCTAVE_VERSION', 'builtin')
    sphere_kernel = zeros(3, 3, 3);
    sphere_kernel(2,2,:) = 1; sphere_kernel(2,:,2) = 1; sphere_kernel(:,2,2) = 1;
    morph_se = strel('arbitrary', sphere_kernel);
else
    morph_se = strel('sphere', 1);
end
morph_min_cc = 10;  % minimum connected component size (voxels)

% Cache for GTV mask loading (avoids repeated disk I/O for same file)
last_rpt_gtv_mat = '';
last_rpt_gtv_mask = [];

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
            % ADC aggregation follows a hierarchy: whole-GTV statistics first,
            % then sub-volume decomposition. This captures both the overall
            % tumor diffusion state and the spatial heterogeneity within it.
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
                
                adc_vec_sub = adc_vec(adc_vec<=adc_thresh);
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
                % Skip k==1: adc_vec and adc_baseline are the same data,
                % so the KS test trivially returns stat=0, p=1.
                if k > 1 && ~isempty(adc_baseline) && numel(adc_vec) >= min_vox_hist && numel(adc_baseline) >= min_vox_hist ...
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
            % IVIM parameters provide biologically specific information beyond ADC:
            %   D (true diffusion): cellularity without perfusion contamination
            %   f (perfusion fraction): microvasculature density
            %   D* (pseudo-diffusion): capillary blood flow velocity
            % D* is the least reliable parameter due to its rapid signal
            % decay (only measurable at very low b-values) and high noise
            % sensitivity. DnCNN denoising and IVIMnet are specifically
            % aimed at improving D* estimation.
            if ~isempty(d_vec)
                % Exclude exact-zero f values that co-occur with failed D*
                % fits (D* == 0 or NaN), preserving genuine zero perfusion.
                % The segmented IVIM fitter returns f=0, D*=0 when the
                % perfusion component cannot be separated from noise.
                % Both f AND D* must be set to NaN for failed fits:
                % leaving D*=0 in the data biases dstar_mean downward,
                % and leaving f=0 biases f_mean downward.
                failed_fit = (f_vec == 0) & (isnan(dstar_vec) | dstar_vec == 0);
                f_vec(failed_fit) = nan;
                dstar_vec(failed_fit) = nan;

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
                % Skip k==1: d_vec and d_baseline are the same data.
                if k > 1 && ~isempty(d_baseline) && numel(d_vec) >= min_vox_hist && numel(d_baseline) >= min_vox_hist ...
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

                    % Apply the same failed-fit filter used in the main
                    % metrics path (lines 297-299) so that repeatability
                    % wCV for f and D* is not biased by spurious zeros.
                    failed_fit_rpt = (f_vec == 0) & (isnan(dstar_vec) | dstar_vec == 0);
                    f_vec(failed_fit_rpt) = nan;
                    dstar_vec(failed_fit_rpt) = nan;

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
                        adc_vec_sub = adc_vec(adc_vec<=adc_thresh);
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

                % --- Spatial repeatability: Dice & Hausdorff between repeat sub-volumes ---
                % Compare threshold-defined sensitive sub-volumes across all
                % pairs of Fx1 repeat scans to assess spatial reproducibility.
                if rp_count >= 2
                    % Determine voxel dimensions for physical Hausdorff distances
                    rpt_vox_dims = data_vectors_gtvp(j,1,1).vox_dims;
                    if isempty(rpt_vox_dims) || ~isnumeric(rpt_vox_dims) || numel(rpt_vox_dims) ~= 3
                        rpt_vox_vol = data_vectors_gtvp(j,1,1).vox_vol;
                        if ~isempty(rpt_vox_vol) && ~isnan(rpt_vox_vol) && rpt_vox_vol > 0
                            side_mm = (rpt_vox_vol * 1000) ^ (1/3);
                            rpt_vox_dims = [side_mm, side_mm, side_mm];
                        else
                            rpt_vox_dims = [1, 1, 1];
                        end
                    end

                    % Collect valid repeat indices and their parameter vectors
                    valid_rpis = [];
                    rpt_vecs = struct('adc', {{}}, 'd', {{}}, 'f', {{}}, 'dstar', {{}});
                    for rpi2 = 1:size(data_vectors_gtvp, 3)
                        switch dwi_type
                            case 1
                                rv_adc = data_vectors_gtvp(j,1,rpi2).adc_vector;
                                rv_d   = data_vectors_gtvp(j,1,rpi2).d_vector;
                                rv_f   = data_vectors_gtvp(j,1,rpi2).f_vector;
                                rv_ds  = data_vectors_gtvp(j,1,rpi2).dstar_vector;
                            case 2
                                rv_adc = data_vectors_gtvp(j,1,rpi2).adc_vector_dncnn;
                                rv_d   = data_vectors_gtvp(j,1,rpi2).d_vector_dncnn;
                                rv_f   = data_vectors_gtvp(j,1,rpi2).f_vector_dncnn;
                                rv_ds  = data_vectors_gtvp(j,1,rpi2).dstar_vector_dncnn;
                            case 3
                                rv_adc = data_vectors_gtvp(j,1,rpi2).adc_vector;
                                rv_d   = data_vectors_gtvp(j,1,rpi2).d_vector_ivimnet;
                                rv_f   = data_vectors_gtvp(j,1,rpi2).f_vector_ivimnet;
                                rv_ds  = data_vectors_gtvp(j,1,rpi2).dstar_vector_ivimnet;
                        end
                        if ~isempty(rv_adc) || ~isempty(rv_d)
                            valid_rpis(end+1) = rpi2; %#ok<AGROW>
                            rpt_vecs.adc{end+1} = rv_adc;
                            rpt_vecs.d{end+1} = rv_d;
                            rpt_vecs.f{end+1} = rv_f;
                            rpt_vecs.dstar{end+1} = rv_ds;
                        end
                    end

                    if numel(valid_rpis) >= 2
                        % Load 3D GTV masks for each valid repeat
                        rpt_masks_3d = cell(numel(valid_rpis), 1);
                        rpt_has_3d = true;
                        for ri = 1:numel(valid_rpis)
                            rpi_idx = valid_rpis(ri);
                            gtv_mat_path = gtv_locations{j, 1, rpi_idx};
                            if ~isempty(gtv_mat_path)
                                path_parts = strsplit(gtv_mat_path, {'/', '\'});
                                gtv_mat_path = fullfile(path_parts{:});
                                if isunix && ~startsWith(gtv_mat_path, filesep) && isempty(path_parts{1})
                                    gtv_mat_path = [filesep gtv_mat_path]; %#ok<AGROW>
                                end
                                if exist(gtv_mat_path, 'file')
                                    if strcmp(gtv_mat_path, last_rpt_gtv_mat)
                                        rpt_masks_3d{ri} = last_rpt_gtv_mask;
                                    else
                                        rpt_masks_3d{ri} = safe_load_mask(gtv_mat_path, 'Stvol3d');
                                        last_rpt_gtv_mat = gtv_mat_path;
                                        last_rpt_gtv_mask = rpt_masks_3d{ri};
                                    end
                                end
                            end
                            if isempty(rpt_masks_3d{ri})
                                rpt_has_3d = false;
                            end
                        end

                        if rpt_has_3d
                            % Accumulate pairwise metrics then average
                            param_names = {'adc', 'd', 'f', 'dstar'};
                            param_thresholds = [adc_thresh, d_thresh, f_thresh, config_struct.dstar_thresh];
                            pair_dice = nan(numel(valid_rpis)*(numel(valid_rpis)-1)/2, 4);
                            pair_hd_max = nan(size(pair_dice));
                            pair_hd95 = nan(size(pair_dice));
                            pi_count = 0;

                            for ri1 = 1:numel(valid_rpis)-1
                                for ri2 = ri1+1:numel(valid_rpis)
                                    pi_count = pi_count + 1;
                                    mask_3d_1 = rpt_masks_3d{ri1};
                                    mask_3d_2 = rpt_masks_3d{ri2};

                                    % Verify compatible 3D mask dimensions
                                    if ~isequal(size(mask_3d_1), size(mask_3d_2))
                                        continue;
                                    end

                                    for pi = 1:4
                                        vec_1 = rpt_vecs.(param_names{pi}){ri1};
                                        vec_2 = rpt_vecs.(param_names{pi}){ri2};
                                        if isempty(vec_1) || isempty(vec_2)
                                            continue;
                                        end

                                        n_gtv_1 = sum(mask_3d_1(:) == 1);
                                        n_gtv_2 = sum(mask_3d_2(:) == 1);
                                        if numel(vec_1) ~= n_gtv_1 || numel(vec_2) ~= n_gtv_2
                                            continue;
                                        end

                                        % Threshold to binary, embed in 3D, morphological cleanup
                                        subvol_3d_1 = false(size(mask_3d_1));
                                        subvol_3d_1(mask_3d_1 == 1) = vec_1 < param_thresholds(pi);
                                        subvol_3d_1 = imclose(imopen(subvol_3d_1, morph_se), morph_se);
                                        subvol_3d_1 = bwareaopen(subvol_3d_1, morph_min_cc);

                                        subvol_3d_2 = false(size(mask_3d_2));
                                        subvol_3d_2(mask_3d_2 == 1) = vec_2 < param_thresholds(pi);
                                        subvol_3d_2 = imclose(imopen(subvol_3d_2, morph_se), morph_se);
                                        subvol_3d_2 = bwareaopen(subvol_3d_2, morph_min_cc);

                                        [d_val, hm_val, h95_val] = compute_dice_hausdorff( ...
                                            subvol_3d_1, subvol_3d_2, rpt_vox_dims);
                                        pair_dice(pi_count, pi) = d_val;
                                        pair_hd_max(pi_count, pi) = hm_val;
                                        pair_hd95(pi_count, pi) = h95_val;
                                    end
                                end
                            end

                            % Average across all repeat pairs
                            if exist('OCTAVE_VERSION', 'builtin')
                                mean_dice = mean(pair_dice(1:pi_count,:), 1, 'omitnan');
                                mean_hd_max = mean(pair_hd_max(1:pi_count,:), 1, 'omitnan');
                                mean_hd95 = mean(pair_hd95(1:pi_count,:), 1, 'omitnan');
                            else
                                mean_dice = nanmean(pair_dice(1:pi_count,:), 1);
                                mean_hd_max = nanmean(pair_hd_max(1:pi_count,:), 1);
                                mean_hd95 = nanmean(pair_hd95(1:pi_count,:), 1);
                            end
                            dice_rpt_adc(j, dwi_type) = mean_dice(1);
                            dice_rpt_d(j, dwi_type)   = mean_dice(2);
                            dice_rpt_f(j, dwi_type)   = mean_dice(3);
                            dice_rpt_dstar(j, dwi_type) = mean_dice(4);
                            hd_max_rpt_adc(j, dwi_type) = mean_hd_max(1);
                            hd_max_rpt_d(j, dwi_type)   = mean_hd_max(2);
                            hd_max_rpt_f(j, dwi_type)   = mean_hd_max(3);
                            hd_max_rpt_dstar(j, dwi_type) = mean_hd_max(4);
                            hd95_rpt_adc(j, dwi_type) = mean_hd95(1);
                            hd95_rpt_d(j, dwi_type)   = mean_hd95(2);
                            hd95_rpt_f(j, dwi_type)   = mean_hd95(3);
                            hd95_rpt_dstar(j, dwi_type) = mean_hd95(4);
                        end
                    end
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

if isfield(config_struct, 'use_checkpoints') && config_struct.use_checkpoints
    fprintf('  [CHECKPOINT] Saving summary_metrics to %s...\n', summary_metrics_file);
    save(summary_metrics_file, 'summary_metrics');
end

end
