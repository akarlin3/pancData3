function config_struct = parse_config(json_path)
% PARSE_CONFIG Reads a JSON config file into a MATLAB struct.
%
% Usage:
%   config_struct = parse_config('config.json')
%
% Inputs:
%   json_path - Path to the config file (e.g., config.json)
%
% Outputs:
%   config_struct - A struct containing all parsed fields with defaults populated
%
% Analytical Notes:
%   This function centralizes all configuration for the pancreatic DWI
%   analysis pipeline.  Every threshold and parameter default encodes a
%   domain-specific decision rooted in the physics of diffusion-weighted
%   MRI and the clinical characteristics of pancreatic tumors.  The
%   isfield-plus-fallback pattern ensures backwards compatibility: older
%   config.json files that lack newer fields continue to work unchanged,
%   which is critical when re-analyzing historical patient cohorts with
%   updated code.
%
%   The threshold values below collectively define a biophysical model of
%   pancreatic tumor tissue: ADC and D thresholds delineate cellularity
%   (restricted diffusion), f and D* thresholds characterize the
%   hypovascular desmoplastic stroma typical of pancreatic adenocarcinoma,
%   and the IVIM b-value threshold separates perfusion-contaminated from
%   diffusion-dominated signal regimes.  Together, these thresholds enable
%   voxel-level classification of tumor biology from DWI data.
%
%
    % Fail early if the config file is missing.  The pipeline cannot
    % proceed without site-specific paths (dataloc, dcm2nii_call) that
    % vary across institutions and cannot be safely defaulted.
    if ~isfile(json_path)

        error('parse_config:fileNotFound', 'Configuration file %s not found. Please copy, rename, and fill out config.example.json.', json_path);
    end
    try
        % Hard limit on config file size to prevent accidental memory
        % exhaustion if a user mistakenly points config_path to a large
        % data file.  A valid config.json should be well under 1 MB;
        % anything larger is almost certainly a misconfigured path.
        MAX_CONFIG_BYTES = 1e6;  % 1 MB hard limit
        file_info = dir(json_path);
        if file_info.bytes > MAX_CONFIG_BYTES
            error('parse_config:fileTooLarge', ...
                'Configuration file %s is too large (%d bytes, limit %d bytes). A valid config.json should be well under 1 MB. Check that config_path does not point to a data file.', ...
                json_path, file_info.bytes, MAX_CONFIG_BYTES);
        end
        
        % Read and parse the JSON config.  jsondecode converts JSON objects
        % to MATLAB structs, arrays to matrices, and strings to char
        % arrays — the native MATLAB types used downstream.  The try-catch
        % around this specific call provides a user-friendly error message
        % when the JSON is malformed (trailing commas, unquoted keys, etc.).
        raw_text = fileread(json_path);
        try
            config_struct = jsondecode(raw_text);
        catch ME
            error('parse_config:invalidJSON', ...
                'Failed to parse %s: %s\nPlease validate your JSON syntax (e.g., https://jsonlint.com/).', ...
                json_path, ME.message);
        end

        % ================================================================
        % Default value assignments for optional configuration fields.
        %
        % Each default is chosen based on domain conventions in pancreatic
        % DWI / IVIM analysis.  The isfield guard ensures that config files
        % from earlier pipeline versions (which lack these fields) still
        % parse correctly — a strict backwards-compatibility requirement
        % for reproducibility of published analyses.
        % ================================================================

        % Assign defaults for missing fields

        % skip_to_reload: When true, bypasses raw DICOM discovery and model
        % fitting (Sections 1-3 of load_dwi_data), jumping straight to
        % loading cached dwi_vectors.mat.  Default false forces the full
        % compute path, which is the safe default for first-time runs.
        if ~isfield(config_struct, 'skip_to_reload')
            config_struct.skip_to_reload = false;
        end

        % skip_tests: When true, bypasses the pre-flight test suite before
        % pipeline execution.  Default false ensures data integrity checks
        % always run, catching regressions before processing patient data.
        if ~isfield(config_struct, 'skip_tests')
            config_struct.skip_tests = false;
        end

        % use_checkpoints: Enables per-patient caching of intermediate
        % results (DICOM conversion, model fits).  Essential for large
        % cohorts where a single patient failure should not require
        % reprocessing the entire dataset.  Default false for simplicity
        % in small-cohort or development runs.
        if ~isfield(config_struct, 'use_checkpoints')
            config_struct.use_checkpoints = false;
        end

        % adc_thresh: ADC threshold (mm^2/s) for defining restricted
        % diffusion sub-volumes.  In pancreatic cancer, ADC values below
        % ~1.0 x 10^-3 mm^2/s indicate highly cellular (viable tumor)
        % tissue.  This threshold separates tumor-like from normal
        % parenchyma in the mono-exponential diffusion model.
        if ~isfield(config_struct, 'adc_thresh')
            config_struct.adc_thresh = 0.001;
        end

        % high_adc_thresh: Upper ADC threshold (mm^2/s) for identifying
        % regions with elevated diffusivity, which may indicate treatment
        % response (cell death / edema).  The 1.15 x 10^-3 mm^2/s value
        % sits above the typical viable-tumor range but below free water,
        % capturing the early cellular breakdown signal seen during RT.
        if ~isfield(config_struct, 'high_adc_thresh')
            config_struct.high_adc_thresh = 0.00115;
        end

        % d_thresh: IVIM true diffusion coefficient (D) threshold (mm^2/s).
        % D represents tissue diffusivity after removing the perfusion
        % contribution.  In IVIM theory, D < 1.0 x 10^-3 mm^2/s indicates
        % restricted diffusion consistent with dense tumor cellularity,
        % similar to ADC but corrected for perfusion contamination at low
        % b-values.
        if ~isfield(config_struct, 'd_thresh')
            config_struct.d_thresh = 0.001;
        end

        % f_thresh: IVIM perfusion fraction (f) threshold (dimensionless,
        % 0-1 range).  f represents the fraction of the DWI signal
        % attributable to microvascular perfusion (capillary blood flow)
        % rather than true Brownian diffusion.  Pancreatic tumors are
        % typically hypovascular, so f < 0.10 (10%) indicates poor
        % perfusion characteristic of desmoplastic pancreatic stroma.
        if ~isfield(config_struct, 'f_thresh')
            config_struct.f_thresh = 0.1;
        end

        % dstar_thresh: IVIM pseudo-diffusion coefficient (D*) threshold
        % (mm^2/s).  D* models the fast signal decay from incoherent
        % capillary blood flow.  Values below 0.01 mm^2/s suggest reduced
        % microvascular flow velocity, consistent with the hypoperfused
        % microenvironment of pancreatic adenocarcinoma.
        if ~isfield(config_struct, 'dstar_thresh')
            config_struct.dstar_thresh = 0.01;
        end

        % ivim_bthr: b-value threshold (s/mm^2) separating the perfusion-
        % sensitive regime (b < bthr) from the diffusion-dominated regime
        % (b >= bthr).  At b = 100 s/mm^2, the fast pseudo-diffusion
        % component (D* ~ 10-100 x 10^-3 mm^2/s) has largely decayed,
        % so signal above this threshold reflects primarily true tissue
        % diffusivity.  This is the standard cutoff in segmented IVIM
        % fitting (Le Bihan et al.).
        if ~isfield(config_struct, 'ivim_bthr')
            config_struct.ivim_bthr = 100;
        end

        % min_vox_hist: Minimum number of voxels within a GTV required to
        % generate reliable histograms and summary statistics.  With fewer
        % than ~100 voxels, histogram-derived metrics (skewness, kurtosis,
        % percentiles) become unreliable due to sampling noise, and the
        % spatial heterogeneity measures lose meaning.  This threshold also
        % prevents small partial-volume contaminated ROIs (e.g., a GTV
        % that is mostly outside the DWI field-of-view) from producing
        % misleading summary statistics that would bias downstream
        % predictive models.
        if ~isfield(config_struct, 'min_vox_hist')
            config_struct.min_vox_hist = 100;
        end

        % adc_max: Upper physiological bound for ADC (mm^2/s).  ADC values
        % above 3.0 x 10^-3 mm^2/s approach free water diffusivity (~3.0
        % at 37C) and indicate non-tissue signal (CSF contamination, cystic
        % regions, or fitting artifacts).  Voxels exceeding this are
        % excluded from tumor characterization.
        if ~isfield(config_struct, 'adc_max')
            config_struct.adc_max = 0.003;
        end

        % td_scan_days: Override for the scan-day schedule used in time-
        % dependent panel construction (build_td_panel).  Empty means use
        % the default schedule [0 5 10 15 20 90] corresponding to
        % fractions 1-5 plus post-RT follow-up.  Sites with non-standard
        % RT schedules can specify their actual scan days here.
        if ~isfield(config_struct, 'td_scan_days')
            config_struct.td_scan_days = [];
        end

        % cause_of_death_column: Column name in the clinical spreadsheet
        % that encodes cause of death for competing-risk survival analysis.
        % Distinguishing cancer-specific death from other causes (e.g.,
        % cardiovascular) is essential for Cause-Specific Hazard modeling
        % in pancreatic cancer where non-cancer mortality is non-negligible.
        % Set to "" (empty string) to disable the warning when no
        % cause-of-death data is available.
        if ~isfield(config_struct, 'cause_of_death_column')
            config_struct.cause_of_death_column = 'CauseOfDeath';
        end

        % patient_ids: Optional cell array of patient ID strings to process.
        % When empty (default), all patients found in the data directory
        % are processed.  When populated, only the specified patients are
        % included - useful for debugging a single case or re-running a
        % subset of the cohort without modifying the data directory.
        if ~isfield(config_struct, 'patient_ids')
            config_struct.patient_ids = {};
        end

        % clear_cache: When true, deletes all cached .mat files from the
        % data directory before pipeline execution (dwi_vectors_*.mat,
        % summary_metrics_*.mat, per-patient checkpoints).  Useful when
        % the cohort, GTV contours, or pipeline code have changed and
        % stale cached data would produce incorrect results.  Default
        % false preserves existing caches for incremental re-runs.
        if ~isfield(config_struct, 'clear_cache')
            config_struct.clear_cache = false;
        end

        % core_method: Algorithm to use for determining the tumor core
        % sub-volume inside the GTV. Default is 'adc_threshold' for
        % backwards compatibility. Other options include 'd_threshold',
        % 'df_intersection', 'otsu', 'gmm', 'kmeans', 'region_growing',
        % 'active_contours', 'percentile', 'spectral', and 'fdm'.
        if ~isfield(config_struct, 'core_method')
            config_struct.core_method = 'adc_threshold';
        end

        % core_percentile: Percentile cutoff for the 'percentile' core
        % method.  Voxels with ADC below the Nth percentile of the
        % patient's own GTV ADC distribution are classified as core.
        % N=25 means the densest quarter of the tumor by cellularity.
        % The result is capped at adc_thresh as a safety floor to avoid
        % labelling normal tissue as core in uniformly necrotic tumors.
        if ~isfield(config_struct, 'core_percentile')
            config_struct.core_percentile = 25;
        end

        % core_n_clusters: Number of tissue sub-populations for the
        % 'spectral' core method.  2 = binary core/non-core separation;
        % 3 = adds a transitional zone between viable tumor and necrosis.
        % Uses spectralcluster (R2019b+) with k-means fallback on older
        % MATLAB or Octave.
        if ~isfield(config_struct, 'core_n_clusters')
            config_struct.core_n_clusters = 2;
        end

        % fdm_parameter: Which diffusion parameter to use for functional
        % diffusion map (fDM) voxel-level change classification.  'adc'
        % is more robust (always available); 'd' removes perfusion
        % contamination but is noisier and may have more fit failures.
        if ~isfield(config_struct, 'fdm_parameter')
            config_struct.fdm_parameter = 'adc';
        end

        % fdm_thresh: Fallback fDM significance threshold (mm^2/s).
        % When repeatability-derived Coefficient of Reproducibility (CoR)
        % is unavailable (no Fx1 repeat scans), voxel-level changes
        % exceeding this value are classified as responding or progressing.
        % 0.4 x 10^-3 mm^2/s is a typical test-retest threshold for
        % abdominal DWI at 1.5-3T.
        if ~isfield(config_struct, 'fdm_thresh')
            config_struct.fdm_thresh = 0.0004;
        end

        % run_compare_cores: When true, the compare_cores step is
        % automatically included in the default pipeline steps so that
        % pairwise Dice/Hausdorff agreement across all 11 core methods
        % is computed without requiring explicit invocation.
        if ~isfield(config_struct, 'run_compare_cores')
            config_struct.run_compare_cores = false;
        end

        % run_all_core_methods: When true, compute_summary_metrics and
        % metrics_dosimetry run all 11 tumor core delineation methods per
        % patient/timepoint, storing per-method sub-volume metrics in
        % summary_metrics.all_core_metrics.<method_name>.  Default false
        % because this multiplies core extraction cost by 11x.
        if ~isfield(config_struct, 'run_all_core_methods')
            config_struct.run_all_core_methods = false;
        end

        % spectral_min_voxels: Minimum number of valid voxels required for
        % spectral clustering.  Unlike histogram-based methods that need
        % ~100 voxels for stable percentile estimates, spectral clustering
        % via the graph Laplacian can produce meaningful 2-class
        % separations with as few as 20 voxels.  Default 20 allows the
        % method to run on small tumors while still preventing degenerate
        % solutions.
        if ~isfield(config_struct, 'spectral_min_voxels')
            config_struct.spectral_min_voxels = 20;
        end

        % store_core_masks: When true and run_all_core_methods is true,
        % also stores the 1D logical core masks for each method in
        % summary_metrics.all_core_metrics.<method>.core_masks.  Enables
        % compare_core_methods to skip re-computation.  Default false to
        % limit memory/disk usage.
        if ~isfield(config_struct, 'store_core_masks')
            config_struct.store_core_masks = false;
        end

        % use_firth_refit: When true, the predictive modeling step refits
        % the final model (and each LOOCV fold) using Firth penalized
        % logistic regression after elastic net feature selection.  Firth's
        % Jeffreys-prior penalty produces finite coefficient estimates even
        % under perfect/quasi-perfect separation, which is common in small
        % pancreatic cancer cohorts.  Default true.
        if ~isfield(config_struct, 'use_firth_refit')
            config_struct.use_firth_refit = true;
        end

        % compute_fine_gray: When true, the survival analysis module fits a
        % Fine-Gray subdistribution hazard model in addition to the default
        % Cause-Specific Hazard Cox model.  The Fine-Gray model estimates
        % cumulative incidence of local failure accounting for competing
        % risks (death from other causes).  Default true.
        if ~isfield(config_struct, 'compute_fine_gray')
            config_struct.compute_fine_gray = true;
        end

        % exclude_motion_volumes: When true, DWI volumes flagged as motion-
        % corrupted by detect_motion_artifacts are excluded from model
        % fitting.  Default false (conservative: log warnings but include
        % all volumes).
        if ~isfield(config_struct, 'exclude_motion_volumes')
            config_struct.exclude_motion_volumes = false;
        end

        % use_texture_features: When true, compute texture features (GLCM,
        % first-order) on ADC maps and include them in the predictive model
        % feature matrix.  Default false to preserve the existing 22-column
        % layout for backward compatibility.
        if ~isfield(config_struct, 'use_texture_features')
            config_struct.use_texture_features = false;
        end

        % texture_3d: When true and the input volume is 3D, GLRLM texture
        % features are computed using all 13 3D directions rather than
        % the 4 in-plane directions from the largest 2D slice.  Default
        % true for volumetric analysis; set false to restrict to 2D
        % (faster, and appropriate when slice thickness >> in-plane
        % resolution makes inter-slice runs physically less meaningful).
        if ~isfield(config_struct, 'texture_3d')
            config_struct.texture_3d = true;
        end

        % texture_quantization_method: Controls how continuous parameter
        % values are discretised into grey levels for GLCM/GLRLM texture
        % features.  IBSI specifies both methods and notes that the choice
        % affects feature values (Section 3.4.1):
        %   'fixed_bin_number' (default): rescales to [1, n_levels] using
        %     the ROI min/max.  Good for cross-patient comparison.
        %   'fixed_bin_width': bins at fixed intervals of (range/n_levels).
        %     Good when absolute values are meaningful (e.g., ADC mm^2/s).
        if ~isfield(config_struct, 'texture_quantization_method')
            config_struct.texture_quantization_method = 'fixed_bin_number';
        end

        % run_imputation_sensitivity: When true, runs imputation sensitivity
        % analysis comparing KNN against LOCF, mean, and linear interpolation.
        % This is computationally expensive (4× LOOCV) so disabled by default.
        if ~isfield(config_struct, 'run_imputation_sensitivity')
            config_struct.run_imputation_sensitivity = false;
        end

        % fit_time_varying_cox: When true, fits stratified and extended Cox
        % models as follow-up when Schoenfeld residual tests detect PH
        % violations.  Includes covariate × log(time) interaction terms.
        if ~isfield(config_struct, 'fit_time_varying_cox')
            config_struct.fit_time_varying_cox = true;
        end

        % export_validation_model: When true, exports the trained elastic net
        % model and preprocessing pipeline to a .mat file for external
        % validation on independent datasets.
        if ~isfield(config_struct, 'export_validation_model')
            config_struct.export_validation_model = false;
        end

        % external_validation_data: Path to folder containing external
        % validation dataset (same directory structure as training data).
        % Empty string means disabled.
        if ~isfield(config_struct, 'external_validation_data')
            config_struct.external_validation_data = '';
        end

        % auxiliary_biomarker_csv: Path to CSV file with non-DWI biomarker
        % data (columns: patient_id, biomarker_name, timepoint, value).
        % Empty string means disabled.
        if ~isfield(config_struct, 'auxiliary_biomarker_csv')
            config_struct.auxiliary_biomarker_csv = '';
        end

        % use_auxiliary_biomarkers: When true, includes auxiliary biomarker
        % features in the predictive model feature matrix.
        if ~isfield(config_struct, 'use_auxiliary_biomarkers')
            config_struct.use_auxiliary_biomarkers = false;
        end

        % use_gpu: When true, offloads computationally intensive operations
        % to a CUDA-capable GPU via gpuArray.  Currently accelerates:
        %   - ADC monoexponential WLS fitting (vectorized matrix ops)
        %   - DnCNN deep learning inference (predict() on GPU)
        % IVIM segmented fitting remains CPU-only because it relies on the
        % read-only IVIMmodelfit dependency.  Default false to ensure the
        % pipeline works on machines without a GPU or the Parallel Computing
        % Toolbox.  When true but no GPU is available, falls back to CPU
        % with a warning.
        if ~isfield(config_struct, 'use_gpu')
            config_struct.use_gpu = false;
        end

        % gpu_device: 1-based index of the CUDA GPU device to use when
        % use_gpu is true.  On multi-GPU workstations, set this to select
        % a specific card (e.g., 2 for the second GPU).  Default 1
        % selects the first available device.
        if ~isfield(config_struct, 'gpu_device')
            config_struct.gpu_device = 1;
        end

        % run_trajectory_plots: generate waterfall, swimmer, and spider
        % plots during the visualization step.  Default true.
        if ~isfield(config_struct, 'run_trajectory_plots')
            config_struct.run_trajectory_plots = true;
        end

        % ================================================================
        % Type validation for critical configuration fields.
        %
        % jsondecode faithfully mirrors JSON types, but user errors in
        % config.json (e.g., quoting a number as a string, or writing a
        % comma-separated string instead of a JSON array) produce MATLAB
        % types that silently misbehave downstream.  For example, MATLAB's
        % '<' operator on a char array compares ASCII code points, not the
        % numeric content of the string, so "adc < adc_thresh" would give
        % wrong results if adc_thresh were accidentally a char.
        %
        % We validate types here, immediately after defaults are assigned,
        % to fail fast with an actionable error message.  Where possible,
        % we coerce recoverable types (e.g., char to double for numeric
        % fields, numeric 0/1 to logical) rather than erroring, but we
        % error if the coercion would be ambiguous or lossy.
        % ================================================================

        % --- Numeric scalar fields (thresholds and counts) ---
        numeric_fields = {'adc_thresh', 'd_thresh', 'f_thresh', ...
            'dstar_thresh', 'ivim_bthr', 'high_adc_thresh', ...
            'adc_max', 'min_vox_hist', 'core_percentile', ...
            'core_n_clusters', 'fdm_thresh', 'spectral_min_voxels', ...
            'gpu_device'};
        for i = 1:numel(numeric_fields)
            fn = numeric_fields{i};
            if isfield(config_struct, fn)
                val = config_struct.(fn);
                if ischar(val) || isstring(val)
                    % User likely quoted a number in JSON, e.g., "0.001".
                    % Attempt to recover by converting to double.
                    converted = str2double(val);
                    if isnan(converted)
                        error('parse_config:invalidType', ...
                            'Configuration field "%s" must be a numeric scalar, but got non-numeric string "%s". Remove the quotes around the value in config.json.', ...
                            fn, char(val));
                    end
                    config_struct.(fn) = converted;
                elseif ~isnumeric(val) || ~isscalar(val)
                    error('parse_config:invalidType', ...
                        'Configuration field "%s" must be a numeric scalar, but got %s (class: %s). Check your config.json — numeric values must not be quoted.', ...
                        fn, mat2str(val), class(val));
                end
            end
        end

        % --- Logical / boolean fields ---
        % JSON booleans become MATLAB logical via jsondecode.  We also
        % accept numeric 0/1 (which is what you get if the user writes
        % the JSON integer 0 or 1 instead of true/false) and coerce to
        % logical.  Strings like "true"/"false" are also coerced.
        logical_fields = {'skip_to_reload', 'skip_tests', 'use_checkpoints', ...
            'clear_cache', 'run_compare_cores', 'run_all_core_methods', ...
            'store_core_masks', 'use_firth_refit', 'compute_fine_gray', ...
            'exclude_motion_volumes', 'use_texture_features', 'texture_3d', ...
            'run_imputation_sensitivity', 'fit_time_varying_cox', ...
            'export_validation_model', 'use_auxiliary_biomarkers', ...
            'use_gpu', 'run_trajectory_plots'};
        for i = 1:numel(logical_fields)
            fn = logical_fields{i};
            if isfield(config_struct, fn)
                val = config_struct.(fn);
                if islogical(val) && isscalar(val)
                    % Already correct type — nothing to do.
                elseif isnumeric(val) && isscalar(val) && (val == 0 || val == 1)
                    % Coerce numeric 0/1 to logical.
                    config_struct.(fn) = logical(val);
                elseif (ischar(val) || isstring(val))
                    % User quoted a boolean, e.g., "true" or "false".
                    val_lower = lower(strtrim(char(val)));
                    if strcmp(val_lower, 'true')
                        config_struct.(fn) = true;
                    elseif strcmp(val_lower, 'false')
                        config_struct.(fn) = false;
                    else
                        error('parse_config:invalidType', ...
                            'Configuration field "%s" must be a logical (true/false), but got string "%s". In config.json, use true or false (unquoted).', ...
                            fn, char(val));
                    end
                else
                    error('parse_config:invalidType', ...
                        'Configuration field "%s" must be a logical (true/false) or numeric 0/1, but got %s (class: %s). In config.json, use true or false (unquoted).', ...
                        fn, mat2str(val), class(val));
                end
            end
        end

        % --- String fields ---
        % These fields must be character vectors (or string scalars, which
        % we coerce to char for consistency with the rest of the pipeline).
        string_fields = {'cause_of_death_column', 'core_method', ...
            'fdm_parameter', 'texture_quantization_method', ...
            'external_validation_data', 'auxiliary_biomarker_csv'};
        for i = 1:numel(string_fields)
            fn = string_fields{i};
            if isfield(config_struct, fn)
                val = config_struct.(fn);
                if isstring(val) && isscalar(val)
                    % Coerce MATLAB string to char for downstream consistency.
                    config_struct.(fn) = char(val);
                elseif ~ischar(val)
                    error('parse_config:invalidType', ...
                        'Configuration field "%s" must be a string, but got %s (class: %s).', ...
                        fn, mat2str(val), class(val));
                end
            end
        end

        % --- Numeric vector / array fields ---
        % td_scan_days must be a numeric vector (or empty).  A common
        % mistake is writing "0,5,10,15,20,90" (a comma-separated string)
        % instead of [0, 5, 10, 15, 20, 90] (a JSON array).
        if isfield(config_struct, 'td_scan_days')
            val = config_struct.td_scan_days;
            if ~isempty(val)
                if ischar(val) || isstring(val)
                    % Attempt to parse comma-separated numeric string.
                    val_str = strtrim(char(val));
                    converted = str2double(strsplit(val_str, ','));
                    if any(isnan(converted))
                        error('parse_config:invalidType', ...
                            'Configuration field "td_scan_days" must be a numeric vector (JSON array of numbers) or empty, but got string "%s". Use a JSON array like [0, 5, 10, 15, 20, 90], not a comma-separated string.', ...
                            val_str);
                    end
                    config_struct.td_scan_days = converted;
                elseif ~isnumeric(val) || ~isvector(val)
                    error('parse_config:invalidType', ...
                        'Configuration field "td_scan_days" must be a numeric vector (JSON array of numbers) or empty, but got %s (class: %s). Use a JSON array like [0, 5, 10, 15, 20, 90], not a comma-separated string.', ...
                        mat2str(val), class(val));
                end
            end
        end

        % --- Cell array fields ---
        % patient_ids should be a cell array of character vectors.
        % jsondecode produces a cell array from a JSON array of strings,
        % but a single-element array may produce a plain char.
        if isfield(config_struct, 'patient_ids')
            val = config_struct.patient_ids;
            if ischar(val)
                % Single string: wrap in a cell.
                config_struct.patient_ids = {val};
            elseif isstring(val)
                % MATLAB string array: convert to cell of char.
                config_struct.patient_ids = cellstr(val);
            elseif ~iscell(val) && ~isempty(val)
                error('parse_config:invalidType', ...
                    'Configuration field "patient_ids" must be a cell array of strings (JSON array of strings) or empty, but got class: %s.', ...
                    class(val));
            end
        end

        fprintf('Successfully loaded configuration from %s\n', json_path);
    catch ME
        % If the error was already wrapped with our ID, rethrow as-is.
        % Otherwise, wrap any unexpected field-access or type error with
        % a descriptive message so callers can distinguish config parsing
        % failures from other errors in their try/catch blocks.
        if strncmp(ME.identifier, 'parse_config:', 13)
            rethrow(ME);
        end
        error('parse_config:invalidJSON', 'Failed to parse JSON configuration file %s: %s', json_path, ME.message);
    end

    % ================================================================
    % Validate core_method against the set of implemented algorithms.
    % ================================================================
    valid_core_methods = {'adc_threshold', 'd_threshold', 'df_intersection', ...
        'otsu', 'gmm', 'kmeans', 'region_growing', 'active_contours', ...
        'percentile', 'spectral', 'fdm'};
    if ~any(strcmpi(config_struct.core_method, valid_core_methods))
        error('parse_config:invalidCoreMethod', ...
            'Unrecognized core_method "%s". Must be one of: %s', ...
            config_struct.core_method, strjoin(valid_core_methods, ', '));
    end

    % ================================================================
    % DWI type validation and mapping to numeric run indices.
    % ================================================================
    if isfield(config_struct, 'dwi_type')
        switch lower(config_struct.dwi_type)
            case 'standard', config_struct.dwi_types_to_run = 1;
            case 'dncnn', config_struct.dwi_types_to_run = 2;
            case 'ivimnet', config_struct.dwi_types_to_run = 3;
            otherwise
                error('parse_config:invalidJSON', ...
                    'Unrecognized dwi_type "%s". Must be one of: Standard, dnCNN, IVIMnet.', ...
                    config_struct.dwi_type);
        end
    else
        config_struct.dwi_types_to_run = 1:3;
    end
end