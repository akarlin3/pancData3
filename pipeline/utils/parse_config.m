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
%   defaults-iteration pattern ensures backwards compatibility: older
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

        % Guard against dir() returning zero or multiple entries.  This
        % can happen when json_path contains glob-like characters (e.g.,
        % brackets in "config[v2].json") or on case-insensitive
        % filesystems.  The same pattern used in safe_load_mask.m.
        if numel(file_info) == 0
            error('parse_config:fileNotFound', ...
                'dir() returned no entries for %s despite isfile() succeeding. Check the path for special characters.', ...
                json_path);
        elseif numel(file_info) > 1
            % Attempt exact name match against the entries returned by dir.
            [~, expected_name, expected_ext] = fileparts(json_path);
            expected_basename = [expected_name, expected_ext];
            exact_idx = find(strcmp({file_info.name}, expected_basename), 1);
            if ~isempty(exact_idx)
                file_info = file_info(exact_idx);
            else
                warning('parse_config:multipleDir', ...
                    'dir() returned %d entries for "%s" and no exact name match was found. Using the first entry ("%s").', ...
                    numel(file_info), json_path, file_info(1).name);
                file_info = file_info(1);
            end
        end

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
        % DWI / IVIM analysis.  Defaults are declared in a single struct
        % and applied via iteration over fieldnames, so adding a new
        % config field requires only a single line here.  The iteration
        % ensures that config files from earlier pipeline versions (which
        % lack newer fields) still parse correctly — a strict backwards-
        % compatibility requirement for reproducibility of published
        % analyses.
        %
        % Threshold rationale (see inline comments for each field):
        %   adc_thresh / d_thresh: 1.0e-3 mm^2/s — restricted diffusion
        %     cutoff separating viable tumor from normal parenchyma.
        %   high_adc_thresh: 1.15e-3 mm^2/s — elevated diffusivity
        %     indicating treatment response (cell death / edema).
        %   f_thresh: 0.10 — perfusion fraction below 10% indicates
        %     hypovascular desmoplastic stroma.
        %   dstar_thresh: 0.01 mm^2/s — reduced microvascular flow.
        %   ivim_bthr: 100 s/mm^2 — standard segmented IVIM cutoff.
        %   adc_max: 3.0e-3 mm^2/s — upper physiological bound (free
        %     water at 37°C).
        %   min_vox_hist: 100 — minimum voxels for reliable histograms.
        %   fdm_thresh: 0.4e-3 mm^2/s — typical abdominal DWI test-retest
        %     threshold.
        %   core_percentile: 25 — densest quarter of the GTV by ADC.
        %   core_n_clusters: 2 — binary core/non-core separation.
        %   spectral_min_voxels: 20 — minimum for graph Laplacian methods.
        %   min_core_voxels: 10 — minimum median core size for retention.
        %   max_core_failure_rate: 1.0 — disabled by default (no pruning).
        % ================================================================

        defaults = struct();

        % --- Boolean flags ---
        defaults.skip_to_reload              = false;
        defaults.skip_tests                  = false;
        defaults.tests_only                  = false;
        defaults.use_checkpoints             = false;
        defaults.clear_cache                 = false;
        defaults.run_compare_cores           = false;
        defaults.run_cross_pipeline_dice     = false;
        defaults.run_core_failure_rates      = false;
        defaults.run_core_method_outcomes    = false;
        defaults.run_all_core_methods        = false;
        defaults.store_core_masks            = false;
        defaults.use_firth_refit             = true;
        defaults.compute_fine_gray           = true;
        defaults.exclude_motion_volumes      = false;
        defaults.use_texture_features        = false;
        defaults.texture_3d                  = true;
        defaults.run_imputation_sensitivity  = false;
        defaults.fit_time_varying_cox        = true;
        defaults.export_validation_model     = false;
        defaults.use_auxiliary_biomarkers     = false;
        defaults.use_gpu                     = false;
        defaults.run_trajectory_plots        = true;
        defaults.run_subvolume_stability     = false;
        defaults.run_dose_response_roc       = false;
        defaults.run_risk_dose_concordance   = false;
        defaults.run_per_method_cor          = false;
        defaults.run_gtv_confounding         = false;
        defaults.run_optimize_threshold      = false;
        defaults.run_baseline_vs_delta       = false;
        % When false, the pipeline ignores nodal GTVs (GTVn) entirely:
        % discovery skips GTVn mask paths, voxel extraction skips nodal
        % data, and downstream metrics treat every patient as primary-
        % only. Default is false because the MSK cohort's nodal contours
        % are largely incidental data and not part of the analysis. Set
        % to true to opt back in for sites that want nodal metrics.
        defaults.process_gtvn                = false;

        % --- Numeric scalar thresholds and counts ---
        defaults.adc_thresh                  = 0.001;
        defaults.high_adc_thresh             = 0.00115;
        defaults.d_thresh                    = 0.001;
        defaults.f_thresh                    = 0.1;
        defaults.dstar_thresh                = 0.01;
        defaults.ivim_bthr                   = 100;
        defaults.min_vox_hist                = 100;
        defaults.adc_max                     = 0.003;
        defaults.core_percentile             = 25;
        defaults.core_n_clusters             = 2;
        defaults.fdm_thresh                  = 0.0004;
        defaults.spectral_min_voxels         = 20;
        defaults.gpu_device                  = 1;
        defaults.max_core_failure_rate       = 1.0;
        defaults.min_core_voxels             = 10;

        % --- String fields ---
        defaults.cause_of_death_column       = 'CauseOfDeath';
        defaults.core_method                 = 'adc_threshold';
        defaults.fdm_parameter               = 'adc';
        defaults.texture_quantization_method = 'fixed_bin_number';
        defaults.external_validation_data    = '';
        defaults.auxiliary_biomarker_csv     = '';

        % --- Numeric vector / array fields ---
        defaults.td_scan_days                = [];

        % --- Cell array fields ---
        defaults.patient_ids                 = {};
        defaults.excluded_core_methods       = {};

        % Apply defaults: for every field in defaults, if the user's
        % config_struct does not already contain it, assign the default.
        default_fields = fieldnames(defaults);
        for i = 1:numel(default_fields)
            fn = default_fields{i};
            if ~isfield(config_struct, fn)
                config_struct.(fn) = defaults.(fn);
            end
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
            'gpu_device', 'max_core_failure_rate', 'min_core_voxels'};
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
        logical_fields = {'skip_to_reload', 'skip_tests', 'tests_only', 'use_checkpoints', ...
            'clear_cache', 'run_compare_cores', 'run_cross_pipeline_dice', 'run_core_failure_rates', 'run_core_method_outcomes', 'run_all_core_methods', ...
            'store_core_masks', 'use_firth_refit', 'compute_fine_gray', ...
            'exclude_motion_volumes', 'use_texture_features', 'texture_3d', ...
            'run_imputation_sensitivity', 'fit_time_varying_cox', ...
            'export_validation_model', 'use_auxiliary_biomarkers', ...
            'use_gpu', 'run_trajectory_plots', ...
            'run_subvolume_stability', 'run_dose_response_roc', ...
            'run_risk_dose_concordance', 'run_per_method_cor', ...
            'run_gtv_confounding', 'run_optimize_threshold', ...
            'run_baseline_vs_delta'};
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
            elseif isnumeric(val) && isempty(val)
                % jsondecode converts JSON [] to double []; normalize to {}
                config_struct.patient_ids = {};
            elseif isnumeric(val) && ~isempty(val)
                % jsondecode converts JSON [10, 20] to numeric array; convert to cell of strings
                config_struct.patient_ids = arrayfun(@(x) num2str(x), val(:)', 'UniformOutput', false);
            elseif ~iscell(val) && ~isempty(val)
                error('parse_config:invalidType', ...
                    'Configuration field "patient_ids" must be a cell array of strings (JSON array of strings) or empty, but got class: %s.', ...
                    class(val));
            end
        end

        % excluded_core_methods: same cell-array normalization as patient_ids.
        if isfield(config_struct, 'excluded_core_methods')
            val = config_struct.excluded_core_methods;
            if ischar(val)
                config_struct.excluded_core_methods = {val};
            elseif isstring(val)
                config_struct.excluded_core_methods = cellstr(val);
            elseif isnumeric(val) && isempty(val)
                config_struct.excluded_core_methods = {};
            elseif ~iscell(val) && ~isempty(val)
                error('parse_config:invalidType', ...
                    'Configuration field "excluded_core_methods" must be a cell array of strings (JSON array of strings) or empty, but got class: %s.', ...
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