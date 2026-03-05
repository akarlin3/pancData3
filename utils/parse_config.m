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
    if ~isfile(json_path)

        error('parse_config:fileNotFound', 'Configuration file %s not found. Please copy, rename, and fill out config.example.json.', json_path);
    end
    try
        raw_json = fileread(json_path);
        config_struct = jsondecode(raw_json);

        % Assign defaults for missing fields
        if ~isfield(config_struct, 'skip_tests')
            config_struct.skip_tests = false;
        end
        if ~isfield(config_struct, 'use_checkpoints')
            config_struct.use_checkpoints = false;
        end
        if ~isfield(config_struct, 'adc_thresh')
            config_struct.adc_thresh = 0.001;   % must be <= high_adc_thresh
        end
        if ~isfield(config_struct, 'high_adc_thresh')
            config_struct.high_adc_thresh = 0.00115;  % must be >= adc_thresh
        end
        % Validate threshold ordering
        if config_struct.adc_thresh > config_struct.high_adc_thresh
            error('parse_config:thresholdOrder', ...
                'adc_thresh (%.6f) must be <= high_adc_thresh (%.6f).', ...
                config_struct.adc_thresh, config_struct.high_adc_thresh);
        end
        if ~isfield(config_struct, 'd_thresh')
            config_struct.d_thresh = 0.001;
        end
        if ~isfield(config_struct, 'f_thresh')
            config_struct.f_thresh = 0.1;
        end
        if ~isfield(config_struct, 'dstar_thresh')
            config_struct.dstar_thresh = 0.01;
        end
        if ~isfield(config_struct, 'ivim_bthr')
            config_struct.ivim_bthr = 100;
        end
        if ~isfield(config_struct, 'min_vox_hist')
            config_struct.min_vox_hist = 100;
        end
        if ~isfield(config_struct, 'adc_max')
            config_struct.adc_max = 0.003;
        end
        if ~isfield(config_struct, 'td_scan_days')
            config_struct.td_scan_days = [];
        end
        if ~isfield(config_struct, 'cause_of_death_column')
            config_struct.cause_of_death_column = 'CauseOfDeath';
        end
        if isfield(config_struct, 'dwi_type')
            switch lower(config_struct.dwi_type)
                case 'standard', config_struct.dwi_types_to_run = 1;
                case 'dncnn', config_struct.dwi_types_to_run = 2;
                case 'ivimnet', config_struct.dwi_types_to_run = 3;
                otherwise
                    error('parse_config:unknownDwiType', ...
                        'Unrecognized dwi_type "%s". Must be one of: Standard, dnCNN, IVIMnet.', ...
                        config_struct.dwi_type);
            end
        else
            config_struct.dwi_types_to_run = 1:3;
        end

        fprintf('Successfully loaded configuration from %s\n', json_path);
    catch ME
        error('parse_config:invalidJSON', 'Failed to parse JSON configuration file: %s', ME.message);
    end
end
