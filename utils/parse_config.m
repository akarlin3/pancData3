function config_struct = parse_config(json_path)
% PARSE_CONFIG Reads a JSON config file into a MATLAB struct.
% 
% Usage:
%   config_struct = parse_config('config.json')

    if ~isfile(json_path)
        error('Configuration file %s not found. Please copy, rename, and fill out config.example.json.', json_path);
    end
    try
        raw_json = fileread(json_path);
        config_struct = jsondecode(raw_json);

        % Assign defaults for missing fields
        if ~isfield(config_struct, 'use_checkpoints')
            config_struct.use_checkpoints = false;
        end
        if ~isfield(config_struct, 'adc_thresh')
            config_struct.adc_thresh = 0.00115;
        end
        if ~isfield(config_struct, 'high_adc_thresh')
            config_struct.high_adc_thresh = 0.001;
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

        fprintf('Successfully loaded configuration from %s\n', json_path);
    catch ME
        error('Failed to parse JSON configuration file: %s', ME.message);
    end
end
