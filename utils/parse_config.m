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
        if ~isfield(config_struct, 'use_checkpoints')
            config_struct.use_checkpoints = false;
        end
        fprintf('Successfully loaded configuration from %s\n', json_path);
    catch ME
        error('Failed to parse JSON configuration file: %s', ME.message);
    end
end
