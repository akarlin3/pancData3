% execute_all_workflows.m
% Automates the /run_data workflow for sequential modeling pipelines

repo_root = pwd;

% --- 1. SET UP ENVIRONMENT AND POOL (max 2 workers) ---
if ~exist('OCTAVE_VERSION', 'builtin')
    p = gcp('nocreate');
    if ~isempty(p)
        delete(p);
    end
    p = parpool('Processes', 2, 'IdleTimeout', Inf);
    addAttachedFiles(p, {fullfile(repo_root, 'core', 'load_dwi_data.m')});
    pctRunOnAll addpath(fullfile(pwd, 'core'));
    pctRunOnAll addpath(fullfile(pwd, 'utils'));
    pctRunOnAll addpath(fullfile(pwd, 'dependencies'));
    pctRunOnAll warning('off', 'all');
else
    addpath(fullfile(pwd, 'core'));
    addpath(fullfile(pwd, 'utils'));
    addpath(fullfile(pwd, 'dependencies'));
    warning('off', 'all');
end

clear global MASTER_OUTPUT_FOLDER;

% Load the configuration JSON explicitly because we are going to modify it
fid = fopen('config.json', 'r');
raw = fread(fid, inf);
str = char(raw');
fclose(fid);
config_struct = jsondecode(str);

% Execute all 9 discrete target modules
steps = {'load', 'sanity', 'visualize', 'metrics_baseline', ...
         'metrics_longitudinal', 'metrics_dosimetry', ...
         'metrics_stats_comparisons', 'metrics_stats_predictive', 'metrics_survival'};

% --- 2. RUN STANDARD PIPELINE ---
disp('====== STARTING STANDARD PIPELINE ======');
config_struct.dwi_type = 'Standard';
config_struct.skip_to_reload = false;
json_str = jsonencode(config_struct);
fid = fopen('config.json', 'w'); fprintf(fid, '%s', json_str); fclose(fid);
run_dwi_pipeline('config.json', steps);

% --- 3. RUN dnCNN PIPELINE ---
disp('====== STARTING dnCNN PIPELINE ======');
config_struct.dwi_type = 'dnCNN';
config_struct.skip_to_reload = true;
json_str = jsonencode(config_struct);
fid = fopen('config.json', 'w'); fprintf(fid, '%s', json_str); fclose(fid);
run_dwi_pipeline('config.json', steps);

% --- 4. RUN IVIMnet PIPELINE ---
disp('====== STARTING IVIMnet PIPELINE ======');
config_struct.dwi_type = 'IVIMnet';
config_struct.skip_to_reload = true;
json_str = jsonencode(config_struct);
fid = fopen('config.json', 'w'); fprintf(fid, '%s', json_str); fclose(fid);
run_dwi_pipeline('config.json', steps);

disp('====== ALL WORKFLOWS COMPLETED ======');
