% execute_all_workflows.m
% Automates the /run_data workflow for sequential modeling pipelines

repo_root = pwd;

% --- 0. START LOGGING ---
logfile = fullfile(repo_root, 'log.out');
diary(logfile);
cleanupObj = onCleanup(@() diary('off'));

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
    warning('on', 'all');  % pctRunOnAll also runs on client — restore client warnings
end

% Reset persistent output folder in run_dwi_pipeline so each full workflow
% sequence starts with a fresh timestamped output directory.
clear run_dwi_pipeline;

% --- 1.5 RUN TEST SUITE BEFORE PIPELINE ---
disp('====== RUNNING TEST SUITE BEFORE PIPELINE ======');
try
    run(fullfile(repo_root, 'tests', 'run_all_tests.m'));
    disp('====== ALL TESTS PASSED — PROCEEDING WITH PIPELINE ======');
catch ME
    fprintf('❌ Test failure: %s\n', ME.message);
    error('PipelineAborted:TestFailure', ...
        'Pipeline aborted: test suite did not pass.');
end

% Load the configuration JSON explicitly because we are going to modify it
fid = fopen('config.json', 'r');
raw = fread(fid, inf);
str = char(raw');
fclose(fid);
config_struct = jsondecode(str);
original_config_str = str;  % preserve original for rollback

% Backup config.json so a mid-run crash does not leave it in a modified state
config_backup = 'config.json.bak';
fid_bak = fopen(config_backup, 'w');
fprintf(fid_bak, '%s', str);
fclose(fid_bak);
restore_config = onCleanup(@() restore_config_file(original_config_str, config_backup));

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

function restore_config_file(original_str, backup_path)
% Restore config.json to its original state when the script exits (whether
% by normal completion or error).  Also removes the backup file.
    try
        fid = fopen('config.json', 'w');
        fprintf(fid, '%s', original_str);
        fclose(fid);
        if exist(backup_path, 'file')
            delete(backup_path);
        end
    catch
        fprintf('⚠️  Could not restore config.json. Backup saved at %s\n', backup_path);
    end
end
