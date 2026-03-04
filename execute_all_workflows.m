% execute_all_workflows.m
% Automates the /run_data workflow for sequential modeling pipelines

repo_root = fileparts(mfilename('fullpath'));
config_file = fullfile(repo_root, 'config.json');

% --- 1. SET UP ENVIRONMENT AND POOL (max 2 workers) ---
if ~exist('OCTAVE_VERSION', 'builtin')
    % Delete any stale parallel jobs before creating a new pool
    try
        allJobs = matlab.internal.parallel.getJobs();
        if ~isempty(allJobs)
            delete(allJobs);
            fprintf('⚙️ Deleted %d stale parallel job(s).\n', numel(allJobs));
        end
    catch
        % Fallback: use the job manager from the default cluster profile
        try
            c = parcluster;
            if ~isempty(c.Jobs)
                delete(c.Jobs);
                fprintf('⚙️ Deleted stale parallel jobs from cluster profile.\n');
            end
        catch
            % No jobs to clean up
        end
    end

    p = gcp('nocreate');
    if ~isempty(p)
        delete(p);
    end
    p = parpool('Processes', 2, 'IdleTimeout', Inf);
    addAttachedFiles(p, {fullfile(repo_root, 'core', 'load_dwi_data.m')});
    pctRunOnAll addpath(fullfile(repo_root, 'core'));
    pctRunOnAll addpath(fullfile(repo_root, 'utils'));
    pctRunOnAll addpath(fullfile(repo_root, 'dependencies'));
    w_state_before_pct = warning;  % save client warning state before pctRunOnAll clobbers it
    pctRunOnAll warning('off', 'all');
    warning(w_state_before_pct);   % restore client warnings (pctRunOnAll also executes on client)
end

% Reset persistent output folder in run_dwi_pipeline so each full workflow
% sequence starts with a fresh timestamped output directory.
clear run_dwi_pipeline;

% Create the timestamped output folder now so the diary can live there.
% run_dwi_pipeline will detect this folder via its persistent variable.
timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
eaw_output_folder = fullfile(repo_root, sprintf('saved_files_%s', timestamp_str));
if ~exist(eaw_output_folder, 'dir'), mkdir(eaw_output_folder); end

% Master diary for execute_all_workflows console output
eaw_diary_file = fullfile(eaw_output_folder, 'execute_all_workflows.log');
diary(eaw_diary_file);

% --- 1.5 RUN TEST SUITE BEFORE PIPELINE ---
disp('====== RUNNING TEST SUITE BEFORE PIPELINE ======');
test_diary_file = fullfile(eaw_output_folder, 'test_suite_output.log');
diary(test_diary_file);
try
    run(fullfile(repo_root, 'tests', 'run_all_tests.m'));
    disp('====== ALL TESTS PASSED — PROCEEDING WITH PIPELINE ======');
    diary(eaw_diary_file);  % switch back to master diary
catch ME
    diary(eaw_diary_file);  % switch back before error
    fprintf('❌ Test failure: %s\n', ME.message);
    error('PipelineAborted:TestFailure', ...
        'Pipeline aborted: test suite did not pass.');
end

% If a backup from a previous crashed run exists, restore it first.
% onCleanup does not fire on kill -9 or MATLAB segfaults, so the .bak
% file is the only recovery mechanism in those cases.
config_backup = fullfile(repo_root, 'config.json.bak');
if exist(config_backup, 'file')
    fprintf('⚠️  Found config.json.bak from a previous crash. Restoring original config.json.\n');
    fid_restore = fopen(config_backup, 'r');
    bak_raw = fread(fid_restore, inf);
    fclose(fid_restore);
    fid_restore_w = fopen(config_file, 'w');
    fwrite(fid_restore_w, bak_raw);
    fclose(fid_restore_w);
    delete(config_backup);
end

% Load the configuration JSON explicitly because we are going to modify it
fid = fopen(config_file, 'r');
raw = fread(fid, inf);
str = char(raw');
fclose(fid);
config_struct = jsondecode(str);
original_config_str = str;  % preserve original for rollback

% Backup config.json so a mid-run crash does not leave it in a modified state
fid_bak = fopen(config_backup, 'w');
fprintf(fid_bak, '%s', str);
fclose(fid_bak);
restore_config = onCleanup(@() restore_config_file(config_file, original_config_str, config_backup));

% Execute all 9 discrete target modules
steps = {'load', 'sanity', 'visualize', 'metrics_baseline', ...
         'metrics_longitudinal', 'metrics_dosimetry', ...
         'metrics_stats_comparisons', 'metrics_stats_predictive', 'metrics_survival'};

% --- 2. RUN STANDARD PIPELINE ---
disp('====== STARTING STANDARD PIPELINE ======');
config_struct.dwi_type = 'Standard';
config_struct.skip_to_reload = false;
json_str = jsonencode(config_struct);
fid = fopen(config_file, 'w'); fprintf(fid, '%s', json_str); fclose(fid);
run_dwi_pipeline(config_file, steps, eaw_output_folder);
diary(eaw_diary_file);  % restart after pipeline run

% --- 3. RUN dnCNN PIPELINE ---
disp('====== STARTING dnCNN PIPELINE ======');
config_struct.dwi_type = 'dnCNN';
config_struct.skip_to_reload = true;
json_str = jsonencode(config_struct);
fid = fopen(config_file, 'w'); fprintf(fid, '%s', json_str); fclose(fid);
run_dwi_pipeline(config_file, steps, eaw_output_folder);
diary(eaw_diary_file);  % restart after pipeline run

% --- 4. RUN IVIMnet PIPELINE ---
disp('====== STARTING IVIMnet PIPELINE ======');
config_struct.dwi_type = 'IVIMnet';
config_struct.skip_to_reload = true;
json_str = jsonencode(config_struct);
fid = fopen(config_file, 'w'); fprintf(fid, '%s', json_str); fclose(fid);
run_dwi_pipeline(config_file, steps, eaw_output_folder);
diary(eaw_diary_file);  % restart after pipeline run

disp('====== ALL WORKFLOWS COMPLETED ======');
diary off;

function restore_config_file(config_path, original_str, backup_path)
% Restore config.json to its original state when the script exits (whether
% by normal completion or error).  Also removes the backup file.
    try
        fid = fopen(config_path, 'w');
        fprintf(fid, '%s', original_str);
        fclose(fid);
        if exist(backup_path, 'file')
            delete(backup_path);
        end
    catch
        fprintf('⚠️  Could not restore config.json. Backup saved at %s\n', backup_path);
    end
end
