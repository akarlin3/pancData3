% Wrapper script to run the Standard Pipeline load phase
repo_root = fileparts(mfilename('fullpath'));
addpath(fullfile(repo_root, 'utils'));
addpath(fullfile(repo_root, 'dependencies'));

% Start pool and strictly force all workers to add the paths
if ~exist('OCTAVE_VERSION', 'builtin')
    p = gcp('nocreate');
    if ~isempty(p)
        delete(p);
    end
    p = parpool('Processes', 2, 'IdleTimeout', Inf);
    addAttachedFiles(p, {fullfile(repo_root, 'core', 'load_dwi_data.m')});
    pctRunOnAll addpath(fullfile(repo_root, 'core'));
    pctRunOnAll addpath(fullfile(repo_root, 'utils'));
    pctRunOnAll addpath(fullfile(repo_root, 'dependencies'));
else
    addpath(fullfile(repo_root, 'core'));
end

clear global MASTER_OUTPUT_FOLDER;
run_dwi_pipeline('config.json', {'load'});
