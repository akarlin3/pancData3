% Wrapper script to run the Standard Pipeline load phase
repo_root = fileparts(mfilename('fullpath'));
addpath(fullfile(repo_root, 'utils'));
addpath(fullfile(repo_root, 'dependencies'));

% Start pool and strictly force all workers to add the paths
p = gcp('nocreate');
if ~isempty(p)
    delete(p);
end
p = parpool('Processes', 2, 'IdleTimeout', Inf);
addAttachedFiles(p, {fullfile(repo_root, 'core', 'load_dwi_data.m')});
pctRunOnAll addpath('C:\Users\karlina1\Desktop\pancData3\core');
pctRunOnAll addpath('C:\Users\karlina1\Desktop\pancData3\utils');
pctRunOnAll addpath('C:\Users\karlina1\Desktop\pancData3\dependencies');

clear global MASTER_OUTPUT_FOLDER;
run_dwi_pipeline('config.json', {'load'});
