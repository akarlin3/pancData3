function add_demo_paths()
% ADD_DEMO_PATHS  Put the pipeline modules (and Octave shims) on the path.
%
% The demo calls the REAL fitter in pipeline/core, which depends on
% pipeline/utils and pipeline/dependencies. On Octave we also need the
% compatibility shims in pipeline/.octave_compat (nanmean, niftiread stubs,
% etc.) that the pipeline relies on. Computing paths relative to THIS file
% keeps the demo runnable from any working directory.
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);                       % repo root (parent of demo/)
    addpath(fullfile(root, 'pipeline', 'core'));
    addpath(fullfile(root, 'pipeline', 'utils'));
    addpath(fullfile(root, 'pipeline', 'dependencies'));
    addpath(here);                                % demo/ itself
    if exist('OCTAVE_VERSION', 'builtin')
        addpath(fullfile(root, 'pipeline', '.octave_compat'));
    end
end
