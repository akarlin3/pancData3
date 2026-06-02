% run_demo.m  —  one-command entry point for the PHI-free synthetic demo.
% SYNTHETIC PHANTOM DATA — not clinical.
%
% Run from the repo root (or anywhere):
%     MATLAB:        run('demo/run_demo.m')
%     Octave (CLI):  octave demo/run_demo.m
%
% This is a thin SCRIPT wrapper so the file is directly runnable (a function
% file is not auto-invoked by `octave demo/run_demo.m`). It locates the demo
% folder relative to itself, puts the pipeline on the path, and calls
% run_synthetic_demo(), which does the real work:
%   synthetic generation -> REAL pipeline fit -> recover-known-truth check
%   -> headline figures, all written to demo/output/.
%
% To customise (SNR, cohort size, seed), call the function directly, e.g.:
%     run_synthetic_demo('snr', 20, 'n_patients', 30)

demo_dir__ = fileparts(mfilename('fullpath'));
addpath(demo_dir__);
add_demo_paths();
run_synthetic_demo();
clear demo_dir__;
