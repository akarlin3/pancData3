% run_demo_tests.m  —  test runner for the synthetic-demo physics + cohort.
% SYNTHETIC PHANTOM DATA — not clinical.
%
%   MATLAB:        run('demo/tests/run_demo_tests.m')
%   Octave (CLI):  octave demo/tests/run_demo_tests.m
%
% Plain assert-based tests (no toolbox dependency) so they run identically on
% MATLAB and Octave. Each test is a function that throws on failure; the
% runner tallies pass/fail and exits non-zero if anything fails.

tests_dir__ = fileparts(mfilename('fullpath'));
demo_dir__  = fileparts(tests_dir__);
addpath(demo_dir__);
addpath(tests_dir__);
add_demo_paths();

tests = { ...
    @test_ivim_signal, ...
    @test_add_rician_noise, ...
    @test_adc_from_signal, ...
    @test_synthetic_ivim };

fprintf('\n=== synthetic-demo test suite ===\n');
n_pass = 0; n_fail = 0; failures = {};
for ti = 1:numel(tests)
    name = func2str(tests{ti});
    try
        tests{ti}();
        fprintf('  PASS  %s\n', name);
        n_pass = n_pass + 1;
    catch err
        fprintf('  FAIL  %s — %s\n', name, err.message);
        n_fail = n_fail + 1;
        failures{end+1} = name; %#ok<SAGROW>
    end
end
fprintf('---------------------------------\n');
fprintf('  %d passed, %d failed\n\n', n_pass, n_fail);

if n_fail > 0
    if exist('OCTAVE_VERSION', 'builtin')
        exit(1);
    else
        error('run_demo_tests:failures', '%d test(s) failed: %s', n_fail, strjoin(failures, ', '));
    end
end
