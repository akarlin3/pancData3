% run_all_tests_octave.m
% Octave-compatible test runner for the pancData3 test suite.
% Discovers classdef test files, parses setup/teardown/test methods,
% and executes them with summary reporting.

testsDir  = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(testsDir);

% Set up paths
addpath(fullfile(repoRoot, '.octave_compat'));
addpath(fullfile(repoRoot, 'core'));
addpath(fullfile(repoRoot, 'utils'));
addpath(fullfile(repoRoot, 'dependencies'));
addpath(repoRoot);
addpath(genpath(testsDir));

% Load Octave Forge packages if available
try; pkg load statistics; catch; end
try; pkg load image; catch; end
try; pkg load io; catch; end

% Suppress common Octave warnings that are harmless during test execution
warning('off', 'Octave:classdef-to-struct');
warning('off', 'Octave:mixed-string-concat');

fprintf('===================================================\n');
fprintf('   Octave CI Test Runner: Initializing Suite       \n');
fprintf('===================================================\n');

% Discover test files (test_*.m and benchmark_*.m) in tests/ and subdirectories
test_files = dir(fullfile(testsDir, 'test_*.m'));
for k = 1:numel(test_files)
    test_files(k).folder = testsDir;
end

sub_dirs = {'benchmarks', 'diagnostics'};
for d = 1:numel(sub_dirs)
    sd = fullfile(testsDir, sub_dirs{d});
    if exist(sd, 'dir')
        extra = dir(fullfile(sd, 'test_*.m'));
        for e = 1:numel(extra); extra(e).folder = sd; end
        test_files = [test_files; extra];
        % Note: benchmark_*.m files are skipped in Octave because their large
        % matrix operations (e.g. 200x2000 correlation) are prohibitively slow
        % without MATLAB's optimised LAPACK/BLAS bindings.
    end
end

% ---- Collect all test entries ----
entry_files      = {};
entry_classes    = {};
entry_tests      = {};
entry_setups     = {};
entry_teardowns  = {};
entry_folders    = {};
nEntries = 0;

for i = 1:numel(test_files)
    fpath = fullfile(test_files(i).folder, test_files(i).name);
    [~, className] = fileparts(test_files(i).name);

    % Check if classdef
    fid = fopen(fpath, 'r');
    if fid == -1; continue; end
    firstLine = fgetl(fid);
    fclose(fid);
    if ~ischar(firstLine) || isempty(strfind(firstLine, 'classdef'))
        continue;
    end

    % Read full file
    fid = fopen(fpath, 'r');
    txt = fread(fid, '*char')';
    fclose(fid);

    % Parse method blocks using external helper
    test_methods     = octave_parse_test_methods(txt, 'Test');
    setup_methods    = octave_parse_test_methods(txt, 'TestMethodSetup');
    teardown_methods = octave_parse_test_methods(txt, 'TestMethodTeardown');

    setupName    = '';
    teardownName = '';
    if ~isempty(setup_methods);    setupName = setup_methods{1};       end
    if ~isempty(teardown_methods); teardownName = teardown_methods{1}; end

    for j = 1:numel(test_methods)
        nEntries = nEntries + 1;
        entry_files{nEntries}     = fpath;
        entry_classes{nEntries}   = className;
        entry_tests{nEntries}     = test_methods{j};
        entry_setups{nEntries}    = setupName;
        entry_teardowns{nEntries} = teardownName;
        entry_folders{nEntries}   = test_files(i).folder;
    end
end

nTests = nEntries;
classNames = unique(entry_classes);
fprintf('Discovered %d test methods across %d class files.\n\n', nTests, numel(classNames));

fprintf('===================================================\n');
fprintf('   Running Tests...                                \n');
fprintf('===================================================\n\n');

% ---- Execute tests ----
nPassed = 0;
nFailed = 0;
fail_names = {};
fail_msgs  = {};

for i = 1:nTests
    className      = entry_classes{i};
    testMethod     = entry_tests{i};
    setupMethod    = entry_setups{i};
    teardownMethod = entry_teardowns{i};
    folder         = entry_folders{i};
    testLabel      = [className '/' testMethod];

    fprintf('  [%3d/%3d] %-55s ', i, nTests, testLabel);

    t0 = tic;
    obj = [];
    try
        % Ensure folder is on path
        addpath(folder);

        % Create instance
        obj = feval(className);

        % Setup
        if ~isempty(setupMethod)
            feval(setupMethod, obj);
        end

        % Run test
        feval(testMethod, obj);

        % Teardown
        if ~isempty(teardownMethod)
            try; feval(teardownMethod, obj); catch; end
        end

        nPassed = nPassed + 1;
        fprintf('PASSED  (%.2fs)\n', toc(t0));
    catch err
        % Teardown on failure
        if ~isempty(obj) && ~isempty(teardownMethod)
            try; feval(teardownMethod, obj); catch; end
        end

        nFailed = nFailed + 1;
        fprintf('FAILED  (%.2fs)\n', toc(t0));
        % Truncate long messages
        msg = err.message;
        if numel(msg) > 200; msg = [msg(1:200) '...']; end
        fprintf('         -> %s\n', msg);
        fail_names{end+1} = testLabel;
        fail_msgs{end+1}  = msg;
    end
end

% ---- Summary ----
fprintf('\n===================================================\n');
fprintf('   Test Execution Completed                        \n');
fprintf('===================================================\n');
fprintf('  Total:  %d\n', nTests);
fprintf('  Passed: %d\n', nPassed);
fprintf('  Failed: %d\n', nFailed);
fprintf('===================================================\n\n');

if nFailed > 0
    fprintf('FAILED TESTS:\n');
    for i = 1:numel(fail_names)
        fprintf('  [%d] %s\n      %s\n', i, fail_names{i}, fail_msgs{i});
    end
    fprintf('\n');
    error('run_all_tests:testsFailed', '%d of %d tests failed.', nFailed, nTests);
else
    fprintf('All %d tests passed successfully!\n', nPassed);
end
