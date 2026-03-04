% run_all_tests.m
% Master test runner for the MATLAB DWI pipeline repository
% Discovers and runs all tests (class-based and function-based) in the 'tests' directory.
% Outputs test results and asserts success for CI pipeline compatibility.

% Define repository root and critical directories
% This file lives in tests/, so repo root is one level up
testsDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(testsDir);
coreDir = fullfile(repoRoot, 'core');
utilsDir = fullfile(repoRoot, 'utils');

% Ensure directories are on the MATLAB path
addpath(coreDir);
addpath(utilsDir);
if exist('OCTAVE_VERSION', 'builtin')
    addpath(fullfile(repoRoot, '.octave_compat'));
else
    % Defensive: remove octave_compat from path if it was added by a
    % previous session or startup.m.  The @table shim inside octave_compat
    % conflicts with MATLAB's built-in table class and causes
    % 'MATLAB:dispatcher:InvalidObjtagReuse' errors.
    oc_dir = fullfile(repoRoot, '.octave_compat');
    if contains(path, oc_dir)
        rmpath(genpath(oc_dir));
    end
end
addpath(repoRoot);

% Save test output to test.out when run standalone (not wrapped by execute_all_workflows)
standalone_diary = strcmp(get(0, 'Diary'), 'off');
if standalone_diary
    diary(fullfile(repoRoot, 'test.out'));
end

% Add test subdirectories to path, excluding transient mock_data dirs that
% may be left over from a previous crashed run.  Including them causes
% spurious "Removed ... from the MATLAB path" warnings when test teardown
% deletes them.
test_paths = genpath(testsDir);
test_paths = strsplit(test_paths, pathsep);
test_paths = test_paths(~contains(test_paths, 'mock_data'));
addpath(strjoin(test_paths, pathsep));

% Suppress figure pop-ups during test execution
set(0, 'DefaultFigureVisible', 'off');

disp('===================================================');
disp('   MATLAB CI Test Runner: Initializing Suite       ');
disp('===================================================');

import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.CodeCoveragePlugin;

% 1. Discover all tests in the tests/ directory (including subdirectories)
suite = TestSuite.fromFolder(testsDir, 'IncludingSubfolders', true);

if isempty(suite)
    error('No tests found in the %s directory.', testsDir);
else
    fprintf('Discovered %d tests.\n\n', numel(suite));
end

% 2. Create a test runner with text output
runner = TestRunner.withTextOutput();

% 3. Add dot-style progress bar plugin
runner.addPlugin(ProgressBarPlugin(numel(suite)));

% 4. Configure code coverage for core and utils directories
% Use Cobertura format or standard coverage plugin that outputs to the console/file.
foldersToCover = {coreDir, utilsDir};

% Add the CodeCoveragePlugin to generate an HTML report or standard coverage.
% Skip coverage when invoked from run_dwi_pipeline's pre-flight check — the
% pre-flight only needs pass/fail, and adding coverage during a nested call
% triggers "Simultaneously generating coverage reports" errors in MATLAB.
is_preflight = strcmp(getenv('PIPELINE_PREFLIGHT_ACTIVE'), '1');
if is_preflight
    disp('Pre-flight mode — skipping code coverage plugin.');
elseif exist('matlab.unittest.plugins.CodeCoveragePlugin', 'class')
    coveragePlugin = CodeCoveragePlugin.forFolder(string(foldersToCover));
    runner.addPlugin(coveragePlugin);
    disp('Code coverage plugin added. Coverage will be generated for /core and /utils.');
else
    disp('CodeCoveragePlugin not available in this MATLAB version.');
end

disp('===================================================');
disp('   Running Tests...                                ');
disp('===================================================');

% 5. Run the test suite
results = runner.run(suite);

disp('===================================================');
disp('   Test Execution Completed                        ');
disp('===================================================');

% 6. Assert success (Throws an error and returns non-zero exit code if any test fails)
% In CI environments running MATLAB with -batch, this will fail the step appropriately.
% assertSuccess(results) was introduced in R2020a, providing backward compatibility
if exist('assertSuccess', 'file') || ismethod(results, 'assertSuccess')
    assertSuccess(results);
else
    if any([results.Failed])
        error('One or more tests failed.');
    end
end

% Restore figure visibility
set(0, 'DefaultFigureVisible', 'on');

if standalone_diary
    diary off;
end

disp('All tests passed successfully!');
