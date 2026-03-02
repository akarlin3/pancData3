% run_all_tests.m
% Master test runner for the MATLAB DWI pipeline repository
% Discovers and runs all tests (class-based and function-based) in the 'tests' directory.
% Outputs test results and asserts success for CI pipeline compatibility.

% Define repository root and critical directories
repoRoot = fileparts(mfilename('fullpath'));
testsDir = fullfile(repoRoot, 'tests');
coreDir = fullfile(repoRoot, 'core');
utilsDir = fullfile(repoRoot, 'utils');

% Ensure directories are on the MATLAB path
addpath(coreDir);
addpath(utilsDir);
addpath(repoRoot);

disp('===================================================');
disp('   MATLAB CI Test Runner: Initializing Suite       ');
disp('===================================================');

import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.CodeCoveragePlugin;

% 1. Discover all tests in the tests/ directory
suite = TestSuite.fromFolder(testsDir);

if isempty(suite)
    error('No tests found in the %s directory.', testsDir);
else
    fprintf('Discovered %d tests.\n\n', numel(suite));
end

% 2. Create a test runner with text output
runner = TestRunner.withTextOutput();

% 3. Configure code coverage for core and utils directories
% Use Cobertura format or standard coverage plugin that outputs to the console/file.
foldersToCover = {coreDir, utilsDir};

% Add the CodeCoveragePlugin to generate an HTML report or standard coverage
if exist('matlab.unittest.plugins.CodeCoveragePlugin', 'class')
    % Create a coverage plugin for the specific folders
    % For CI, producing Cobertura or just a console report is best.
    % In modern MATLAB (R2017b+), we can use:
    coveragePlugin = CodeCoveragePlugin.forFolder(foldersToCover);
    runner.addPlugin(coveragePlugin);
    disp('Code coverage plugin added. Coverage will be generated for /core and /utils.');
else
    disp('CodeCoveragePlugin not available in this MATLAB version.');
end

disp('===================================================');
disp('   Running Tests...                                ');
disp('===================================================');

% 4. Run the test suite
results = runner.run(suite);

disp('===================================================');
disp('   Test Execution Completed                        ');
disp('===================================================');

% 5. Assert success (Throws an error and returns non-zero exit code if any test fails)
% In CI environments running MATLAB with -batch, this will fail the step appropriately.
% assertSuccess(results) was introduced in R2020a, providing backward compatibility
if exist('assertSuccess', 'file') || ismethod(results, 'assertSuccess')
    assertSuccess(results);
else
    if any([results.Failed])
        error('One or more tests failed.');
    end
end

disp('All tests passed successfully!');
