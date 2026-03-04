% run_all_tests.m
% Master test runner for the MATLAB DWI pipeline repository
% Discovers and runs all tests (class-based and function-based) in the 'tests' directory.
% Outputs test results and asserts success for CI pipeline compatibility.
%
% When the Parallel Computing Toolbox is available (and not in preflight
% mode), tests that are known to be parallel-safe are executed via
% runInParallel for faster wall-clock time.  All remaining tests run
% sequentially as before.

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

% --- Parallel execution configuration ---
% Tests listed here are safe to run concurrently: they do not open diary
% files, do not call core modules that open diary files, and do not share
% mutable filesystem state.  When adding a new test, include it here ONLY
% if it satisfies all three criteria.
parallel_safe_classes = { ...
    'test_escape_shell_arg', ...
    'test_init_scan_structs', ...
    'test_perform_statistical_test', ...
    'test_source_code_standards', ...
    'test_statistical_methods', ...
    'test_build_td_panel', ...
    'test_scale_td_panel', ...
    'test_compute_summary_metrics', ...
    'test_text_progress_bar', ...
    'test_IVIMmodelfit', ...
    'test_fit_adc_mono', ...
    'test_fit_models', ...
    'test_dvh', ...
    'test_apply_dir_mask_propagation', ...
    'test_calculate_subvolume_metrics', ...
    'test_landmark_cindex', ...
    'test_landmark_cindex_mock', ...
    'test_perf_knn', ...
    'test_fix_verify', ...
    'test_octave', ...
    'benchmark_filter_collinear', ...
    'benchmark_make_grouped_folds', ...
    'benchmark_metrics_opt', ...
    'benchmark_scale_td', ...
    'test_opt', ...
    'test_opt2', ...
    'test_perf', ...
    'test_accumarray', ...
    'test_make_grouped_folds' ...
};

% 1. Discover all tests in the tests/ directory (including subdirectories)
suite = TestSuite.fromFolder(testsDir, 'IncludingSubfolders', true);

if isempty(suite)
    error('No tests found in the %s directory.', testsDir);
end

% 2. Partition suite into parallel-safe and serial groups
is_parallel = false(1, numel(suite));
for k = 1:numel(suite)
    % Extract the class name from the test name (format: 'ClassName/MethodName')
    tokens = strsplit(suite(k).Name, '/');
    className = tokens{1};
    % Handle subfolder-qualified names (e.g. 'benchmarks.test_opt')
    parts = strsplit(className, '.');
    className = parts{end};
    is_parallel(k) = ismember(className, parallel_safe_classes);
end
parallel_suite = suite(is_parallel);
serial_suite   = suite(~is_parallel);

fprintf('Discovered %d tests (%d parallel-safe, %d serial).\n\n', ...
    numel(suite), numel(parallel_suite), numel(serial_suite));

% --- Waitbar setup (GUI environments only) ---
has_display = false;
if ~exist('OCTAVE_VERSION', 'builtin')
    try
        has_display = usejava('desktop') || ...
                      (~isempty(getenv('DISPLAY')) && usejava('jvm'));
    catch
        has_display = false;
    end
end
hWaitbar = [];
if has_display
    hWaitbar = waitbar(0, 'Initializing test suite...', ...
        'Name', 'Test Suite Progress', ...
        'Visible', 'on');
end

% 3. Check whether parallel execution is available
is_preflight = strcmp(getenv('PIPELINE_PREFLIGHT_ACTIVE'), '1');

can_run_parallel = false;
if ~is_preflight && ~exist('OCTAVE_VERSION', 'builtin') && ~isempty(parallel_suite)
    has_pct = license('test', 'Distrib_Computing_Toolbox');
    % Check if runInParallel method exists (R2018a+)
    has_method = false;
    try
        m = ?matlab.unittest.TestRunner;
        method_names = {m.MethodList.Name};
        has_method = ismember('runInParallel', method_names);
    catch
        has_method = false;
    end
    can_run_parallel = has_pct && has_method;
end

% 4. Configure code coverage for core and utils directories
foldersToCover = {coreDir, utilsDir};

disp('===================================================');
disp('   Running Tests...                                ');
disp('===================================================');

% 5. Run the test suite
if can_run_parallel
    % --- Phase 1: parallel-safe tests via runInParallel ---
    fprintf('Running %d parallel-safe tests with runInParallel...\n', numel(parallel_suite));
    if ~isempty(hWaitbar) && isvalid(hWaitbar)
        waitbar(0, hWaitbar, sprintf('Running %d parallel-safe tests...', numel(parallel_suite)));
    end
    par_runner = TestRunner.withTextOutput();
    parallel_results = par_runner.runInParallel(parallel_suite);

    parallel_done = numel(parallel_suite);
    if ~isempty(hWaitbar) && isvalid(hWaitbar)
        waitbar(parallel_done / numel(suite), hWaitbar, ...
            sprintf('Parallel phase complete (%d/%d). Starting serial tests...', ...
            parallel_done, numel(suite)));
    end

    % --- Phase 2: serial tests sequentially ---
    fprintf('\nRunning %d serial tests sequentially...\n', numel(serial_suite));
    ser_runner = TestRunner.withTextOutput();
    ser_runner.addPlugin(ProgressBarPlugin(numel(serial_suite)));
    if ~isempty(hWaitbar) && isvalid(hWaitbar)
        ser_runner.addPlugin(WaitbarProgressPlugin(hWaitbar, numel(suite), parallel_done));
    end

    % Add coverage plugin only to the serial runner — CodeCoveragePlugin is
    % not compatible with runInParallel.  Serial tests exercise all core
    % modules, so coverage remains meaningful.
    if exist('matlab.unittest.plugins.CodeCoveragePlugin', 'class')
        coveragePlugin = CodeCoveragePlugin.forFolder(string(foldersToCover));
        ser_runner.addPlugin(coveragePlugin);
        disp('Code coverage plugin added (serial tests only).');
    end

    serial_results = ser_runner.run(serial_suite);

    % Merge results from both phases
    results = [parallel_results, serial_results];
else
    % --- Fallback: fully sequential execution (original behavior) ---
    if ~isempty(parallel_suite) && ~is_preflight
        disp('Parallel execution not available; running all tests sequentially.');
    end

    runner = TestRunner.withTextOutput();
    runner.addPlugin(ProgressBarPlugin(numel(suite)));
    if ~isempty(hWaitbar) && isvalid(hWaitbar)
        runner.addPlugin(WaitbarProgressPlugin(hWaitbar, numel(suite), 0));
    end

    if is_preflight
        disp('Pre-flight mode — skipping code coverage plugin.');
    elseif exist('matlab.unittest.plugins.CodeCoveragePlugin', 'class')
        coveragePlugin = CodeCoveragePlugin.forFolder(string(foldersToCover));
        runner.addPlugin(coveragePlugin);
        disp('Code coverage plugin added. Coverage will be generated for /core and /utils.');
    else
        disp('CodeCoveragePlugin not available in this MATLAB version.');
    end

    results = runner.run(suite);
end

disp('===================================================');
disp('   Test Execution Completed                        ');
disp('===================================================');

% 6. Write failure summary file if any tests failed
failureSummaryFile = fullfile(repoRoot, 'failure_summary.out');
failedIdx = find([results.Failed]);
if ~isempty(failedIdx)
    fid = fopen(failureSummaryFile, 'w');
    if fid ~= -1
        fprintf(fid, 'Test Failure Summary — %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fprintf(fid, '===================================================\n');
        fprintf(fid, '%d of %d tests failed.\n\n', numel(failedIdx), numel(results));
        for fi = 1:numel(failedIdx)
            idx = failedIdx(fi);
            fprintf(fid, '--- %s ---\n', results(idx).Name);
            if ~isempty(results(idx).Details) && isfield(results(idx).Details, 'DiagnosticRecord')
                diagRec = results(idx).Details.DiagnosticRecord;
                for di = 1:numel(diagRec)
                    report = diagRec(di).Report;
                    if ~isempty(report)
                        fprintf(fid, '%s\n', report);
                    end
                end
            end
            fprintf(fid, '\n');
        end
        fclose(fid);
        fprintf('Failure summary written to: %s\n', failureSummaryFile);
    end
elseif exist(failureSummaryFile, 'file')
    % Clean up stale summary from a previous failing run
    delete(failureSummaryFile);
end

% 7. Assert success (Throws an error and returns non-zero exit code if any test fails)
% In CI environments running MATLAB with -batch, this will fail the step appropriately.
% assertSuccess(results) was introduced in R2020a, providing backward compatibility
if exist('assertSuccess', 'file') || ismethod(results, 'assertSuccess')
    assertSuccess(results);
else
    if any([results.Failed])
        error('One or more tests failed.');
    end
end

% Close waitbar
if ~isempty(hWaitbar) && isvalid(hWaitbar)
    close(hWaitbar);
end

% Restore figure visibility
set(0, 'DefaultFigureVisible', 'on');

if standalone_diary
    diary off;
end

disp('All tests passed successfully!');
