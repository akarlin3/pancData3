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

% 2. Partition suite into parallel-safe and serial groups.
%    Tests in parallel_safe_classes run via runInParallel (Phase 1).
%    All others run sequentially (Phase 2) to avoid diary/filesystem conflicts.
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

% --- Progress GUI setup (GUI environments only) ---
hGUI = [];
if ProgressGUI.isDisplayAvailable()
    hGUI = ProgressGUI('Running Tests', numel(suite));
end

% 3. Check whether parallel execution is available.
%    Parallel mode is disabled during preflight (quick sanity runs invoked by
%    the pipeline before a full data run) and in Octave (no PCT support).
is_preflight = strcmp(getenv('PIPELINE_PREFLIGHT_ACTIVE'), '1');

can_run_parallel = false;
if ~is_preflight && ~exist('OCTAVE_VERSION', 'builtin') && ~isempty(parallel_suite)
    has_pct = license('test', 'Distrib_Computing_Toolbox');
    % Check if runInParallel method exists (R2018a+) via metaclass introspection
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
    if ~isempty(hGUI) && hGUI.isValid()
        counts = struct('completed', 0, 'total', numel(suite), 'passed', 0, 'failed', 0);
        hGUI.update(0, counts, sprintf('Running %d parallel-safe tests...', numel(parallel_suite)), 'running');
    end
    par_runner = TestRunner.withTextOutput();
    parallel_results = par_runner.runInParallel(parallel_suite);

    parallel_done = numel(parallel_suite);
    if ~isempty(hGUI) && hGUI.isValid()
        par_passed = sum(~[parallel_results.Failed]);
        par_failed = sum([parallel_results.Failed]);
        counts = struct('completed', parallel_done, 'total', numel(suite), ...
                        'passed', par_passed, 'failed', par_failed);
        hGUI.update(parallel_done / numel(suite), counts, ...
            'Parallel phase complete. Starting serial...', 'running');
    end

    % --- Phase 2: serial tests sequentially ---
    fprintf('\nRunning %d serial tests sequentially...\n', numel(serial_suite));
    ser_runner = TestRunner.withTextOutput();
    ser_runner.addPlugin(ProgressBarPlugin(numel(serial_suite)));
    if ~isempty(hGUI) && hGUI.isValid()
        ser_runner.addPlugin(WaitbarProgressPlugin(hGUI, numel(suite), parallel_done));
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
    if ~isempty(hGUI) && hGUI.isValid()
        runner.addPlugin(WaitbarProgressPlugin(hGUI, numel(suite), 0));
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

% 6. Print a detailed results table to the console
failedIdx = find([results.Failed]);
passedIdx = find(~[results.Failed]);
durations = [results.Duration];

fprintf('\n');
fprintf('===================================================\n');
fprintf('   Results by Test Class\n');
fprintf('===================================================\n');

% Group results by class name for a compact per-class summary.
% Each result's Name has the format 'ClassName/MethodName'; we extract the class part.
classNames = cell(1, numel(results));
for ri = 1:numel(results)
    tokens = strsplit(results(ri).Name, '/');
    classNames{ri} = tokens{1};
end
[uniqueClasses, ~, classIdx] = unique(classNames, 'stable');
for ci = 1:numel(uniqueClasses)
    mask = (classIdx(:)' == ci);
    nClass = sum(mask);
    nFailed = sum([results(mask).Failed]);
    classDur = sum([results(mask).Duration]);
    if nFailed == 0
        fprintf('  ✅ %s — %d tests passed (%.1fs)\n', uniqueClasses{ci}, nClass, classDur);
    else
        fprintf('  ❌ %s — %d/%d failed (%.1fs)\n', uniqueClasses{ci}, nFailed, nClass, classDur);
        % List the individual failures within this class
        classMask = find(mask);
        for fi = classMask
            if results(fi).Failed
                tokens = strsplit(results(fi).Name, '/');
                methodName = tokens{end};
                fprintf('       ↳ %s\n', methodName);
            end
        end
    end
end

fprintf('\n===================================================\n');
fprintf('   Summary: %d passed, %d failed, %d total (%.1fs)\n', ...
    numel(passedIdx), numel(failedIdx), numel(results), sum(durations));
fprintf('===================================================\n');

% 7. Write failure summary file if any tests failed
failureSummaryFile = fullfile(repoRoot, 'failure_summary.out');
if ~isempty(failedIdx)
    fprintf('\n  Failed tests:\n');
    for fi = 1:numel(failedIdx)
        fprintf('    ❌ %s\n', results(failedIdx(fi)).Name);
    end
    fprintf('\n');

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

% Close progress GUI before assertSuccess (which may throw on failure)
if ~isempty(hGUI) && hGUI.isValid()
    hGUI.close();
end

% 8. Assert success (Throws an error and returns non-zero exit code if any test fails)
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

% Clean up any saved_files_* folders created by tests in the repo root.
% Only do this when running standalone — when called from
% execute_all_workflows, the parent orchestrator's output folder is also a
% saved_files_* directory and must not be deleted.
if standalone_diary
    stray_dirs = dir(fullfile(repoRoot, 'saved_files_*'));
    for k = 1:numel(stray_dirs)
        if stray_dirs(k).isdir
            rmdir(fullfile(repoRoot, stray_dirs(k).name), 's');
            fprintf('Cleaned up stray test artifact: %s\n', stray_dirs(k).name);
        end
    end
end

% Shut down parallel pool when running standalone (not under execute_all_workflows)
if standalone_diary && ~exist('OCTAVE_VERSION', 'builtin')
    p = gcp('nocreate');
    if ~isempty(p)
        delete(p);
        disp('Parallel pool shut down.');
    end
end

if standalone_diary
    diary off;
end

% 7. Report results — clean summary for interactive use, error for CI/batch
num_failed = nnz([results.Failed]);
if num_failed == 0
    disp('All tests passed successfully!');
else
    % Print a concise failure summary to the console
    fprintf('\n❌ %d of %d tests FAILED:\n', num_failed, numel(results));
    failedIdx = find([results.Failed]);
    for fi = 1:numel(failedIdx)
        fprintf('    - %s\n', results(failedIdx(fi)).Name);
    end
    fprintf('\n');

    % In batch/CI mode, throw an error for non-zero exit code.
    % In interactive mode, just display the summary — no error wall.
    % batchStartupOptionUsed() requires R2019a+; fall back to checking
    % whether the desktop is running.
    if exist('batchStartupOptionUsed', 'builtin') || exist('batchStartupOptionUsed', 'file')
        is_batch = batchStartupOptionUsed();
    else
        is_batch = ~usejava('desktop');
    end
    if is_batch
        error('pancData3:testFailure', '%d of %d tests failed.', num_failed, numel(results));
    end
end

