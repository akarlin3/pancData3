classdef run_test_local < matlab.unittest.TestCase
% RUN_TEST_LOCAL  Diagnostic harness that runs test_compute_summary_metrics
%   in isolation using a bare TestRunner (no plugins, no coverage).
%
%   Useful for quick local verification of compute_summary_metrics without
%   running the full test suite. Discovers and executes all methods in
%   test_compute_summary_metrics, then asserts that none failed.

    methods (Test)
        function testComputeSummaryMetricsHarness(testCase)
        %TESTCOMPUTESUMMARYMETRICSHARNESS Discover and run all tests in
        %   test_compute_summary_metrics, then verify zero failures.
            import matlab.unittest.TestSuite;
            % Navigate up from diagnostics/ -> tests/ -> repo root
            repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fullfile(repo_root, 'tests'));
            % Build suite from the target test class and run with no plugins
            suite = TestSuite.fromClass(?test_compute_summary_metrics);
            runner = matlab.unittest.TestRunner.withNoPlugins();
            results = runner.run(suite);
            % Verify all tests in the target class passed
            testCase.verifyTrue(all(~[results.Failed]), ...
                'test_compute_summary_metrics had failures');
        end
    end
end
