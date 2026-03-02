classdef run_test_local < matlab.unittest.TestCase
    % Quick local harness that runs test_compute_summary_metrics

    methods (Test)
        function testComputeSummaryMetricsHarness(testCase)
            import matlab.unittest.TestSuite;
            suite = TestSuite.fromClass(?test_compute_summary_metrics);
            runner = matlab.unittest.TestRunner.withNoPlugins();
            results = runner.run(suite);
            testCase.verifyTrue(all(~[results.Failed]), ...
                'test_compute_summary_metrics had failures');
        end
    end
end
