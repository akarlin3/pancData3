classdef test_perform_statistical_test < matlab.unittest.TestCase
    % TEST_PERFORM_STATISTICAL_TEST Unit tests for perform_statistical_test.
    %
    % Tests the NaN-safe Wilcoxon rank-sum wrapper used throughout the
    % metrics_stats_comparisons pipeline.  Validates correct p-value
    % computation, graceful NaN handling, and edge-case guards (single
    % group, insufficient sample size, tiny groups).

    methods(TestMethodSetup)
        function setupPath(testCase)
            % Add the utils/ directory so perform_statistical_test is on path.
            folder = fullfile(fileparts(mfilename('fullpath')), '..', 'utils');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(folder));
        end
    end

    methods(Test)
        function testRanksumBasic(testCase)
            % Two completely separated groups (1-5 vs 10-14) should
            % produce a significant p-value (p <= 0.05) from Wilcoxon
            % rank-sum, confirming the function detects true differences.
            data = [1; 2; 3; 4; 5; 10; 11; 12; 13; 14];
            groups = [1; 1; 1; 1; 1; 2; 2; 2; 2; 2];

            p = perform_statistical_test(data, groups, 'ranksum');

            testCase.verifyTrue(~isnan(p));
            testCase.verifyLessThanOrEqual(p, 0.05);
        end

        function testRanksumIdentical(testCase)
            % Two groups with identical distributions ([1..5] vs [1..5])
            % should yield a very high p-value (> 0.9), confirming no
            % false positives for indistinguishable groups.
            data = [1; 2; 3; 4; 5; 1; 2; 3; 4; 5];
            groups = [1; 1; 1; 1; 1; 2; 2; 2; 2; 2];

            p = perform_statistical_test(data, groups, 'ranksum');

            testCase.verifyTrue(~isnan(p));
            testCase.verifyGreaterThan(p, 0.9);
        end

        function testRanksumInsufficientData(testCase)
            % A single observation cannot support a rank-sum test.
            % The function should return NaN rather than error.
            data = [1];
            groups = [1];

            p = perform_statistical_test(data, groups, 'ranksum');

            testCase.verifyTrue(isnan(p));
        end

        function testRanksumOneGroup(testCase)
            % All observations belong to group 1 (no group 2).  A two-
            % sample test is impossible, so the function should return NaN.
            data = [1; 2; 3; 4];
            groups = [1; 1; 1; 1];

            p = perform_statistical_test(data, groups, 'ranksum');

            testCase.verifyTrue(isnan(p));
        end

        function testRanksumWithNaN(testCase)
            % One NaN value in group 1 should be stripped before testing.
            % With 5 valid observations per group the test should still
            % produce a valid (non-NaN) p-value.
            data = [1; 2; 3; 4; 5; NaN; 10; 11; 12; 13; 14];
            groups = [1; 1; 1; 1; 1; 1; 2; 2; 2; 2; 2];

            p = perform_statistical_test(data, groups, 'ranksum');

            testCase.verifyTrue(~isnan(p));
        end
        function testRanksumTinyGroups(testCase)
            % Groups with fewer than 5 observations each should return NaN
            % to avoid unreliable p-values from very small samples.
            data = [1; 2; 3; 4; 10; 11; 12; 13];
            groups = [1; 1; 1; 1; 2; 2; 2; 2];   % n1=4, n2=4

            p = perform_statistical_test(data, groups, 'ranksum');

            testCase.verifyTrue(isnan(p), ...
                'Rank-sum with fewer than 5 per group should return NaN');
        end

    end
end
