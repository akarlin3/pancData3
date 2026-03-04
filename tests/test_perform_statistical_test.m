classdef test_perform_statistical_test < matlab.unittest.TestCase
    % TEST_PERFORM_STATISTICAL_TEST Unit tests for robust rank-sum tests
    methods(TestMethodSetup)
        function setupPath(testCase)
            folder = fullfile(fileparts(mfilename('fullpath')), '..', 'utils');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(folder));
        end
    end

    methods(Test)
        function testRanksumBasic(testCase)
            data = [1; 2; 3; 10; 11; 12];
            groups = [1; 1; 1; 2; 2; 2];

            p = perform_statistical_test(data, groups, 'ranksum');

            % Check p-value for ranksum
            % We expect a small p-value for completely separated distributions
            testCase.verifyTrue(~isnan(p));
            testCase.verifyLessThanOrEqual(p, 0.1);
        end

        function testRanksumIdentical(testCase)
            data = [1; 2; 3; 1; 2; 3];
            groups = [1; 1; 1; 2; 2; 2];

            p = perform_statistical_test(data, groups, 'ranksum');

            testCase.verifyTrue(~isnan(p));
            testCase.verifyGreaterThan(p, 0.9);
        end

        function testRanksumInsufficientData(testCase)
            data = [1];
            groups = [1];

            p = perform_statistical_test(data, groups, 'ranksum');

            testCase.verifyTrue(isnan(p));
        end

        function testRanksumOneGroup(testCase)
            data = [1; 2; 3; 4];
            groups = [1; 1; 1; 1];

            p = perform_statistical_test(data, groups, 'ranksum');

            testCase.verifyTrue(isnan(p));
        end

        function testRanksumWithNaN(testCase)
            data = [1; 2; 3; NaN; 10; 11; 12];
            groups = [1; 1; 1; 1; 2; 2; 2];

            p = perform_statistical_test(data, groups, 'ranksum');

            testCase.verifyTrue(~isnan(p));
        end
        function testRanksumTinyGroups(testCase)
            % Groups with fewer than 3 observations each should return NaN
            data = [1; 2; 10; 11];
            groups = [1; 1; 2; 2];   % n1=2, n2=2

            p = perform_statistical_test(data, groups, 'ranksum');

            testCase.verifyTrue(isnan(p), ...
                'Rank-sum with fewer than 3 per group should return NaN');
        end

    end
end
