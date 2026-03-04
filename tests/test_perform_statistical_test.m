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
            data = [1; 2; 3; 4; 5; 10; 11; 12; 13; 14];
            groups = [1; 1; 1; 1; 1; 2; 2; 2; 2; 2];

            p = perform_statistical_test(data, groups, 'ranksum');

            % Check p-value for ranksum
            % We expect a small p-value for completely separated distributions
            testCase.verifyTrue(~isnan(p));
            testCase.verifyLessThanOrEqual(p, 0.05);
        end

        function testRanksumIdentical(testCase)
            data = [1; 2; 3; 4; 5; 1; 2; 3; 4; 5];
            groups = [1; 1; 1; 1; 1; 2; 2; 2; 2; 2];

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
            data = [1; 2; 3; 4; 5; NaN; 10; 11; 12; 13; 14];
            groups = [1; 1; 1; 1; 1; 1; 2; 2; 2; 2; 2];

            p = perform_statistical_test(data, groups, 'ranksum');

            testCase.verifyTrue(~isnan(p));
        end
        function testRanksumTinyGroups(testCase)
            % Groups with fewer than 5 observations each should return NaN
            data = [1; 2; 3; 4; 10; 11; 12; 13];
            groups = [1; 1; 1; 1; 2; 2; 2; 2];   % n1=4, n2=4

            p = perform_statistical_test(data, groups, 'ranksum');

            testCase.verifyTrue(isnan(p), ...
                'Rank-sum with fewer than 5 per group should return NaN');
        end

    end
end
