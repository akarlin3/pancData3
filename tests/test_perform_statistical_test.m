classdef test_perform_statistical_test < matlab.unittest.TestCase
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
            testCase.verifyLessThan(p, 0.1);
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
    end
end
