classdef test_detect_baseline_outliers < matlab.unittest.TestCase
% TEST_DETECT_BASELINE_OUTLIERS — Unit tests for detect_baseline_outliers.m
%
% Validates outcome-blinded IQR outlier detection:
%   - Normal data (no outliers)
%   - Extreme outlier detection
%   - Output dimensions
%   - Zero-IQR and few-samples edge cases
%   - Multi-metric aggregation

    properties
        origPath
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            testCase.origPath = path;
            baseDir = fullfile(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods (TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.origPath);
        end
    end

    methods (Test)
        function testNoOutliers(testCase)
            % Normal distribution data tightly clustered; no outliers expected.
            rng(42);
            n = 50;
            col = randn(n, 1) * 0.1 + 5;  % tight cluster around 5
            lf = zeros(n, 1);
            lf(1:15) = 1;
            evalc('[is_outlier, n_out] = detect_baseline_outliers({col}, {''metric1''}, lf);');
            testCase.verifyEqual(n_out, 0, ...
                'Tightly clustered data should have no outliers.');
        end

        function testExtremeOutlierDetected(testCase)
            % Insert a value far beyond 3*IQR; it should be flagged.
            rng(42);
            n = 50;
            col = randn(n, 1) * 0.1 + 5;
            col(1) = 5 + 100 * iqr(col);  % extreme outlier
            lf = zeros(n, 1);
            lf(1:15) = 1;
            evalc('[is_outlier, n_out] = detect_baseline_outliers({col}, {''metric1''}, lf);');
            testCase.verifyTrue(is_outlier(1), ...
                'Extreme value should be flagged as outlier.');
            testCase.verifyGreaterThan(n_out, 0, ...
                'At least one outlier should be detected.');
        end

        function testOutputDimensions(testCase)
            rng(42);
            n = 30;
            col = randn(n, 1);
            lf = zeros(n, 1);
            evalc('[is_outlier, ~] = detect_baseline_outliers({col}, {''m1''}, lf);');
            testCase.verifyEqual(size(is_outlier), size(lf), ...
                'is_outlier should have the same size as lf.');
            testCase.verifyTrue(islogical(is_outlier), ...
                'is_outlier should be a logical vector.');
        end

        function testZeroIQRSkipped(testCase)
            % Constant column has IQR=0; should not flag anything.
            n = 20;
            col = ones(n, 1) * 3;
            lf = zeros(n, 1);
            evalc('[is_outlier, n_out] = detect_baseline_outliers({col}, {''const''}, lf);');
            testCase.verifyEqual(n_out, 0, ...
                'Constant column (IQR=0) should not flag any outliers.');
        end

        function testMultipleMetrics(testCase)
            % Outlier in one metric should flag the patient overall.
            rng(42);
            n = 50;
            col1 = randn(n, 1) * 0.1 + 5;
            col2 = randn(n, 1) * 0.1 + 10;
            col2(3) = 10 + 100 * iqr(col2);  % outlier in metric 2 only
            lf = zeros(n, 1);
            evalc('[is_outlier, ~] = detect_baseline_outliers({col1, col2}, {''m1'', ''m2''}, lf);');
            testCase.verifyTrue(is_outlier(3), ...
                'Patient 3 should be flagged due to outlier in metric 2.');
        end

        function testFewSamplesSkipped(testCase)
            % Column with <3 non-NaN values should be skipped.
            n = 10;
            col = NaN(n, 1);
            col(1) = 5;
            col(2) = 100;  % would be outlier if column were analysed
            lf = zeros(n, 1);
            evalc('[is_outlier, n_out] = detect_baseline_outliers({col}, {''sparse''}, lf);');
            testCase.verifyEqual(n_out, 0, ...
                'Column with <3 non-NaN values should be skipped entirely.');
        end
    end
end
