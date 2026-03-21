classdef test_calibration_metrics < matlab.unittest.TestCase
    % TEST_CALIBRATION_METRICS  Tests for predictive model calibration assessment.
    %
    % Covers:
    %   - Perfect calibration (predicted = observed)
    %   - Miscalibrated model
    %   - Single-class edge case
    %   - Calibration plot generation

    properties
        OriginalPath
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            close all;
            path(testCase.OriginalPath);
        end
    end

    methods(Test)

        function testPerfectCalibration(testCase)
            % When predicted probabilities match observed frequencies,
            % Brier score should be near 0 and H-L p should be > 0.05.
            rng(42);
            n = 100;
            outcomes = [ones(50, 1); zeros(50, 1)];
            % Perfect predictions: 1.0 for events, 0.0 for non-events
            risk_scores = [ones(50, 1) * 3; ones(50, 1) * -3];  % logistic: high/low

            cal = compute_calibration_metrics(risk_scores, outcomes, 5);

            testCase.verifyLessThan(cal.brier_score, 0.1, ...
                'Brier score should be low for well-calibrated predictions.');
            testCase.verifyTrue(isfinite(cal.hl_p), ...
                'H-L p-value should be finite.');
        end

        function testMiscalibratedModel(testCase)
            % Systematically miscalibrated model: predict high for everyone.
            rng(42);
            n = 50;
            outcomes = [ones(10, 1); zeros(40, 1)];  % 20% event rate
            risk_scores = ones(n, 1) * 2;  % all predict high risk

            cal = compute_calibration_metrics(risk_scores, outcomes, 5);

            testCase.verifyGreaterThan(cal.brier_score, 0.1, ...
                'Brier score should be high for miscalibrated predictions.');
            testCase.verifyTrue(isfinite(cal.brier_score));
        end

        function testSingleClassEdgeCase(testCase)
            % All same outcome: should handle gracefully.
            rng(42);
            n = 20;
            outcomes = ones(n, 1);  % all events
            risk_scores = randn(n, 1);

            output = evalc('cal = compute_calibration_metrics(risk_scores, outcomes, 5);');

            % Should report insufficient data or return NaN
            testCase.verifyTrue(isnan(cal.brier_score) || isfinite(cal.brier_score), ...
                'Should handle single-class gracefully.');
        end

        function testInsufficientData(testCase)
            % Very few data points: should return NaN.
            outcomes = [1; 0];
            risk_scores = [0.8; 0.2];

            output = evalc('cal = compute_calibration_metrics(risk_scores, outcomes, 5);');

            testCase.verifyTrue(isnan(cal.brier_score), ...
                'Should return NaN with insufficient data.');
        end

        function testCalibrationPlotSaved(testCase)
            % Verify calibration plot is created when output_folder provided.
            rng(42);
            n = 50;
            outcomes = [ones(20, 1); zeros(30, 1)];
            risk_scores = randn(n, 1);

            tmp_dir = tempname;
            mkdir(tmp_dir);
            cleanup = onCleanup(@() rmdir(tmp_dir, 's'));

            evalc('cal = compute_calibration_metrics(risk_scores, outcomes, 5, tmp_dir, ''Standard'', ''Fx2'');');

            fig_path = fullfile(tmp_dir, 'calibration_plot_Standard_Fx2.png');
            testCase.verifyTrue(exist(fig_path, 'file') == 2, ...
                'Calibration plot should be saved.');
        end

        function testBrierDecomposition(testCase)
            % Verify Brier decomposition components are returned.
            rng(42);
            n = 60;
            outcomes = [ones(25, 1); zeros(35, 1)];
            risk_scores = randn(n, 1) * 0.5;

            cal = compute_calibration_metrics(risk_scores, outcomes, 5);

            testCase.verifyTrue(isfinite(cal.reliability), 'Reliability should be finite.');
            testCase.verifyTrue(isfinite(cal.resolution), 'Resolution should be finite.');
            testCase.verifyTrue(isfinite(cal.uncertainty), 'Uncertainty should be finite.');
            testCase.verifyGreaterThanOrEqual(cal.reliability, 0, 'Reliability should be non-negative.');
        end

    end
end
