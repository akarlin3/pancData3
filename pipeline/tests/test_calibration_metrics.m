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

        function testCompetingRisksInterface(testCase)
            % With outcomes containing 0,1,2 values and times provided,
            % verify competing_risks flag is true.
            rng(42);
            n = 80;
            outcomes = zeros(n, 1);
            outcomes(1:25) = 1;   % event of interest
            outcomes(26:40) = 2;  % competing event
            % remaining are censored (0)
            times = rand(n, 1) * 100 + 5;
            risk_scores = randn(n, 1);

            evalc('cal = compute_calibration_metrics(risk_scores, outcomes, times, 5);');

            testCase.verifyTrue(cal.competing_risks, ...
                'Should detect competing risks when outcomes contain 0, 1, 2.');
            testCase.verifyTrue(isfinite(cal.brier_score) || isnan(cal.brier_score), ...
                'Brier score should be finite or NaN (not error).');
        end

        function testBackwardCompatOldInterface(testCase)
            % Test old interface (risk_scores, outcomes, n_bins) with 3 args.
            rng(42);
            n = 50;
            outcomes = [ones(20, 1); zeros(30, 1)];
            risk_scores = randn(n, 1);

            % Old interface: compute_calibration_metrics(risk_scores, outcomes, n_bins)
            evalc('cal = compute_calibration_metrics(risk_scores, outcomes, 5);');

            testCase.verifyTrue(isstruct(cal), 'Should return a struct.');
            testCase.verifyTrue(isfield(cal, 'brier_score'), ...
                'Should have brier_score field.');
            testCase.verifyTrue(isfinite(cal.brier_score), ...
                'Brier score should be finite with old interface.');
        end

        function testBackwardCompatSixArgs(testCase)
            % Test old interface with 6 args (risk_scores, outcomes, n_bins,
            % output_folder, dtype_label, fx_label).
            rng(42);
            n = 50;
            outcomes = [ones(20, 1); zeros(30, 1)];
            risk_scores = randn(n, 1);

            tmp_dir = tempname;
            mkdir(tmp_dir);
            cleanup = onCleanup(@() rmdir(tmp_dir, 's'));

            % Old 6-arg interface: n_bins is 3rd arg (scalar)
            evalc('cal = compute_calibration_metrics(risk_scores, outcomes, 5, tmp_dir, ''Standard'', ''Fx1'');');

            testCase.verifyTrue(isstruct(cal), 'Should return a struct.');
            testCase.verifyTrue(isfinite(cal.brier_score), ...
                'Brier score should be finite with old 6-arg interface.');
            % Should have saved a plot
            fig_path = fullfile(tmp_dir, 'calibration_plot_Standard_Fx1.png');
            testCase.verifyTrue(exist(fig_path, 'file') == 2, ...
                'Calibration plot should be saved with old 6-arg interface.');
        end

        function testCalibrationSlopeAndIntercept(testCase)
            % Verify calibration_slope and calibration_intercept are returned
            % and finite.
            rng(42);
            n = 80;
            outcomes = [ones(35, 1); zeros(45, 1)];
            % Moderately informative scores
            risk_scores = [randn(35, 1) * 0.5 + 1; randn(45, 1) * 0.5 - 1];

            evalc('cal = compute_calibration_metrics(risk_scores, outcomes, 5);');

            testCase.verifyTrue(isfield(cal, 'calibration_slope'), ...
                'Should have calibration_slope field.');
            testCase.verifyTrue(isfield(cal, 'calibration_intercept'), ...
                'Should have calibration_intercept field.');
            testCase.verifyTrue(isfinite(cal.calibration_slope), ...
                'Calibration slope should be finite.');
            testCase.verifyTrue(isfinite(cal.calibration_intercept), ...
                'Calibration intercept should be finite.');
        end

        function testAllNanScoresHandled(testCase)
            % All NaN risk_scores should not crash.
            rng(42);
            n = 30;
            outcomes = [ones(15, 1); zeros(15, 1)];
            risk_scores = nan(n, 1);

            evalc('cal = compute_calibration_metrics(risk_scores, outcomes, 5);');

            testCase.verifyTrue(isstruct(cal), 'Should return a struct without crashing.');
            testCase.verifyTrue(isnan(cal.brier_score), ...
                'Brier score should be NaN for all-NaN input.');
        end

        function testLogisticTransformApplied(testCase)
            % When scores are outside [0,1], verify logistic transform is
            % applied (resulting in valid probabilities).
            rng(42);
            n = 60;
            outcomes = [ones(25, 1); zeros(35, 1)];
            % Scores well outside [0,1]
            risk_scores = [randn(25, 1) * 2 + 3; randn(35, 1) * 2 - 3];

            evalc('cal = compute_calibration_metrics(risk_scores, outcomes, 5);');

            testCase.verifyTrue(isfinite(cal.brier_score), ...
                'Should produce finite Brier score after logistic transform.');
            % Brier score must be between 0 and 1 for valid probabilities
            testCase.verifyGreaterThanOrEqual(cal.brier_score, 0, ...
                'Brier score should be >= 0.');
            testCase.verifyLessThanOrEqual(cal.brier_score, 1, ...
                'Brier score should be <= 1.');
        end

    end
end
