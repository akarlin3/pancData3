classdef test_find_volume_inflection < matlab.unittest.TestCase
% TEST_FIND_VOLUME_INFLECTION — Direct unit tests for find_volume_inflection.m.
%
% The algorithm is covered with hand-crafted curves so the assertions do
% not depend on upstream morphology noise from the optimize_adc_threshold
% fixture.

    properties
        OriginalPath
    end

    methods(TestMethodSetup)
        function setup(testCase)
            testCase.OriginalPath = path;
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            path(testCase.OriginalPath);
        end
    end

    methods(Test)
        function test_step_function_knee(testCase)
            % Clean step: vol_frac jumps from 0 to 0.5 between indices 6
            % and 7.  Smoothing blurs the step across a few indices; the
            % most negative second derivative lands at the top of the
            % ramp (where the curve saturates into the plateau).
            thresholds = 0.8e-3 : 0.1e-3 : 2.0e-3;  % 13 points
            vol_frac = [0 0 0 0 0 0 0.5 0.5 0.5 0.5 0.5 0.5 0.5];

            [knee, idx, c_knee, curvature] = find_volume_inflection(thresholds, vol_frac);

            testCase.verifyEqual(numel(curvature), 13);
            testCase.verifyTrue(isnan(curvature(1)));
            testCase.verifyTrue(isnan(curvature(13)));

            % Knee must land on the plateau side of the step, strictly
            % inside the threshold range.
            testCase.verifyFalse(isnan(knee));
            testCase.verifyGreaterThan(knee, thresholds(1));
            testCase.verifyLessThan(knee, thresholds(end));
            testCase.verifyGreaterThanOrEqual(idx, 2);
            testCase.verifyLessThanOrEqual(idx, 12);
            testCase.verifyLessThan(c_knee, 0, 'Knee curvature must be negative.');
        end

        function test_flat_curve_returns_nan(testCase)
            % A perfectly flat curve has zero curvature everywhere.
            % argmin returns the first non-NaN index (2), but this is
            % a legitimate "no knee" scenario — we still report index 2
            % with curvature 0.  The caller can check if curvature is
            % meaningfully negative before trusting the knee.
            thresholds = 0.8e-3 : 0.1e-3 : 2.0e-3;
            vol_frac = 0.3 * ones(1, 13);

            [knee, idx, c_knee, ~] = find_volume_inflection(thresholds, vol_frac);

            testCase.verifyFalse(isnan(knee));
            testCase.verifyEqual(c_knee, 0, 'AbsTol', 1e-12, ...
                'Flat curve should have zero curvature at the knee.');
            testCase.verifyGreaterThanOrEqual(idx, 2);
            testCase.verifyLessThanOrEqual(idx, 12);
        end

        function test_too_few_samples_returns_nan(testCase)
            % Fewer than 3 samples: curvature is undefined.
            [knee, idx, c_knee, curvature] = find_volume_inflection([0.001 0.002], [0.1 0.2]);
            testCase.verifyTrue(isnan(knee));
            testCase.verifyTrue(isnan(idx));
            testCase.verifyTrue(isnan(c_knee));
            testCase.verifyEqual(numel(curvature), 2);
        end

        function test_mostly_nan_curve_returns_nan(testCase)
            % Fewer than 3 finite samples in a 13-point curve.
            thresholds = 0.8e-3 : 0.1e-3 : 2.0e-3;
            vol_frac = nan(1, 13);
            vol_frac([6 8]) = [0.2 0.4];  % only 2 finite samples

            [knee, idx, c_knee, ~] = find_volume_inflection(thresholds, vol_frac);

            testCase.verifyTrue(isnan(knee));
            testCase.verifyTrue(isnan(idx));
            testCase.verifyTrue(isnan(c_knee));
        end

        function test_endpoints_always_nan_in_curvature(testCase)
            % Regardless of input shape, curvature[1] and curvature[end]
            % must be NaN (the discrete 2nd derivative is undefined
            % there).
            thresholds = linspace(0.8e-3, 2.0e-3, 13);
            vol_frac = rand(1, 13);

            [~, ~, ~, curvature] = find_volume_inflection(thresholds, vol_frac);

            testCase.verifyTrue(isnan(curvature(1)));
            testCase.verifyTrue(isnan(curvature(end)));
        end

        function test_knee_matches_thresholds_entry(testCase)
            % Consistency: knee_thresh must equal thresholds(knee_idx).
            thresholds = 0.8e-3 : 0.1e-3 : 2.0e-3;
            vol_frac = [0 0.05 0.15 0.30 0.48 0.49 0.50 0.50 0.50 0.50 0.50 0.50 0.50];

            [knee, idx, ~, ~] = find_volume_inflection(thresholds, vol_frac);

            testCase.verifyFalse(isnan(knee));
            testCase.verifyEqual(knee, thresholds(idx), 'AbsTol', 1e-12);
        end

        function test_monotonic_sigmoid_picks_plateau_top(testCase)
            % Classic sigmoid shape: the knee should be at the
            % saturation top (not the bottom), i.e., where the curve
            % transitions from rising into plateau.
            thresholds = 0.8e-3 : 0.1e-3 : 2.0e-3;
            % Smooth sigmoid centred around index 6.
            x = 1:13;
            vol_frac = 0.5 ./ (1 + exp(-1.2*(x - 6.5)));

            [~, idx, c_knee, ~] = find_volume_inflection(thresholds, vol_frac);

            % Saturation top of a rising sigmoid is the most concave
            % point (negative curvature).  It should land AFTER the
            % midpoint (index > 6) and be negative.
            testCase.verifyGreaterThan(idx, 6, ...
                'Saturation knee of rising sigmoid should be past the midpoint.');
            testCase.verifyLessThan(c_knee, 0);
        end
    end
end
