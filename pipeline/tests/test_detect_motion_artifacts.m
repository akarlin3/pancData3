classdef test_detect_motion_artifacts < matlab.unittest.TestCase
    % TEST_DETECT_MOTION_ARTIFACTS  Tests for DWI motion artifact detection.
    %
    % Covers:
    %   - Clean synthetic DWI passes
    %   - Volume with injected signal dropout is flagged
    %   - Single-slice input edge case
    %   - Output structure validation

    properties
        OriginalPath
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.OriginalPath);
        end
    end

    methods(Test)

        function testCleanDWIPasses(testCase)
            % Clean synthetic DWI with exponential decay should pass.
            rng(42);
            sz = [20, 20, 5];
            b_values = [0, 50, 100, 200, 400, 800];
            n_b = length(b_values);

            % Generate synthetic DWI: S = S0 * exp(-b * ADC)
            S0 = 1000 * ones(sz);
            adc = 1e-3;  % typical ADC
            dwi_4d = zeros([sz, n_b]);
            for bi = 1:n_b
                dwi_4d(:,:,:,bi) = S0 .* exp(-b_values(bi) * adc) + randn(sz) * 10;
            end
            dwi_4d = max(dwi_4d, 0);
            mask = true(sz);

            motion = detect_motion_artifacts(dwi_4d, b_values, mask);

            testCase.verifyEqual(motion.n_flagged, 0, ...
                'Clean DWI should have no flagged volumes.');
            testCase.verifyEqual(length(motion.per_volume), n_b, ...
                'Should have per-volume metrics for each b-value.');
        end

        function testSignalDropoutFlagged(testCase)
            % Volume with injected signal dropout should be flagged.
            rng(42);
            sz = [20, 20, 5];
            b_values = [0, 100, 400, 800];
            n_b = length(b_values);

            S0 = 500 * ones(sz);
            dwi_4d = zeros([sz, n_b]);
            for bi = 1:n_b
                dwi_4d(:,:,:,bi) = S0 .* exp(-b_values(bi) * 1e-3) + randn(sz) * 5;
            end
            dwi_4d = max(dwi_4d, 1);

            % Inject severe dropout in volume 3
            dwi_4d(:,:,:,3) = 0.5;  % near-zero signal

            mask = true(sz);

            motion = detect_motion_artifacts(dwi_4d, b_values, mask);

            testCase.verifyTrue(motion.flagged(3), ...
                'Volume with signal dropout should be flagged.');
            testCase.verifyGreaterThan(motion.n_flagged, 0, ...
                'At least one volume should be flagged.');
        end

        function testSingleSliceInput(testCase)
            % Single-slice (2D) input should work.
            rng(42);
            sz = [20, 20, 1];
            b_values = [0, 200, 800];
            n_b = length(b_values);

            dwi_4d = zeros([sz, n_b]);
            for bi = 1:n_b
                dwi_4d(:,:,:,bi) = 500 * exp(-b_values(bi) * 1e-3) + randn(sz) * 5;
            end
            dwi_4d = max(dwi_4d, 1);
            mask = true(sz);

            motion = detect_motion_artifacts(dwi_4d, b_values, mask);

            testCase.verifyEqual(length(motion.per_volume), n_b);
            testCase.verifyTrue(isfinite(motion.noise_floor));
        end

        function testOutputStructure(testCase)
            % Verify output struct has all expected fields.
            sz = [10, 10, 3];
            dwi_4d = rand([sz, 4]) * 100;
            b_values = [0, 100, 400, 800];
            mask = true(sz);

            motion = detect_motion_artifacts(dwi_4d, b_values, mask);

            testCase.verifyTrue(isfield(motion, 'per_volume'));
            testCase.verifyTrue(isfield(motion, 'flagged'));
            testCase.verifyTrue(isfield(motion, 'n_flagged'));
            testCase.verifyTrue(isfield(motion, 'noise_floor'));
            testCase.verifyEqual(length(motion.flagged), 4);
        end

        function testEmptyMaskReturnsEarly(testCase)
            % Empty mask should return without flagging anything.
            dwi_4d = rand(10, 10, 3, 4) * 100;
            b_values = [0, 100, 400, 800];
            mask = false(10, 10, 3);

            motion = detect_motion_artifacts(dwi_4d, b_values, mask);

            testCase.verifyEqual(motion.n_flagged, 0);
        end

    end
end
