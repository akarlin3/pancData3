classdef test_fit_models < matlab.unittest.TestCase
    % TEST_FIT_MODELS Unit tests for the fit_models wrapper function.
    %
    % Validates:
    %   - 3D-to-1D flattening and reconstruction roundtrip
    %   - Padding to even voxel count (odd mask sum)
    %   - Insufficient b-values above threshold -> graceful NaN return
    %   - Zero output from IVIM fit replaced with NaN
    %   - Non-positive signal values handled correctly
    %   - Empty mask (no true voxels) returns all-NaN maps

    methods(TestMethodSetup)
        function addPaths(testCase)
            % Add core/ (fit_models.m) and dependencies/ (IVIMmodelfit.m,
            % fit_adc_mono.m) to the MATLAB path.
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'core'));
            addpath(fullfile(repoRoot, 'dependencies'));
        end
    end

    methods(Test)

        function testOutputDimensionsMatchInput(testCase)
            % Output maps must have [Ny, Nx, Nz] dimensions matching the
            % spatial dimensions of the input DWI volume, regardless of
            % the number of b-values.
            bvals = [0; 50; 100; 200; 400; 800];
            Ny = 4; Nx = 3; Nz = 2;
            % Ground-truth IVIM parameters for synthetic biexponential signal
            S0 = 100; true_D = 1e-3; true_f = 0.15; true_Dstar = 15e-3;

            % Generate biexponential (IVIM) signal: S = S0 * [f*exp(-b*D*) + (1-f)*exp(-b*D)]
            sig = S0 * (true_f * exp(-bvals * true_Dstar) + ...
                        (1 - true_f) * exp(-bvals * true_D));
            dwi = repmat(reshape(sig, [1, 1, 1, numel(bvals)]), [Ny, Nx, Nz, 1]);

            mask = true(Ny, Nx, Nz);
            opts.bthr = 100;

            [d_map, f_map, dstar_map, adc_map] = fit_models(dwi, bvals, mask, opts);

            testCase.verifyEqual(size(d_map), [Ny, Nx, Nz]);
            testCase.verifyEqual(size(f_map), [Ny, Nx, Nz]);
            testCase.verifyEqual(size(dstar_map), [Ny, Nx, Nz]);
            testCase.verifyEqual(size(adc_map), [Ny, Nx, Nz]);
        end

        function testFlattenReconstructRoundtrip(testCase)
            % Verify that the 3D->1D->3D reconstruction preserves spatial
            % correspondence: values at masked voxels should be non-NaN and
            % values at unmasked voxels should be NaN.
            bvals = [0; 50; 100; 200; 400; 800];
            S0 = 100; true_D = 1e-3; true_f = 0.15; true_Dstar = 15e-3;
            sig = S0 * (true_f * exp(-bvals * true_Dstar) + ...
                        (1 - true_f) * exp(-bvals * true_D));

            dwi = repmat(reshape(sig, [1, 1, 1, numel(bvals)]), [4, 4, 2, 1]);

            % Create a partial mask (checkerboard pattern)
            mask = false(4, 4, 2);
            mask(1:2:end, 1:2:end, :) = true;

            opts.bthr = 100;
            [d_map, ~, ~, adc_map] = fit_models(dwi, bvals, mask, opts);

            % Masked voxels should have finite values
            testCase.verifyTrue(all(isfinite(d_map(mask))), ...
                'Masked voxels should have finite D values.');
            testCase.verifyTrue(all(isfinite(adc_map(mask))), ...
                'Masked voxels should have finite ADC values.');

            % Unmasked voxels should be NaN
            testCase.verifyTrue(all(isnan(d_map(~mask))), ...
                'Unmasked voxels should be NaN in D map.');
            testCase.verifyTrue(all(isnan(adc_map(~mask))), ...
                'Unmasked voxels should be NaN in ADC map.');
        end

        function testOddMaskSumPadding(testCase)
            % When the mask has an odd number of true voxels, fit_models
            % pads to even. Verify no crash and correct output.
            bvals = [0; 50; 100; 200; 400; 800];
            S0 = 100; true_D = 1e-3; true_f = 0.15; true_Dstar = 15e-3;
            sig = S0 * (true_f * exp(-bvals * true_Dstar) + ...
                        (1 - true_f) * exp(-bvals * true_D));

            dwi = repmat(reshape(sig, [1, 1, 1, numel(bvals)]), [3, 3, 1, 1]);

            % Odd number of masked voxels = 3
            mask = false(3, 3, 1);
            mask(1, 1) = true;
            mask(2, 2) = true;
            mask(3, 3) = true;
            testCase.verifyEqual(sum(mask(:)), 3, 'Sanity: mask should have 3 true voxels.');

            opts.bthr = 100;
            [d_map, f_map, dstar_map, adc_map] = fit_models(dwi, bvals, mask, opts);

            % Should produce valid output without crash
            % Note: MATLAB drops trailing singleton dims, so size is [3, 3]
            testCase.verifyEqual(size(d_map), [3, 3]);
            testCase.verifyTrue(all(isfinite(d_map(mask))), ...
                'Odd mask voxels should still produce finite D.');
        end

        function testInsufficientBvaluesReturnsNaN(testCase)
            % When fewer than 2 b-values are above the threshold, IVIM
            % maps should be all-NaN, but ADC should still be computed.
            bvals = [0; 50; 100; 150];  % Only 1 value >= 200 threshold
            S0 = 100;
            dwi = repmat(S0, [2, 2, 1, 4]);

            mask = true(2, 2, 1);
            opts.bthr = 200;  % High threshold -> insufficient b-values

            [d_map, f_map, dstar_map, adc_map] = fit_models(dwi, bvals, mask, opts);

            % IVIM maps should be all NaN
            testCase.verifyTrue(all(isnan(d_map(:))), ...
                'D map should be NaN when insufficient b-values for IVIM.');
            testCase.verifyTrue(all(isnan(f_map(:))), ...
                'f map should be NaN when insufficient b-values for IVIM.');
            testCase.verifyTrue(all(isnan(dstar_map(:))), ...
                'D* map should be NaN when insufficient b-values for IVIM.');
        end

        function testEmptyMaskReturnsNaN(testCase)
            % When the mask has no true voxels, all outputs should be NaN.
            bvals = [0; 50; 100; 200; 400; 800];
            S0 = 100;
            dwi = repmat(S0, [3, 3, 1, numel(bvals)]);
            mask = false(3, 3, 1);  % Empty mask

            opts.bthr = 100;
            [d_map, f_map, dstar_map, adc_map] = fit_models(dwi, bvals, mask, opts);

            testCase.verifyTrue(all(isnan(d_map(:))), 'D map should be all-NaN for empty mask.');
            testCase.verifyTrue(all(isnan(f_map(:))), 'f map should be all-NaN for empty mask.');
            testCase.verifyTrue(all(isnan(dstar_map(:))), 'D* map should be all-NaN for empty mask.');
            testCase.verifyTrue(all(isnan(adc_map(:))), 'ADC map should be all-NaN for empty mask.');
        end

        function testNonPositiveSignalADCHandled(testCase)
            % Voxels with zero or negative signal at any b-value should
            % produce NaN in ADC (log of non-positive is undefined).
            bvals = [0; 200; 800];
            dwi = zeros(2, 2, 1, 3);

            % Voxel (1,1): valid positive signal
            dwi(1, 1, 1, :) = [100, 80, 50];
            % Voxel (1,2): zero signal at b=800
            dwi(1, 2, 1, :) = [100, 80, 0];
            % Voxel (2,1): negative signal (noise)
            dwi(2, 1, 1, :) = [100, -5, 50];
            % Voxel (2,2): valid positive signal
            dwi(2, 2, 1, :) = [100, 90, 60];

            mask = true(2, 2, 1);
            opts.bthr = 100;

            [~, ~, ~, adc_map] = fit_models(dwi, bvals, mask, opts);

            % Valid voxels should have finite ADC
            testCase.verifyTrue(isfinite(adc_map(1, 1, 1)), ...
                'Voxel with all-positive signal should have finite ADC.');
            testCase.verifyTrue(isfinite(adc_map(2, 2, 1)), ...
                'Voxel with all-positive signal should have finite ADC.');

            % Non-positive signal voxels should be NaN
            testCase.verifyTrue(isnan(adc_map(1, 2, 1)), ...
                'Voxel with zero signal should have NaN ADC.');
            testCase.verifyTrue(isnan(adc_map(2, 1, 1)), ...
                'Voxel with negative signal should have NaN ADC.');
        end

        function testKnownADCValue(testCase)
            % Verify that fit_models produces correct ADC for a known
            % mono-exponential decay: S = S0 * exp(-b * ADC).
            % With true_adc = 1.5e-3 mm^2/s, the fitted ADC should match
            % within 5% relative tolerance.
            bvals = [0; 200; 400; 800];
            true_adc = 1.5e-3;   % mm^2/s — typical tissue ADC
            S0 = 100;
            sig = S0 * exp(-bvals * true_adc);

            dwi = reshape(sig, [1, 1, 1, 4]);
            mask = true(1, 1, 1);
            opts.bthr = 100;

            [~, ~, ~, adc_map] = fit_models(dwi, bvals, mask, opts);

            testCase.verifyEqual(adc_map(1), true_adc, 'RelTol', 0.05, ...
                'ADC should match known ground truth within 5%%.');
        end

        function testBiexponentialIVIMRecovery(testCase)
            % Test that IVIM segmented fit recovers D from a biexponential
            % signal. D should be close to the tissue diffusion ground truth.
            bvals = [0; 30; 50; 100; 150; 400; 800];
            S0 = 100; true_D = 1.2e-3; true_f = 0.12; true_Dstar = 15e-3;
            sig = S0 * (true_f * exp(-bvals * true_Dstar) + ...
                        (1 - true_f) * exp(-bvals * true_D));

            dwi = repmat(reshape(sig, [1, 1, 1, numel(bvals)]), [2, 2, 1, 1]);
            mask = true(2, 2, 1);
            opts.bthr = 100;

            [d_map, f_map, dstar_map, ~] = fit_models(dwi, bvals, mask, opts);

            % D should be in physically plausible range for tissue
            testCase.verifyTrue(all(d_map(mask) > 0), ...
                'D values should be positive.');
            testCase.verifyTrue(all(d_map(mask) < 5e-3), ...
                'D values should be < 5e-3 mm^2/s for tissue.');

            % f should be in [0, 1] range
            valid_f = f_map(mask & ~isnan(f_map));
            if ~isempty(valid_f)
                testCase.verifyTrue(all(valid_f >= 0 & valid_f <= 1), ...
                    'Perfusion fraction f should be in [0, 1].');
            end

            % D* should be positive
            valid_dstar = dstar_map(mask & ~isnan(dstar_map));
            if ~isempty(valid_dstar)
                testCase.verifyTrue(all(valid_dstar > 0), ...
                    'D* should be positive.');
            end
        end

        function testZeroFitReplacedWithNaN(testCase)
            % When the IVIM segmented fitter returns D=0 (non-positive slope),
            % fit_models should replace all three IVIM outputs with NaN.
            % We test this implicitly: constant signal (no decay) should
            % yield NaN or zero D, which gets converted to NaN.
            bvals = [0; 50; 100; 200; 400; 800];
            S0 = 100;
            % Constant signal across b-values → slope = 0 → D = 0 → NaN
            dwi = repmat(S0, [2, 2, 1, numel(bvals)]);
            mask = true(2, 2, 1);
            opts.bthr = 100;

            [d_map, f_map, dstar_map, ~] = fit_models(dwi, bvals, mask, opts);

            % With constant signal, IVIM should produce NaN (zero→NaN conversion)
            testCase.verifyTrue(all(isnan(d_map(mask))), ...
                'Constant signal should produce NaN D (zero-fit replacement).');
            testCase.verifyTrue(all(isnan(f_map(mask))), ...
                'Constant signal should produce NaN f (zero-fit replacement).');
            testCase.verifyTrue(all(isnan(dstar_map(mask))), ...
                'Constant signal should produce NaN D* (zero-fit replacement).');
        end

        function testNegativeADCReplacedWithNaN(testCase)
            % Signal that increases with b-value (physically impossible)
            % should produce negative ADC, which must be set to NaN.
            bvals = [0; 200; 800];
            dwi = zeros(1, 1, 1, 3);
            % Reversed decay: signal INCREASES with b → negative ADC
            dwi(1,1,1,:) = [50, 80, 120];
            mask = true(1, 1, 1);
            opts.bthr = 100;

            [~, ~, ~, adc_map] = fit_models(dwi, bvals, mask, opts);

            testCase.verifyTrue(isnan(adc_map(1)), ...
                'Negative ADC (increasing signal) should be replaced with NaN.');
        end

        function testNonZeroFirstBvalueErrors(testCase)
            % ADC computation requires bvalues(1)==0 as S0 reference.
            % A non-zero first b-value should error.
            bvals = [50; 200; 800];  % No b=0!
            S0 = 100;
            dwi = repmat(S0, [2, 2, 1, 3]);
            mask = true(2, 2, 1);
            opts.bthr = 100;

            testCase.verifyError(@() fit_models(dwi, bvals, mask, opts), ...
                'fit_models:noB0', ...
                'Should error when first b-value is not 0.');
        end

        function testADCComputedIndependentlyOfIVIM(testCase)
            % ADC should be computed even when IVIM is skipped (insufficient
            % b-values above threshold).
            bvals = [0; 50; 100; 150];
            S0 = 100; true_adc = 1.0e-3;
            sig = S0 * exp(-bvals * true_adc);
            dwi = reshape(sig, [1, 1, 1, 4]);
            mask = true(1, 1, 1);
            opts.bthr = 200;  % Too high → IVIM skipped

            [d_map, f_map, dstar_map, adc_map] = fit_models(dwi, bvals, mask, opts);

            % IVIM should be NaN
            testCase.verifyTrue(isnan(d_map(1)), 'D should be NaN when IVIM skipped.');
            % ADC should still be computed
            testCase.verifyTrue(isfinite(adc_map(1)), ...
                'ADC should be finite even when IVIM is skipped.');
            testCase.verifyEqual(adc_map(1), true_adc, 'RelTol', 0.1, ...
                'ADC should be close to ground truth.');
        end

        function testMultipleVoxelsDifferentSignals(testCase)
            % Verify that voxels with different signal profiles produce
            % different ADC values (not all identical).
            bvals = [0; 200; 400; 800];
            S0 = 100;
            adc1 = 0.5e-3;
            adc2 = 2.0e-3;

            dwi = zeros(1, 2, 1, 4);
            dwi(1,1,1,:) = S0 * exp(-bvals * adc1);
            dwi(1,2,1,:) = S0 * exp(-bvals * adc2);
            mask = true(1, 2, 1);
            opts.bthr = 100;

            [~, ~, ~, adc_map] = fit_models(dwi, bvals, mask, opts);

            testCase.verifyTrue(adc_map(1,1,1) ~= adc_map(1,2,1), ...
                'Different signal profiles should produce different ADC values.');
            testCase.verifyEqual(adc_map(1,1,1), adc1, 'RelTol', 0.05);
            testCase.verifyEqual(adc_map(1,2,1), adc2, 'RelTol', 0.05);
        end

        function testSingleBvalueAboveThreshold(testCase)
            % Only 1 b-value above threshold → IVIM skipped, but ADC works.
            bvals = [0; 50; 200];  % Only 1 value >= 100
            S0 = 100;
            dwi = repmat(S0, [2, 2, 1, 3]);
            mask = true(2, 2, 1);
            opts.bthr = 100;

            [d_map, ~, ~, ~] = fit_models(dwi, bvals, mask, opts);

            testCase.verifyTrue(all(isnan(d_map(:))), ...
                'D should be all NaN with only 1 b-value above threshold.');
        end

    end
end
