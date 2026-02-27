classdef test_fit_adc_mono < matlab.unittest.TestCase
    % TEST_FIT_ADC_MONO Unit test for the monoexponential ADC fitting function.

    methods(TestMethodSetup)
        function addDependenciesToPath(testCase)
            % Ensure the dependencies folder is on the path
            import matlab.unittest.fixtures.PathFixture
            % Get the root of the repository relative to this file
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            depPath = fullfile(repoRoot, 'dependencies');

            % Use PathFixture to add path temporarily for the test
            testCase.applyFixture(PathFixture(depPath));
        end
    end

    methods(Test)
        function testIdealSignal(testCase)
            % Test with a single voxel and perfect exponential decay
            bvals = [0; 500; 1000];
            true_adc = 1.5e-3;
            S0 = 100;

            % Generate signal: S = S0 * exp(-b * ADC)
            signal = S0 * exp(-bvals * true_adc);

            % Reshape to [Ny, Nx, Nz, Nb] -> [1, 1, 1, 3]
            dwi = reshape(signal, [1, 1, 1, 3]);

            % Run fit
            adc_map = fit_adc_mono(dwi, bvals);

            % Verify
            testCase.verifyEqual(adc_map, true_adc, 'AbsTol', 1e-6, ...
                'Fitted ADC should match the ground truth for ideal signal.');
        end

        function testVaryingADC(testCase)
            % Test with 2x2x1 grid with different ADCs
            bvals = [0; 800];
            S0 = 100;

            % Define true ADCs
            % Voxel (1,1): 0.5e-3
            % Voxel (1,2): 1.0e-3
            % Voxel (2,1): 1.5e-3
            % Voxel (2,2): 2.0e-3
            true_adcs = [0.5e-3, 1.0e-3; 1.5e-3, 2.0e-3];

            dwi = zeros(2, 2, 1, 2);
            for y = 1:2
                for x = 1:2
                    adc = true_adcs(y, x);
                    s_vec = S0 * exp(-bvals * adc);
                    dwi(y, x, 1, :) = reshape(s_vec, [1, 1, 1, 2]);
                end
            end

            % Run fit
            adc_map = fit_adc_mono(dwi, bvals);

            % Verify spatial correspondence
            testCase.verifyEqual(adc_map, true_adcs, 'AbsTol', 1e-6, ...
                'Fitted ADC map should spatially match the ground truth map.');
        end

        function testDimensions(testCase)
            % Verify output dimensions are [Ny, Nx, Nz]
            Ny = 4; Nx = 3; Nz = 2; Nb = 3;
            dwi = rand(Ny, Nx, Nz, Nb) + 10; % Ensure positive signal
            bvals = [0; 100; 200];

            adc_map = fit_adc_mono(dwi, bvals);

            expectedSize = [Ny, Nx, Nz];
            testCase.verifyEqual(size(adc_map), expectedSize, ...
                'Output ADC map must have dimensions [Ny, Nx, Nz].');
        end

        function testNoiseRobustness(testCase)
            % Test stability with small noise
            % Note: Monoexponential fitting on log signal is sensitive to noise at low SNR,
            % but with high SNR and small noise, it should be stable.

            rng(42); % Deterministic noise
            bvals = [0; 50; 400; 800];
            true_adc = 1.2e-3;
            S0 = 1000;

            signal_ideal = S0 * exp(-bvals * true_adc);

            % Add small Gaussian noise (SNR ~ 100)
            noise = 10 * randn(size(signal_ideal));
            signal_noisy = signal_ideal + noise;

            % Ensure positive signal for log
            signal_noisy = abs(signal_noisy);

            dwi = reshape(signal_noisy, [1, 1, 1, 4]);

            adc_map = fit_adc_mono(dwi, bvals);

            % Allow 15% error due to noise and OLS bias on log-transformed data
            rel_error = abs(adc_map - true_adc) / true_adc;

            testCase.verifyLessThan(rel_error, 0.15, ...
                sprintf('Fitted ADC (%.4e) should be within 15%% of truth (%.4e) with moderate noise.', adc_map, true_adc));
        end
    end
end
