classdef test_fit_adc_mono < matlab.unittest.TestCase
    % TEST_FIT_ADC_MONO Unit tests for the monoexponential ADC fitting function.
    %
    % Validates fit_adc_mono.m (in dependencies/) which fits the mono-
    % exponential diffusion model: S(b) = S0 * exp(-b * ADC) using ordinary
    % least squares on log-transformed signal data.
    %
    % Tests cover:
    %   - Ideal (noiseless) signal recovery for a single voxel
    %   - Spatially varying ADC map (2x2x1 grid)
    %   - Output dimension verification ([Ny, Nx, Nz])
    %   - Noise robustness at SNR ~100

    methods(TestMethodSetup)
        function addDependenciesToPath(testCase)
            % Add the dependencies/ folder so fit_adc_mono.m is accessible.
            % Uses PathFixture in MATLAB (auto-cleaned) or plain addpath in Octave.
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            depPath = fullfile(repoRoot, 'dependencies');
            if exist('OCTAVE_VERSION', 'builtin')
                addpath(depPath);
            else
                import matlab.unittest.fixtures.PathFixture
                testCase.applyFixture(PathFixture(depPath));
            end
        end
    end

    methods(Test)
        function testIdealSignal(testCase)
            % Verify exact ADC recovery from a noiseless monoexponential signal.
            % With perfect data and 3 b-values, log-linear OLS should recover
            % the true ADC to machine precision (AbsTol = 1e-6 mm^2/s).
            bvals = [0; 500; 1000];
            true_adc = 1.5e-3; % mm^2/s, typical for soft tissue
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
            % Verify spatial fidelity: a 2x2x1 grid where each voxel has a
            % different ADC. The fitted map should match the ground truth at
            % each spatial location, ensuring the fitter indexes voxels correctly.
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
            % Verify the output ADC map is 3D [Ny, Nx, Nz] (one scalar per voxel),
            % not 4D like the input. Guards against shape/squeeze bugs.
            Ny = 4; Nx = 3; Nz = 2; Nb = 3;
            dwi = rand(Ny, Nx, Nz, Nb) + 10; % Offset by 10 to ensure positive signal for log
            bvals = [0; 100; 200];

            adc_map = fit_adc_mono(dwi, bvals);

            expectedSize = [Ny, Nx, Nz];
            testCase.verifyEqual(size(adc_map), expectedSize, ...
                'Output ADC map must have dimensions [Ny, Nx, Nz].');
        end

        function testNoiseRobustness(testCase)
            % Verify that the log-linear ADC fit remains accurate under moderate
            % Gaussian noise (SNR = S0/sigma = 1000/10 = 100). The fitted ADC
            % should be within 15% of the ground truth. Note: OLS on log-
            % transformed data has a known noise-floor bias at low SNR, but at
            % SNR=100 this effect is negligible.
            rng(42); % Fixed seed for deterministic noise
            bvals = [0; 50; 400; 800];
            true_adc = 1.2e-3;
            S0 = 1000;

            signal_ideal = S0 * exp(-bvals * true_adc);

            % Add small Gaussian noise (SNR ~ 100)
            noise = 10 * randn(size(signal_ideal));
            signal_noisy = signal_ideal + noise;

            % Take absolute value to simulate magnitude MRI signal (Rician-like);
            % required because log-linear fit needs positive values.
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
