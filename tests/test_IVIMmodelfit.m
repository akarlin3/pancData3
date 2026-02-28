classdef test_IVIMmodelfit < matlab.unittest.TestCase
    % TEST_IVIMMODELFIT Unit test for the biexponential IVIM fitting function.

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
        function testIdealSignalSeg(testCase)
            % Test with a single voxel and perfect biexponential decay
            bvals = [0; 10; 20; 50; 100; 200; 400; 800];
            true_D = 1.0e-3;
            true_f = 0.2;
            true_Dstar = 15e-3;
            S0 = 100;

            % Generate signal: S = S0 * (f * exp(-b * D*) + (1-f) * exp(-b * D))
            signal = S0 * (true_f * exp(-bvals * true_Dstar) + (1 - true_f) * exp(-bvals * true_D));

            % Reshape to [Ny, Nx, Nz, Nb] -> [1, 1, 1, 8]
            dwi = repmat(reshape(signal, [1, 1, 1, length(bvals)]), [3, 2, 2, 1]);

            % Options
            opts.bthr = 200;
            opts.dispprog = false;

            % Run fit
            maps = IVIMmodelfit(dwi, bvals, 'seg', true(size(dwi,1), size(dwi,2), size(dwi,3)), opts);

            % maps is 4D: [Ny, Nx, Nz, nPars]
            % pars = {'D','S0','f','Dstar'}
            fitted_D = maps(1, 1, 1, 1);
            fitted_S0 = maps(1, 1, 1, 2);
            fitted_f = maps(1, 1, 1, 3);
            fitted_Dstar = maps(1, 1, 1, 4);

            % Verify
            testCase.verifyEqual(fitted_D, true_D, 'RelTol', 0.05, ...
                'Fitted D should match the ground truth for ideal signal.');
            testCase.verifyEqual(fitted_f, true_f, 'RelTol', 0.10, ...
                'Fitted f should match the ground truth for ideal signal.');
            testCase.verifyEqual(fitted_Dstar, true_Dstar, 'RelTol', 0.15, ...
                'Fitted Dstar should match the ground truth for ideal signal.');
            testCase.verifyEqual(fitted_S0, S0, 'RelTol', 0.05, ...
                'Fitted S0 should match the ground truth for ideal signal.');
        end

        function testVaryingParameters(testCase)
            % Test with 2x2x1 grid with different parameters
            bvals = [0; 10; 20; 50; 100; 200; 400; 800];
            S0 = 100;

            % Define true parameters
            true_D = [0.8e-3, 1.2e-3; 1.0e-3, 1.5e-3];
            true_f = [0.1, 0.2; 0.3, 0.15];
            true_Dstar = [10e-3, 20e-3; 15e-3, 25e-3];

            dwi = zeros(3, 2, 2, 8);
            for y = 1:2
                for x = 1:2
                    D = true_D(y, x);
                    f = true_f(y, x);
                    Dstar = true_Dstar(y, x);
                    s_vec = S0 * (f * exp(-bvals * Dstar) + (1 - f) * exp(-bvals * D));
                    dwi(y, x, 1, :) = reshape(s_vec, [1, 1, 1, 8]);
                end
            end

            % Run fit
            opts.bthr = 200;
            opts.dispprog = false;
            maps = IVIMmodelfit(dwi, bvals, 'seg', true(size(dwi,1), size(dwi,2), size(dwi,3)), opts);

            fitted_D = maps(1:2,1:2,1,1);
            fitted_S0 = maps(1:2,1:2,1,2);
            fitted_f = maps(1:2,1:2,1,3);
            fitted_Dstar = maps(1:2,1:2,1,4);

            % Verify spatial correspondence
            testCase.verifyEqual(fitted_D, true_D, 'RelTol', 0.10, ...
                'Fitted D map should spatially match the ground truth map.');
            testCase.verifyEqual(fitted_f, true_f, 'RelTol', 0.20, ...
                'Fitted f map should spatially match the ground truth map.');
            testCase.verifyEqual(fitted_Dstar, true_Dstar, 'RelTol', 0.50, ...
                'Fitted Dstar map should spatially match the ground truth map.');
        end

        function testDimensions(testCase)
            % Verify output dimensions are [Ny, Nx, Nz, 4]
            Ny = 4; Nx = 3; Nz = 2; Nb = 8;
            bvals = [0; 10; 20; 50; 100; 200; 400; 800];

            % Use somewhat realistic signals to avoid fit failure / warnings
            S0 = 100; true_D = 1e-3; true_f = 0.2; true_Dstar = 15e-3;
            s_vec = S0 * (true_f * exp(-bvals * true_Dstar) + (1 - true_f) * exp(-bvals * true_D));

            dwi = zeros(Ny, Nx, Nz, Nb);
            for z = 1:Nz
                for y = 1:Ny
                    for x = 1:Nx
                        dwi(y, x, z, :) = reshape(s_vec, [1, 1, 1, Nb]);
                    end
                end
            end

            opts.bthr = 200;
            opts.dispprog = false;
            maps = IVIMmodelfit(dwi, bvals, 'seg', true(size(dwi,1), size(dwi,2), size(dwi,3)), opts);

            expectedSize = [Ny, Nx, Nz, 4];
            testCase.verifyEqual(size(maps), expectedSize, ...
                'Output IVIM maps must have dimensions [Ny, Nx, Nz, 4].');
        end

        function testNoiseRobustness(testCase)
            % Test stability with small noise
            rng(42); % Deterministic noise
            bvals = [0; 10; 20; 50; 100; 200; 400; 800; 1000];
            true_D = 1.2e-3;
            true_f = 0.15;
            true_Dstar = 20e-3;
            S0 = 1000;

            signal_ideal = S0 * (true_f * exp(-bvals * true_Dstar) + (1 - true_f) * exp(-bvals * true_D));

            % Add small Gaussian noise (SNR ~ 100)
            noise = 10 * randn(size(signal_ideal));
            signal_noisy = signal_ideal + noise;
            signal_noisy = abs(signal_noisy);

            dwi = repmat(reshape(signal_noisy, [1, 1, 1, length(bvals)]), [3, 2, 2, 1]);

            opts.bthr = 200;
            opts.dispprog = false;
            maps = IVIMmodelfit(dwi, bvals, 'seg', true(size(dwi,1), size(dwi,2), size(dwi,3)), opts);

            fitted_D = maps(1, 1, 1, 1);
            fitted_S0 = maps(1, 1, 1, 2);
            fitted_f = maps(1, 1, 1, 3);
            fitted_Dstar = maps(1, 1, 1, 4);

            % Allow 15% error due to noise
            rel_error_D = abs(fitted_D - true_D) / true_D;
            rel_error_f = abs(fitted_f - true_f) / true_f;
            % Dstar is notoriously sensitive to noise, so we allow a bit more margin
            rel_error_Dstar = abs(fitted_Dstar - true_Dstar) / true_Dstar;

            testCase.verifyLessThan(rel_error_D, 0.15, ...
                sprintf('Fitted D (%.4e) should be within 15%% of truth (%.4e) with moderate noise.', fitted_D, true_D));
            testCase.verifyLessThan(rel_error_f, 0.20, ...
                sprintf('Fitted f (%.4e) should be within 20%% of truth (%.4e) with moderate noise.', fitted_f, true_f));
            testCase.verifyLessThan(rel_error_Dstar, 0.50, ...
                sprintf('Fitted Dstar (%.4e) should be within 50%% of truth (%.4e) with moderate noise.', fitted_Dstar, true_Dstar));
        end

        function testMasking(testCase)
            % Test with a mask to verify it only fits in the masked region
            bvals = [0; 10; 20; 50; 100; 200; 400; 800];
            S0 = 100; true_D = 1e-3; true_f = 0.2; true_Dstar = 15e-3;
            s_vec = S0 * (true_f * exp(-bvals * true_Dstar) + (1 - true_f) * exp(-bvals * true_D));

            dwi = zeros(3, 2, 2, 8);
            for y = 1:2
                for x = 1:2
                    dwi(y, x, 1, :) = reshape(s_vec, [1, 1, 1, 8]);
                end
            end

            % Only fit the first voxel
            mask = false(3, 2, 2);
            mask(1, 1) = true;

            opts.bthr = 200;
            opts.dispprog = false;
            maps = IVIMmodelfit(dwi, bvals, 'seg', mask, opts);

            % Expected maps: Voxel 1,1 should be fitted, rest should be NaN
            testCase.verifyFalse(isnan(maps(1, 1, 1, 1)), 'Masked voxel should be fitted.');
            testCase.verifyTrue(isnan(maps(1, 2, 1, 1)), 'Unmasked voxel should be NaN.');
            testCase.verifyTrue(isnan(maps(2, 1, 1, 1)), 'Unmasked voxel should be NaN.');
            testCase.verifyTrue(isnan(maps(2, 2, 1, 1)), 'Unmasked voxel should be NaN.');
        end
    end
end
