classdef test_IVIMmodelfit < matlab.unittest.TestCase
    % TEST_IVIMMODELFIT Unit tests for the biexponential IVIM fitting function.
    %
    % Validates IVIMmodelfit.m (in dependencies/) which fits the IVIM
    % biexponential model: S(b) = S0 * [f*exp(-b*D*) + (1-f)*exp(-b*D)]
    % using the segmented fitting approach ('seg' method).
    %
    % Tests cover:
    %   - Ideal (noiseless) signal recovery for a single voxel
    %   - Spatially varying parameter maps (2x2x1 grid)
    %   - Output dimension verification ([Ny, Nx, Nz, 4])
    %   - Noise robustness at SNR ~100
    %   - Mask-based selective fitting (only masked voxels fitted)

    methods(TestMethodSetup)
        function addDependenciesToPath(testCase)
            % Add the dependencies/ folder to the MATLAB path so that
            % IVIMmodelfit.m is accessible. Uses PathFixture in MATLAB
            % (auto-cleaned on teardown) or plain addpath in Octave.
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
        function testIdealSignalSeg(testCase)
            % Verify that the segmented IVIM fit recovers ground-truth parameters
            % from a noiseless biexponential signal. This is the most basic
            % sanity check: if fitting fails on ideal data, something is
            % fundamentally broken.
            %
            % Ground truth: D=1.0e-3, f=0.2, D*=15e-3 (typical pancreas values)
            % Expected tolerances: D within 5%, f within 10%, D* within 15%
            % (D* is inherently less stable due to fast-decaying component)
            bvals = [0; 10; 20; 50; 100; 200; 400; 800];
            true_D = 1.0e-3;
            true_f = 0.2;
            true_Dstar = 15e-3;
            S0 = 100;

            % Generate signal: S = S0 * (f * exp(-b * D*) + (1-f) * exp(-b * D))
            signal = S0 * (true_f * exp(-bvals * true_Dstar) + (1 - true_f) * exp(-bvals * true_D));

            % Reshape signal to 4D array [Ny, Nx, Nz, Nb] and replicate to a
            % 3x2x2 volume so the fitter has enough voxels to process.
            dwi = repmat(reshape(signal, [1, 1, 1, length(bvals)]), [3, 2, 2, 1]);

            % bthr=200: b-values above this threshold are used for the
            % monoexponential D fit in the first segmented step.
            opts.bthr = 200;
            opts.dispprog = false;

            % Run fit
            maps = IVIMmodelfit(dwi, bvals, 'seg', true(size(dwi,1), size(dwi,2), size(dwi,3)), opts);

            % Output maps is 4D: [Ny, Nx, Nz, nPars] where nPars=4
            % Parameter order: maps(:,:,:,1)=D, maps(:,:,:,2)=S0,
            %                  maps(:,:,:,3)=f, maps(:,:,:,4)=D*
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
            % Verify spatial fidelity: each voxel in a 2x2x1 grid has different
            % IVIM parameters. The fitted maps should recover the correct spatial
            % pattern, ensuring the fitter does not mix up voxel locations.
            bvals = [0; 10; 20; 50; 100; 200; 400; 800];
            S0 = 100;

            % Define spatially varying ground-truth parameters across a 2x2 grid
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

            % Verify each voxel's fitted parameters match the correct ground truth.
            % D* has the loosest tolerance (50%) because the perfusion component
            % decays very quickly and is sensitive to the segmented fit threshold.
            testCase.verifyEqual(fitted_D, true_D, 'RelTol', 0.10, ...
                'Fitted D map should spatially match the ground truth map.');
            testCase.verifyEqual(fitted_f, true_f, 'RelTol', 0.20, ...
                'Fitted f map should spatially match the ground truth map.');
            testCase.verifyEqual(fitted_Dstar, true_Dstar, 'RelTol', 0.50, ...
                'Fitted Dstar map should spatially match the ground truth map.');
        end

        function testDimensions(testCase)
            % Verify that the output parameter maps have the expected shape
            % [Ny, Nx, Nz, 4] regardless of input volume size. This guards
            % against reshape or indexing bugs in the fitting loop.
            Ny = 4; Nx = 3; Nz = 2; Nb = 8;
            bvals = [0; 10; 20; 50; 100; 200; 400; 800];

            % Use realistic biexponential signals to avoid NaN or degenerate fits
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
            % Verify that the segmented IVIM fit remains stable under moderate
            % Gaussian noise (SNR ~100). D and f should be within 15-20% of
            % truth; D* is allowed 50% error because the perfusion pseudo-
            % diffusion coefficient is notoriously noise-sensitive in IVIM.
            rng(42); % Fixed seed for deterministic noise
            bvals = [0; 10; 20; 50; 100; 200; 400; 800; 1000];
            true_D = 1.2e-3;
            true_f = 0.15;
            true_Dstar = 20e-3;
            S0 = 1000;

            signal_ideal = S0 * (true_f * exp(-bvals * true_Dstar) + (1 - true_f) * exp(-bvals * true_D));

            % Add Gaussian noise with sigma=10 (SNR = S0/sigma = 1000/10 = 100)
            noise = 10 * randn(size(signal_ideal));
            signal_noisy = signal_ideal + noise;
            signal_noisy = abs(signal_noisy); % Magnitude signal (Rician-like)

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
            % Verify that IVIMmodelfit respects the binary mask argument:
            % voxels inside the mask should be fitted (non-NaN), while voxels
            % outside the mask should remain NaN. This is critical for
            % restricting fitting to the GTV region in clinical use.
            bvals = [0; 10; 20; 50; 100; 200; 400; 800];
            S0 = 100; true_D = 1e-3; true_f = 0.2; true_Dstar = 15e-3;
            s_vec = S0 * (true_f * exp(-bvals * true_Dstar) + (1 - true_f) * exp(-bvals * true_D));

            dwi = zeros(3, 2, 2, 8);
            for y = 1:2
                for x = 1:2
                    dwi(y, x, 1, :) = reshape(s_vec, [1, 1, 1, 8]);
                end
            end

            % Create a mask selecting only the single voxel at (1,1,1)
            mask = false(3, 2, 2);
            mask(1, 1) = true;  % Only this voxel should be fitted

            opts.bthr = 200;
            opts.dispprog = false;
            maps = IVIMmodelfit(dwi, bvals, 'seg', mask, opts);

            % The masked voxel (1,1,1) should have a valid D value; all others NaN
            testCase.verifyFalse(isnan(maps(1, 1, 1, 1)), 'Masked voxel should be fitted.');
            testCase.verifyTrue(isnan(maps(1, 2, 1, 1)), 'Unmasked voxel should be NaN.');
            testCase.verifyTrue(isnan(maps(2, 1, 1, 1)), 'Unmasked voxel should be NaN.');
            testCase.verifyTrue(isnan(maps(2, 2, 1, 1)), 'Unmasked voxel should be NaN.');
        end
    end
end
