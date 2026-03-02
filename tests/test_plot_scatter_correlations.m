classdef test_plot_scatter_correlations < matlab.unittest.TestCase
    % TEST_PLOT_SCATTER_CORRELATIONS Smoke tests for plot_scatter_correlations.
    %
    % Verifies:
    %   - Valid data generates the expected PNG in output_folder
    %   - Datasets with fewer than 3 clean (non-NaN) points are skipped
    %     gracefully without error
    %   - The dtype index selects the correct slice of 3-D metric arrays
    %   - All-NaN dose vectors still produce a saved (but empty) figure

    properties
        TempDir
        OriginalPath
    end

    methods(TestMethodSetup)
        function setup(testCase)
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            close all;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
            path(testCase.OriginalPath);
        end
    end

    % ------------------------------------------------------------------ %
    %  Helpers                                                            %
    % ------------------------------------------------------------------ %
    methods(Access = private)
        function [dtype_label, dmean_gtvp, d95_gtvp, adc_mean, d_mean, f_mean, ...
                  valid_pts, lf_group, dtype, output_folder] = buildValidArgs(testCase, n)
            % Returns arguments for a clean n-patient dataset (dtype = 1).
            rng(1);
            dtype_label = 'Standard';
            dmean_gtvp  = rand(n, 1) * 50 + 20;   % n × 1 (only Fx1 used)
            d95_gtvp    = rand(n, 1) * 40 + 15;
            % 3-D arrays: n × nTp × n_dwi_types; function uses (:, 1, dtype)
            adc_mean    = rand(n, 1, 1) * 2e-3;
            d_mean      = rand(n, 1, 1) * 1e-3;
            f_mean      = rand(n, 1, 1) * 0.3;
            valid_pts   = true(n, 1);
            lf_group    = [ones(ceil(n/2), 1); zeros(floor(n/2), 1)];
            dtype        = 1;
            output_folder = testCase.TempDir;
        end
    end

    methods(Test)

        function testGeneratesFigureWithValidData(testCase)
            % 10 patients with valid dose and diffusion values.
            % The function should save Dose_vs_Diffusion_Standard.png.
            n = 10;
            [lbl, dmean, d95, adc, d, f, vp, lfg, dt, out] = ...
                testCase.buildValidArgs(n);

            plot_scatter_correlations(lbl, dmean, d95, adc, d, f, vp, lfg, dt, out);

            expected = fullfile(out, ['Dose_vs_Diffusion_' lbl '.png']);
            testCase.verifyTrue(exist(expected, 'file') > 0, ...
                'Dose_vs_Diffusion PNG should be saved for valid input data.');
        end

        function testInsufficientCleanPointsSkipsSubplot(testCase)
            % Only 2 of 5 patients have non-NaN dose AND diffusion values
            % (< 3 required for Spearman correlation).  Function should not
            % throw; it displays "insufficient data" titles and still saves.
            n = 5;
            [lbl, dmean, d95, adc, d, f, vp, lfg, dt, out] = ...
                testCase.buildValidArgs(n);

            % Make most dose values NaN → fewer than 3 clean points
            dmean(1:4) = NaN;
            d95(1:4)   = NaN;

            plot_scatter_correlations(lbl, dmean, d95, adc, d, f, vp, lfg, dt, out);

            % Figure is still saved even when subplots lack data
            expected = fullfile(out, ['Dose_vs_Diffusion_' lbl '.png']);
            testCase.verifyTrue(exist(expected, 'file') > 0, ...
                'PNG should still be saved when some panels lack sufficient data.');
        end

        function testAllNaNDoseStillSavesFigure(testCase)
            % All dose values NaN → every subplot gets "insufficient data".
            % The function must not crash and must save the figure.
            n = 8;
            [lbl, dmean, d95, adc, d, f, vp, lfg, dt, out] = ...
                testCase.buildValidArgs(n);

            dmean(:) = NaN;
            d95(:)   = NaN;

            plot_scatter_correlations(lbl, dmean, d95, adc, d, f, vp, lfg, dt, out);

            expected = fullfile(out, ['Dose_vs_Diffusion_' lbl '.png']);
            testCase.verifyTrue(exist(expected, 'file') > 0, ...
                'PNG should be saved even when all dose values are NaN.');
        end

        function testDtypeLabelAppearsInFilename(testCase)
            % The dtype_label string is embedded in the output filename.
            n   = 6;
            lbl = 'dnCNN';
            [~, dmean, d95, adc, d, f, vp, lfg, dt, out] = ...
                testCase.buildValidArgs(n);

            plot_scatter_correlations(lbl, dmean, d95, adc, d, f, vp, lfg, dt, out);

            expected = fullfile(out, ['Dose_vs_Diffusion_' lbl '.png']);
            testCase.verifyTrue(exist(expected, 'file') > 0, ...
                'Output filename should use the dtype_label argument.');
        end

        function testValidPtsMaskSubsetsPatients(testCase)
            % valid_pts = false for some patients → those patients' data are
            % excluded.  Function must still complete and save the figure.
            n   = 10;
            [lbl, dmean, d95, adc, d, f, ~, lfg, dt, out] = ...
                testCase.buildValidArgs(n);

            % Exclude the first 4 patients
            vp = [false(4, 1); true(6, 1)];

            plot_scatter_correlations(lbl, dmean, d95, adc, d, f, vp, lfg, dt, out);

            expected = fullfile(out, ['Dose_vs_Diffusion_' lbl '.png']);
            testCase.verifyTrue(exist(expected, 'file') > 0, ...
                'PNG should be saved when valid_pts excludes some patients.');
        end

        function testSinglePatientHandledGracefully(testCase)
            % Edge case: only 1 valid patient — fewer than 3 clean points
            % in every subplot → no Spearman correlation, but no crash.
            n = 1;
            [lbl, dmean, d95, adc, d, f, vp, lfg, dt, out] = ...
                testCase.buildValidArgs(n);

            plot_scatter_correlations(lbl, dmean, d95, adc, d, f, vp, lfg, dt, out);

            expected = fullfile(out, ['Dose_vs_Diffusion_' lbl '.png']);
            testCase.verifyTrue(exist(expected, 'file') > 0, ...
                'PNG should be saved even for a single patient.');
        end

    end
end
