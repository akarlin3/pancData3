classdef test_visualize_results < matlab.unittest.TestCase
% TEST_VISUALIZE_RESULTS  Smoke tests for core/visualize_results.m.
%
% visualize_results depends on full pipeline data structures and calls
% plot_parameter_maps, plot_feature_distributions, and plot_scatter_correlations.
% These tests exercise the function signature and early code paths using
% minimal synthetic data to verify the function does not crash on valid input.

    properties
        OriginalPath
        TempDir
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            if exist('OCTAVE_VERSION', 'builtin')
                addpath(fullfile(baseDir, '.octave_compat'));
            end
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            diary off;
            close all;
            path(testCase.OriginalPath);
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)

        function testFunctionAcceptsFourArgs(testCase)
            % Verify visualize_results accepts (data_vectors, summary_metrics,
            % calculated_results, config_struct) without immediately crashing
            % on the argument count. The function will likely error on missing
            % data or missing sub-functions, but the signature should be valid.
            rng(42);
            nPat = 3;
            nTp = 2;

            % Minimal data_vectors_gtvp struct array
            data_vectors_gtvp = struct();
            for p = 1:nPat
                data_vectors_gtvp(p).id = sprintf('P%03d', p);
                data_vectors_gtvp(p).adc = rand(50, nTp);
            end

            % Minimal summary_metrics
            sm = struct();
            sm.id_list = arrayfun(@(p) sprintf('P%03d', p), 1:nPat, 'UniformOutput', false);
            sm.mrn_list = sm.id_list;
            sm.lf = [0; 1; 0];
            sm.adc_mean = rand(nPat, nTp, 1) * 0.002;
            sm.d_mean = rand(nPat, nTp, 1) * 0.002;
            sm.f_mean = rand(nPat, nTp, 1) * 0.2;
            sm.dstar_mean = rand(nPat, nTp, 1) * 0.02;
            sm.d95_gtvp = rand(nPat, nTp) * 50;
            sm.dmean_gtvp = rand(nPat, nTp) * 45;

            calc = struct();

            cfg = struct();
            cfg.dataloc = testCase.TempDir;
            cfg.output_folder = testCase.TempDir;
            cfg.dwi_types_to_run = 1;

            % The function will try to call plot_parameter_maps which needs
            % real DICOM data. We expect it to error or warn during execution,
            % but it should at least parse its inputs without a nargin error.
            try
                evalc('visualize_results(data_vectors_gtvp, sm, calc, cfg);');
            catch ME
                % Expected to fail on missing DICOM data or sub-functions,
                % but NOT on argument count or field access of inputs.
                testCase.verifyFalse(contains(ME.message, 'Not enough input arguments'), ...
                    'Should not fail on argument count.');
                testCase.verifyFalse(contains(ME.identifier, 'MATLAB:narginchk'), ...
                    'Should not fail on narginchk.');
            end
        end

        function testOutputFolderCreated(testCase)
            % Verify that the function creates the output folder if it
            % doesn't exist (the mkdir call near the top of visualize_results).
            rng(42);
            nPat = 2;
            nTp = 2;

            data_vectors_gtvp = struct();
            for p = 1:nPat
                data_vectors_gtvp(p).id = sprintf('P%03d', p);
            end

            sm = struct();
            sm.id_list = {'P001', 'P002'};
            sm.mrn_list = sm.id_list;
            sm.lf = [0; 1];
            sm.adc_mean = rand(nPat, nTp, 1) * 0.002;
            sm.d_mean = rand(nPat, nTp, 1) * 0.002;
            sm.f_mean = rand(nPat, nTp, 1) * 0.2;
            sm.dstar_mean = rand(nPat, nTp, 1) * 0.02;
            sm.d95_gtvp = rand(nPat, nTp) * 50;
            sm.dmean_gtvp = rand(nPat, nTp) * 45;

            out_dir = fullfile(testCase.TempDir, 'new_output');
            cfg = struct();
            cfg.dataloc = testCase.TempDir;
            cfg.output_folder = out_dir;
            cfg.dwi_types_to_run = 1;

            try
                evalc('visualize_results(data_vectors_gtvp, sm, struct(), cfg);');
            catch
                % Expected to fail downstream
            end

            % The output folder should have been created by the mkdir call.
            testCase.verifyTrue(exist(out_dir, 'dir') == 7, ...
                'Output folder should be created by visualize_results.');
        end

    end
end
