classdef test_plot_dose_vs_delta < matlab.unittest.TestCase
    % TEST_PLOT_DOSE_VS_DELTA  Tests for plot_dose_vs_delta.m

    properties
        TempDir
        BaselineResults
        DosimetryResults
        ConfigStruct
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'core'));
            addpath(fullfile(repoRoot, 'utils'));
            if exist('OCTAVE_VERSION', 'builtin')
                addpath(fullfile(repoRoot, '.octave_compat'));
            end

            set(0, 'DefaultFigureVisible', 'off');

            [testCase.BaselineResults, testCase.DosimetryResults, testCase.ConfigStruct] = ...
                buildDoseVsDeltaData(testCase.TempDir);
        end
    end

    methods(TestMethodTeardown)
        function cleanup(testCase)
            diary off;
            close all;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testPngFilesCreated(testCase)
            % Verify that PNG files are created with the expected names.
            plot_dose_vs_delta(testCase.BaselineResults, testCase.DosimetryResults, ...
                testCase.ConfigStruct);

            dwi_types = testCase.ConfigStruct.dwi_types;
            for dt = 1:numel(dwi_types)
                for fx = [2, 3]
                    fname = sprintf('dose_vs_deltaD_Fx%d_%s.png', fx, dwi_types{dt});
                    full_path = fullfile(testCase.ConfigStruct.output_folder, fname);
                    testCase.verifyTrue(exist(full_path, 'file') == 2, ...
                        sprintf('Expected PNG not found: %s', fname));
                end
            end
        end

        function testMissingDosimetryHandled(testCase)
            % Should not error when dosimetry_results is an empty struct.
            empty_dosi = struct();
            testCase.verifyWarningFree(@() plot_dose_vs_delta( ...
                testCase.BaselineResults, empty_dosi, testCase.ConfigStruct));

            % PNG should still be created (left panel only)
            fname = sprintf('dose_vs_deltaD_Fx2_%s.png', testCase.ConfigStruct.dwi_types{1});
            full_path = fullfile(testCase.ConfigStruct.output_folder, fname);
            testCase.verifyTrue(exist(full_path, 'file') == 2, ...
                'PNG should still be created with missing dosimetry.');
        end

        function testAllNanDpctHandled(testCase)
            % Should not error when D_pct is all NaN.
            br = testCase.BaselineResults;
            br.D_pct(:) = NaN;
            testCase.verifyWarningFree(@() plot_dose_vs_delta( ...
                br, testCase.DosimetryResults, testCase.ConfigStruct));
        end

        function testMissingRequiredFieldHandled(testCase)
            % Should not error when required baseline fields are missing.
            br = struct();
            testCase.verifyWarningFree(@() plot_dose_vs_delta( ...
                br, testCase.DosimetryResults, testCase.ConfigStruct));
        end
    end
end


function [baseline, dosi, cfg] = buildDoseVsDeltaData(tempDir)
% Build 10-patient synthetic dataset: 3 LC, 3 LF, 2 CR, 2 NaN. 3 Tp, 2 DWI.
    rng(7);
    n_pts = 10;
    n_tp = 3;
    n_dwi = 2;

    D_pct = 10 * randn(n_pts, n_tp, n_dwi);
    m_d95_gtvp = 40 + 5 * randn(n_pts, n_tp);
    m_lf = [0; 0; 0; 1; 1; 1; 2; 2; NaN; NaN];

    baseline = struct();
    baseline.D_pct = D_pct;
    baseline.m_d95_gtvp = m_d95_gtvp;
    baseline.m_lf = m_lf;
    baseline.m_id_list = arrayfun(@(x) sprintf('P%02d', x), 1:n_pts, 'UniformOutput', false);

    dosi = struct();
    dosi.d95_d_sub = 38 + 6 * randn(n_pts, n_tp);

    cfg = struct();
    cfg.output_folder = tempDir;
    cfg.dwi_types = {'Standard', 'dnCNN'};
end
