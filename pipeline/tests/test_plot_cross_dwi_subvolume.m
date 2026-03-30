classdef test_plot_cross_dwi_subvolume < matlab.unittest.TestCase
% TEST_PLOT_CROSS_DWI_SUBVOLUME  Tests for plot_cross_dwi_subvolume_comparison.
%
% Validates:
%   - Early return when fewer than 2 DWI types available
%   - Function accepts correct arguments without nargin error
%   - Graceful handling when checkpoint files are missing

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

        function testSkipsWhenNoCheckpoints(testCase)
            % When no DWI type checkpoint files exist, the function should
            % print a skip message and return without error.
            rng(42);
            nPat = 3;
            sm = struct();
            sm.id_list = {'P001', 'P002', 'P003'};
            sm.adc_sub_vol_pc = rand(nPat, 2, 3);

            cfg = struct();
            cfg.dataloc = testCase.TempDir;
            cfg.output_folder = testCase.TempDir;
            cfg.adc_thresh = 0.001;

            % Should not error — just prints skip message
            output = evalc('plot_cross_dwi_subvolume_comparison(sm, cfg);');
            testCase.verifyTrue(contains(output, 'skipped') || contains(output, 'need'), ...
                'Should report skipping when < 2 DWI types available.');
        end

        function testFunctionAcceptsTwoArgs(testCase)
            % Verify the function signature accepts (summary_metrics, config_struct).
            sm = struct();
            sm.id_list = {'P001'};
            sm.adc_sub_vol_pc = rand(1, 2, 3);

            cfg = struct();
            cfg.dataloc = testCase.TempDir;
            cfg.output_folder = testCase.TempDir;

            % Should not fail on argument count
            try
                evalc('plot_cross_dwi_subvolume_comparison(sm, cfg);');
            catch ME
                testCase.verifyFalse(contains(ME.message, 'Not enough input arguments'), ...
                    'Should not fail on argument count.');
            end
        end

        function testNoFiguresWhenInsufficient(testCase)
            % When fewer than 2 DWI types have data, no figures should be generated.
            fig_count_before = numel(findall(0, 'Type', 'figure'));

            sm = struct();
            sm.id_list = {'P001', 'P002'};
            sm.adc_sub_vol_pc = rand(2, 2, 3);

            cfg = struct();
            cfg.dataloc = testCase.TempDir;
            cfg.output_folder = testCase.TempDir;

            evalc('plot_cross_dwi_subvolume_comparison(sm, cfg);');

            fig_count_after = numel(findall(0, 'Type', 'figure'));
            testCase.verifyEqual(fig_count_after, fig_count_before, ...
                'No figures should be created when insufficient DWI types available.');
        end

    end
end
