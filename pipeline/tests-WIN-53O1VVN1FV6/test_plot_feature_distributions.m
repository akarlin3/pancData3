classdef test_plot_feature_distributions < matlab.unittest.TestCase
    % TEST_PLOT_FEATURE_DISTRIBUTIONS Unit tests for plot_feature_distributions.
    %
    % plot_feature_distributions is the multi-parameter wrapper that creates
    % a 4-panel histogram figure and a 4-panel boxplot figure (one panel each
    % for ADC, D, f, D*), saved as PNG files. These tests verify that the
    % output files are created under both normal and edge-case conditions.

    properties
        TempDir  % Temporary directory for output PNG files
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            % Create temp dir for output files
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            % Add core/utils to path
            baseDir = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'dependencies'));

            % Ensure no figures are visible during testing to prevent popups
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function cleanup(testCase)
            % Close any left-over figures
            close all;

            % Remove temp dir
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testHappyPath(testCase)
            % Verify that plot_feature_distributions produces the expected
            % histogram and boxplot PNG files for a standard 5-patient cohort
            % with mixed LC/LF groups and realistic parameter ranges.

            % Setup inputs
            dtype_label = 'Standard';
            dtype = 1;                         % DWI type index (1 = Standard)
            valid_pts = [1, 2, 3, 4, 5];       % Patient indices to include
            lf_group = [0, 1, 0, 1, 0];        % 0 = LC (local control), 1 = LF (local failure)
            output_folder = testCase.TempDir;

            % Create mock data arrays: [patients x timepoints x dtypes]
            % We only need timepoint 1, dtype 1
            adc_mean = zeros(5, 1, 1);
            d_mean = zeros(5, 1, 1);
            f_mean = zeros(5, 1, 1);
            dstar_mean = zeros(5, 1, 1);

            % Populate with deterministic random values in clinically plausible ranges
            rng(42);
            adc_mean(:, 1, 1) = rand(5, 1) * 2e-3;    % ADC: 0-2e-3 mm^2/s
            d_mean(:, 1, 1) = rand(5, 1) * 2e-3;      % D: 0-2e-3 mm^2/s
            f_mean(:, 1, 1) = rand(5, 1) * 0.4;        % f: 0-0.4 (perfusion fraction)
            dstar_mean(:, 1, 1) = rand(5, 1) * 0.1;    % D*: 0-0.1 mm^2/s

            % Call the function
            plot_feature_distributions(dtype_label, adc_mean, d_mean, f_mean, dstar_mean, valid_pts, lf_group, dtype, output_folder);

            % Verify that both output PNG files were created
            hist_file = fullfile(output_folder, ['Feature_Histograms_' dtype_label '.png']);
            box_file = fullfile(output_folder, ['Feature_BoxPlots_' dtype_label '.png']);

            testCase.verifyTrue(exist(hist_file, 'file') > 0, ...
                'Expected Feature_Histograms_Standard.png to be created');
            testCase.verifyTrue(exist(box_file, 'file') > 0, ...
                'Expected Feature_BoxPlots_Standard.png to be created');
        end

        function testEmptyValidPts(testCase)
            % Edge case: when valid_pts is empty (no patients to plot), the
            % function should still create the PNG files (with empty axes)
            % rather than crashing.

            dtype_label = 'Standard';
            dtype = 1;
            valid_pts = [];
            lf_group = [];
            output_folder = testCase.TempDir;

            % Zero-row arrays simulate no patients
            adc_mean = zeros(0, 1, 1);
            d_mean = zeros(0, 1, 1);
            f_mean = zeros(0, 1, 1);
            dstar_mean = zeros(0, 1, 1);

            % Should run without crashing
            plot_feature_distributions(dtype_label, adc_mean, d_mean, f_mean, dstar_mean, valid_pts, lf_group, dtype, output_folder);

            hist_file = fullfile(output_folder, ['Feature_Histograms_' dtype_label '.png']);
            box_file = fullfile(output_folder, ['Feature_BoxPlots_' dtype_label '.png']);

            testCase.verifyTrue(exist(hist_file, 'file') > 0, ...
                'Expected Feature_Histograms_Standard.png to be created even with no data');
            testCase.verifyTrue(exist(box_file, 'file') > 0, ...
                'Expected Feature_BoxPlots_Standard.png to be created even with no data');
        end
    end
end
