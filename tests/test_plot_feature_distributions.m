classdef test_plot_feature_distributions < matlab.unittest.TestCase
    % TEST_PLOT_FEATURE_DISTRIBUTIONS Unit tests for plot_feature_distributions
    properties
        TempDir
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            % Create temp dir for output files
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);

            % Add core/utils to path
            addpath(fullfile(pwd, 'core'));
            addpath(fullfile(pwd, 'utils'));
            addpath(fullfile(pwd, 'dependencies'));

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
            % Should run without error and generate the two PNG files

            % Setup inputs
            dtype_label = 'Standard';
            dtype = 1;
            valid_pts = [1, 2, 3, 4, 5];
            lf_group = [0, 1, 0, 1, 0]; % 0 = LC, 1 = LF
            output_folder = testCase.TempDir;

            % Create mock data arrays: [patients x timepoints x dtypes]
            % We only need timepoint 1, dtype 1
            adc_mean = zeros(5, 1, 1);
            d_mean = zeros(5, 1, 1);
            f_mean = zeros(5, 1, 1);
            dstar_mean = zeros(5, 1, 1);

            % Populate with some random values
            rng(42); % Make deterministic
            adc_mean(:, 1, 1) = rand(5, 1) * 2e-3;
            d_mean(:, 1, 1) = rand(5, 1) * 2e-3;
            f_mean(:, 1, 1) = rand(5, 1) * 0.4;
            dstar_mean(:, 1, 1) = rand(5, 1) * 0.1;

            % Call the function
            plot_feature_distributions(dtype_label, adc_mean, d_mean, f_mean, dstar_mean, valid_pts, lf_group, dtype, output_folder);

            % Verify outputs
            hist_file = fullfile(output_folder, ['Feature_Histograms_' dtype_label '.png']);
            box_file = fullfile(output_folder, ['Feature_BoxPlots_' dtype_label '.png']);

            testCase.verifyTrue(exist(hist_file, 'file') > 0, ...
                'Expected Feature_Histograms_Standard.png to be created');
            testCase.verifyTrue(exist(box_file, 'file') > 0, ...
                'Expected Feature_BoxPlots_Standard.png to be created');
        end

        function testEmptyValidPts(testCase)
            % Test when valid_pts is empty

            dtype_label = 'Standard';
            dtype = 1;
            valid_pts = [];
            lf_group = [];
            output_folder = testCase.TempDir;

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
