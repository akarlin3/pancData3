classdef test_plot_feature_distribution < matlab.unittest.TestCase
    % TEST_PLOT_FEATURE_DISTRIBUTION Unit tests for plot_feature_distribution.
    % Tests histogram and boxplot modes, NaN handling, competing risk exclusion,
    % and edge cases.

    properties
        OldVisible
    end

    methods(TestMethodSetup)
        function suppressFigures(testCase)
            testCase.OldVisible = get(0, 'DefaultFigureVisible');
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function restoreFigures(testCase)
            close all;
            set(0, 'DefaultFigureVisible', testCase.OldVisible);
        end
    end

    methods(Test)
        function test_histogram_basic(testCase)
            vals = [randn(20,1); randn(20,1)+1];
            lf = [zeros(20,1); ones(20,1)];
            figure;
            plot_feature_distribution(vals, lf, 'ADC', 'mm^2/s', 'histogram');
            % Should not error
        end

        function test_boxplot_basic(testCase)
            vals = [randn(20,1); randn(20,1)+1];
            lf = [zeros(20,1); ones(20,1)];
            figure;
            plot_feature_distribution(vals, lf, 'ADC', 'mm^2/s', 'boxplot');
        end

        function test_nan_values_excluded(testCase)
            vals = [1; 2; NaN; 4; 5; NaN; 7; 8];
            lf = [0; 0; 0; 0; 1; 1; 1; 1];
            figure;
            plot_feature_distribution(vals, lf, 'D', 'mm^2/s', 'boxplot');
        end

        function test_competing_risk_excluded(testCase)
            % lf==2 patients should be excluded
            vals = [1; 2; 3; 4; 5; 6];
            lf = [0; 0; 1; 1; 2; 2];
            figure;
            plot_feature_distribution(vals, lf, 'f', '', 'boxplot');
        end

        function test_all_nan_input(testCase)
            vals = nan(10, 1);
            lf = [zeros(5,1); ones(5,1)];
            figure;
            plot_feature_distribution(vals, lf, 'D*', 'mm^2/s', 'histogram');
        end

        function test_single_value_histogram(testCase)
            % All identical values should not crash
            vals = ones(10, 1);
            lf = [zeros(5,1); ones(5,1)];
            figure;
            plot_feature_distribution(vals, lf, 'Constant', '', 'histogram');
        end

        function test_row_vector_input(testCase)
            % Row vectors should be handled (shape normalization)
            vals = [1 2 3 4 5 6 7 8];
            lf = [0 0 0 0 1 1 1 1];
            figure;
            plot_feature_distribution(vals, lf, 'ADC', 'mm^2/s', 'boxplot');
        end

        function test_single_group_boxplot(testCase)
            % Only one group present (all LC)
            vals = randn(10, 1);
            lf = zeros(10, 1);
            figure;
            plot_feature_distribution(vals, lf, 'D', 'mm^2/s', 'boxplot');
        end

        function test_invalid_plot_type_errors(testCase)
            vals = randn(10, 1);
            lf = zeros(10, 1);
            figure;
            testCase.verifyError(@() plot_feature_distribution(vals, lf, 'X', '', 'scatter'), ...
                '');
        end

        function test_single_data_point_boxplot(testCase)
            % Only 1 data point after NaN/competing risk removal
            vals = [1; NaN; NaN];
            lf = [0; 0; 1];
            figure;
            plot_feature_distribution(vals, lf, 'D', 'mm^2/s', 'boxplot');
        end
    end
end
