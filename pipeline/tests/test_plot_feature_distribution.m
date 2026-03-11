classdef test_plot_feature_distribution < matlab.unittest.TestCase
    % TEST_PLOT_FEATURE_DISTRIBUTION Unit tests for plot_feature_distribution.
    % Tests histogram and boxplot modes, NaN handling, competing risk exclusion,
    % and edge cases.
    %
    % plot_feature_distribution renders a single feature's distribution split
    % by local failure (LF) group. These tests verify that the function runs
    % without error across various input scenarios; visual correctness is not
    % asserted (smoke tests only).

    properties
        OldVisible  % Stores the original DefaultFigureVisible setting for restoration
    end

    methods(TestMethodSetup)
        function suppressFigures(testCase)
            % Hide all figures during tests to prevent GUI popups.
            testCase.OldVisible = get(0, 'DefaultFigureVisible');
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function restoreFigures(testCase)
            % Close any figures created during the test and restore visibility.
            close all;
            set(0, 'DefaultFigureVisible', testCase.OldVisible);
        end
    end

    methods(Test)
        function test_histogram_basic(testCase)
            % Verify histogram mode runs without error on two well-separated
            % groups (LC=0, LF=1) with normally distributed values.
            vals = [randn(20,1); randn(20,1)+1];
            lf = [zeros(20,1); ones(20,1)];
            figure;
            plot_feature_distribution(vals, lf, 'ADC', 'mm^2/s', 'histogram');
            % Should not error
        end

        function test_boxplot_basic(testCase)
            % Verify boxplot mode runs without error on standard two-group input.
            vals = [randn(20,1); randn(20,1)+1];
            lf = [zeros(20,1); ones(20,1)];
            figure;
            plot_feature_distribution(vals, lf, 'ADC', 'mm^2/s', 'boxplot');
        end

        function test_nan_values_excluded(testCase)
            % NaN entries in vals should be silently excluded from the plot
            % without causing errors. Two NaNs are placed (one per group).
            vals = [1; 2; NaN; 4; 5; NaN; 7; 8];
            lf = [0; 0; 0; 0; 1; 1; 1; 1];
            figure;
            plot_feature_distribution(vals, lf, 'D', 'mm^2/s', 'boxplot');
        end

        function test_competing_risk_excluded(testCase)
            % Patients with lf==2 (competing risk / non-cancer death) should
            % be excluded from the visualization, leaving only LC (0) and LF (1).
            vals = [1; 2; 3; 4; 5; 6];
            lf = [0; 0; 1; 1; 2; 2];
            figure;
            plot_feature_distribution(vals, lf, 'f', '', 'boxplot');
        end

        function test_all_nan_input(testCase)
            % When every value is NaN, the function should handle gracefully
            % (empty plot or skip) rather than crash.
            vals = nan(10, 1);
            lf = [zeros(5,1); ones(5,1)];
            figure;
            plot_feature_distribution(vals, lf, 'D*', 'mm^2/s', 'histogram');
        end

        function test_single_value_histogram(testCase)
            % All identical values (zero variance) should not crash the
            % histogram binning logic.
            vals = ones(10, 1);
            lf = [zeros(5,1); ones(5,1)];
            figure;
            plot_feature_distribution(vals, lf, 'Constant', '', 'histogram');
        end

        function test_row_vector_input(testCase)
            % Row vectors (1xN) should be handled via internal shape
            % normalization; the function expects column vectors but should
            % not crash on row input.
            vals = [1 2 3 4 5 6 7 8];
            lf = [0 0 0 0 1 1 1 1];
            figure;
            plot_feature_distribution(vals, lf, 'ADC', 'mm^2/s', 'boxplot');
        end

        function test_single_group_boxplot(testCase)
            % Only one group present (all LC, lf==0). The boxplot should
            % still render without error even though there is no LF group.
            vals = randn(10, 1);
            lf = zeros(10, 1);
            figure;
            plot_feature_distribution(vals, lf, 'D', 'mm^2/s', 'boxplot');
        end

        function test_invalid_plot_type_errors(testCase)
            % An unsupported plot type ('scatter') should raise an error.
            % The empty string '' matches any error identifier.
            vals = randn(10, 1);
            lf = zeros(10, 1);
            figure;
            testCase.verifyError(@() plot_feature_distribution(vals, lf, 'X', '', 'scatter'), ...
                '');
        end

        function test_single_data_point_boxplot(testCase)
            % Edge case: only 1 valid data point remains after NaN removal
            % and competing risk exclusion. The boxplot should not crash.
            vals = [1; NaN; NaN];
            lf = [0; 0; 1];
            figure;
            plot_feature_distribution(vals, lf, 'D', 'mm^2/s', 'boxplot');
        end
    end
end
