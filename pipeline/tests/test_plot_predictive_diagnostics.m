classdef test_plot_predictive_diagnostics < matlab.unittest.TestCase
% TEST_PLOT_PREDICTIVE_DIAGNOSTICS — Unit tests for plot_predictive_diagnostics.m
%
% Validates ROC curve, sanity check panels, and 2D scatter plot generation:
%   - ROC plot saved when valid data exists
%   - Sanity check figures saved per feature
%   - 2D scatter plots generated for 2+ selected features
%   - Graceful handling of insufficient data (no valid ROC points)
%   - No crash on single-feature input (skip 2D scatter)
%   - Output files created in the correct folder

    properties
        OutputFolder
        OldFigVis
    end

    methods (TestMethodSetup)
        function addPaths(testCase)
            testDir = fileparts(mfilename('fullpath'));
            addpath(fullfile(testDir, '..', 'utils'));
            addpath(fullfile(testDir, '..', 'dependencies'));
        end

        function setupTempFolder(testCase)
            testCase.OutputFolder = fullfile(tempdir, ['test_ppd_' char(datetime('now','Format','yyyyMMddHHmmss'))]);
            mkdir(testCase.OutputFolder);
            % Suppress figure display during tests
            testCase.OldFigVis = get(0, 'DefaultFigureVisible');
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods (TestMethodTeardown)
        function cleanupTempFolder(testCase)
            set(0, 'DefaultFigureVisible', testCase.OldFigVis);
            if exist(testCase.OutputFolder, 'dir')
                rmdir(testCase.OutputFolder, 's');
            end
        end
    end

    methods (Test)
        function test_roc_plot_saved(testCase)
            % ROC figure should be saved when valid risk scores exist
            args = testCase.buildDefaultArgs(2);

            plot_predictive_diagnostics(args{:});

            roc_files = dir(fullfile(testCase.OutputFolder, 'ROC_OOF_*.png'));
            testCase.verifyGreaterThanOrEqual(numel(roc_files), 1, ...
                'ROC plot should be saved.');
        end

        function test_sanity_check_figures_saved(testCase)
            % One sanity check figure per selected feature
            n_features = 3;
            args = testCase.buildDefaultArgs(n_features);

            plot_predictive_diagnostics(args{:});

            sanity_files = dir(fullfile(testCase.OutputFolder, 'Sanity_Checks_*.png'));
            testCase.verifyEqual(numel(sanity_files), n_features, ...
                'Should produce one sanity check figure per feature.');
        end

        function test_2d_scatter_with_two_features(testCase)
            % 2D scatter plots should be generated for pairs of features
            args = testCase.buildDefaultArgs(2);

            plot_predictive_diagnostics(args{:});

            scatter_files = dir(fullfile(testCase.OutputFolder, '2D_Space_*.png'));
            testCase.verifyGreaterThanOrEqual(numel(scatter_files), 1, ...
                'Should produce 2D scatter plot for 2+ features.');
        end

        function test_single_feature_no_scatter(testCase)
            % Single feature should skip 2D scatter without error
            args = testCase.buildDefaultArgs(1);

            plot_predictive_diagnostics(args{:});

            scatter_files = dir(fullfile(testCase.OutputFolder, '2D_Space_*.png'));
            testCase.verifyEqual(numel(scatter_files), 0, ...
                'Single feature should not produce 2D scatter plots.');
        end

        function test_no_crash_all_nan_risk(testCase)
            % All-NaN risk scores should not crash (insufficient data path)
            args = testCase.buildDefaultArgs(1);
            % Replace risk scores with all NaN
            args{12} = nan(size(args{12}));

            % Should not error
            plot_predictive_diagnostics(args{:});
        end

        function test_competing_risk_exclusion(testCase)
            % Patients with lf_group==2 should be excluded from ROC
            args = testCase.buildDefaultArgs(2);

            plot_predictive_diagnostics(args{:});

            % Just verify no crash — internal exclusion logic
            roc_files = dir(fullfile(testCase.OutputFolder, 'ROC_OOF_*.png'));
            testCase.verifyGreaterThanOrEqual(numel(roc_files), 1);
        end

        function test_dstar_filename_sanitized(testCase)
            % D* feature name should have asterisk replaced with 'star' in filenames
            n = 30;
            rng(42);
            n_sig = 1;
            selected_indices = [4];  % D* index
            sig_names = {'D*'};
            sig_data = {randn(n, 3)};
            sig_is_abs = false;
            sig_is_pct = true;
            sig_disp = {'\Delta D*'};
            sig_units = {'%'};
            sig_col_idx = [3];
            sig_abs_data = {randn(n, 1)};
            sig_pct_data = {randn(n, 3)};

            valid_pts = true(n, 1);
            lf_group = [ones(10,1); zeros(15,1); 2*ones(5,1)];
            risk = randn(n, 1);
            m_gtv_vol = abs(randn(n, 3));
            adc_sd = abs(randn(n, 3, 3));
            ADC_abs = abs(randn(n, 3)) * 0.001 + 0.001;

            args = {selected_indices, n_sig, sig_data, sig_names, sig_is_abs, ...
                sig_is_pct, sig_disp, sig_units, sig_col_idx, ...
                sig_abs_data, sig_pct_data, ...
                risk, lf_group, valid_pts, ...
                m_gtv_vol, adc_sd, ADC_abs, ...
                3, 'Fx3', 'Standard', 1, testCase.OutputFolder, false};

            plot_predictive_diagnostics(args{:});

            % Verify no file contains literal '*' in name
            all_files = dir(fullfile(testCase.OutputFolder, '*.png'));
            for i = 1:numel(all_files)
                testCase.verifyFalse(contains(all_files(i).name, '*'), ...
                    'Filenames should not contain asterisk.');
            end
        end
    end

    methods (Access = private)
        function args = buildDefaultArgs(testCase, n_sig)
            % Build a consistent set of arguments for plot_predictive_diagnostics
            n = 30;
            rng(42);

            selected_indices = 1:n_sig;
            sig_names = arrayfun(@(i) sprintf('F%d', i), 1:n_sig, 'UniformOutput', false);
            sig_data = cell(1, n_sig);
            sig_abs_data = cell(1, n_sig);
            sig_pct_data = cell(1, n_sig);
            sig_is_abs = false(1, n_sig);
            sig_is_pct = true(1, n_sig);
            sig_disp = arrayfun(@(i) sprintf('\\Delta F%d', i), 1:n_sig, 'UniformOutput', false);
            sig_units = repmat({'%'}, 1, n_sig);
            sig_col_idx = 3 * ones(1, n_sig);  % target fraction column

            for i = 1:n_sig
                sig_data{i} = randn(n, 3);
                sig_abs_data{i} = randn(n, 1);
                sig_pct_data{i} = randn(n, 3);
            end

            valid_pts = true(n, 1);
            lf_group = [ones(10,1); zeros(15,1); 2*ones(5,1)];
            risk = randn(n, 1);
            m_gtv_vol = abs(randn(n, 3));
            adc_sd = abs(randn(n, 3, 3));
            ADC_abs = abs(randn(n, 3)) * 0.001 + 0.001;

            args = {selected_indices, n_sig, sig_data, sig_names, sig_is_abs, ...
                sig_is_pct, sig_disp, sig_units, sig_col_idx, ...
                sig_abs_data, sig_pct_data, ...
                risk, lf_group, valid_pts, ...
                m_gtv_vol, adc_sd, ADC_abs, ...
                3, 'Fx3', 'Standard', 1, testCase.OutputFolder, false};
        end
    end
end
